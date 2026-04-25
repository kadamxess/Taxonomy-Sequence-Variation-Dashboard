# =====================================================
# IMPORTS
# =====================================================
import pandas as pd
import numpy as np
import plotly.express as px
import dash
from tqdm import tqdm
from dash import dcc, html
from Bio.Align import PairwiseAligner
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import warnings

warnings.filterwarnings("ignore")

# =====================================================
# CONFIGURATION
# =====================================================
FILE_PATH = r"C:\Users\it_ka\Desktop\kadam\BOLD_NEW\Code to filter and merge arthropoda data\clean_merged_output.xlsx"

IDENTITY_THRESHOLD = 0.97
MIN_SEQ_LEN = 50

MAX_SPECIES_HEATMAP = 3
TOP_SPECIES_BAR = 20
TOP_OTU = 10

# =====================================================
# VISUAL CONFIG (EDIT HERE ONLY)
# =====================================================
FIG_HEIGHT_SMALL = 450
FIG_HEIGHT_MEDIUM = 700
FIG_HEIGHT_LARGE = 700
FIG_WIDTH_NARROW = 900
COLOR_BLUE = px.colors.sequential.Blues
COLOR_RED = px.colors.sequential.Reds
COLOR_GREEN = px.colors.sequential.Viridis

# =====================================================
# LOAD & VALIDATE DATA
# =====================================================
df = pd.read_excel(FILE_PATH)
df.columns = df.columns.str.strip().str.lower()

REQUIRED_COLS = [
    "phylum", "class_name", "order_name",
    "family_name", "genus_name", "species",
    "country", "nucleotides"
]

missing = set(REQUIRED_COLS) - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

TAXON_COLS = [
    "phylum", "class_name", "order_name",
    "family_name", "genus_name", "species"
]

df[TAXON_COLS] = df[TAXON_COLS].fillna("Unclassified")
df["country"] = df["country"].fillna("Unknown").astype(str)
df["nucleotides"] = df["nucleotides"].fillna("").astype(str)

df = df[df["nucleotides"].str.len() >= MIN_SEQ_LEN]

# =====================================================
# PAIRWISE IDENTITY FUNCTION
# =====================================================
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = 0
aligner.extend_gap_score = 0

def identity(a, b):
    return aligner.score(a, b) / max(len(a), len(b))

# =====================================================
# SUNBURST FIGURE
# =====================================================
sunburst_df = (
    df.groupby(TAXON_COLS)
    .size()
    .reset_index(name="count")
)

sunburst_fig = px.sunburst(
    sunburst_df,
    path=TAXON_COLS,
    values="count",
    color="count",
    color_continuous_scale=COLOR_BLUE,
    title="Taxonomic Composition"
)

sunburst_fig.update_layout(
    height=FIG_HEIGHT_LARGE,
    title=dict(x=0.5, font=dict(size=24)),
    margin=dict(t=60, l=20, r=20, b=20)
)

# =====================================================
# SPECIES COUNT BAR (TOP + OTHERS)
# =====================================================
count_df = (
    df.groupby("species")
    .size()
    .reset_index(name="count")
    .sort_values("count", ascending=False)
)

top_df = count_df.head(TOP_SPECIES_BAR)
others_count = count_df.iloc[TOP_SPECIES_BAR:]["count"].sum()

if others_count > 0:
    top_df = pd.concat([
        top_df,
        pd.DataFrame([{"species": "Others", "count": others_count}])
    ])

count_fig = px.bar(
    top_df,
    x="species",
    y="count",
    color="count",
    title="Sequences per Species",
    color_continuous_scale=COLOR_GREEN
)

count_fig.update_layout(
    width=FIG_WIDTH_NARROW,   # ⬅ HERE
    height=FIG_HEIGHT_SMALL,
    xaxis_tickangle=45,
    margin=dict(t=60, b=80, l=60, r=40)
)

# =====================================================
# COUNTRY DISTRIBUTION MAP
# =====================================================
map_df = (
    df[df["country"] != "Unknown"]
    .groupby("country")
    .size()
    .reset_index(name="count")
)

map_fig = px.choropleth(
    map_df,
    locations="country",
    locationmode="country names",
    color="count",
    color_continuous_scale=COLOR_GREEN,
    title="Country-wise Sequence Distribution"
)

map_fig.update_layout(height=FIG_HEIGHT_MEDIUM)

# =====================================================
# INTRA-SPECIES VARIATION + HEATMAPS
# =====================================================
variation_records = []
heatmaps = []

for i, (species, g) in enumerate(df.groupby("species")):
    if len(g) < 2:
        continue

    seqs = g["nucleotides"].tolist()
    n = len(seqs)

    m = np.zeros((n, n))
    for x in range(n):
        for y in range(n):
            m[x, y] = identity(seqs[x], seqs[y])

    min_id = np.min(m[m > 0])

    variation_records.append({
        "species": species,
        "n_sequences": n,
        "min_identity": round(min_id, 4)
    })

    if i < MAX_SPECIES_HEATMAP:
        heatmaps.append(
            dcc.Graph(
                figure=px.imshow(
                    m,
                    title=f"Intra-species Identity: {species}",
                    color_continuous_scale=COLOR_GREEN
                ).update_layout(
                    height=400,
                    margin=dict(t=50, b=40, l=40, r=40)
                )
            )
        )

variation_df = pd.DataFrame(variation_records)

# =====================================================
# IDENTITY ALERTS
# =====================================================
alert_df = variation_df[
    variation_df["min_identity"] < IDENTITY_THRESHOLD
]
# convert to percentage 0–100
alert_df["min_identity_pct"] = (alert_df["min_identity"] * 100).round(2)
alert_fig = px.bar(
    alert_df.sort_values("min_identity_pct"),
    x="species",
    y="min_identity_pct",
    color="min_identity_pct",
    color_continuous_scale=COLOR_RED,
    title=f"Species Below {int(IDENTITY_THRESHOLD*100)}% Identity",
    text="min_identity_pct"
)
alert_fig.update_traces(width=0.4)   # value between 0 and 1
alert_fig.update_layout(
    width=FIG_WIDTH_NARROW,   # ⬅ HERE
    height=FIG_HEIGHT_SMALL,
    xaxis_tickangle=45,
    yaxis_title="Min identity (%)",  
)

# =====================================================
# OTU RE-CLUSTERING LOGIC
# =====================================================
def suggest_otus(seqs, threshold):
    if len(seqs) < 2:
        return 1

    n = len(seqs)
    m = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            m[i, j] = identity(seqs[i], seqs[j])
            m[j, i] = m[i, j]

    dist = 1 - m
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")
    clusters = fcluster(Z, t=(1 - threshold), criterion="distance")

    return len(set(clusters))

otu_records = []

for _, row in alert_df.iterrows():
    seqs = df[df["species"] == row["species"]]["nucleotides"].tolist()
    otu_records.append({
        "species": row["species"],
        "suggested_otus": suggest_otus(seqs, IDENTITY_THRESHOLD)
    })

otu_df = pd.DataFrame(otu_records)

otu_fig = px.bar(
    otu_df.sort_values("suggested_otus", ascending=False).head(TOP_OTU),
    x="species",
    y="suggested_otus",
    color="suggested_otus",
    color_continuous_scale=COLOR_RED,
    title="OTU Re-clustering Suggestions",
    text="suggested_otus"
)
otu_fig.update_traces(width=0.4)   # value between 0 and 1
otu_fig.update_layout(
    width=FIG_WIDTH_NARROW, 
    height=FIG_HEIGHT_SMALL,
    xaxis_tickangle=45
)

# =====================================================
# DASH APP
# =====================================================
app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("Taxonomy & Sequence Variation Dashboard"),

    dcc.Graph(figure=sunburst_fig),
    dcc.Graph(figure=count_fig),
    dcc.Graph(figure=map_fig),

    html.H2("Intra-species Genetic Variation"),
    html.Div(heatmaps),

    html.H2("Identity Threshold Alerts"),
    dcc.Graph(figure=alert_fig),

    html.H2("OTU Re-clustering Recommendations"),
    dcc.Graph(figure=otu_fig)
], style={"padding": "90px"})

# =====================================================
# RUN
# =====================================================
if __name__ == "__main__":
    app.run(debug=True)
