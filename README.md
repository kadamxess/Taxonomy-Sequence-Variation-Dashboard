#  Taxonomy & Sequence Variation Dashboard

##  Overview

This project is an **interactive bioinformatics dashboard** built using **Dash and Plotly** to analyze DNA sequence datasets (e.g., BOLD COI data).

It integrates **taxonomy visualization, sequence quality analysis, intra-species variation, and OTU re-clustering** into a single analytical pipeline.

The dashboard is designed for:

* Biodiversity informatics
* DNA barcoding studies
* Sequence quality validation
* Species delimitation analysis

---

##  Key Features

###  Taxonomic Analysis

* Hierarchical **sunburst visualization** from phylum → species
* Distribution of sequences across taxa

###  Species-Level Insights

* Top species by sequence count (with “Others” aggregation)
* Helps identify dataset bias or overrepresentation

###  Geographic Distribution

* Country-wise sequence distribution using a **choropleth map**

###  Genetic Variation Analysis

* Pairwise **sequence identity calculation**
* Intra-species variation heatmaps
* Detects sequence divergence within species

###  Quality Alerts

* Flags species with identity below threshold (default: **97%**)
* Useful for:

  * Misidentification detection
  * Data quality issues

###  OTU Re-clustering

* Suggests number of **Operational Taxonomic Units (OTUs)**
* Uses hierarchical clustering based on sequence similarity

---

##  Input Requirements

The input file must be an **Excel file (.xlsx)** containing:

### Required Columns

* `phylum`
* `class_name`
* `order_name`
* `family_name`
* `genus_name`
* `species`
* `country`
* `nucleotides`

### Notes

* Missing taxonomy values are replaced with `"Unclassified"`
* Missing country values are replaced with `"Unknown"`
* Sequences shorter than **50 bp** are removed

---

##  Methodology

### 1. Sequence Filtering

* Removes sequences below minimum length (`MIN_SEQ_LEN = 50`)

### 2. Pairwise Identity Calculation

* Uses **global alignment** via `Bio.Align.PairwiseAligner`
* Identity defined as:

```
identity = alignment_score / max(sequence_length)
```

### 3. Variation Analysis

* Computes identity matrix per species
* Extracts **minimum intra-species identity**

### 4. Threshold-Based Alerts

* Species flagged if:

```
min_identity < 0.97
```

### 5. OTU Clustering

* Converts identity → distance:

```
distance = 1 - identity
```

* Applies hierarchical clustering (`average linkage`)
* Suggests OTU count based on threshold

---

##  Dashboard Components

| Section       | Description                       |
| ------------- | --------------------------------- |
| Sunburst Plot | Taxonomic hierarchy visualization |
| Bar Chart     | Top species by sequence count     |
| Map           | Country-wise distribution         |
| Heatmaps      | Intra-species similarity matrices |
| Alert Chart   | Species below identity threshold  |
| OTU Chart     | Suggested OTU counts              |

---

##  Installation

### 1. Clone Repository

```bash
git clone https://github.com/your-username/sequence-dashboard.git
cd sequence-dashboard
```

### 2. Install Dependencies

```bash
pip install pandas numpy plotly dash tqdm biopython scipy openpyxl
```

---

##  Usage

### 1. Update File Path

Edit in script:

```python
FILE_PATH = "path_to_your_excel_file.xlsx"
```

### 2. Run the Application

```bash
python app.py
```

### 3. Open Dashboard

Navigate to:

```
http://127.0.0.1:8050/
```

---

##  Configuration

You can easily tune parameters:

```python
IDENTITY_THRESHOLD = 0.97   # species similarity cutoff
MIN_SEQ_LEN = 50            # minimum sequence length
MAX_SPECIES_HEATMAP = 3     # number of heatmaps
TOP_SPECIES_BAR = 20        # top species displayed
TOP_OTU = 10                # OTU suggestions shown
```

---

##  Project Structure

```
project/
│── app.py
│── data/
│   └── clean_merged_output.xlsx
│── README.md
```

---

##  Performance Notes

* Pairwise alignment is **O(n²)** per species
* Large datasets may be slow for:

  * High sequence counts per species
* Recommended:

  * Pre-filter large datasets
  * Limit sequences per species if needed

---

##  Use Cases

* DNA barcoding validation
* Species delimitation studies
* Biodiversity dataset exploration
* Detecting mislabeled sequences
* Preprocessing for phylogenetic analysis

---

##  Author

**Kadam Xess**
Bioinformatics & Clinical Research

---

##  License

This project is licensed under the MIT License.

---

##  Future Improvements

* Add FASTA input support
* Parallelize alignment for speed
* Export reports (PDF/CSV)
* Integrate phylogenetic tree visualization

---
