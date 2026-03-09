
# Project Documentation: Epigenetic Clock Analysis
Repository contains the code I did for my systematic review (data analysis of technical methods used for epigenetic clock development)
## Requirements

**Python Version:** 3.9.7

**Required Packages:**

```
pandas
matplotlib
seaborn
numpy
collections (Counter, defaultdict)
```

---

## Data Requirements

The notebooks and scripts expect the following data files:
- `supplementary_table_2.xlsx` - Publication list (Excel)
- `supplementary_table_3.xlsx` - GEO dataset counts by year (Excel)
- `supplementary_table_1.xlsx` - Epigenetic clock data (Excel)

---

## File Descriptions

### 1. **analysis.py**
**Purpose:** Analyze and visualize publication data for epigenetic clock systematic review

**Main Classes:**
- `PublicationsSummary`: Summarizes publications metadata
  - Counts human vs. non-human studies
  - Lists all species studied
  - Creates timeline visualization of epigenetic clocks and GEO datasets over years 

- `MethylationMeasurement`: Analyzes methylation measurement methods
  - Categorizes and counts different methylation techniques (array, pyrosequencing, sequencing, etc.) according to animal model (human vs non-human models)

**Key Outputs:**
- Bar visualizations
- Counts of animal models used
- Breakdown of methylation methods by type for the 2 main types of animal models used 

**analysis.ipynb:**
- Workflow notebook that reads analysis.py

---

### 2. **comput_methods_analysis.py**
**Purpose:** Analyze computational methods used in epigenetic clock development

**Main Classes:**
- `ImputationMethods`: Analyzes imputation strategies
  - Counts studies using imputation vs. not using imputation
  - Separates analysis by human and non-human studies

- `FeatureSelectionMethods`: Categorizes feature selection approaches
  - Categorizes methods (penalized regression, tree-based, correlation-based, etc.)
  - Generates bar charts showing frequency of each category
  - Separate analysis for human and non-human clocks

- `AgePredictionMethods`: Analyzes modeling approaches
  - Categorizes prediction models (linear, neural networks, SVMs, tree-based, etc.)
  - Counts usage frequency of each model type

**Key Outputs:**
- Bar charts of feature selection frequencies
- Model category distributions
- Imputation method proportions
- PDF exports with titles specifying "human" or "non_human"

**comput_methods.ipynb:**
- Workflow notebook that reads comput_methods_analysis.py


### 3. **try.notebook**
- Jupyter notebook file generating the Imputation bar plot
---
