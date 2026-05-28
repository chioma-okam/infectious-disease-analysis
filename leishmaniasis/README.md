# Infectious Disease Analysis: Transcriptomics & GSEA Pipeline

A professional computational biology workflow dedicated to the differential expression and functional enrichment analysis of infectious disease transcriptomic datasets. The current pipeline focuses on downstream analysis of Visceral Leishmaniasis (VL) compared to asymptomatic cases using RNA-Seq data.

## 🧬 Project Overview
This repository contains a robust, end-to-end downstream bioinformatics pipeline built in R. It processes top-table differential expression data to extract biological insights regarding host immune responses, specifically contrasting active symptomatic disease states with asymptomatic infections.

* **Current Dataset Focus:** Visceral Leishmaniasis vs Asymptomatic (GSE77528)
* **Primary Objective:** Identify key regulatory genes, enriched Gene Ontology (GO) terms, and disrupted KEGG pathways to understand the mechanisms of disease progression and host tolerance.

---

## 🛠️ Key Technical Features
The pipeline is fully vectorized, modular, and leverages standard Bioconductor packages for genomic data science:

* **Data Preprocessing & Filtering:** Automated handling of missing gene symbols, removing duplicates, and ensuring strict data-type compliance for quantitative variables (e.g., Log2 Fold Change).
* **Gene Set Enrichment Analysis (GSEA):** Implements `clusterProfiler` to execute un-biased functional enrichment.
* **Advanced Visualizations:** Generates publication-ready figures using `enrichplot`, including:
    * **Volcano Plots** for statistical significance vs. effect size mapping.
    * **Ridge Plots** for visualizing density distributions of enriched categories.
    * **Dot Plots** to clearly rank the most statistically significant biological processes.
    * **Pathway Mapping:** Annotation and extraction of core enrichment genes from key immune signaling pathways.

---

## 🧰 Tech Stack & Libraries
This project is built using the **R Programming Language** and the following core bioinformatics libraries:
* `clusterProfiler` - For statistical analysis of functional profiles and gene cluster comparisons.
* `org.Hs.eg.db` - Genome-wide annotation for Human database.
* `enrichplot` - For visualization of functional enrichment results.
* `readxl` & `readr` - Robust data import workflows for `.xlsx` and `.tsv` files.
* `writexl` - Clean data export pipelines.

---

## 📁 Repository Structure
```text
├── leishmaniasis/
│   ├── Leishmaniasis_Complete_GSEA_Pipeline.R   # Main analysis and plotting script
│   └── README.md                                # Project documentation
