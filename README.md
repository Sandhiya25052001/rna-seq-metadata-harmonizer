# rna-seq-metadata-harmonizer
RNA-seq Metadata Harmonizer – A reproducible pipeline for curating, harmonizing, and integrating RNA-seq study metadata. This project automates the extraction of sample information from multiple GEO datasets, assigns experimental conditions (e.g., control vs treated), and generates a unified metadata table for downstream analysis.
# 🧬 RNA-seq Metadata Harmonization Tool

## 📌 Overview
Different RNA-seq studies often provide metadata in inconsistent formats, making it difficult to perform cross-study analysis.  
This project provides a **Python-based metadata harmonization tool** that standardizes sample annotations across multiple studies.  

It:
- ✅ Reads metadata from GEO series matrix files (`.txt.gz`)  
- ✅ Extracts sample IDs and treatment/control information  
- ✅ Handles inconsistent annotation formats (e.g., "control", "normal", "tumor")  
- ✅ Performs quality checks (missing values, sample counts)  
- ✅ Produces harmonized metadata in **CSV format**  
- ✅ Includes summary plots of sample distributions  

---


