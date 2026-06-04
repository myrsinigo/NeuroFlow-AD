# NeuroFlow-AD

## Reproducible Single-Cell Alzheimer's Disease Analysis Pipeline

NeuroFlow-AD is a reproducible neuroinformatics workflow for the analysis of single-cell and single-nucleus transcriptomics data from Alzheimer's disease studies.

The pipeline combines Nextflow, Docker, and Scanpy to perform quality control, normalization, dimensionality reduction, clustering, marker-gene identification, and visualization in a fully reproducible environment.

This project was developed as a portfolio project at the intersection of bioinformatics, neuroscience, and scientific data engineering.

---

## Project Overview

Single-cell transcriptomics has become a key technology for understanding cellular heterogeneity in neurodegenerative diseases such as Alzheimer's disease.

This workflow enables:

* Quality control of single-cell datasets
* Gene expression normalization
* Highly variable gene selection
* Dimensionality reduction (PCA/UMAP)
* Leiden clustering
* Marker-gene identification
* Reproducible execution using Docker and Nextflow

---

## Dataset

This project was tested using data from the Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD).
(link for downloading the data: https://datasets.cellxgene.cziscience.com/92b37feb-aa2c-40d7-bd90-0a9b5ddb3b27.h5ad) 

SEA-AD is a large-scale human brain atlas containing over 1.3 million cells and provides an important resource for studying cellular changes associated with Alzheimer's disease.

For computational efficiency, a 5,000-cell subset was generated and analyzed using the following workflow.



## Workflow

<p align="center">
  <img src="docs/workflow.png" width="900" alt="Single-cell RNA-seq workflow">
</p>

---

## Technologies

### Bioinformatics

* Scanpy
* AnnData
* Single-cell RNA-seq analysis

### Programming

* Python

### Workflow Orchestration

* Nextflow

### Reproducibility

* Docker
* Conda

---

## Results

Dataset subset analyzed:

* 5,000 cells
* 35,483 genes

Pipeline outputs:

* Processed AnnData object
* UMAP visualization
  
<img width="1794" height="429" alt="umap_clusters_qc" src="https://github.com/user-attachments/assets/dbb113b6-cd75-40f7-ac64-2173a29629fa" />

* QC reports
  
<img width="1508" height="490" alt="violin_qc" src="https://github.com/user-attachments/assets/57699894-5fda-42b9-8c97-9f093fac880e" />

* Marker gene tables

Example analysis identified:

* 23 transcriptionally distinct clusters

---

## Repository Structure

```text
.
├── main.nf
├── nextflow.config
├── Dockerfile
├── environment.yml
├── bin/
│   └── run_scanpy_pipeline.py
├── data/
├── results/
└── README.md
```

---

## Running the Pipeline

Build Docker image:

```bash
docker build -t ad-singlecell:latest .
```

Run workflow:

```bash
nextflow run main.nf \
  --input data/sea_ad_subset_5000.h5ad \
  --outdir results \
  -profile docker
```

---

## Future Development

Planned extensions include:

* Cell-type annotation
* Differential expression analysis
* Alzheimer's vs control comparisons
* Automated report generation
* Cloud-native execution on AWS

