# TYK2 mediates neuroinflammation in dementia with Alzheimer’s disease

This repository contains code to reproduce the findings of [TYK2 mediates neuroinflammation in dementia with Alzheimer’s disease](https://www.biorxiv.org/content/10.1101/2024.06.04.595773v1). It uses an updated version of the [Drug Repurposing In Alzheimer's Disease (DRIAD)](https://www.nature.com/articles/s41467-021-21330-0) prediction pipeline, designed to identify candidates for drug repurposing in Alzheimer's disease.

* Primary DRIAD repository: https://github.com/labsyspharm/DRIAD
* DRIAD as a webapp: https://labsyspharm.shinyapps.io/DRIAD/

## Reproducing the analysis

Fully reproducing the analysis requires access to the [AMP-AD consortium](https://adknowledgeportal.synapse.org/Explore/Programs/DetailsPage?Program=AMP-AD) dataset. Because of data sharing restrictions we do not make this dataset available in this repository directly. However, we include a small synthetic example dataset in the `data` directory that can be used to test the pipeline.

### Step 1: Install dependencies

For directly downloading the necessary datasets the non-CRAN [synapser R package](https://github.com/Sage-Bionetworks/synapser) and a free [Synapse user account](https://www.synapse.org) are required. Installing all dependencies should take around 10min.

```r
install.packages("synapser", repos=c("http://ran.synapse.org", "https://cloud.r-project.org"))
```

Additionally the following R packages are required:

```r
install.packages(
  c("tidyverse", "data.table", "powerjoin", "here", "qs", "remotes", "BiocManager")
)
remotes::install_github("ArtemSokolov/synExtra")
remotes::install_github("labsyspharm/DRIAD")
remote::install_github("labsyspharm/ordinalRidge")
BiocManager::install("tximport")
```

### Step 2: Download the data

The pipeline requires two datasets:

1. The AMP-AD consortium RNA-seq dataset containing RNA-seq data from Alzheimer's disease patients at different stages.
2. Drug perturbation RNA-seq data for the compounds that are to be investigated for repurposing.

#### Downloading the AMP-AD dataset

We use a subset of the ROSMAP dataset. It can be downloaded using the following steps:

```bash
# Compile a list of Synapse IDs. Creates 'selected_rosmap_samples.csv'
Rscript download/01_select_rosmap_samples.R

# Downloads the data
bash job_scripts/download_synapse.sh selected_rosmap_samples.csv 2
```

#### Download the drug perturbation datasets

The drug perturbation gene signatures are located on Synapse at [syn20820056](https://www.synapse.org/#!Synapse:syn20820056) and [syn20820060](https://www.synapse.org/#!Synapse:syn20820060).

They are automatically downloaded in the `prepare_rosmap_tasks.R` script.

#### Download RNA-seq data from TDP-43 sorted neurons

This data was used to validate the quantification of TDP-43-loss induced cryptic exons ([Liu EY et al. Loss of Nuclear TDP-43 Is Associated with Decondensation of LINE Retrotransposons. Cell Rep 2019](https://www.ncbi.nlm.nih.gov/pubmed/31042469)). It is available from GEO at [GSE126542](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126542).

```bash
bash job_scripts/GSE126542_sra_download.sh
```

### Step 3: Quantify gene expression

RNA-seq samples from the ROSMAP dataset are quantified using Salmon. Some of the samples are only available as bam files, which first must be converted to fastq using `job_scripts/bam_to_fastq.sh`.

Salmon quantification is performed using `job_scripts/salmon.sh` using a custom index with additional transcripts containing TDP-43 pathology-induced cryptic exons. The index is created using the `job_scripts/make_salmon_index.sh` script.

Once Salmon quantification is complete, the samples are combined into a single count matrix using `driad/prep_rosmap_counts.R`.

### Step 4: Run the DRIAD pipeline

DRIAD background gene set and drug gene set tasks are created using `driad/prepare_rosmap_tasks.R`. Note that this step requires a HPC cluster with a job scheduler.

### Step 5: Analyze the results

The results of the DRIAD pipeline can be analyzed using the `driad/driad_plots.R` script. This script generates plots showing the performance of the drug repurposing candidates in TDP-43+ and TDP-43- patient subgroups.

## Example dataset

A small synthetic example dataset is available at [syn63663018](https://www.synapse.org/#!Synapse:syn63663018). It was generated using `download/synthetic_data.Rmd`. It can be used in lieu of the AMP-AD dataset patient sample RNA-seq data to test the pipeline.

## Funding

We gratefully acknowledge support by NIA grant R01 AG058063, RF1 AG078297, and NCI U54 CA225088.
