# UPS-TF Correlation Analysis

The Ubiquitin-Proteasome System (UPS) regulates proteins, but generally there is very limited measurements of the proteome. To understand the UPS in cancer, this project intends to find whether mutated genes in the 
UPS may modulate transcription factor (TF) activity, since TF activity can be inferred by the RNA expression of their target genes. 

## Installation

Install the dependencies using conda:

```bash
$ conda env create -f environment.yaml 
```

Then activate the environment:

```bash
$ source activate UPS
```

Next, you need to install rabit (see http://rabit.dfci.harvard.edu/download/).

## Data

Please download tumor purity estimates, gene expression and mutation calls from the Genomic Data Commons (GDC) website (https://gdc.cancer.gov/about-data/publications/pancanatlas):
* Mutations (mc3.v0.2.8.PUBLIC.maf.gz)
* Gene expression (EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv)
* Tumor purity (TCGA_mastercalls.abs_tables_JSedit.fixed.txt)

Next, filter the mutations (MAF file) by following instructions here (https://www.dropbox.com/sh/wglgggbgketh982/AABJEqQ2QdCEruy9c6UXBdjba?dl=0).

Now, update the configuration file (config.yaml) to point towards the mutation, gene expression and tumor purity files.

Lastly, you will also need to download the cistrome interaction and background file for Rabit.

I downloaded copy number from GDC on 4/4/2019.

## Commands

The first step is to prepare/format the mutation and gene expression data.

```bash
$ snakemake preprocess
```
