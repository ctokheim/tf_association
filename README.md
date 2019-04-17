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

Next, you need to install rabit (see http://rabit.dfci.harvard.edu/download/). Then update the configuration file (config.yaml) to point to where you installed rabit.

## Data

Please download tumor purity estimates, gene expression and mutation calls from the Genomic Data Commons (GDC) website (https://gdc.cancer.gov/about-data/publications/pancanatlas):
* Mutations (mc3.v0.2.8.PUBLIC.maf.gz)
* Gene expression (EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv)
* Tumor purity (TCGA_mastercalls.abs_tables_JSedit.fixed.txt)

Next, filter the mutations (MAF file) by following instructions [here](https://www.dropbox.com/sh/wglgggbgketh982/AABJEqQ2QdCEruy9c6UXBdjba?dl=0).

Now, update the configuration file (config.yaml) to point towards the mutation, gene expression and tumor purity files.

Lastly, you will also need to download the cistrome interaction and background file for [Rabit](https://www.dropbox.com/sh/l64jxw8ucwuiov6/AADkuapIpTvk9vniqqWkBBONa?dl=0). Update the config file to point to the rabit data files.

## Command

Once you have everything set up, you should be able to run the entire pipeline with one snakemake 
command.

```bash
$ snakemake -p -k --cores=10 --config output=output
```

The above command will save results into the directory specified by the "output" parameter. Additionally, you can specify the number of cores to use with the `--cores` parameter. One thing to note is that the Statsmodels package tends to use all cores that are available on a machine for each thread. To ensure each job uses only one core, you will need to set the "OMP_NUM_THREADS" variable:

```bash
$ OMP_NUM_THREADS=1 snakemake -p -k --cores=10 --config output=output
```

If you want to run on a computer cluster, please see the snakemake documentation (https://snakemake.readthedocs.io/en/stable/).
