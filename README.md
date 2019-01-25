% Whole-genome genotyping in experimental pedigrees from outbred founders utilizing low coverage individual based sequencing

Carlborg's Lab (Yanjun Zan, Thibaut Payen, Örjan Carlborg) - BMC - Uppsala University

## Problem solved

Short introduction

## Installation

- If you use conda ( https://conda.io/projects/conda/en/latest/user-guide/install/index.html ) the requirements.txt in the project will help you deploy it, if not you should have:
  - Python3
  - bcftools
  - Snakemake
- In any case some things are considered installed on the computer here, if it is not the case you must install them via conda or accessible by the user in an other way:
  - Java (for TIGER)
  - Perl (for TIGER)
  - R

- Download git clone Stripes
- Command to install the environnement


```
git clone https://github.com/CarlborgGenomics/Stripes
cd Stripes
conda create --name Stripes --file requirements.txt
source activate Stripes
```


## Usage

The first step is to give to the pipeline the folders containing your data in the configuration file named config.yaml.

This pipeline need **at least**:

- a pedigree file,
- a VCF of the founders,
- a table containing informations about the scaffolds (see example),
- and the way to obtain BAM files for the F2 (either directy the BAMs or FASTQs + the reference genome)

Once the configuration file is set, launch the pipeline with:
```
snakemake -k -s scripts/Snakefile
```

NB: Creating the graphe can take some time (more than 10 minutes for 900 individuals is still in the reasonnable range)

## Cite

If you are using this pipeline in a scientific publication please cite:

Zan, Yanjun, Thibaut Payen, Mette Lillie, Christa Honaker, Paul B Siegel, and Örjan Carlborg. “Genotyping by Low-Coverage Whole-Genome Sequencing in Intercross Pedigrees from Outbred Founders: A Cost Efficient Approach.” BioRxiv, January 1, 2018, 421768. <https://doi.org/10.1101/421768>.


As the inputation is made by TIGER also cite:

Rowan, B. A., Patel, V., Weigel, D., & Schneeberger, K. (2015). Rapid and inexpensive whole-genome genotyping-by-sequencing for crossover localization and fine-scale genetic mapping. G3: Genes, Genomes, Genetics, g3-114 (<https://doi.org/10.1534/g3.114.016501> )


## Pipeline steps and files description

### 1 fastq2vcf

TODO

### 2 vcf2individualGenotypes

TODO

### 3 individualGenotypes2breaks

copied from <https://github.com/betharowan/TIGER_Scripts-for-distribution> with one modification to take the scaffold name correctly into account. The README contained in this folder explain the different steps for the inputation that have been automated and linked in the general Snakefile.

### 4 figures

To obtain the figures in the paper see Stripes_downstream

## Comparison with existing methods

Why use Stripes ?

A lot of methods exists if you don't have a pedigree or for outbred populations, like using TIGER directy for example
