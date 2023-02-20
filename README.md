# **Jornada_Creosote_Tarbush_Fungal_Metabarcoding**

pipeline for metabarcoding of soil README

This pipeline is divided into 7 sections. Each section has a filename and a step. 

## Hardware Configuration
Ran using New Mexico State University's Discovery HPC Cluster

Operating System: CentOS 7

Scheduler: Slurm 21.08.4

## Build Conda Environments: 

File name: Conda_environment_qc

Step: Build the 'qc' Conda environment. Download and use the 'qc.conda.yml' file for this step.

File name: Conda_environment_qiime2_2021_11

Step: Build the 'qiime2_2021_11' Conda environment. Download and use the 'qiime2_2021.11-custom.yml' file for this step.

## Section 1

File name: 1_Initial_quality_reports

Step: Generate quality reports of raw paired-end sequences.

## Section 2

File name: 2_Extract_ITS2_region

Step: Use QIIME2 to extract the ITS2 region from the paired-end sequences.

## Section 3

File name: 3_Post_quality_reports

Step: Generate quality reports of post-trimmed paired-end sequences.

## Section 4 

File name: 4_Taxonomic_tables

Step: Use RStudio to process paired-end sequences and generate a taxonomic table.

## Section 5

File name: 5_Format_ASV_table_for_FUNGuild

Step: Use RStudio to format the 'otu.tax.table.csv' for input in FUNGuild.

## Section 6

File name: 6_Assignments_on_FUNGuild

Step: Use command line to make functional assignments using FUNGuild.

## Section 7

File name: 7_Format_table_for_stat_analyses

Step: Use RStudio to format the output table from FUNGuild for statistical analyses.