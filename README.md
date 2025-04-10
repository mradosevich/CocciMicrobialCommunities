# CocciMicrobialCommunities
Repository for code and data used to analyze the soil microbial communities associated with Coccidioides immitis

Raw 16S and ITS2 sequencing files for this project can be found on the NCBI Sequence Read Archives under Accession Numbers PRJNA1201319 (16S) and PRJNA1201328 (ITS2). 

This repository contains the following data files: 

1. feature-table outputs from qiime2 (.tsv tables with columns as samples and rows as fungal or bacterial taxa, with values corresponding to sequence read counts)
2. a .csv file with all sample metadata used in analyses
3. two .csv files matching sample barcodes in the metadata file to sequencing library filenames in the feature-tables
4. taxonomy tables for ITS2 data
5. taxonomy tables for 16S data

And the following code scripts: 
1. Shell scripts for taking the raw fasta files through the QIIME2 pipeline using the command-line interface
2. An R script for preprocessing feature-table.tsv files to create taxonomy tables
3. An R script for generating alpha diversity metrics
4. An R script for generating beta diversity metrics
5. R scripts for analyzing fungal and bacterial co-occurrence patterns with Coccidioides immitis
