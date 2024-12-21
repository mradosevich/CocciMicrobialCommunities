# CocciMicrobialCommunities
Repository for code and data used to analyze the soil microbial communities associated with Coccidioides immitis

Raw 16S and ITS2 sequencing files for this project can be found on the NCBI Sequence Read Archives under accession numbers #### and ####. 

This repository contains the following data files: 

1. feature-table outputs from qiime2 (.tsv tables with columns as samples and rows as fungal or bacterial taxa, with values corresponding to sequence read counts)
2. a .csv file with all sample metadata used in analyses
3. two .csv files matching sample barcodes in the metadata file to sequencing library filenames in the feature-tables

And the following code scripts: 
1. A shell script for taking the raw fasta files through the QIIME2 pipeline using the command-line interface
2. An R script for analyzing ITS2 alpha diversity within each sample
3. An R script for analyzing 16S alpha diversity within each sample
4. An R script for analyzing ITS2 beta diversity across the sample set
5. An R script for analyzing 16S beta diversity across the sample set
6. An R script for analyzing fungal co-occurrence patterns with Coccidioides immitis
7.  An R script for analyzing bacterial co-occurrence patterns with Coccidioides immitis
