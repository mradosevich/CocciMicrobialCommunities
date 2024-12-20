# CocciMicrobialCommunities
Repository for code and data used to analyze the soil microbial communities associated with Coccidioides immitis

Raw 16S and ITS2 sequencing files for this project can be found on the NCBI Sequence Read Archives under accession numbers #### and ####. 

This repository contains the following data folders and files: 

1. A folder with quality control outputs for the ITS2 sequencing data
2. A folder with quality control outputs for the 16S sequencing data
3. A folder with all ITS2 QIIME2 objects
4. A folder with all 16S QIIME2 objects
5. tsv table with all ITS2 sequences matched to fungal taxa using the UNITE database
6. tsv table with all 16S sequences matched to bacterial taxa using the SILVA database

And the following code scripts: 
1. A shell script for taking the raw fasta files through the QIIME2 pipeline using the command-line interface
2. An R script for analyzing ITS2 alpha diversity within each sample
3. An R script for analyzing 16S alpha diversity within each sample
4. An R script for analyzing ITS2 beta diversity across the sample set
5. An R script for analyzing 16S beta diversity across the sample set
6. An R script for analyzing fungal co-occurrence patterns with Coccidioides immitis
7.  An R script for analyzing bacterial co-occurrence patterns with Coccidioides immitis
