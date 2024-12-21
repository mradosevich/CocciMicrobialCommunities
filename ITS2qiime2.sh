##Start QIIME2 environment:
module load qiime2/2021.8

##CasavaOneEightSingleLanePerSampleDirFmt - the format of data that comes off the illumina platform. 

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format CasavaOneEightSingleLanePerSampleDirFmt --input-path /fastq --output-path demux.qza

##to visualize the data
qiime demux summarize --i-data demux.qza --p-n 10000 --o-visualization demux.qzv

##cutadapt to trim primers
qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza --p-cores 16 --p-front-f 'AACTTTYRRCAAYGGATCWCT' --p-front-r 'AGCCTCCGCTTATTGATATGCTTAART' --p-match-adapter-wildcards --p-match-read-wildcards --p-discard-untrimmed --o-trimmed-sequences demux.trimmed.qza --verbose

##visualize trimmed data
qiime demux summarize --i-data demux.trimmed.qza --p-n 10000 --o-visualization demux.trimmed.qzv


##Make an OTU table


##To denoise and merge
qiime dada2 denoise-paired --i-demultiplexed-seqs demux.trimmed.qza --p-trunc-len-f 269 --p-trunc-len-r 207 --p-chimera-method consensus --p-pooling-method independent --p-n-reads-learn 1000000 --p-n-threads 0 --o-table demux.table.qza --o-representative-sequences demux.rep-seqs.qza --o-denoising-stats demux.stats.qza --verbose


##Turning the stats into qzv objects for visualization:
qiime feature-table summarize --i-table demux.table.qza --o-visualization demux.table.qzv

qiime feature-table tabulate-seqs --i-data demux.rep-seqs.qza --o-visualization demux.rep-seqs.qzv

qiime metadata tabulate --m-input-file demux.stats.qza --o-visualization demux.stats.qzv

##Taxonomy classification/assignment:
##import UNITE database and convert fasta and taxonomy files to .qza files
##use UNITE version fungi-only, RefS

qiime tools import --type 'FeatureData[Sequence]' --input-path sh_qiime_release_25.07.2023/sh_refs_qiime_ver9_99_25.07.2023.fasta --output-path UNITE.reference.99.25072023.qza

qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path sh_qiime_release_25.07.2023/sh_taxonomy_qiime_ver9_99_25.07.2023.txt --output-path UNITE.taxonomy.99.25072023.qza

##Train naive bayes classifier

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads UNITE.reference.99.25072023.qza --i-reference-taxonomy UNITE.taxonomy.99.25072023.qza --o-classifier classifier.qza

#run sklearn naive bayes classifier

qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads demux.rep-seqs.qza --o-classification taxonomy-output.qza 

##To visualize:

qiime metadata tabulate --m-input-file taxonomy-output.qza --o-visualization taxonomy-output.qzv

qiime taxa collapse --i-table demux.table.qza --i-taxonomy taxonomy-output.qza --p-level 7 --o-collapsed-table table-collapsed.qza

##export feature table to biom file

qiime tools export --input-path table-collapsed.qza --output-path exported-feature-table

##convert biom file to tsv

cd exported-feature-table

biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
