# open the tsv file in R and make taxa tables (full taxonomy and at each taxonomic level)
# load libraries
library(tidyverse)
library(reshape2)

# import metadata
# read in csv
meta <- read.csv("samplemetadata.csv")
# import sampleid to match metadata to sequence filenames (for 16S use 16S_sampleid.csv)
sampleid <- read.csv("ITS2_sampleid.csv")
# note: sampleid file includes eight PCR control samples. NEGCTRL is water, POSCTRL is Neurospora lab stock (for ITS2) or Zymobiomics Microbial Community Standard (for 16S). 

meta <- full_join(meta, sampleid, by = join_by(barcode == sampleid))

# load feature-tables (for 16S use directory qiime2_data/16S)
# "\t" specifies tab-delimited, skip = 1 means skip first line of tsv file
tbl1 <- read.table("qiime2_data/ITS2/seqrun1/feature-table.tsv", sep = "\t", skip = 1,
                     comment.char = "")
tbl2 <- read.table("qiime2_data/ITS2/seqrun2/feature-table.tsv", sep = "\t", skip = 1,
                     comment.char = "")

# transpose matrix, adjust sample column title, convert characters to numeric, convert to dataframe
tbl1 <- t(tbl1)
colnames(tbl1) <- tbl1[1,]
colnames(tbl1)[1] <- "sample"
tbl1 <- tbl1[-1,]
tbl1 <- as.data.frame(tbl1)

tbl2 <- t(tbl2)
colnames(tbl2) <- tbl2[1,]
colnames(tbl2)[1] <- "sample"
tbl2 <- tbl2[-1,]
tbl2 <- as.data.frame(tbl2)

tbl <- full_join(tbl1, tbl2)

tbl[is.na(tbl)] <- 0

sample <- tbl %>% pull(sample)

tbl = sapply(tbl, as.numeric)
tbl <- as.data.frame(tbl)
tbl$sample <- sample

# melt otu data into long format
tbl.melt <- melt(tbl, id=c("sample"))
# remove taxa level prefix from otu names
tbl.melt$variable <- gsub("[a-z]__", "", tbl.melt$variable)
# make kingdom column
function.kingdom <- function(x) substr(x, 1, unlist(gregexpr(";", x))[1] -
                                        1)
tbl.melt$kingdom <- unlist(lapply(tbl.melt$variable, function.kingdom))
# make phylum column
function.phylum <- function(x) substr(x, unlist(gregexpr(";", x))[1] + 1,
                                     unlist(gregexpr(";", x))[2] - 1)
tbl.melt$phylum <- unlist(lapply(tbl.melt$variable, function.phylum))
# make class column
function.class <- function(x) substr(x, unlist(gregexpr(";", x))[2] + 1,
                                    unlist(gregexpr(";", x))[3] - 1)
tbl.melt$class <- unlist(lapply(tbl.melt$variable, function.class))

# make order column
function.order <- function(x) substr(x, unlist(gregexpr(";", x))[3] + 1,
                                    unlist(gregexpr(";", x))[4] - 1)
tbl.melt$order <- unlist(lapply(tbl.melt$variable, function.order))
# make family column
function.family <- function(x) substr(x, unlist(gregexpr(";", x))[4] + 1,
                                     unlist(gregexpr(";", x))[5] - 1)
tbl.melt$family <- unlist(lapply(tbl.melt$variable, function.family))
# make genus column
function.genus <- function(x) substr(x, unlist(gregexpr(";", x))[5] + 1,
                                    unlist(gregexpr(";", x))[6] - 1)
tbl.melt$genus <- unlist(lapply(tbl.melt$variable, function.genus))
# make species column
function.species<- function(x) substr(x, unlist(gregexpr(";", x))[6] + 1,
                                     nchar(x))
tbl.melt$species <- unlist(lapply(tbl.melt$variable, function.species))
# replace __ wih Unknown 
tbl.melt[tbl.melt == "__"] <- "Unknown"
# append "unspecified" with next highest identified taxa level
unspecified.phylum <- function (x) paste("unspecified_", substr(x, 1,
                                                               unlist(gregexpr(";", x))[1]-1), sep = "")
tbl.melt[tbl.melt$phylum=="unspecified",][,5:10] <-
  unlist(lapply(tbl.melt[tbl.melt$phylum=="unspecified",]$variable,
                unspecified.phylum))
unspecified.class <- function (x) paste("unspecified_", substr(x,
                                                              unlist(gregexpr(";", x))[1]+1, unlist(gregexpr(";", x))[2]-1), sep = "")
tbl.melt[tbl.melt$class=="unspecified",][,6:10] <-
  unlist(lapply(tbl.melt[tbl.melt$class=="unspecified",]$variable,
                unspecified.class))
unspecified.order <- function (x) paste("unspecified_", substr(x,
                                                              unlist(gregexpr(";", x))[2]+1, unlist(gregexpr(";", x))[3]-1), sep = "")
tbl.melt[tbl.melt$order=="unspecified",][,7:10] <-
  unlist(lapply(tbl.melt[tbl.melt$order=="unspecified",]$variable,
                unspecified.order))
unspecified.family <- function (x) paste("unspecified_", substr(x,
                                                               unlist(gregexpr(";", x))[3]+1, unlist(gregexpr(";", x))[4]-1), sep = "")
tbl.melt[tbl.melt$family=="unspecified",][,8:10] <-
  unlist(lapply(tbl.melt[tbl.melt$family=="unspecified",]$variable,
                unspecified.family))
unspecified.genus <- function (x) paste("unspecified_", substr(x,
                                                              unlist(gregexpr(";", x))[4]+1, unlist(gregexpr(";", x))[5]-1), sep = "")
tbl.melt[tbl.melt$genus=="unspecified",][,9:10] <-
  unlist(lapply(tbl.melt[tbl.melt$genus=="unspecified",]$variable,
                unspecified.genus))

unspecified.species <- function (x) paste("unspecified_", substr(x,
                                                                unlist(gregexpr(";", x))[5]+1, unlist(gregexpr(";", x))[6]-1), sep = "")
tbl.melt[tbl.melt$species=="unspecified",][,10] <-
  unlist(lapply(tbl.melt[tbl.melt$species=="unspecified",]$variable,
                unspecified.species))

# cast molten data at each taxonomic level
tbl.cast.phylum <- dcast(tbl.melt, sample ~ phylum, sum)
tbl.cast.class <- dcast(tbl.melt, sample ~ class, sum)
tbl.cast.order <- dcast(tbl.melt, sample ~ order, sum)
tbl.cast.family <- dcast(tbl.melt, sample ~ family, sum)
tbl.cast.genus <- dcast(tbl.melt, sample ~ genus, sum)
tbl.cast.species <- dcast(tbl.melt, sample ~ species, sum)

# cast full taxonomy table for reference
tbl.cast.full <- dcast(tbl.melt, phylum + class + order + family + genus + species ~ sample, sum)
# note different table structure for this file, where rows are taxa and columns are samples
tbl.full <- tbl.cast.full[,7:dim(tbl.cast.full)[2]]
tbl.full <- as.data.frame(t(tbl.full))
tbl.full$sample <- rownames(tbl.full)
tbl.full <- full_join(tbl.full, sampleid, by = join_by(sample == Library_Name))
tbl.full <- tbl.full %>% filter(sampleid )

# join taxonomic tables with metadata, filter out PCR controls in preparation for ecological analyses
tbl.combined.phylum <- meta %>% full_join(tbl.cast.phylum, by = c("Library_Name" = "sample")) %>% 
  filter(site != 'NA')
tbl.combined.class <- meta %>% full_join(tbl.cast.class, by = c("Library_Name" = "sample")) %>% 
  filter(site != 'NA')
tbl.combined.order <- meta %>% full_join(tbl.cast.order, by = c("Library_Name" = "sample")) %>% 
  filter(site != 'NA')
tbl.combined.family <- meta %>% full_join(tbl.cast.family, by = c("Library_Name" = "sample")) %>% 
  filter(site != 'NA')
tbl.combined.genus <- meta %>% full_join(tbl.cast.genus, by = c("Library_Name" = "sample")) %>% 
  filter(site != 'NA')
tbl.combined.species <- meta %>% full_join(tbl.cast.species, by = c("Library_Name" = "sample")) %>% 
  filter(site != 'NA')

# write to csv
write.csv(tbl.combined.phylum, "ITS2/ITS2.phylum.csv", row.names = FALSE)
write.csv(tbl.combined.class, "ITS2/ITS2.class.csv", row.names = FALSE)
write.csv(tbl.combined.order, "ITS2/ITS2.order.csv", row.names = FALSE)
write.csv(tbl.combined.family, "ITS2/ITS2.family.csv", row.names = FALSE)
write.csv(tbl.combined.genus, "ITS2/ITS2.genus.csv", row.names = FALSE)
write.csv(tbl.combined.species, "ITS2/ITS2.species.csv", row.names = FALSE)
write.csv(tbl.cast.full, "ITS2/ITS2.phylum.csv")

# make a full taxonomy table
tbl.full <- full_join(sampleid, tbl, by = c("Library_Name" = "sample"))

tbl.full = tbl.full[,-1]

barcode <- tbl.full %>% pull(sampleid)
tbl.full = sapply(tbl.full, as.numeric)
tbl.full = as.data.frame(tbl.full)
tbl.full$sampleid <- barcode

tbl.full = t(tbl.full)
colnames(tbl.full) = tbl.full[1,]
tbl.full = tbl.full[-1,]
tbl.full = as.data.frame(tbl.full)
taxanames <- rownames(tbl.full)

tbl.full = select(tbl.full, -c(NEGCTRL, POSCTRL))
tbl.full$taxa <- taxanames
tbl.full = tbl.full %>% relocate(taxa)

# remove taxa level prefix from otu names
tbl.full$taxa = gsub("[a-z]__", "", tbl.full$taxa)

# make kingdom column
function.kingdom = function(x) substr(x, 1, unlist(gregexpr(";", x))[1] -
                                        1)
tbl.full$kingdom = unlist(lapply(tbl.full$taxa, function.kingdom))
# make phylum column
function.phylum = function(x) substr(x, unlist(gregexpr(";", x))[1] + 1,
                                     unlist(gregexpr(";", x))[2] - 1)
tbl.full$phylum = unlist(lapply(tbl.full$taxa, function.phylum))
# make class column
function.class = function(x) substr(x, unlist(gregexpr(";", x))[2] + 1,
                                    unlist(gregexpr(";", x))[3] - 1)
tbl.full$class = unlist(lapply(tbl.full$taxa, function.class))

# make order column
function.order = function(x) substr(x, unlist(gregexpr(";", x))[3] + 1,
                                    unlist(gregexpr(";", x))[4] - 1)
tbl.full$order = unlist(lapply(tbl.full$taxa, function.order))
# make family column
function.family = function(x) substr(x, unlist(gregexpr(";", x))[4] + 1,
                                     unlist(gregexpr(";", x))[5] - 1)
tbl.full$family = unlist(lapply(tbl.full$taxa, function.family))
# make genus column
function.genus = function(x) substr(x, unlist(gregexpr(";", x))[5] + 1,
                                    unlist(gregexpr(";", x))[6] - 1)
tbl.full$genus = unlist(lapply(tbl.full$taxa, function.genus))
# make species column
function.species= function(x) substr(x, unlist(gregexpr(";", x))[6] + 1,
                                     nchar(x))
tbl.full$species = unlist(lapply(tbl.full$taxa, function.species))
# replace __ wih Unknown
tbl.full[tbl.full == "__"] = "unspecified"

#exported csv here
tbl.full = tbl.full %>% relocate(taxa, kingdom, phylum, class, order, family, genus, species) 

# append "unspecified" with next highest identified taxa level
unspecified.phylum = function (x) paste("unspecified_", substr(x, 1,
                                                               unlist(gregexpr(";", x))[1]-1), sep = "")
tbl.full[tbl.full$phylum=="unspecified",][,3:8] =
  unlist(lapply(tbl.full[tbl.full$phylum=="unspecified",]$taxa,
                unspecified.phylum))
unspecified.class = function (x) paste("unspecified_", substr(x,
                                                              unlist(gregexpr(";", x))[1]+1, unlist(gregexpr(";", x))[2]-1), sep = "")
tbl.full[tbl.full$class=="unspecified",][,4:8] =
  unlist(lapply(tbl.full[tbl.full$class=="unspecified",]$taxa,
                unspecified.class))
unspecified.order = function (x) paste("unspecified_", substr(x,
                                                              unlist(gregexpr(";", x))[2]+1, unlist(gregexpr(";", x))[3]-1), sep = "")
tbl.full[tbl.full$order=="unspecified",][,5:8] =
  unlist(lapply(tbl.full[tbl.full$order=="unspecified",]$taxa,
                unspecified.order))
unspecified.family = function (x) paste("unspecified_", substr(x,
                                                               unlist(gregexpr(";", x))[3]+1, unlist(gregexpr(";", x))[4]-1), sep = "")
tbl.full[tbl.full$family=="unspecified",][,6:8] =
  unlist(lapply(tbl.full[tbl.full$family=="unspecified",]$taxa,
                unspecified.family))
unspecified.genus = function (x) paste("unspecified_", substr(x,
                                                              unlist(gregexpr(";", x))[4]+1, unlist(gregexpr(";", x))[5]-1), sep = "")
tbl.full[tbl.full$genus=="unspecified",][,7:8] =
  unlist(lapply(tbl.full[tbl.full$genus=="unspecified",]$taxa,
                unspecified.genus))

unspecified.species = function (x) paste("unspecified_", substr(x,
                                                                unlist(gregexpr(";", x))[5]+1, unlist(gregexpr(";", x))[6]-1), sep = "")
tbl.full[tbl.full$species=="unspecified",][,8] =
  unlist(lapply(tbl.full[tbl.full$species=="unspecified",]$taxa,
                unspecified.species))
tbl.full<- tbl.full %>% select(-taxa)

write.csv(tbl.full, "ITS2/ITS2.fulltaxa.csv", row.names = FALSE)
