### PCoA and PERMANOVA
# load libraries

library(vegan)
library(tidyverse)

# load data
d = read.csv("ITS2/ITS2.species.csv")
# remove unidentified and unspecified taxa
d = d[,-grep("Unknown", colnames(d))]
d = d[,-grep("unidentified", colnames(d))]
# remove rows containing no taxa
d = d[which(rowSums(d[,39:dim(d)[2]]) != 0),]

# for 16S, load file 16S.family.csv instead
#d = read.csv("16S/16S.family.csv")
#d = d[,-grep("Unknown", colnames(d))]
#d = d[,-grep("uncultured", colnames(d))]
#d = d[,colSums(d[,39:dim(d)[2]]) > 0]

# separate out meta data
d.meta = d[,1:38]
d.meta$treatment = as.character(d.meta$treatment)
d.meta$site = as.character(d.meta$site)
d.meta$exclosure = as.character(d.meta$exclosure)
d.meta$positive = as.character(d.meta$positive)
d.meta$lib <- substr(d.meta$Library_Name, 1, 7)

# make data matrix
d.matrix = as.matrix(d[,39:dim(d)[2]])

# squareroot and wisconsin 2x transform data
d.matrix.transformed = sqrt(d.matrix)
d.matrix.transformed = wisconsin(d.matrix.transformed)

# next, create dissimilarity matrix with transformed data (bray-curtis)
d.dist = vegdist(d.matrix.transformed, method="bray")

######## Permutational multivariate analysis of variance (PerMANOVA) #########
myPermaNova = adonis2(d.dist ~ site + treatment + exclosure + positive +
                        positive/treatment/exclosure/site, d.meta, parallel = 8, method = "bray",
                      permutations = 1000)

pvals = myPermaNova$`Pr(>F)`[1:7]
p.adjust(pvals, method = "fdr")

adj <- p.adjust(pvals, method = "fdr")
adj <- append(adj, c(NA, NA))
myPermaNova$`P value` <- adj
