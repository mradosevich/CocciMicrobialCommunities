############ 16S Co-occurrence analysis ###############
library(tidyverse)
library(ecospat)

d <- read.csv("16S/16S.family.csv", sep = ',', check.names=FALSE)
d.meta = d %>% select(1:50)

d <- d[,52:dim(d)[2]]
d = d[,-grep("unspecified", colnames(d))]
rownames <- d.meta[,13]
d$samples <- rownames

Coccidioides_immitis <- as.data.frame(d.meta$positive) %>%
  rename('Coccidioides_immitis' = 'd.meta$positive')

Coccidioides_immitis$samples <- rownames


d <- full_join(d, Coccidioides_immitis, by = 'samples')
rownames(d) <- rownames
d <- d %>% select(!'samples')
data <- data.frame(ifelse(d == 0, 0, 1))

# using ecospat
# algorithm is fixed-equiprobable (sim2 equivalent), colSums are fixed
outpath <- getwd()
ecospat.Cscore(data = data, nperm = 1000, outpath = outpath, verbose = TRUE)
