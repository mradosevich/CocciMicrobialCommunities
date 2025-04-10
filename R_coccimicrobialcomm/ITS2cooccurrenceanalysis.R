############ Co-occurrence analysis ###############
library(tidyverse)
library(ecospat)

d <- read.csv("ITS2/ITS2.fulltaxa.csv")
d <- d %>% filter(rowSums(d[,8:dim(d)[2]])>0)
d <- d %>% filter(!stringr::str_detect(species,'unspecified'))
d <- d %>% filter(!stringr::str_detect(species,'unidentified'))
rownames(d) <- d[,7]
d <- d[,7:dim(d)[2]]
colnames <- colnames(d)

d.meta = read.csv("CPNM_MCA2.csv") %>% select(1:50)
rownames <- d.meta[,13]
Coccidioides_immitis <- as.data.frame(d.meta$positive) %>%
  rename('Coccidioides_immitis' = 'd.meta$positive')
Coccidioides_immitis <- as.data.frame(t(Coccidioides_immitis))
colnames(Coccidioides_immitis) <- rownames
Coccidioides_immitis$species <- c("Coccidioides_immitis")

d2 <- data.frame(species = d$species, ifelse(d[,-1] == 0, 0, 1))

d <- full_join(d, Coccidioides_immitis)
speciesData <- data.frame(species = d$species, ifelse(d[,-1] == 0, 0, 1))

# using ecospat

rownames(speciesData) <- speciesData[,1]
data <- speciesData[,-1]
data <- as.data.frame(t(data))
outpath <- getwd()
ecospat.Cscore(data = data, nperm = 1000, outpath = outpath, verbose = TRUE)
