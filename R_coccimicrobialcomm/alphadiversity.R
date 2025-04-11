### alpha diversity ###
# load libraries

library(phyloseq)
library(tidyverse)

# for ITS2 
#### load data and metadata #####
d = read.csv("ITS2/ITS2.fulltaxa.csv")
d <- d %>% filter(rowSums(d[,8:dim(d)[2]])>0)
d <- d %>% filter(!stringr::str_detect(species,'unspecified'))
d <- d %>% filter(!stringr::str_detect(species,'unidentified'))

d.meta = read.csv("samplemetadata.csv")
rownames(d.meta) = d.meta[,13]
d.meta$treatment=as.factor(d.meta$treatment)
d.meta = sample_data(d.meta)

# phyloseq 
OTU = d[,8:dim(d)[2]]
TAX = as(d[1:7], "matrix")

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

d.physeq = phyloseq(OTU, TAX, d.meta)
d.physeq

# for ITS2 or 16S 
###### estimate richness ######
est = estimate_richness(d.physeq, measures = c("Observed", "Shannon"))
est = est %>% mutate(
  evenness = Shannon/log(Observed)
)
est = est %>%  
  rename(richness = Observed)

######### Example of generating summary statistics for treatment variable #######
measures <- colnames(est)

est = est %>% mutate(barcode = rownames(est))

d.meta = as.data.frame(d.meta)
est <- full_join(est, d.meta, by = "barcode")

treat <- est %>% 
  group_by(treatment) %>%
  reframe(means = c(mean(richness), mean(Shannon), mean(evenness)),
          Q1 = c(quantile(richness, probs = .25), quantile(Shannon, probs = .25), quantile(evenness, probs = .25)),
          Q3 = c(quantile(richness, probs = .75), quantile(Shannon, probs = .75), quantile(evenness, probs = .75)))

treat$measure <- rep(measures, length(unique(treat$treatment)))
treat <- pivot_wider(treat, names_from = measure, values_from = c(means,Q1,Q3))

##### Example of running unpaired Wilcox test for treatment variable ##### 
burrow <- est %>% filter(treatment == 1)
surface <- est %>% filter(treatment == 0)

wilcox.test(burrow$richness, surface$richness, paired = F)
wilcox.test(burrow$evenness, surface$evenness, paired = F)
wilcox.test(burrow$Shannon, surface$Shannon, paired = F)

######### Matching paired samples for Wilcoxon sign-rank test #######
n <- est %>% group_by(cluster) %>% summarise(n = n())

est2 <- est %>% group_by(cluster) %>%
  mutate(
    avg_richness = mean(richness),
    avg_Shannon = mean(Shannon),
    avg_evenness = mean(evenness),
    tot_pos = sum(positive)) %>% 
  select(c("avg_richness", 
           "avg_Shannon",
           "avg_evenness",
           "cluster_pair",
           "cluster",
           "treatment",
           "tot_pos",
           "site",
           "exclosure")) %>%
  full_join(.,n, by = 'cluster')
est3 <- est2 %>% 
  distinct() %>%
  mutate(prop_pos = tot_pos/n)

######## Example of running wilcox test for paired samples ########
bur <- est3 %>% filter(treatment==1)
sur <- est3 %>% filter(treatment==0)
est4 <- full_join(bur, sur, by = "cluster_pair")

wilcox.test(est4$avg_richness.x, est4$avg_richness.y, paired=T)
wilcox.test(est4$avg_evenness.x, est4$avg_evenness.y, paired=T)
wilcox.test(est4$avg_Shannon.x, est4$avg_Shannon.y, paired=T)

########## code for importing and cleaning 16S data #########
# load data and metadata
d <- read.csv("16S/16S.family.csv")
rownames(d) <- d[,13]
d <- d %>% filter(rowSums(d[,39:dim(d)[2]])>0)
d <- d[,-grep("uncultured", colnames(d))]
d <- d[,-grep("Unknown", colnames(d))]
d <- d[,colSums(d[,39:dim(d)[2]]) > 0]

d.meta = d[,1:38]
rownames(d.meta) = d.meta[,13]
d.meta$treatment=as.factor(d.meta$treatment)
d.meta = sample_data(d.meta)

d <- d[,39:dim(d)[2]]
dt <- t(d)
tax <- rownames(dt)
tax <- as.data.frame(tax) %>%
  rename('Family' = 'tax')

# phyloseq analysis
OTU = dt
TAX = tax

setdiff(rownames(OTU),TAX$Family)

rownames(TAX) <- TAX$Family
TAX =  as.matrix(TAX)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

d.physeq = phyloseq(OTU, TAX, d.meta)
d.physeq

# use code above for ITS2 for remainder of analyses