### alpha diversity ###
# load libraries

library(phyloseq)
library(ggplot2)
library(tidyverse)
library(vegan)

# FOR ITS2
# load data and metadata
d = read.csv("ITS2/ITS2.fulltaxa.csv")
d <- d %>% filter(rowSums(d[,8:dim(d)[2]])>0)
d <- d %>% filter(!stringr::str_detect(species,'unspecified'))
d <- d %>% filter(!stringr::str_detect(species,'unidentified'))

d.meta = read.csv("CPNM_MCA2.csv") %>% select(1:50)
rownames(d.meta) = d.meta[,13]

d.meta$treatment=as.factor(d.meta$treatment)
d.meta = sample_data(d.meta)

# phyloseq analysis
OTU = d[,8:dim(d)[2]]
TAX = as(d[1:7], "matrix")

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

d.physeq = phyloseq(OTU, TAX, d.meta)
d.physeq

########## FOR 16S #########
# load data and metadata
#d <- read.csv("16S/16S.family.csv")
#rownames(d) <- d[,13]
#d <- d %>% filter(rowSums(d[,52:dim(d)[2]])>0)
#d <- d[,-grep("unspecified", colnames(d))]
#d <- d[,-grep("uncultured", colnames(d))]
#d <- d[,-grep("Unknown", colnames(d))]
#d <- d[,colSums(d[,52:dim(d)[2]]) > 0]

#d.meta = d[,1:50]
#rownames(d.meta) = d.meta[,13]
#d.meta$treatment=as.factor(d.meta$treatment)
#d.meta = sample_data(d.meta)

#d <- d[,52:dim(d)[2]]
#dt <- t(d)
#tax <- rownames(dt)
#tax <- as.data.frame(tax) %>%
#  rename('Family' = 'tax')

# phyloseq analysis
#OTU = dt
#TAX = tax

#setdiff(rownames(OTU),TAX$Family)

#rownames(TAX) <- TAX$Family
#TAX =  as.matrix(TAX)

#OTU = otu_table(OTU, taxa_are_rows = TRUE)
#TAX = tax_table(TAX)

#d.physeq = phyloseq(OTU, TAX, d.meta)
#d.physeq

##############################################
###### FOR ITS2 or 16S ########
#estimate richness
est = estimate_richness(d.physeq, measures = c("Observed", "Shannon"))
est = est %>% mutate(
  evenness = Shannon/log(Observed)
)
est = est %>%  
  rename(richness = Observed)

# estimator
measures <- colnames(est)
meanIQR <- function(est) {
  est %>%
   summarise(
    mean = sapply(est, FUN = mean),
    Q1 = lapply(est, quantile, prob=c(.25)),
    Q3 = lapply(est, quantile, prob=c(.75))
   )
}

#result <- meanIQR(est[1:3])
#rownames(result) <- measures
############# Wilcox Tests #######

est = est %>% mutate(barcode = rownames(est))

d.meta = as.data.frame(d.meta)
est <- full_join(est, d.meta, by = "barcode")

burrow <- est %>% filter(treatment == 1)
surface <- est %>% filter(treatment == 0)
positive <- est %>% filter(positive == 1)
negative <- est %>% filter(positive == 0)
exclosure <- est %>% filter(exclosure == 0)
nonexclosure <- est %>% filter(exclosure ==1)
north <- est %>% filter(site == 0)
south <- est %>% filter(site == 1)

wilcox.test(burrow$richness, surface$richness, paired = F)
wilcox.test(burrow$evenness, surface$evenness, paired = F)
wilcox.test(burrow$Shannon, surface$Shannon, paired = F)

wilcox.test(positive$richness, negative$richness, paired = F)
wilcox.test(positive$evenness, negative$evenness, paired = F)
wilcox.test(positive$Shannon, negative$Shannon, paired = F)

wilcox.test(exclosure$richness, nonexclosure$richness, paired = F)
wilcox.test(exclosure$evenness, nonexclosure$evenness, paired = F)
wilcox.test(exclosure$Shannon, nonexclosure$Shannon, paired = F)

wilcox.test(north$richness, south$richness, paired = F)
wilcox.test(north$evenness, south$evenness, paired = F)
wilcox.test(north$Shannon, south$Shannon, paired = F)

burrow$exclosure=as.factor(burrow$exclosure)
surface$exclosure=as.factor(surface$exclosure)

posbur <- positive %>% filter(treatment == 1)
negbur <- negative %>% filter(treatment == 1)

wilcox.test(posbur$richness, negbur$richness, paired = F)
wilcox.test(posbur$evenness, negbur$evenness, paired = F)
wilcox.test(posbur$Shannon, negbur$Shannon, paired = F)

posburex <- posbur %>% filter(exclosure == 0)
posburnonex <- posbur %>% filter(exclosure == 1)

wilcox.test(posburex$richness, posburnonex$richness, paired = F)
wilcox.test(posburex$evenness, posburnonex$evenness, paired = F)
wilcox.test(posburex$Shannon,  posburnonex$Shannon, paired = F)

possur <- positive %>% filter(treatment == 0)
negsur <- negative %>% filter(treatment == 0)

wilcox.test(possur$richness, negsur$richness, paired = F)
wilcox.test(possur$evenness, negsur$evenness, paired = F)
wilcox.test(possur$Shannon, negsur$Shannon, paired = F)

possurex <- possur %>% filter(exclosure == 0)
possurnonex <- possur %>% filter(exclosure == 1)

wilcox.test(possurex$richness, possurnonex$richness, paired = F)
wilcox.test(possurex$evenness, possurnonex$evenness, paired = F)
wilcox.test(possurex$Shannon, possurnonex$Shannon, paired = F)

negburex <- negbur %>% filter(exclosure == 0)
negburnonex <- negbur %>% filter(exclosure == 1)

wilcox.test(negburex$richness, negburnonex$richness, paired = F)
wilcox.test(negburex$evenness, negburnonex$evenness, paired = F)
wilcox.test(negburex$Shannon, negburnonex$Shannon, paired = F)

negsurex <- negsur %>% filter(exclosure == 0)
negsurnonex <- negsur %>% filter(exclosure == 1)

wilcox.test(negsurex$richness, negsurnonex$richness, paired = F)
wilcox.test(negsurex$evenness, negsurnonex$evenness, paired = F)
wilcox.test(negsurex$Shannon, negsurnonex$Shannon, paired = F)

wilcox.test(posburex$richness, negburex$richness, paired = F)
wilcox.test(posburex$evenness, negburex$evenness, paired = F)
wilcox.test(posburex$Shannon, negburex$Shannon, paired = F)

wilcox.test(posburnonex$richness, negburnonex$richness, paired = F)
wilcox.test(posburnonex$evenness, negburnonex$evenness, paired = F)
wilcox.test(posburnonex$Shannon, negburnonex$Shannon, paired = F)

surexcl <- surface %>% filter(exclosure == 0)
surnonexcl <- surface %>% filter(exclosure == 1)
burexcl <- burrow %>% filter(exclosure == 0)
burnonexcl <- burrow %>% filter(exclosure == 1)

wilcox.test(surexcl$Shannon, surnonexcl$Shannon, paired = F)
wilcox.test(burexcl$Shannon, burnonexcl$Shannon, paired = F)
wilcox.test(surexcl$Shannon, burexcl$Shannon, paired = F)
wilcox.test(surnonexcl$evenness, burnonexcl$evenness, paired = F)

mean(surexcl$richness)
mean(surnonexcl$richness)
mean(burexcl$richness)
mean(burnonexcl$richness)

Npasture <- est %>% filter(site == 0)
result16 <- meanIQR(Npasture[1:3])
rownames(result16) <- measures

Spasture <- est %>% filter(site == 1)
result17 <- meanIQR(Spasture[1:3])
rownames(result17) <- measures

negburrow <- burrow %>% filter(positive == 0)
posburrow <- burrow %>% filter(positive == 1)

wilcox.test(negburrow$Observed, posburrow$Observed, paired = F)

exclosure <- burrow %>% filter(exclosure == 0)
nonexclosure <- burrow %>% filter(exclosure == 1)
exclneg <- exclosure %>% filter(positive==0)
exclpos <- exclosure %>% filter(positive==1)
wilcox.test(exclneg$Observed, exclpos$Observed, paired = F)

nonexclneg <- nonexclosure %>% filter(positive==0)
nonexclpos <- nonexclosure %>% filter(positive==1)
wilcox.test(nonexclneg$Observed, nonexclpos$Observed, paired = F)

######### matching paired samples for wilcoxon sign-rank test #######
clusters <- est %>% group_by(cluster) %>% 
  summarise(
    mean(richness),
    mean(Shannon),
    sd(richness),
    sd(Shannon))

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
write.csv(est3, file = "ITS2_est_alpha.csv", row.names = FALSE)

# wilcox test for paired samples 
bur <- est3 %>% filter(treatment==1)
sur <- est3 %>% filter(treatment==0)
est4 <- full_join(bur, sur, by = "cluster_pair")

est5 <- est4 %>% mutate(richness_diff = abs(avg_richness.x - avg_richness.y),
                        shannon_diff = abs(avg_Shannon.x - avg_Shannon.y))

mean(est5$richness_diff)
mean(est5$shannon_diff)

sd(est5$richness_diff)
sd(est5$shannon_diff)

wilcox.test(est4$avg_richness.x, est4$avg_richness.y, paired=T)
wilcox.test(est4$avg_evenness.x, est4$avg_evenness.y, paired=T)
wilcox.test(est4$avg_Shannon.x, est4$avg_Shannon.y, paired=T)
