# Load Packages -----------------------------------------------------------
## CRAN packages
install.packages("pacman")
pacman::p_load(
  dplyr, 
  seqinr,
  ape,
  phangorn,
  ggplot2,
  ggtree,
  stringr
  )

## Other packages
# devtools::install_github("lucymli/EpiGenR")
library(EpiGenR)

# date of last tip
latestdate <- "20230810"
## fasta file name
fn <- paste0("./longlabs_seq_AFGPAK_",latestdate,".fasta")

# Load data --------------------------------------------------------------
## sequence data with labels in "SEQUENCENAME_SAMPLETYPE_ADMIN0NAME_SAMPLEDECIMALDATE_CLUSTER_REGION" format
seqs <- read.FASTA(fn, type = "DNA")

# Check root-to-tip ------------------------------------------------------

## convert to phydat
seqs_phydat <- phyDat(seqs, type = "DNA")

## run a model test - v.slow
# mt <- modelTest(seqs_phydat)

dna_dist <- dist.ml(seqs_phydat, model="JC69")
NJ_tree <- NJ(dna_dist)

## Names and colours for regions used
getPalette <- c("#f4766c", "#ffd712", "#ef9b95", "#0a93d3", "#e50767", "#71e244", "#008121", "#ff80d9", "#858585", "#8999ff", "#0000b0",  "#a10052", "#922d92", "#c87dff", "#ff8a00", "#00da94")
names(getPalette) <- c("CENTRAL-CORRIDOR-AF", "CENTRE-AF", "CENTRE-PK", "EAST-PK", "ENDEMIC-ZONE", "GB", "KARACHI", "KP", "NORTH-AF", "NORTH-CORRIDOR-AF", "NORTH-CORRIDOR-PK", "SINDH", "SOUTH-CORRIDOR-AF", "SOUTH-CORRIDOR-PK", "SOUTH-PUNJAB", "WEST-AF")

tipdata <- do.call(rbind.data.frame, str_split(NJ_tree$tip.label, "_"))
colnames(tipdata) <- c("seqname", "sample_type", "adm0_name", "sample_date_num", "cluster", "region")
tipdata$node <- as.numeric(rownames(tipdata))


## Run IQtree from the path - use iqtree1
system( paste0( 'iqtree -nt AUTO -redo -m HKY -s ', fn ), intern=FALSE)

## Load in the generated tree
IQtr <- read.tree(paste0(fn, ".treefile"))

## Root to give a consistent molecular clock
dates <- sapply( strsplit( IQtr$tip.label, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )})
IQ_root <- rtt(IQtr, dates)


## Produce a data frame with time and divergence for each sequence
mc.df <- EpiGenR::root2tip.divergence(tr=IQ_root,tip.dates = dates)

## Linear regression
lm.res = lm(divergence~time,data=mc.df)

## Plot root-to-tip regression
ggplot(mc.df)+
  geom_point(aes(x=time, y=divergence))+
  geom_abline(intercept = lm.res$coefficients[1], slope = lm.res$coefficients[2], color = "red", size = 1.5)+
  theme_classic()+
  annotate("text", label = paste0("No. subst/site/year = ",round(lm.res$coefficients[2],3)), x= 2014.5, y = 0.16)+
  annotate("text", label = paste0("R2 = ",round(summary(lm.res)$r.squared,2)), x= 2014.5, y = 0.156)+
  labs(x = "Year", y = "Root-to-tip distance")

ggsave(paste0("./rtt",latestdate,".png"), width = 5, height = 5)

## Plot ML tree
tipdata <- do.call(rbind.data.frame, str_split(IQ_root$tip.label, "_"))
colnames(tipdata) <- c("seqname", "sample_type", "adm0_name", "sample_date_num", "cluster", "region")
tipdata$node <- as.numeric(rownames(tipdata))

ggtree(IQ_root)  %<+% tipdata +
  geom_tippoint(aes(fill = region), shape = 21, size =2, color = "black") +
  scale_fill_manual(values = getPalette, name = "Location")+
  guides(fill = guide_legend(override.aes = list(size = 5, shape = 22), ncol=2))+
  theme(legend.position = c(.25, .8),
        legend.box.background = element_rect(colour = "black"))

ggsave(paste0("./ML_root_location_",latestdate,".png"), width = 10, height = 7)
ggsave(paste0("./ML_root_location_",latestdate,".svg"), width = 10, height = 7)
