# Plot orphan lineages

# relies on the downstreamtree file and tree file from the LTL analysis
load("downstreamtree_all.rdata")
# also need to load the sequence data
seqdata <- read.FASTA("longlabs_seq_20230810.fasta")


# Functions ---------------------------------------------------------------
# give mat_meta an ancestral region column
region_min_anc <- function(name_string, dmat){
  ancestral_dates <- mat_meta[which(mat_meta$date < mat_meta[mat_meta$fullname == name_string,]$date),]$fullname
  locn <- dmat[name_string,ancestral_dates]
  mat_meta[mat_meta$fullname ==  names(which.min(locn[locn>0])),]$region
}

## distance matrix with names
min_anc <- function(name_string, dmat){
  ancestral_dates <- mat_meta[which(mat_meta$date < mat_meta[mat_meta$fullname == name_string,]$date),]$fullname
  # name these with the tipnames
  locn <- dmat[name_string,ancestral_dates]
  min(locn[locn>0])
}

# mean distance to other sequences
mean_min_dist <- function(z, dmat){
  # list of tips in subtree
  tips <-  str_subset(unlist(strsplit(z, ", ")), pattern = "collapse_*", negate = T)
  tips <- str_subset(tips, pattern = "^NA", negate = T)
  ## date of every tip in the tree
  mat_meta <- colsplit(rownames(dmat), "\\_", names = c("ID", "source", "country", "date", "clade", "region"))
  mat_meta$fullname <- rownames(dmat)
  # apply over the subtree tips
  
  dists <- sapply(tips, min_anc, dmat = dmat)
  mean(dists[is.finite(dists)])
}

# median distance to other sequences
median_min_dist <- function(z, dmat){
  # list of tips in subtree
  tips <-  str_subset(unlist(strsplit(z, ", ")), pattern = "collapse_*", negate = T)
  tips <- str_subset(tips, pattern = "^NA", negate = T)
  ## date of every tip in the tree
  mat_meta <- colsplit(rownames(dmat), "\\_", names = c("ID", "source", "country", "date", "clade", "region"))
  mat_meta$fullname <- rownames(dmat)
  # apply over the subtree tips
  
  dists <- sapply(tips, min_anc, dmat = dmat)
  median(dists[is.finite(dists)])
}


# count orphan detections
count_orphans <- function(z, dmat){
  tips <-  str_subset(unlist(strsplit(z, ", ")), pattern = "collapse_*", negate = T)
  tips <- str_subset(tips, pattern = "^NA", negate = T)
  ## date of every tip in the tree
  mat_meta <- colsplit(rownames(dmat), "\\_", names = c("ID", "source", "country", "date", "clade", "region"))
  mat_meta$fullname <- rownames(dmat)
  # apply over the subtree tips
  dists <- sapply(tips, min_anc, dmat = dmat)
  sum(dists[is.finite(dists)]>0.015)
}

# SD of minimum distances
sd_min_dist <- function(z, dmat){
  # list of tips in subtree
  tips <-  str_subset(unlist(strsplit(z, ", ")), pattern = "collapse_*", negate = T)
  tips <- str_subset(tips, pattern = "^NA", negate = T)
  ## date of every tip in the tree
  mat_meta <- colsplit(rownames(dmat), "\\_", names = c("ID", "source", "country", "date", "clade", "region"))
  mat_meta$fullname <- rownames(dmat)
  # apply over the subtree tips
  
  dists <- sapply(tips, min_anc, dmat = dmat)
  sd(dists[is.finite(dists)])
}

# number of tips
n_min_dist <- function(z){
  # list of tips in subtree
  tips <-  str_subset(unlist(strsplit(z, ", ")), pattern = "collapse_*", negate = T)
  tips <- str_subset(tips, pattern = "^NA", negate = T)
  length(tips)
}

# Run and plot ------------------------------------------------------------
dna_dist_polio <- dist.dna(seqdata, as.matrix = T)


mat_meta <- colsplit(rownames(dna_dist_polio), "\\_", names = c("ID", "source", "country", "date", "clade", "region"))
mat_meta$fullname <- rownames(dna_dist_polio)

mat_meta$anc_region <- sapply(mat_meta$fullname, region_min_anc, dmat = dna_dist_polio)
mat_meta$internal <- mat_meta$region == mat_meta$anc_region


# ancestor in same or different region
mean_detect_all_r <- mat_meta %>%
  mutate(date_group =  lubridate::floor_date(date_decimal(date), "6 months")) %>%
  group_by(date_group, internal) %>%
  summarise(mean = mean_min_dist(fullname, dmat = dna_dist_polio),
            median = median_min_dist(fullname, dmat = dna_dist_polio),
            s = sd_min_dist(fullname, dmat = dna_dist_polio),
            n = n_min_dist(fullname),
            error = 1.96*(s/sqrt(n)),
            lower = ifelse(mean-error>=0,mean-error,0),
            upper = mean+error,
            max_div = max(sapply(fullname, min_anc, dmat = dna_dist_polio)),
            min_div = min(sapply(fullname, min_anc, dmat = dna_dist_polio)),
            count = count_orphans(fullname, dmat = dna_dist_polio))

## current names of dists are the tip names not the ancestor names
#
## relabel internal 
mean_detect_all_r$internal <- as.factor(mean_detect_all_r$internal)

levels(mean_detect_all_r$internal) <- list(Different = FALSE, Matching = TRUE)

dist_mean_all_r <- ggplot(mean_detect_all_r) +
  geom_col(aes(x = date_group, y = count, group = internal, fill = internal), show.legend = T, just = 0)+
  theme_bw()+
  labs(x = "Date",
       y = "Number of orphan detections (>1.5% divergence)",
       fill = "Region of ancestor and \norphan detection")+
  coord_cartesian(expand = F)+
  theme(legend.position = c(0.85,0.8),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



dist_mean_all_r_range <- ggplot(mean_detect_all_r) +
  geom_rect(aes(ymax = max_div*100, ymin = min_div*100, xmin  = date_group, xmax = date_group %m+% months(6), group = internal, , fill = internal), alpha = 0.5, show.legend = F)+
  geom_segment(aes(y = median*100, yend = median*100, x = date_group, xend = date_group %m+% months(6), color = internal, group = internal), show.legend = F)+
  coord_cartesian(ylim = c(0,6), xlim = c(as.POSIXct("2012-01-01"), as.POSIXct("2024-01-01")), expand = F)+
  geom_hline(yintercept = 1.5, color = "black") +
  geom_hline(yintercept = 3, color = "black", linetype = 2)+
  labs(x = "Date",
       y = "Difference to previous sequences (% divergence)")+
  theme_bw()+
  facet_wrap(~internal, nrow = 2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# region plots internal and external - just subset and plot separately
require(cowplot)

orphan_detections <- plot_grid(dist_mean_all_r, dist_mean_all_r_range, ncol = 1, rel_heights = c(0.33, 0.66))


ggsave(orphan_detections, file = "n_orphans_median.svg", width = 6, height = 8,  dpi = 700)
