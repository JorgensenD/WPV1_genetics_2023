# Local transmission lineages



#install_github("YuLab-SMU/ggtree")
pacman::p_load(ggtree,
               treeio,
               tidyr,
               parallel,
               doParallel,
               foreach,
               lubridate,
               stringr,
               dplyr,
               ggplot2,
               egg
)

## functions
source("./github_code/3_Local_Transmsision_Lineages/LTL_functions.R")

# Tree input --------------------------------------------------------------

## Import tree in two formats
trees <- read.beast("./all_20230810/ds_trees.trees")
tres  <- read.nexus("./all_20230810/ds_trees.trees")

cons_tree <- read.beast("./all_20230810/MJ_MCC_CA_DS.trees")

mrsd=2023.605
mrsd=as.Date(date_decimal(mrsd))

# Sometimes get a 50/50 split in the BEAST calculation of a node location - example of a manual fix for this
# tree@data$state[tree@data$state=="SOUTH-CORRIDOR-AF+SOUTH-CORRIDOR-PK"] <- "SOUTH-CORRIDOR-PK"
n_cores <- detectCores()
cluster <- makeCluster(n_cores-3)
registerDoParallel(cluster)

return_ds <- foreach(g = 1:length(trees), .packages = c("ggplot2","gginnards","plyr","dplyr",
                                           "lubridate","ape","grid","gridExtra",
                                           "phylobase","stringr","treeio",
                                           "tidytree","egg","igraph", "phytools", "ggtree")) %dopar% {
  source("./github_code/3_Local_Transmsision_Lineages/LTL_functions.R")
  tree <- trees[[g]]
  tre <- tres[[g]]

# node heights  (from tip back in beast tree) 
heights <- distinct(data.frame(height = max(c(nodeHeights(tree@phylo)))-c(nodeHeights(tree@phylo)), node = as.character(c(tree@phylo$edge))))
lengths <- data.frame(node = as.character(tree@phylo$edge[,2]), length = tree@phylo$edge.length)

tree@data <- left_join(tree@data, heights, by = "node") %>%
              left_join(lengths, by = "node")

## extract the ancestral node location and number for each node and add this to the dataframe 
tree@data$ancestor <- mapply(ancestor, tree@data$node)

ancestor_loc <- tree@data[,c("node","state")]
names(ancestor_loc) <- c("ancestor","ancestor_loc")
tree@data <- as_tibble(join(tree@data,ancestor_loc, by="ancestor"))

#add branch length to node age to get the branching point
tree@data$branching <- as.numeric(tree@data$height) + as.numeric(tree@data$length)#/2
tree@data$branching <- date_decimal(decimal_date(mrsd)-tree@data$branching)



# Plotting tree -----------------------------------------------------------

getPalette <- c("#f4766c", "#ffd712",
                "#ef9b95", "#0a93d3",
                "#e50767", "#71e244",
                "#008121", "#ff80d9",
                "#858585", "#8999ff",
                "#0000b0", "#a10052",
                "#922d92", "#c87dff",
                "#ff8a00", "#00da94")
names(getPalette) <- c("CENTRAL-CORRIDOR-AF", "CENTRE-AFG",
                       "CENTRE-PAK", "EAST-PAK",
                       "ENDEMIC-ZONE", "GB",
                       "KARACHI", "KP",
                       "NORTH-AFG", "NORTH-CORRIDOR-AF",
                       "NORTH-CORRIDOR-PK", "SINDH",
                       "SOUTH-CORRIDOR-AF", "SOUTH-CORRIDOR-PK",
                       "SOUTH-PUNJAB", "WEST-AFG")
newpal <- getPalette # save a copy
  
## need to plot tree in full first as the subtree functions just crop the original plot
tr <- ggtree(tree, mrsd=mrsd, aes(colour=state)) +
  scale_y_reverse()+
  scale_color_manual(values = getPalette, drop=FALSE, na.value="white")+
  geom_tippoint()+
  theme_tree2() +
  theme(legend.position = "none")+
  geom_rootpoint()

## same tree but less complex plot
tr2 <- ggtree(tree, mrsd=mrsd, aes(colour=state)) +
  scale_y_reverse()+
  geom_tippoint()+
  theme_tree2() +
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_rootpoint()

# Local Transmisison Lineages ---------------------------------------------
## extract subtree at any node where the ancestor location and current location do not match 
## collapse any downstream subtrees and label with their location 

## Find nodes where ancestor doesn't match the current location
treelvls <- tree@data[order(tree@data$branching),]
treelvls <- unique(treelvls$state)
subtree <- tree@data[tree@data$state!=tree@data$ancestor_loc & !is.na(tree@data$state)& !is.na(tree@data$ancestor),]
subtree$state <- factor(subtree$state, levels=treelvls)
subtree <- subtree[order(subtree$state, subtree$branching),] # Sort this by location
subtree$state <- as.character(subtree$state)

# tree split and plot -----------------------------------------------------

# prepare to run
tstcl <- ggplot_build(tr) ## need to build to get the coordinates of each node
tipsh <- c(21,22)
names(tipsh) <- c("AFP","ES")
treeroot <- rootnode(tre)

subtreelist <- list()
treedata <- list()

## apply subtree custom fuction
## Cannot deal with plotting the root with this function so seperately done after
treedata <- lapply(1:nrow(subtree),subtrefunc_noplot, subtree=subtree, tree=tree)

downstreamtree <- data.frame(treedata)
downstreamtree <- data.table::transpose(downstreamtree)

## convert the colours in the list into locations
downstreamtree[,1] <- unlist(lapply(downstreamtree[,1], locatefunc)) 
downstreamtree[,2] <- unlist(lapply(downstreamtree[,2], locatefunc))

#number of importations per area
downstreamtree$V1 <- factor(downstreamtree$V1, levels=sort(unique(downstreamtree$V2)))
downstreamtree$V3 <- as.Date(downstreamtree$V3)
downstreamtree$V4 <- as.Date(downstreamtree$V4)
downstreamtree$diff <- as.numeric(downstreamtree$V4-downstreamtree$V3)
downstreamtree$V5 <- as.numeric(downstreamtree$V5)
downstreamtree$V6 <- as.numeric(downstreamtree$V6)
# Plot split trees --------------------------------------------------------
downstreamtree$ID <- 2:(nrow(downstreamtree)+1)
downstreamtree[downstreamtree$ID==nrow(downstreamtree)+1,]$ID <- 1
downstreamtree <- downstreamtree[order(downstreamtree$ID),]
rownames(downstreamtree) <- NULL
return(downstreamtree)
print(g)
}

stopCluster(cl = cluster)

ds_copy_safe <- return_ds
rolling_type_clusters <- function(date_end, downstreamtree, roll_period){
  order_regi <- levels(as.factor(downstreamtree$V2))
  day_include <- seq(date_end-roll_period, date_end, by = "days")
  suppressMessages({
    matchdata <-  downstreamtree  %>%
      drop_na() %>%
      rowwise() %>%
      mutate(match = ifelse(any(between(day_include, V3, V4)), 1, 0)) %>%
      filter(match == 1)
    
    matchdata <- matchdata %>%
      group_by(V2, Category) %>%
      tally() %>% 
      mutate(freq = prop.table(n)) %>%
      right_join(expand.grid(V2 = order_regi), by = "V2") %>%
      mutate(date = date_end)
  })
  
  return(matchdata[order(match(names(matchdata), order_regi))])
}

## Classify all of these clusters and then count them in each time period


for(g in 1:length(trees)){
return_ds[[g]]$onwards <- do.call(rbind.data.frame, lapply(return_ds[[g]]$V8, function(x) str_detect(x, "collapse")))
return_ds[[g]]$onwards <- return_ds[[g]]$onwards[,1]


return_ds[[g]] <- return_ds[[g]] %>%
  mutate(Category = case_when(diff<=183 & V5<5 & onwards == F ~ "dead end",
                              diff>=365 & V5>10 ~ "persistent",
                              onwards == T ~ "export",
                              T ~ "other")) 
}


 
cluster <- makeCluster(n_cores-3)
registerDoParallel(cluster)

return_cluster_types <- foreach(g = 1:length(trees), .packages = c("tidyr", "dplyr")) %dopar% {
alldays <- seq(as.Date("2012-01-01"), as.Date("2023-08-09"), by = 15)

cluster_types <- lapply(alldays,rolling_type_clusters, downstreamtree=return_ds[[g]], roll_period = 15)
return(do.call(rbind, cluster_types))

}
stopCluster(cl = cluster)

#save(return_ds, return_cluster_types, file = "sample300_downstreamtree.RData")


all_clustertypes <- do.call(rbind, return_cluster_types)
xmax = decimal_date(mrsd)

all_clustertypes_summary <- all_clustertypes %>%
  group_by(V2, Category, date) %>%
  dplyr::summarise(median = median(n, na.rm=T), lower = quantile(n, 0.025, na.rm = T), upper = quantile(n, 0.975, na.rm = T))

all_clustertypes_summary$Category <- relevel(as.factor(all_clustertypes_summary$Category), 'other')  
all_clustertypes_summary[all_clustertypes_summary$V2 == "ENDEMIC-ZONE",]$V2 <- "CENTRAL-CORRIDOR-PK"

all_clustertypes_summary[which(all_clustertypes_summary$date>as.Date(date_decimal(xmax-0.5)) & all_clustertypes_summary$Category=='dead end'),]$Category <- 'other'

all_clustertypes_summary %<>%
  mutate(region = recode(V2, "CENTRAL-CORRIDOR-AF"="CENTRAL CORRIDOR AF", "CENTRE-AFG"="CENTRAL AF", "CENTRE-PAK"="CENTRAL PK", "EAST-PAK"="EAST PK",
                         "CENTRAL-CORRIDOR-PK"="CENTRAL CORRIDOR PK", "GB"="GB", "KARACHI"="KARACHI", "KP"="KP",
                         "NORTH-AFG"="NORTH AF", "NORTH-CORRIDOR-AF"="NORTH CORRIDOR AF", "NORTH-CORRIDOR-PK"="NORTH CORRIDOR PK", "SINDH"="SINDH",
                         "SOUTH-CORRIDOR-AF"="SOUTH CORRIDOR AF", "SOUTH-CORRIDOR-PK"="SOUTH CORRIDOR PK", "SOUTH-PUNJAB"="SOUTH PUNJAB", "WEST-AFG"="WEST AF"))

all_clustertypes_summary %<>%
  ungroup()%>%
  select(-V2) %>%
  select(region, everything())

save(all_clustertypes_summary, file = "./Plot_data/Figure4.RData")

bounds <- ggplot() +
  geom_line(data = all_clustertypes_summary, aes(x = date, y = median, color = Category, group = Category))+
  geom_ribbon(data = all_clustertypes_summary, aes(x = date, ymin = lower, ymax = upper, fill = Category, group = Category), alpha = 0.4)+
  geom_vline(xintercept = as.Date(date_decimal(xmax-1)))+
  #scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values = c("gray80", "orange", "deeppink", "purple4"), na.translate = F)+
  scale_fill_manual(values = c("gray80", "orange", "deeppink", "purple4"), na.translate = F)+
  theme_bw()+
  scale_x_date(limits = c(as.Date("2012-01-01"), ymd(mrsd)), expand = c(0,0))+
  labs(x = "Date",
       y = "Number of lineages") +
  facet_wrap(~region, ncol = 4) +
  theme(legend.direction = "horizontal", legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5), color = guide_legend(title.position = "top", title.hjust = 0.5))

# plot the barplot

barchart <- ggplot() +
  geom_bar(data = all_clustertypes_summary, aes(x = date, y = median, color = Category, fill = Category, group = Category), stat = "identity")+
  geom_vline(xintercept = as.Date(date_decimal(xmax-1)))+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values = c("gray80", "orange", "deeppink", "purple4"), na.translate = F)+
  scale_fill_manual(values = c("gray80", "orange", "deeppink", "purple4"), na.translate = F)+
  theme_bw()+
  scale_x_date(limits = c(as.Date("2012-01-01"), ymd(mrsd)), expand = c(0,0))+
  labs(x = "Date",
       y = "Number of lineages") +
  facet_wrap(~region, ncol = 4) +
  theme(legend.direction = "horizontal", legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5), color = guide_legend(title.position = "top", title.hjust = 0.5))

ggsave("./classification_new.svg",barchart, height = 6, width = 8, units = "in", dpi = 700)
ggsave("./classification_new.png",barchart, height = 6, width = 8, units = "in", dpi = 700)

ggsave("./classification_err_new.svg",bounds, height = 6, width = 8, units = "in", dpi = 700)
ggsave("./classification_err_new.png",bounds, height = 6, width = 8, units = "in", dpi = 700)

# number of clusters and number of dead ends for each tree
n_clusters <- vector()
n_deadend <- vector()
for(g in 1:length(trees)){
  n_clusters[g] <- nrow(return_ds[[g]])
  n_deadend[g] <- length(which(return_ds[[g]]$Category == "dead end"))
}
median(n_deadend)
quantile(n_deadend, 0.025, na.rm = T)
quantile(n_deadend, 0.975, na.rm = T)

# number persistent by region
return_ds[[1]]$V2 <- as.factor(return_ds[[1]]$V2)
return_ds[[1]]$Category <- as.factor(return_ds[[1]]$Category)
levels(return_ds[[1]]$Category) <- c("dead end","export","other","persistent")
levels(return_ds[[1]]$V2) <-  levels(return_ds[[1]]$V1)
n_lineages <- return_ds[[1]] %>%
  group_by(V2, Category, .drop = F) %>%
  summarise(count = n())

for(g in 2:length(trees)){
# dataframe not a vector
return_ds[[g]]$V2 <- as.factor(return_ds[[g]]$V2)
return_ds[[g]]$Category <- as.factor(return_ds[[g]]$Category)
levels(return_ds[[g]]$Category) <- c("dead end","export","other","persistent")
levels(return_ds[[g]]$V2) <-  levels(return_ds[[g]]$V1)
temp_df <- return_ds[[g]] %>%
    group_by(V2, Category, .drop = F) %>%
    summarise(count = n())
n_lineages <- cbind(n_lineages, temp_df$count)
}

n_lineages$median <- apply(n_lineages[,3:302], 1, median)
n_lineages$low <- apply(n_lineages[,3:302], 1, quantile, prob = 0.025, na.rm = T )
n_lineages$high <- apply(n_lineages[,3:302], 1, quantile, prob = 0.975, na.rm = T )
n_lineages[n_lineages$Category == "persistent",]$median

# time to detection
for(g in 1:length(trees)){
return_ds[[g]] <- return_ds[[g]][return_ds[[g]]$V5>0,]
detection <- do.call(rbind.data.frame, lapply(return_ds[[g]]$V8, function(x) min(as.numeric(gsub("_", "",str_extract_all(x, "\\_\\d+\\.?\\d*", simplify = T)))))) ## more complex as some dates don't have decimals
names(detection) <- "detection"
return_ds[[g]]$detection <- as.Date(date_decimal(detection$detection))
return_ds[[g]]$ttd <- return_ds[[g]]$detection - return_ds[[g]]$V3
}

do.call(rbind, return_ds) %>%
  summarise(mean = mean(ttd, rm.na=T),
            n = nrow(.),
            s = sd(ttd),
            error = 1.96*(s/sqrt(n)),
            lower = ifelse(mean-error>=0,mean-error,0),
            upper = mean+error)



save(all_clustertypes, file = "all_clustertypes.RData")


