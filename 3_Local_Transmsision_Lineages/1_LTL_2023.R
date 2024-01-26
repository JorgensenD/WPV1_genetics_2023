# Local transmission lineages

pacman::p_load(
  ggplot2,
  gginnards,
  plyr,
  dplyr,
  lubridate,
  ape,
  grid,
  gridExtra,
  phylobase,
  stringr,
  treeio,
  tidytree,
  egg,
  igraph
  )

#install_github("YuLab-SMU/ggtree")
library(ggtree)


## functions
source("LTL_functions.R")

# Tree input --------------------------------------------------------------

## Import tree in two formats
tree <- read.beast("../MJ_MCC_CA_DS.trees")
tre  <- read.nexus("../MJ_MCC_CA_DS.trees")

mrsd=2023.605
mrsd=as.Date(date_decimal(mrsd))

# Sometimes get a 50/50 split in the BEAST calculation of a node location - example of a manual fix for this
# tree@data$state[tree@data$state=="SOUTH-CORRIDOR-AF+SOUTH-CORRIDOR-PK"] <- "SOUTH-CORRIDOR-PK"


## extract the ancestral node location and number for each node and add this to the dataframe 
tree@data$ancestor <- mapply(ancestor, tree@data$node)

ancestor_loc <- tree@data[,c("node","state")]
names(ancestor_loc) <- c("ancestor","ancestor_loc")
tree@data <- as_tibble(join(tree@data,ancestor_loc, by="ancestor"))

#add branch length to node age to get the branching point
tree@data$branching <- as.numeric(tree@data$height) + as.numeric(tree@data$length)/2
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

## Grey vertical stripes on final plot
low_season <- matrix(c("2006-01-01", "2007-01-01",
                       "2008-01-01", "2009-01-01",
                       "2010-01-01", "2011-01-01",
                       "2012-01-01", "2013-01-01",
                       "2014-01-01", "2015-01-01",
                       "2016-01-01", "2017-01-01",
                       "2018-01-01", "2019-01-01",
                       "2020-01-01", "2021-01-01",
                       "2022-01-01", "2023-01-01"),ncol=2,byrow=TRUE)

colnames(low_season) <- c("Start","End")
low_season <- as.data.frame(low_season)
low_season$Start <- as.Date(low_season$Start)
low_season$End <- as.Date(low_season$End)

# prepare to run
tstcl <- ggplot_build(tr) ## need to build to get the coordinates of each node
tipsh <- c(21,22)
names(tipsh) <- c("AFP","ES")
treeroot <- rootnode(tre)

subtreelist <- list()
treedata <- list()

## apply subtree custom fuction
## Cannot deal with plotting the root with this function so seperately done after
treedat <- lapply(1:nrow(subtree),subtrefunc, subtree=subtree, tree=tree)

## split the lists, first 9 elements to the ggplot and the rest to my table
subtreelist <- list()
treedata <- list()
for(i in 1:length(treedat)){
  subtreelist[[i]] <- append_layers(treedat[[i]][[1]], geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=decimal_date(Start), xmax=decimal_date(End), ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")
  class(subtreelist[[i]]) <- class(treedat[[i]][[1]])
  treedata[[i]] <- treedat[[i]][[2]]
}

downstreamtree <- data.frame(treedata)
downstreamtree <- data.table::transpose(downstreamtree)


xmin <- decimal_date(min(tree@data$branching))
xmax <- decimal_date(mrsd)
## parent list of any node in the tree gives the root as the last entry
treroot <- last(ancestors(phylo4(tre), 1))
nodenum <- treroot
tiplist <- tidytree::offspring(tre,treroot)
# find all nodes to collapse
tocollapse <- unlist(lapply(as.list(tiplist), off))
# find all of the tips which are ancestors of these
collapse <- as.vector(unlist(descendants(phylo4(tre), tocollapse, type = c("tips"))))
droptree <- custom.drop.tip(tre, tip=collapse, tree=tree)
tipnames <- droptree$tip.label
tipdates <- date_decimal(sapply( strsplit( tipnames, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )}))
tipdate <- max(tipdates, na.rm = T)

# Same analysis for the root ----------------------------------------------

tipshapes <-     sapply( strsplit( tipnames, '\\_' ), function(x){  head(x,2)[2]})
tipshapes[which(tipshapes %!in% c("AFP", "ES"))] <- NA

edges <- as.data.frame(tre[["edge"]])
rootlng <- 1
droptree$root.edge <- tre[["edge.length"]][rootlng]

length(tipshapes) <- length(droptree$tip.label)+droptree$Nnode
roottips <- length(tipshapes[!is.na(tipshapes)])

rootdate <- as.character(as.Date(tree@data[tree@data$node==treroot,]$branching))

treelist <- list()
treelist[[1]] <- ggplot()
treelist[[2]] <- ggtree(droptree, mrsd=tipdate,colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour) +
  scale_y_reverse()+
  coord_cartesian(xlim=c(xmin-1.5, xmax+.2), clip = 'off')+ #clip off v. important to allow the plots to overlap
  geom_tippoint(aes(subset = !grepl("collapse", label), shape=tipshapes), fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour, size=2) +
  scale_shape_manual(values=tipsh)+
  geom_tiplab(aes(subset = grepl("collapse", label),label = "\u2B9E", color = c(sapply( strsplit( droptree$tip.label, '\\_' ), function(x){ tail(x,1)[1] }), rep("NA", Nnode2(droptree)-Ntip(droptree)))), show.legend = F, offset = -.02, vjust = 0.4) +
  scale_color_manual(values = getPalette, drop=FALSE, na.value="white") +
  geom_rootedge(colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour)+
  theme_tree2() +
  geom_rootpoint(shape=23, fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$parent,]$colour, position = position_nudge(x=-droptree$root.edge), size=2)+
  theme_bw()+
  theme(legend.position = "none",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill = "transparent"),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_rect(fill = "transparent"),
        plot.title = element_blank(),
        legend.spacing = unit(0,"cm"),
        plot.margin=unit(c(-0.095,0,-0.095,0), "lines"))
treelist[[2]] <- append_layers(treelist[[2]], geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=decimal_date(Start), xmax=decimal_date(End), ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")

treelist <- append(treelist, subtreelist)

treaxis <- ggtree(tree, mrsd=mrsd, aes(colour=state)) +
  coord_cartesian(xlim=c(xmin-1.5, xmax+.2), clip = 'off')+ 
  scale_y_reverse()+
  scale_color_manual(values = getPalette, drop=FALSE, na.value="white")+
  geom_tippoint()+
  theme_tree2() +
  theme(legend.position = "none")+
  geom_rootpoint()


gTable <- ggplot_gtable(ggplot_build(treaxis))
grid.newpage()
groblst <- list()
groblst[[1]] <- as_ggplot(gTable$grobs[[7]])
treetst <- append(treelist, groblst)
tipls <-  as.character(table(factor(tipshapes, levels=c("AFP", "ES"))))


# Downstream analysis -----------------------------------------------------

downstreamtree[nrow(downstreamtree)+1,] <- c(NA, tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==treroot,]$colour, rootdate ,max(tipdates),nTips(droptree),tipls, toString(tipnames))

## convert the colours in the list into locations
downstreamtree[,1] <- unlist(lapply(downstreamtree[,1], locatefunc)) 
downstreamtree[,2] <- unlist(lapply(downstreamtree[,2], locatefunc))

#number of importations per area
imp <- as.data.frame(table(downstreamtree$V2))

sum1 <- ggplot(downstreamtree, aes(factor(downstreamtree$V2), fill=downstreamtree$V2)) +
  geom_bar(stat="count", position = "dodge") + 
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  theme(legend.position = "none",
        #axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  xlab("Region")+
  ylab(paste("Number of importations"))


#number of exportations per region
downstreamtree$V1 <- factor(downstreamtree$V1, levels=sort(unique(downstreamtree$V2)))
dta <- as.data.frame(table(downstreamtree$V1))
sum5 <- ggplot(dta, aes(x=dta$Var1,y=dta$Freq, fill=dta$Var1)) +
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  theme(legend.position = "none",
        #axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  xlab("Region")+
  ylab(paste("Number of exportations"))

#proportion of these importations which come from each of the other locations
sum2 <- ggplot(subset(downstreamtree, !is.na(downstreamtree$V1)), aes(factor(V2), fill=V1)) +
  geom_bar(stat="count", position = "fill", na.rm = TRUE) + 
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  theme(legend.position = "none",
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.title.x=element_blank())+
  guides(fill = guide_legend(title="Importation from",title.position="top", title.hjust = 0.5))+
  xlab("Region")+
  ylab(paste("Proportion of importations \n from each region")
  )

#stripchart of the survival time of lineages in each area
downstreamtree$V3 <- as.Date(downstreamtree$V3)
downstreamtree$V4 <- as.Date(downstreamtree$V4)
downstreamtree$diff <- as.numeric(downstreamtree$V4-downstreamtree$V3)
#root
rootlin <- as.numeric(as.Date(max(tipdates))-as.Date(min(tree@data$branching)))


sum3 <- ggplot(downstreamtree, aes(x=diff)) + 
  geom_histogram(aes(fill=V2), binwidth = 20)+
  facet_grid(downstreamtree$V2 ~.)+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  scale_x_continuous(expand = c(0,0))+
  guides(fill = guide_legend(title="Region",title.position="top"))+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  ylab("Count")+
  xlab("Lineage survival (days)")

#similar representation of the number of points per tree in each area
downstreamtree$V5 <- as.numeric(downstreamtree$V5)

sum4 <- ggplot(downstreamtree, aes(x=V5)) + 
  geom_histogram(aes(fill=V2), binwidth = 2)+
  facet_grid(downstreamtree$V2 ~.)+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  guides(fill = guide_legend(title="Region",title.position="top", title.hjust = 0.5, nrow = 2))+
  theme(legend.position = "bottom",
        legend.justification="center",
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  ylab("Count")+
  xlab("Number of tips")
downstreamtree$V6 <- as.numeric(downstreamtree$V6)


# AFP only
AFP_tips <- ggplot(downstreamtree, aes(x=V6)) + 
  geom_histogram(aes(fill=V2), binwidth = 2)+
  facet_grid(downstreamtree$V2 ~.)+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  #scale_y_log10()+
  #annotation_logticks(sides="l")+
  guides(fill = guide_legend(title="Region",title.position="top", title.hjust = 0.5, nrow = 2))+
  theme(legend.position = "bottom",
        legend.justification="center",
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  ylab("Count")+
  xlab("Number of AFP tips")

grid_arrange_shared_legend(sum4, AFP_tips, sum3, arrangeGrob( sum1,sum5,sum2, ncol=1, nrow=3, heights = c(1/4, 1/4, 1/2)), ncol=4, nrow=1, position = "bottom")


# Plot split trees --------------------------------------------------------
downstreamtree$ID <- 2:(nrow(downstreamtree)+1)
downstreamtree[downstreamtree$ID==nrow(downstreamtree)+1,]$ID <- 1
downstreamtree <- downstreamtree[order(downstreamtree$ID),]
rownames(downstreamtree) <- NULL

require(egg)
tiplist <- as.list(as.numeric(downstreamtree$V5))

alltips <- unlist(tiplist)

##rescale with tanh to ensure the plots are not dominated by large or small subtrees but there is still some scaling with size.
alltips <- 1-2/((exp(1)^(.1*alltips)) +1)
alltips[length(alltips)+1] <- 0.4
spacing <- alltips/sum(alltips)

clusternames <- c("East Pakistan"="EAST-PAK", "South Corridor Afghanistan"="SOUTH-CORRIDOR-AF", "South Corridor Pakistan"="SOUTH-CORRIDOR-PK", "North Corridor Afghanistan"="NORTH-CORRIDOR_ AF",
                  "North Corridor Pakistan"="NORTH-CORRIDOR-PK", "West Afghanistan"="WEST-AFG", "Central Pakistan"="CENTRE-PAK", "Central Corridor Afghanistan"="CENTRAL-CORRIDOR-AF", "Endemic Zone Pakistan"="ENDEMIC-ZONE", 
                  "Karachi"="KARACHI", "Sindh"="SINDH"  , "South Punjab"="SOUTH-PUNJAB", "North Afghanistan"="NORTH-AFG", "Khyber Pakhtunkhwa"="KP", "Gilgit Baltistan"="GB", "Central Afghanistan"="CENTRE-AFG")
cluster_shrt <- unique(downstreamtree$V2)



## EAST PAK
plttrees <- treetst

plttrees[[1]] <- ggplot() + 
  coord_cartesian(xlim=c(xmin-1.5, xmax+.2), ylim=c(0,2), clip = 'off')+ 
  theme(legend.position = "none",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill = "transparent"),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_rect(fill = "transparent"),
        plot.title = element_blank(),
        legend.spacing = unit(0,"cm"),
        plot.margin=unit(c(-0.15,0,-0.15,0), "lines"))
plttrees[[1]] <- append_layers(plttrees[[1]], geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=decimal_date(Start), xmax=decimal_date(End), ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")

cluster_block <- list()
for (i in 1:length(cluster_shrt)){
  plttrees_tips <- c(0.4,alltips)
  plttrees_spacing <- plttrees_tips[c(1,downstreamtree[downstreamtree$V2==cluster_shrt[i],]$ID+1,length(plttrees_tips))]#/sum(plttrees_tips[c(1,downstreamtree[downstreamtree$V2==cluster_shrt[i],]$ID+1,length(plttrees_tips))])
  arrange <- ggarrange(plots = plttrees[c(1,downstreamtree[downstreamtree$V2==cluster_shrt[i],]$ID+1,length(treetst))], heights = plttrees_spacing, ncol=1,padding = unit(0, "line"))
  cluster_block[[i]] <- annotate_figure(arrange, top = grobTree( rectGrob(gp=gpar(fill="grey"), height = unit(3, "npc")), textGrob(names(which(clusternames==cluster_shrt[i])), gp=gpar(fontsize=14, col="black"), vjust = 0.8)))
}

