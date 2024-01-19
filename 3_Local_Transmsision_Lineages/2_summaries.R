# plots based on the split trees - require tree splitting code to be run first

## make a full tree with tip shapes
tipshapes_n <- read.table(text = as.character(tre$tip.label), sep = "_")$V2
length(tipshapes_n) <- length(tre$tip.label)+tre$Nnode
tipsh <- c(21,22)
names(tipsh) <- c("AFP","ES")

earliest_date = as.Date("2012-01-07")

tr2 <- ggtree(tree, mrsd=mrsd, aes(colour=state), lineend = "square") +
  scale_y_reverse()+
  scale_color_manual(values = getPalette, drop=FALSE, na.value="white", name="Regions")+
  geom_tippoint(aes( fill=state, shape = tipshapes_n), colour="black")+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white", name="Regions")+
  scale_shape_manual(values = tipsh, name = "Tip shape")+
  coord_cartesian(xlim= c(xmin, xmax+.2), clip = 'off')+
  geom_vline(xintercept = decimal_date(earliest_date), size = 1, color = "darkorange")+
  theme_tree2() +
  theme(legend.position=c(.17, .45),
        axis.text.x=element_text(size=11),
        legend.key.width=unit(.6, "cm"),
        legend.spacing.x = unit(.3, 'cm'),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", 
                                         colour ="black"),
        legend.margin= margin(6,6,6,10) ,
        plot.margin=unit(c(0.1,0.1,-0.5,0.4), "cm"))+
  guides(shape = guide_legend(order = 1),
         colour = guide_legend(override.aes = list(shape = NA, linewidth=5)))
#geom_rootpoint()+
#geom_vline(xintercept = 2018, colour="red", linetype=2, size=1)


tr2 <- append_layers(tr2, geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=decimal_date(Start), xmax=decimal_date(End), ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")


png("./2023_tree.png",width = 7, height = 7 , units = 'in', res = 200)
tr2
dev.off()


## circos plot of movements during and after the outbreak period

# order palette how I want the layout
# then reassign the N/S
pal2 <- getPalette[c(16, 9, 2, 1, 5, 10, 11, 8, 6, 4, 15, 12, 7, 3, 14 ,13)]


grouping = c("North", "North", "North", "North", "North", "North", "North", "North", "North", "North", "South", "South", "South", "South", "South", "South")
names(grouping) <- names(pal2)
levels(grouping)= c("North", "South")


#during outbreak 1
downstream_outbreak_1 <- downstreamtree %>%
  filter(V3>="2013-01-01" & V3<"2015-01-01")

outbreak_1_cons <- graph.edgelist(as.matrix(downstream_outbreak_1[c("V1", "V2")]), directed = FALSE)
outbreak_1_adj <-get.adjacency(outbreak_1_cons)

# during outbreak 2
downstream_outbreak_2 <- downstreamtree %>%
  filter(V3>="2019-01-01" & V3<"2021-01-01")

outbreak_2_cons <- graph.edgelist(as.matrix(downstream_outbreak_2[c("V1", "V2")]), directed = FALSE)
outbreak_2_adj <-get.adjacency(outbreak_2_cons)

# rewrite gap function as we have less than n sectors plotted
small_gap = 1
big_gap = 20
outbreak2_gap = (14*small_gap)+(big_gap*2)
outbreak1_gap = 12*small_gap

percent = outbreak1_gap/outbreak2_gap

blank.degree = (360 - outbreak2_gap) * (1 - percent)
out = (blank.degree - outbreak1_gap)/2

circos.clear()

circos.par(start.degree = 125)

png(paste0(directory_interest, "outbreak_mvmt_13-14.png"), units="in", width=5, height=5, res=350)
chordDiagram(downstream_outbreak_1[,c("V1", "V2")], grid.col = pal2, transparency = 0.6, group = grouping, big.gap = out, small.gap = 1)
dev.off()
circos.clear()

circos.par(start.degree = 125)

svg(paste0(directory_interest, "outbreak_mvmt_13-14.svg"), width=5, height=5)
chordDiagram(downstream_outbreak_1[,c("V1", "V2")], grid.col = getPalette, transparency = 0.6,  group = grouping, big.gap = out, small.gap = 1)
dev.off()
circos.clear()


circos.par(start.degree = 135)

png(paste0(directory_interest, "outbreak_mvmt_19-20.png"), units="in", width=5, height=5, res=350)
chordDiagram(downstream_outbreak_2[,c("V1", "V2")], grid.col = getPalette, transparency = 0.6, group = grouping, big.gap = 20, small.gap = 1)
dev.off()
circos.clear()

circos.par(start.degree = 135)

svg(paste0(directory_interest, "outbreak_mvmt_19-20.svg"), width=5, height=5)
chordDiagram(downstream_outbreak_2[,c("V1", "V2")], grid.col = getPalette, transparency = 0.6,  group = grouping, big.gap = 20, small.gap = 1)
dev.off()
circos.clear()



#since outbreak
downstream_recent <- downstreamtree %>%
  filter(V3>="2021-01-01")

png(paste0(directory_interest, "outbreak_mvmt_recent.png"), units="in", width=5, height=5, res=350)
chordDiagram(downstream_recent[,c("V1", "V2")], grid.col = getPalette, transparency = 0.6,  group = grouping, big.gap = 20)
dev.off()
circos.clear()

svg(paste0(directory_interest, "outbreak_mvmt_recent.svg"), width=5, height=5)
chordDiagram(downstream_recent[,c("V1", "V2")], grid.col = getPalette, transparency = 0.6,  group = grouping, big.gap = 20)
dev.off()
circos.clear()




# plot outbreaks as matrices


## matrices not based on the adjacency matrix
downstream_outbreak_1$V2 <- factor(downstream_outbreak_1$V2, levels = levels(downstream_outbreak_1$V1))

plot_outbreak1 <- data.frame(table(downstream_outbreak_1[,c("V1", "V2")])) %>%
  ggplot(aes(x = V2, y = V1, fill=Freq, label=ifelse(Freq !=0, Freq, ""))) +
  geom_raster() +
  scale_fill_gradient(low="lightblue", high="red", limits = c(1,55), na.value = "white") +
  geom_text() +
  theme_bw() +
  labs(x="Imported to",
       y="Exported from", 
       title = "2013-2014",
       fill = "Number of\nmovements") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

downstream_outbreak_2$V2 <- factor(downstream_outbreak_2$V2, levels = levels(downstream_outbreak_2$V1))

plot_outbreak2 <- data.frame(table(downstream_outbreak_2[,c("V1", "V2")])) %>%
  ggplot(aes(x = V2, y = V1, fill=Freq, label = ifelse(Freq !=0, Freq, ""))) +
  geom_raster() +
  scale_fill_gradient(low="lightblue", high="red", limits = c(1,55), na.value = "white") +
  geom_text() +
  theme_bw() +
  labs(x="Imported to",
       y="Exported from", 
       title = "2019-2020",
       fill = "Number of\nmovements") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(panel.grid = element_blank(), legend.position = "none")


plot_grid(plot_grid(plot_outbreak1+ theme(legend.position="none"), 
                    plot_outbreak2+ theme(legend.position="none"), labels = c('A', 'B')),
          get_legend(plot_outbreak1), ncol = 1, rel_heights = c(.9,.1))

ggsave("movement_matrix.png", height = 6, width = 10)


save(downstreamtree, file = "downstreamtree_all.rdata")




