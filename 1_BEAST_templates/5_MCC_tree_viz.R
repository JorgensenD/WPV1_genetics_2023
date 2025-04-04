## MCC tree viz
require(ggplot2)
require(ggtree)
require(treeio)


# load the MCC tree
MCC_MJ_tree <- read.beast(paste0(directory_interest, "MJ_AFP_MCC.trees"))

getPalette <- c("#f4766c", "#ffd712", "#ef9b95", "#0a93d3", "#e50767", "#71e244", "#008121", "#ff80d9", "#858585", "#8999ff", "#0000b0",  "#a10052", "#922d92", "#c87dff", "#ff8a00", "#00da94")
names(getPalette) <- c("CENTRAL-CORRIDOR-AF", "CENTRE-AF", "CENTRE-PK", "EAST-PK", "ENDEMIC-ZONE", "GB", "KARACHI", "KP", "NORTH-AF", "NORTH-CORRIDOR-AF", "NORTH-CORRIDOR-PK", "SINDH", "SOUTH-CORRIDOR-AF", "SOUTH-CORRIDOR-PK", "SOUTH-PUNJAB", "WEST-AF")

mrsd=2023.605
mrsd=as.Date(date_decimal(mrsd))

ggtree(MCC_MJ_tree, aes(color = state), mrsd = mrsd) +
  geom_tippoint(aes(fill = state), shape = 21, size =2, color = "black") +
  scale_fill_manual(values = getPalette, name = "Location")+
  scale_color_manual(values = getPalette, name = "Location")+
  theme_tree2() +
 # guides(fill = guide_legend(override.aes = list(size = 5, shape = 22), ncol=2))+
  theme(legend.position = c(.25, .8),
        legend.box.background = element_rect(colour = "black"))

ggsave(paste0("./MCC_root_location_",latestdate,".png"), width = 10, height = 7)
ggsave(paste0("./MCC_root_location_",latestdate,".svg"), width = 10, height = 7)
