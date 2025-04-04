# Figure 3
pacman::p_load(
  ggplot2,
  dplyr,
  lubridate,
  circlize
)

getPalette <- c("#f4766c", "#ffd712",
                "#ef9b95", "#0a93d3",
                "#e50767", "#71e244",
                "#008121", "#ff80d9",
                "#858585", "#8999ff",
                "#0000b0", "#a10052",
                "#922d92", "#c87dff",
                "#ff8a00", "#00da94")
names(getPalette) <- c("CENTRAL-CORRIDOR-AF", "CENTRAL-AF",
                       "CENTRAL-PK", "EAST-PK",
                       "CENTRAL-CORRIDOR-PK", "GB",
                       "KARACHI", "KP",
                       "NORTH-AF", "NORTH-CORRIDOR-AF",
                       "NORTH-CORRIDOR-PK", "SINDH",
                       "SOUTH-CORRIDOR-AF", "SOUTH-CORRIDOR-PK",
                       "SOUTH-PUNJAB", "WEST-AF")


# 3A ----
mvmts <- read.csv("./Plot_data/Figure3A.csv")

mvmts$Location <- factor(mvmts$Location)

## Imports
mvmts %>%
  ggplot(aes(y = reorder(Location, IN))) +
  geom_col(aes(x = IN, fill = Location), show.legend = F) +
  geom_errorbarh(aes(xmax = IN_U, xmin = IN_L), height = 0.4)+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white") +
  theme_bw() +
  scale_x_continuous(expand = expansion(add = c(0, 2)), position = "top") +
  xlab("Importations") +
  ylab(NULL) +
  theme(panel.grid.major.y = element_blank())

## Exports
mvmts %>%
  ggplot(aes(y = reorder(Location, OUT))) +
  geom_col(aes(x = OUT, fill = Location), show.legend = F) +
  geom_errorbarh(aes(xmax = OUT_U, xmin = OUT_L), height = 0.4)+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white") +
  theme_bw() +
  scale_x_continuous(expand = expansion(add = c(0, 2)), position = "top") +
  xlab("Exportations") +
  ylab(NULL) +
  theme(panel.grid.major.y = element_blank())

# 3B + C
# rewrite gap function as we have less than n sectors plotted
small_gap = 1
big_gap = 20
outbreak2_gap = (14*small_gap)+(big_gap*2)
outbreak1_gap = 12*small_gap

percent = outbreak1_gap/outbreak2_gap

blank.degree = (360 - outbreak2_gap) * (1 - percent)
out = (blank.degree - outbreak1_gap)/2
circos.clear()

# palette order
pal2 <- getPalette[c(16, 9, 2, 1, 5, 10, 11, 8, 6, 4, 15, 12, 7, 3, 14 ,13)]

grouping = c("North", "North", "North", "North", "North", "North", "North", "North", "North", "North", "South", "South", "South", "South", "South", "South")
names(grouping) <- names(pal2)
levels(grouping)= c("North", "South")


load("./Plot_data/Figure3C.rdata")
circos.par(start.degree = 125)
chordDiagram(downstream_outbreak_1, grid.col = pal2, transparency = 0.6, group = grouping, big.gap = out, small.gap = 1)
circos.clear()

load("./Plot_data/Figure3D.rdata")
circos.par(start.degree = 135)
chordDiagram(downstream_outbreak_2, grid.col = getPalette, transparency = 0.6, group = grouping, big.gap = 20, small.gap = 1)
circos.clear()



