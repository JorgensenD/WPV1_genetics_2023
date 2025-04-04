# Figure 2
pacman::p_load(
  ggtree,
  treeio,
  ggplot2,
  gginnards,
  dplyr,
  lubridate
)

# Region color palette ----
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

# Annual stripes
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
low_season %<>% as_tibble() %>%
  mutate(across(everything(), as.Date),
         across(everything(), decimal_date) 
  )

mrsd = as.Date("2023-08-09")
earliest_date = as.Date("2012-01-07")


# A ----
# Consensus tree tip data cannot be provided due to sensitive metadata.
load("./Plot_data/minimal_tree.RData")

# shapes
tipsh <- c(21,22)
names(tipsh) <- c("AFP","ES")
tipshapes_n <- read.table(text = as.character(cons_tree@phylo$tip.label), sep = "_")$V1
length(tipshapes_n) <- length(cons_tree@phylo$tip.label)+cons_tree@phylo$Nnode


# tree plot
tr2 <- ggtree(cons_tree, mrsd=mrsd, aes(colour=state), lineend = "square") +
  scale_y_reverse() + # top to bottom
  scale_color_manual(values = getPalette, drop=FALSE, na.value="white", name="Regions") +
  geom_tippoint(aes( fill=state, shape = tipshapes_n), colour="black") + # shape and colour of the sampled tips
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white", name="Regions") +
  scale_shape_manual(values = tipsh, name = "Tip shape") +
  coord_cartesian(xlim= c(2005.9, 2023.7+.2), clip = 'off') +
  geom_vline(xintercept = decimal_date(earliest_date), linewidth = 1, color = "darkorange") + # date of earliest sample
  theme_tree2() +
  theme(legend.position=c(.17, .45),
        axis.text.x=element_text(size=11),
        legend.key.width=unit(.6, "cm"),
        legend.spacing.x = unit(.3, 'cm'),
        legend.background = element_rect(fill="white",
                                         linewidth=1, linetype="solid", 
                                         colour ="black"),
        legend.margin= margin(6,6,6,10) ,
        plot.margin=unit(c(0.1,0.4,0.1,0.2), "cm"))+
  guides(shape = guide_legend(order = 1),
         colour = guide_legend(override.aes = list(shape = NA, linewidth=5)))

# add background
append_layers(tr2,   
              geom_rect(data = low_season, 
                        aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf), 
                        fill = "lightgrey", 
                        inherit.aes = F),
              position = "bottom")

# B ----
# Lineages through time
load("./Plot_data/Figure2B.rdata")


ltt <- ggplot(coord_df, aes(x = time2, y = N))+
  geom_line(stat = "summary", fun = mean)+
  stat_summary(fun = mean,
               fun.min = function(N) mean(N) - 2*sd(N), 
               fun.max = function(N) mean(N) + 2*sd(N), 
               geom = "ribbon", 
               alpha = 0.3,
               fill = "blue")+
  scale_y_continuous("Lineages", breaks = seq(0,300, by = 100), expand = c(0,0)) +
  coord_cartesian(xlim= c(2006, 2023.7+.2), ylim = c(0,350))+
  geom_vline(xintercept = 2012, size = 1, color = "darkorange")+
  theme_cowplot() +
  labs(x = "Date") +
  theme(plot.margin=unit(c(-0.24,0.1,0.1,0.4), "cm"))

ltt <- append_layers(ltt, geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=decimal_date(Start), xmax=decimal_date(End), ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")


ltt_no <- ltt +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(),
        axis.title.x = element_blank())


# C ----
# Effective viral population size - calculated from AFP data only
sky <- ggplot(bs_AFP) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), alpha = 0.3, fill = "blue")+
  geom_line(aes(x = time, y = mean), color = "black") +
  scale_y_continuous("Ne(t)", breaks = seq(0,120, by = 40), expand = c(0,0)) +
  coord_cartesian(xlim= c(2006, 2023.7+.2), ylim = c(0,170))+
  geom_vline(xintercept = 2012, size = 1, color = "darkorange")+
  theme_cowplot() +
  labs(x = "Date") +
  theme(plot.margin=unit(c(-0.24,0.1,0.1,0.4), "cm"))

sky <- append_layers(sky, geom_rect(data=low_season,inherit.aes=FALSE, aes(xmin=decimal_date(Start), xmax=decimal_date(End), ymin=-Inf, ymax=+Inf), fill='lightgrey'), position = "bottom")

