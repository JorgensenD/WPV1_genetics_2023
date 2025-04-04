# Figure 4
pacman::p_load(
  ggplot2
)

# this summary data is generated with 1_LTL_2023_multitree
load("./Plot_data/Figure4.RData")

# Bounds over 300 samples - in supplement
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



# plot just median
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
