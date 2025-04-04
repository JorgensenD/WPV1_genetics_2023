pacman::p_load(
  circlize
)

transitions <- read.csv("supported_out.csv", header = T)
colnames(transitions) <- c("From", "To", "Bayes_Factor", "Posterior_Probability", "Median_Transitions", "Lower_HPD", "Upper_HPD")

## Circos plots - whole period
png(paste0(directory_interest, "Supported_mvmt_circ_",latestdate,".png"), units="in", width=5, height=5, res=350)
chordDiagram(transitions[,c("From", "To", "Median_Transitions")], grid.col = getPalette, transparency = 0.6)
dev.off()
circos.clear()

svg(paste0(directory_interest, "Supported_mvmt_circ.svg"), width=5, height=5)
chordDiagram(transitions[,c("From", "To", "Median_Transitions")], grid.col = newpal, transparency = 0.6)
dev.off()
circos.clear()



## Present as barcharts - total movements calculated separately in MJ-DTA and stored in movement_bounds.csv

transitions$To <- factor(transitions$To, levels = names(getPalette))
transitions$From <- factor(transitions$From, levels = names(getPalette))

mvmts <- read.csv("movement_bounds.csv")

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


ggsave("import_bar.png", width = 8, height = 5, units = "in")
ggsave("import_bar.svg", width = 8, height = 5, units = "in")

## Exports
mvmts %>%
  ggplot(aes(y = reorder(Location, OUT))) +
  geom_col(aes(x = OUT, fill = Location), show.legend = F) +
  geom_errorbarh(aes(xmax = OUT_U, xmin = OUT_L), height = 0.4)+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white") +
  theme_bw() +
  scale_x_continuous(expand = expansion(add = c(0, 2)), position = "top") +
  xlab("Importations") +
  ylab(NULL) +
  theme(panel.grid.major.y = element_blank())


ggsave("export_bar.png", width = 8, height = 5, units = "in")
ggsave("export_bar.svg", width = 8, height = 5, units = "in")
