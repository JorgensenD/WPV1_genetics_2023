# smoothed bayesian skyline - mlesky with a set of trees from beast.

## CRAN packages
pacman::p_load(
  ggplot2,
  gginnards,
  cowplot
)

## Other packages
#devtools::install_github('emvolz-phylodynamics/mlesky')
library(mlesky)


## Read consensus tree
pv_trees <- ape::read.nexus("../MJ_MCC_CA_DS.trees")

# Skygrowth model ---------------------------------------------------------

# extract dates from tree
datelabs <- as.numeric(sapply(strsplit(pv_trees$tip.label, "\\_"), "[", 4))
names(datelabs) <- pv_trees$tip.label

## find the optimal res value
res=optim_res_aic(pv_trees,ncpu=6) 

skygrowth <- mlskygrid(pv_trees, sampleTimes = datelabs, tau = 20, ncpu = 6,res=50,model=3)
skygrowthboot <- parboot(skygrowth, dd=F)
mlesky:::.neplot(skygrowthboot, logy = F, ggplot = T)

## Alternate skyline method
# skykappa <- mlskygrid(pv_trees, sampleTimes = datelabs, tau = 20, ncpu = 6,res=50,model=1)
# skykappaboot <- parboot(skykappa, dd=F)
# mlesky:::.neplot(skykappaboot, logy = F, ggplot = T) 



# Background grey vertical stripes ----------------------------------------


low_season <- matrix(c("2006-01-01","2007-01-01",
                       "2008-01-01","2009-01-01",
                       "2010-01-01","2011-01-01",
                       "2012-01-01","2013-01-01",
                       "2014-01-01","2015-01-01",
                       "2016-01-01","2017-01-01",
                       "2018-01-01","2019-01-01",
                       "2020-01-01", "2021-01-01",
                       "2022-01-01", "2023-01-01"),ncol=2,byrow=TRUE)

colnames(low_season) <- c("Start","End")
low_season <- as.data.frame(low_season)
low_season$Start <- as.Date(low_season$Start)
low_season$End <- as.Date(low_season$End)



# Plotting data -----------------------------------------------------------

plotdata <- data.frame(t = as.Date(date_decimal(skygrowthboot$time)), # time
                       nelb = skygrowthboot$ne_ci[,1],                # lower ci bound on Ne
                       nemed = skygrowthboot$ne_ci[,2],               # mean Ne
                       neub = skygrowthboot$ne_ci[,3],                # upper ci bound on Ne
                       glb = skygrowthboot$growthrate_ci[,1],         # lower ci bound on growth rate
                       gmed = skygrowthboot$growthrate_ci[,2],        # mean growth rate
                       gub = skygrowthboot$growthrate_ci[,3]).        # upper ci bound on growth rate


# Plot skyline ------------------------------------------------------------

ggsky <- ggplot(plotdata)+
  geom_line(aes(x = t, y = nemed))+
  geom_ribbon(aes(ymin = nelb, ymax = neub, x = t), fill = "blue", alpha = 0.2) +
  scale_y_continuous("Ne(t)", breaks = seq(0,100, by = 25), expand = c(0,0)) +
  theme_cowplot()+  
  coord_cartesian(xlim= c(as.Date(date_decimal(xmin)), as.Date(date_decimal(xmax+.2))))+
  geom_vline(xintercept = earliest_date, size = 1, color = "darkorange")+
  theme(plot.margin=unit(c(-0.2,0.1,0.1,0.4), "cm"),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

skyplot <- append_layers(ggsky, 
                         geom_rect(data=low_season,
                                   inherit.aes=FALSE, 
                                   aes(xmin=Start, xmax=End, ymin=-Inf, ymax=+Inf), 
                                   fill='lightgrey'), 
                         position = "bottom")


# Plot growth rate --------------------------------------------------------

gggrowth <- ggplot(plotdata)+
  geom_line(aes(x = t, y = (1/10*gmed)+1 ))+
  geom_ribbon(aes(ymin = (1/10*glb)+1, ymax = (1/10*gub)+1, x = t), fill = "blue", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype =2, color = "red") +
  scale_y_continuous("R(t)", breaks = seq(0.7, 1.3, by = 0.1), expand = c(0,0)) +
  coord_cartesian(xlim= c(as.Date(date_decimal(xmin)), as.Date(date_decimal(xmax+.2))))+
  geom_vline(xintercept = earliest_date, size = 1, color = "darkorange")+
  theme_cowplot() +
  theme(plot.margin=unit(c(-0.2,0.1,0.1,0.4), "cm"),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

growthplot <- append_layers(gggrowth, 
                            geom_rect(data=low_season,
                                      inherit.aes=FALSE, 
                                      aes(xmin=Start, xmax=End, ymin=-Inf, ymax=+Inf), 
                                      fill='lightgrey'), 
                            position = "bottom")


# Combine -----------------------------------------------------------------

plot_grid(skyplot, growthplot, ncol = 1, align = "v")  
save(plotdata, file = "R_estimate.rdata")