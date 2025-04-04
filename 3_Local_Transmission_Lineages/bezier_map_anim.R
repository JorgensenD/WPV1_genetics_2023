# Animated movement

pacman::p_load(
  sf,
  ggplot2,
  arcgislayers,
  parallel,
  treeio,
  lubridate,
  plyr,
  dplyr,
  bezier,
  matlib,
  foreach,
  gtools,
  magick
)

adm2_Pk_Af <- st_read(dsn = "./all_20230810/tipswap/adm2_Pk_Af.geojson", layer = "adm2_Pk_Af.geojson")

furl_0 <- "https://services.arcgis.com/5T5nSi527N4F7luB/arcgis/rest/services/POLIO_ADMINISTRATIVE_BOUNDARIES/FeatureServer/6" 
# Open connection
country_fl_0 <- arc_open(furl_0)
# Subset with SQL
adm0_Pk_Af <- arc_select(country_fl_0,
                         where = "ISO_2_CODE = 'AF' or ISO_2_CODE = 'PK'")


adm0_Pk_Af <- adm0_Pk_Af[as.Date(adm0_Pk_Af$ENDDATE)>=Sys.Date(),]


#### color palette ----
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

#### curve plotting function ----
plotbezier <- function(loc1,loc2,ID,height,mapdata,name_column){
  #get x and y of each polygon and add to a dataframe
  coords <- list()
  print(paste0(loc1,loc2))
  for(i in 1:nrow(mapdata)){
    # sf tibble so convert to df for easier manipulation
    name <- as.data.frame(mapdata)[i,name_column]
    #get midpoint of each column
    coord <- c(st_coordinates(st_centroid(mapdata[i,])))
    #add name of the region to the coordinates
    coords[[name]] <- coord
  }
  if(loc2==loc1){
    # control points for internal transmission (4 to give loop)
    controlx1 <- coords[[loc1]][1]+sample(c(-1,1),1)*height
    controly1 <- coords[[loc1]][2]
    controlx2 <- coords[[loc1]][1]
    controly2 <- coords[[loc1]][2]+sample(c(-1,1),1)*height
    x <- c(coords[[loc1]][1], controlx1, controlx2, coords[[loc2]][1])
    y <- c(coords[[loc1]][2], controly1, controly2, coords[[loc2]][2])
  } else{
    # control point calculations for lines linking two locations
    m <- (coords[[loc2]][2] - coords[[loc1]][2])/(coords[[loc2]][1] - coords[[loc1]][1])
    b <-  coords[[loc2]][2] - m*coords[[loc2]][1]
    #print the equation
    #if (b>0){bprint <- paste("+",b, sep="")} else {bprint <- paste(b)}
    #paste("y=",m,"x",bprint, sep="")
    #intercept +- height to give parallel line (ifelse so different for each direction)
    ifelse(coords[[loc2]][1]-coords[[loc1]][1]<=0, bplus <- b - height, bplus <- b + height) 
    #if (bplus>0){bplusprint <- paste("+",bplus, sep="")} else {bplusprint <- paste(bplus)}
    #paste("parallel line: y=",m,"x",bplusprint, sep="")
    #find the midpoint between the locations
    midx <- (coords[[loc2]][1] + coords[[loc1]][1])/2
    midy <- (coords[[loc2]][2] + coords[[loc1]][2])/2
    midm <- -1/m
    midb <- midy-midm*midx
    #if (midb>0){midbprint <- paste("+",midb, sep="")} else {midbprint <- paste(midb)}
    #paste("perpendicular line: y=",midm,"x",midbprint, sep="")
    # find intercept of perpendicular and parallel lines
    A <- matrix(c(1,1,-midm,-m),2,2)
    B <- c(midb,bplus)
    #showEqn(A,b)
    controlyx <- solve(A,B)
    # standardise the height of the arc so it doesnt depend on the angle of the line
    d <- sqrt((controlyx[2]-midx)^2+(controlyx[1]-midy)^2)
    t <- height/d
    controlx <- ((1-t)*midx+t*controlyx[2])
    controly <- ((1-t)*midy+t*controlyx[1])
    
    #set up 3 control points for the curve (start, middle and end)
    x <- c(coords[[loc1]][1], controlx, coords[[loc2]][1])
    y <- c(coords[[loc1]][2], controly, coords[[loc2]][2])
  }
  coordsp <- cbind(x,y)
  curve <- pointsOnBezier(n=20,method='evenly_spaced',p=coordsp)
  
  curvegg <- data.frame(curve$points)
  
  curvegg$ID <- as.numeric(ID)
  curvegg$START <- loc1
  Curve <- list(curvegg)
  names(Curve)[[1]] <- paste(loc1, loc2, sep="_")
  return(Curve)
}

#### load tree ----
directory_interest <- "./all_20230810"
tree <- read.beast(paste0(directory_interest, "/MJ_MCC_CA_DS.trees"))
mrsd=as.Date(date_decimal(2023.605))

tree@data$state[tree@data$state=="SOUTH-CORRIDOR-AF+SOUTH-CORRIDOR-PK"] <- "SOUTH-CORRIDOR-PK"

             
## apply a lookup function to each row of the data slot
ancestor <- function(node){
  if(!is.na(node)){
    return(tree@phylo[["edge"]][tree@phylo[["edge"]][,2]==node,][1])
  }
}
tree@data$ancestor <- mapply(ancestor, tree@data$node)

ancestor_loc <- tree@data[,c("node","state")]
names(ancestor_loc) <- c("ancestor","ancestor_loc")
tree@data <- as_tibble(join(tree@data,ancestor_loc, by="ancestor"))

#add branch length to node age to get the branching point
tree@data$branching <- as.numeric(tree@data$height) + as.numeric(tree@data$length)/2
tree@data$branching <- date_decimal(decimal_date(mrsd)-tree@data$branching)

region_col="epiblock" #name of the column containing the region names you want to use

# match map to the old names in the tree
model_map <- adm2_Pk_Af %>%
  group_by(epiblock) %>%
  mutate(epiblock = case_when(epiblock == "CENTRAL-AF" ~ "CENTRE-AFG",
                              epiblock == "CENTRAL-PK" ~ "CENTRE-PAK",
                              epiblock == "EAST-PK" ~ "EAST-PAK",
                              epiblock == "CENTRAL-CORRIDOR-PK" ~ "ENDEMIC-ZONE",
                              epiblock == "NORTH-AF" ~ "NORTH-AFG",
                              epiblock == "WEST-AF" ~ "WEST-AFG",
                              TRUE~epiblock)) %>%
  dplyr::summarize(geometry = st_union(geometry))




locationmatrix <- unique(na.omit(tree[,c("ancestor_loc","state")]))
ID <- c(1:nrow(locationmatrix))
locationmatrix <- cbind(locationmatrix, ID)
# these are named characters for some reason
locationmatrix$ancestor_loc <- as.character(locationmatrix$ancestor_loc)
locationmatrix$state <- as.character(locationmatrix$state)

## generate. curves
ptm <- proc.time()
out <- mapply(plotbezier,locationmatrix[,1], locationmatrix[,2], locationmatrix[,3],  MoreArgs = list(height=2, mapdata=model_map, name_column=1))
proc.time() - ptm

plot_bg <- ggplot() +
  geom_sf(data = model_map, aes(fill = epiblock), color = "white", linewidth = 0.2, show.legend = F, alpha = 0.2)+
  geom_sf(data = adm0_Pk_Af, fill = NA, linewidth = 0.7, color = "black")+
  scale_fill_manual(values = getPalette)+
  theme_void()

# centroids
coords <- list()
for(i in 1:nrow(model_map)){
  # sf tibble so convert to df for easier manipulation
  name <- as.data.frame(model_map)[i,1]
  #get midpoint of each column
  coord <- c(st_coordinates(st_centroid(model_map[i,])))
  #add name of the region to the coordinates
  coords[[name]] <- coord
}
centroids <- data.frame(do.call(rbind,coords))
centroids$group <- rownames(centroids)



linedata <- as.data.frame(tree@data)
linedata <- linedata[!is.na(linedata$ancestor_loc),]
linedata$branching <- as.Date(linedata$branching)
frames <- seq(from=min(as.Date(linedata$branching)), to=max(as.Date(mrsd)),by='days')
linedata$drawing <- 0  


foreach(j=1:length(frames), .packages=c("ggplot2", "dplyr", "plyr", "lubridate", "RColorBrewer")) %do% {
  ## add 1 to the drawing parameter each day the the date is less than the plotting date
  linedata[linedata$branching <= frames[j],]$drawing <- difftime(frames[j] ,linedata[linedata$branching <= frames[j],]$branching , units = c("days"))+1
  ##CAN DO 200/j * j if we want it to be variable how many days the lines take to travel and use the 200 point beziers
  #counts over 30 do not contribute to the plots
  current_lines <- linedata[linedata$drawing > 0 & linedata$drawing <= 30,]
  if(nrow(current_lines)>0){
    current_lines$ref <- paste(current_lines$ancestor_loc, ".", current_lines$ancestor_loc, "_", current_lines$state, sep="")
    ## look up the line in the list of lines and reference drawing in the number of points of the line to draw...
    current_lines <- current_lines[order(-current_lines$drawing),]
    rownames(current_lines) <- NULL
    pltlsts <- list()
    for (k in 1:nrow(current_lines)){
      if(current_lines[k,]$drawing <=10){
        line <- out[[current_lines[k,]$ref]][1:current_lines[k,]$drawing,]
      } else if (current_lines[k,]$drawing >10 & current_lines[k,]$drawing <=20) {
        line <- out[[current_lines[k,]$ref]][(current_lines[k,]$drawing-10):current_lines[k,]$drawing,]
      } else {
        line <- out[[current_lines[k,]$ref]][(current_lines[k,]$drawing-10):20,]
      }
      line$ID2 <- as.numeric(k)
      pltlsts[[k]] <- line
    }
    pltlst <- do.call(rbind, pltlsts)
    pltlst$size <- 1.3
    bg_list <- pltlst
    bg_list$START <- NA
    bg_list$size <- 1.5
    
    bg_list$ID2 <- bg_list$ID2-0.00001
    tst <- rbind(pltlst, bg_list)
    tst <- tst[order(tst$ID2),]
    tst$ID2 <- ordered(factor(tst$ID2))
    tst$size <- factor(tst$size)
  } else {pltlst <- data.frame()}
  
  bg_month <- as.Date(floor_date(frames[j], "month"))
  if(nrow(pltlst)>0){
    plot <-  plot_bg +
      #geom_path(data=pltlst, aes(y=X2, x=X1, group=ID2), size=2.2, colour="white", lineend="round")+
      geom_path(data=tst, aes(y=X2, x=X1, group=ID2, colour=factor(START), linewidth =size), lineend="round", show.legend = F)+
      geom_point(data = centroids, aes(y=X2, x=X1, group = group, fill = group), shape = 21, color = "white", size = 5, show.legend = F) +
      scale_color_manual(values = getPalette, drop=FALSE, na.value="white")+
      scale_linewidth_manual(values=c(1.3,2))+
      annotate("text", x = Inf, y = Inf, label = format(bg_month, "%Y-%m"), hjust = 2, vjust = 35, color = "black", fontface = "bold")
    
  } else {
    plot <- plot_bg+
      geom_point(data = centroids, aes(y=X2, x=X1, group = group, fill = group), shape = 21, color = "white", size = 5, show.legend = F) +
      annotate("text", x = Inf, y = Inf, label =  format(bg_month, "%Y-%m"), hjust = 2, vjust = 35, color = "black", fontface = "bold")
  }
  # plots will be daily
  ggsave(paste0("PkAf_Anim_", j, ".png"), plot = print(plot), height = 5.3, width = 5, units = "in", device="png", dpi=80, path = paste0(directory_interest, "/anim/"), bg = "white")
   print(plot)
  gc()
  rm(plot)
}



# List all images
image_files <- mixedsort(list.files(path = paste0(directory_interest, "/anim/"), pattern = "*.png", full.names = TRUE))

# Read and combine images into a GIF
img_list <- image_read(image_files[seq(1000,length(image_files),2)])  # Read all images
gif <- image_animate(img_list, fps = 10)  # Create GIF 

# Save GIF
image_write(gif, paste0(directory_interest, "/anim/animated_map_dropframe_lores.gif"))


 