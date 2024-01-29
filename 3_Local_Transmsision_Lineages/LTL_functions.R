# Tree splitting functions

#' @title Find the ancestor of a node from the tree data
#' @param node node number to return ancestor for - also requires the tree to exist as "tree"
#' @return ancestral node number
#' @export
ancestor <- function(node){
  if(!is.na(node)){
    return(tree@phylo[["edge"]][tree@phylo[["edge"]][,2]==node,][1])
  }
}

#' @title colour for the offspring of a node - used in conjuntion with the tidytree offspring function

off <- function(offspring){
  if(tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==offspring,]$parent,]$colour &
     tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==offspring,]$colour!=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==offspring,]$parent,]$colour
  ){  
    return(offspring) 
  }
}


#' @title simple not in function
'%!in%' <- Negate('%in%')


#' @title Tip dropping function extending the ape::drop.tip function to collapse downstream transmission trees
#' @param phy phylogenetic tree of class phylo
#' @param tip tip numbers to drop
#' @param tree tre to plot from
#' @param trim.internal T/F to remove internal subtrees
#' @param  subtree T/F is a subtree of a larger tree
#' @param root.edge length of edge to add to root
#' @param rooted T/F is phylogeny rooted
#' @param collapse.singles a logical specifying whether to delete the internal nodes of degree 2.
#' @param interactive if TRUE the user is asked to select the tips or the node by clicking on the tree which must be plotted. Inhereted from ape and not tested.
#' @return trimmed subtree for plotting.
#' @export

custom.drop.tip <- function (phy, tip, tree, trim.internal = TRUE, subtree = TRUE, root.edge = 0, 
                             rooted = is.rooted(phy), collapse.singles = FALSE, interactive = FALSE) 
{
  if (!inherits(phy, "phylo")) 
    stop("object \"phy\" is not of class \"phylo\"")
  Ntip <- length(phy$tip.label)  # total tips
  if (interactive) {
    cat("Left-click close to the tips you want to drop; right-click when finished...\n")
    xy <- locator()
    nToDrop <- length(xy$x)
    tip <- integer(nToDrop)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    for (i in 1:nToDrop) {
      d <- sqrt((xy$x[i] - lastPP$xx)^2 + (xy$y[i] - lastPP$yy)^2)
      tip[i] <- which.min(d)
    }
  }
  else {
    if (is.character(tip)) 
      tip <- which(phy$tip.label %in% tip) # convert to tip numbers
  }
  out.of.range <- tip > Ntip
  if (any(out.of.range)) {
    warning("some tip numbers were larger than the number of tips: they were ignored")
    tip <- tip[!out.of.range]
  }
  if (!length(tip)) 
    return(phy)
  if (length(tip) == Ntip) {
    if (Nnode(phy) < 3 || trim.internal) {
      warning("drop all tips of the tree: returning NULL")
      return(NULL)
    }
  }
  og_tree <- phy
  wbl <- !is.null(phy$edge.length)
  if (length(tip) == Ntip - 1 && trim.internal) { ## only works if there is only one remaining
    i <- which(phy$edge[, 2] == (1:Ntip)[-tip])   ## drop edges ending with the dropped tips
    res <- list(edge = matrix(2:1, 1, 2), tip.label = phy$tip.label[phy$edge[i, 
                                                                             2]], Nnode = 1L)
    class(res) <- "phylo"
    if (wbl) 
      res$edge.length <- phy$edge.length[i]
    if (!is.null(phy$node.label)) 
      res$node.label <- phy$node.label[phy$edge[i, 1] - 
                                         Ntip]
    return(res)
  }
  if (!rooted && subtree) {
    phy <- root(phy, (1:Ntip)[-tip][1])
    root.edge <- 0
  }
  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- dim(phy$edge)[1]
  if (subtree) {
    trim.internal <- TRUE
    tr <- reorder(phy, "postorder")
    N <- .C(node_depth, as.integer(Ntip), as.integer(tr$edge[, 
                                                             1]), as.integer(tr$edge[, 2]), as.integer(Nedge), 
            double(Ntip + Nnode), 1L)[[5]]
  }
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !logical(Nedge)
  keep[match(tip, edge2)] <- FALSE
  if (trim.internal) {
    ints <- edge2 > Ntip
    repeat {
      sel <- !(edge2 %in% edge1[keep]) & ints & keep
      if (!sum(sel)) 
        break
      keep[sel] <- FALSE
    }
    if (subtree) {
      subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
      keep[subt] <- TRUE
    }
    if (root.edge && wbl) {
      degree <- tabulate(edge1[keep])
      if (degree[ROOT] == 1) {
        j <- integer(0)
        repeat {
          i <- which(edge1 == NEWROOT & keep)
          j <- c(i, j)
          NEWROOT <- edge2[i]
          degree <- tabulate(edge1[keep])
          if (degree[NEWROOT] > 1) 
            break
        }
        keep[j] <- FALSE
        if (length(j) > root.edge) 
          j <- 1:root.edge
        NewRootEdge <- sum(phy$edge.length[j])
        if (length(j) < root.edge && !is.null(phy$root.edge)) 
          NewRootEdge <- NewRootEdge + phy$root.edge
        phy$root.edge <- NewRootEdge
      }
    }
  }
  if (!root.edge) 
    phy$root.edge <- NULL
  phy$edge <- phy$edge[keep, ]
  if (wbl) 
    phy$edge.length <- phy$edge.length[keep]
  TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
  oldNo.ofNewTips <- phy$edge[TERMS, 2]
  if (subtree) {
    i <- which(tip %in% oldNo.ofNewTips)
    if (length(i)) {
      # phy$tip.label[tip[i]] <- "[1_tip]"
      phy$tip.label[tip[i]] <- paste("collapse",sapply( strsplit( phy$tip.label[tip[i]], '\\_' ), function(x){ tail(x,1)[1] }), sep = "_")
      tip <- tip[-i]
    }
  }
  n <- length(oldNo.ofNewTips)
  phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
  if (length(tip)) 
    phy$tip.label <- phy$tip.label[-tip]
  if (subtree || !trim.internal) {
    node2tip <- oldNo.ofNewTips[oldNo.ofNewTips > Ntip]
    length(node2tip)
    new.tip.label <- if (!length(node2tip)) {
      character(0)
    } else if (subtree) {
      
      paste("collapse", lapply(node2tip, function(x) tree@data[tree@data$node==MRCA(tree, og_tree$tip.label[offspring(og_tree, x, type = "tips")]),]$state), sep = "_")
      
    } else {
      if (is.null(phy$node.label)) 
        rep("NA", length(node2tip))
      else phy$node.label[node2tip - Ntip]
    }
    phy$tip.label <- c(phy$tip.label, new.tip.label)
  }
  if(any(grepl("collapse", phy$tip.label))){
    # half any collapsed branches
    edgenos <- which(phy$edge[,2] %in% which(grepl("collapse", phy$tip.label)))
    phy$edge.length[edgenos] <- phy$edge.length[edgenos]/2
  }
  phy$Nnode <- dim(phy$edge)[1] - n + 1L
  newNb <- integer(Ntip + Nnode)
  newNb[NEWROOT] <- n + 1L
  sndcol <- phy$edge[, 2] > n
  newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
  phy$edge[, 1] <- newNb[phy$edge[, 1]]
  storage.mode(phy$edge) <- "integer"
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[which(newNb > 0) - 
                                       Ntip]
  if (collapse.singles) 
    phy <- collapse.singles(phy)
  phy
}





#' @title Iterative plotting of the subtrees i can't deal with the root
#' @param i index of subtree to plot
#' @param subtree list of subtrees from which I is pulled
#' @param tree phylogenetic tree to be used
#' @return list of plotted subtrees
#' @export

subtrefunc <- function(i,subtree, tree){ 
  nodenum <<- as.numeric(subtree[i,]$node)
  offspring <- NULL
  xmin <- decimal_date(min(tree@data$branching))
  xmax <- decimal_date(mrsd)
  if(nrow(tstcl[["data"]][[2]][tstcl[["data"]][[2]]$parent==nodenum,])>0){
    #droptips from different locations
    offspring <- tidytree::offspring(tre,subtree[i,]$node)  # find all offspring of the node of interest
    
    ## get tips
    tips <- offspring[offspring<=nTips(tre)]
    parentlist <- list()
    for(k in 1:length(tips)){
      parent <- tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tips[k],]$parent
      ancest <- parent
      while(treeroot%!in%parent){
        ancest <- tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==ancest,]$parent
        parent <- c(parent, ancest)
      }
      parentlist[[k]] <- parent
    }
    tocollapse <- unlist(lapply(as.list(offspring), off))#unlist the offspring of interest (root of subclades)
    
    
    ## use parentlist instead of working out the tips
    
    
    collapsetips <- list()
    for(j in 1:length(tips)){
      if(any(parentlist[[j]]%in%tocollapse | tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tips[j],]$colour!=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tips[j],]$parent,]$colour)){
        collapsetips[[j]] <- tips[j]  
      }
    } 
    collapsetips <- unlist(collapsetips)
    if(!is.null(collapsetips)){
      collapsenames <- list()
      for(r in 1:length(collapsetips)){
        ## get names for these tips
        collapsenames[[r]] <- tre[["tip.label"]][collapsetips[r]]
      }
      collapsenames <- unlist(collapsenames)
      #collapsenames <- collapsenames[!is.na(collapsenames)]
      droptree <-  custom.drop.tip(extract.clade(tre,nodenum), tip=collapsenames, tree=tree, trim.internal = T, collapse.singles = F, subtree = T)
      
    }else{
      droptree <- extract.clade(tre,nodenum)
    }
    # blank plot if no tree
    if(is.null(droptree)){
      subtreelist[[i]] <- ggplot()+
        theme_bw() +
        coord_cartesian(xlim=c(xmin-1.5, xmax+.2), ylim=c(0,2))+ #clip off v. important to allow the plots to overlap
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
      tipshapes <- NA
      tipnames <- NA
      alltips <- 0
      edges <- as.data.frame(tre[["edge"]])
      rootdate <- as.character(as.Date(subtree[i,]$branching))
      #as.character(as.Date(tree@data[tree@data$node==as.numeric(edges[edges[,2]==nodenum,][1]),]$branching))
      tipdateout <- NA
    }else if (nTips(droptree)>1){   ###  ADD A LESS THAN 5 TO THIS ONE
      edges <- as.data.frame(tre[["edge"]])
      ## can we get this from the branching date somehow? 
      rootdate <- decimal_date(as.Date(subtree[i,]$branching))
      current_root <- decimal_date(mrsd) - as.numeric(tree@data[tree@data$node==nodenum,]$height)
      droptree$root.edge <- current_root - rootdate
      rootdate <- as.character(date_decimal(rootdate))
      ## extract mrsd from the tip dates
      tipnames <- droptree$tip.label
      tipdates <- date_decimal(sapply( strsplit( tipnames, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )}))
      tipdate <- max(tipdates, na.rm = T)
      
      tipshapes <-     sapply( strsplit( tipnames, '\\_' ), function(x){  head(x,2)[2]})
      tipshapes[which(tipshapes %!in% c("AFP", "ES"))] <- NA
      ## extract AFP and ES and set all the other tips to NA
      
      length(tipshapes) <- length(droptree$tip.label)+droptree$Nnode
      alltips <- length(tipshapes[!is.na(tipshapes)])
      subtreelist[[i]] <- ggtree(droptree, mrsd=tipdate,colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour) +
        scale_y_reverse()+
        coord_cartesian(xlim=c(xmin-1.5, xmax+.2), clip = 'off')+ #clip off v. important to allow the plots to overlap
        ###ADD YLIMS A BIT WIDER THAN THE TREEE
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
      tipdateout <- as.character(tipdate)
      
    }else{
      tipnames <- droptree[["tip.label"]][1]
      # extract the date from this string
      tipdate <- date_decimal(sapply( strsplit( tipnames, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )}))
      tipshapes <- read.table(text = as.character(tipnames), sep = "_")$V2
      length(tipshapes) <- length(droptree$tip.label)+droptree$Nnode
      #tipdate <- as.POSIXct(max(tipdates))
      edges <- as.data.frame(tre[["edge"]])
      rootdate <- as.character(as.Date(subtree[i,]$branching))
      
      alltips <- 1
      subtreelist[[i]] <- ggplot() +
        geom_line(aes(x=c(decimal_date(as.Date(rootdate)),decimal_date(tipdate)), y=c(1,1)), colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour)+
        geom_point(aes(x=decimal_date(as.Date(rootdate)),y=1),shape=23, fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$parent,]$colour, size=2)+
        geom_point(aes(x=decimal_date(tipdate),y=1, shape=tipshapes), fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour, size=2)+
        scale_shape_manual(values=tipsh)+
        coord_cartesian(xlim= c(xmin-1.5, xmax+.2), clip = 'off')+
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
      tipdateout <- as.character(tipdate)
    }
    #nudge the point by the length of the root edge
  }else{
    # if a single tip is left need to plot it as a point on an empty ggplot axis then add a line left of it the length of the distance to the mrca and then another point in the colour of the mrca
    # therefore need the date of the node from it's name
    nodename <- tre[["tip.label"]][nodenum]
    tipnames <- nodename
    # extract the date from this string
    tipdate <- date_decimal(sapply( strsplit( tipnames, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )}))
    tipshapes <- read.table(text = as.character(nodename), sep = "_")$V2
    tipshapes <- tipshapes[!is.na(tipshapes)]
    #tipdate <- as.POSIXct(max(tipdates))
    
    edges <- as.data.frame(tre[["edge"]])
    ## can we get this from the branching date somehow? 
    rootdate <- decimal_date(as.Date(subtree[i,]$branching))
    current_root <- decimal_date(mrsd) - as.numeric(tree@data[tree@data$node==nodenum,]$height)
    root.edge <- current_root - rootdate
    rootdate <- as.character(date_decimal(rootdate))
    
    alltips <- 1
    subtreelist[[i]] <- ggplot() +
      geom_line(aes(x=c(decimal_date(as.Date(rootdate)),decimal_date(tipdate)), y=c(1,1)), colour=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour)+
      geom_point(aes(x=decimal_date(as.Date(rootdate)),y=1), shape=23, fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$parent,]$colour, size=2)+
      geom_point(aes(x=decimal_date(tipdate),y=1, shape=tipshapes), fill=tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour, size=2)+
      scale_shape_manual(values=tipsh)+
      coord_cartesian(xlim= c(xmin-1.5, xmax+.2), clip = 'off')+
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
    
    tipdateout <- as.character(tipdate)
    #need to set the x axis based on tr
    #from the lowest date in tree@data$branching to the MRSD
  }
  # if(length(tipnames)>10){
  #   locns <- data.frame(table(read.table(text = as.character(droptree$tip.label), sep = "_")$V4))
  #   main <- locns[which.max(locns$Freq),]$Var1
  #   freq <- round((max(locns$Freq)/sum(locns$Freq))*100, 2)
  #   my_grob = grid.text(paste(main," (",freq, "%)", sep=""), x=0.1,  y=0.5, gp=gpar(col="black", fontsize=12))
  #   subtreelist[[i]] <- subtreelist[[i]]+annotation_custom(my_grob) + geom_tippoint(aes(subset= (node %in% ggtree::nodeid(droptree,subset(droptree$tip.label, grepl(paste(main), droptree$tip.label)))),shape=tipshapes), fill="white", size=2)
  # }
  
  subtreelist[[i]] <- print(subtreelist[[i]])
  ancestor_loc <-  tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$parent,]$colour
  location <- tstcl[["data"]][[2]][tstcl[["data"]][[2]]$node==nodenum,]$colour
  tipshap <- as.character(table(factor(tipshapes, levels=c("AFP", "ENV"))))
  out <- c(ancestor_loc,location,rootdate,tipdateout,alltips,tipshap, toString(tipnames))
  print(i)
  
  return(list(subtreelist[[i]],out))
}

#'@title convert colours in downstreamtree to true location names
#'@return name from palette
#'@export

locatefunc <- function(x){
  if(!is.na(x)){
    return(names(getPalette[getPalette==x]))}
  else(return(NA))}


#'@title Find all identical sequences and return names
#'@param sequences aligned sequences in DNAbin format
#'@return names of identical sequences as a vector (returns both names in a pair of identical sequences)
#'@export
find_all_identical_sequence_names <- function(sequences) {
  sequence_table <- new.env(hash = TRUE)
  identical_sequence_names <- list()
  
  for (i in seq_along(sequences)) {
    seq_string <- paste(sequences[[i]], collapse = "")
    seq_hash <- digest::digest(seq_string)
    
    if (exists(seq_hash, sequence_table)) {
      # Add the name of the sequence to the existing list
      identical_sequence_names[[seq_hash]] <- c(identical_sequence_names[[seq_hash]], names(sequences)[i])
    } else {
      # Create a new entry in the list with the name of this sequence
      identical_sequence_names[[seq_hash]] <- names(sequences)[i]
      assign(seq_hash, TRUE, envir = sequence_table)
    }
  }
  
  # Filter out the entries that have only one sequence (i.e., no duplicates)
  identical_sequence_names <- identical_sequence_names[sapply(identical_sequence_names, function(x) length(x) > 1)]
  identical_sequence_names <- unlist(identical_sequence_names)
  return(unname(identical_sequence_names))
}

