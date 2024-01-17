# Data subsets and starting trees - want these in their own folders
pacman::p_load(
treedater,
dplyr
)

# Tree generating function ------------------------------------------------
## 
#' @title Make starting trees for BEAST analysis
#' @param tree IQtree maximum likelihood phylogeny used to generate startign trees
#' @param treeoutfn name and path of output starting tree file (multiple trees will be bundled in a single nwk file) 
#' @param ncpu number of cores to use when dating phylogenies
#' @param ntres number of starting trees to output
#' @return a newick file containing ntres starting trees
#' @references Derived from sarscov2Rutils (DOI: 10.5281/zenodo.5726953)
#' @export
make_starting_trees = function(tree , treeoutfn = 'startTrees.nwk', ncpu = 1, ntres = 1){
  tr <- tree
  tr <- di2multi( tr, tol = 1e-5 ) 
  dates <- sapply( strsplit( tr$tip.label, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )})
  names(dates) <- tr$tip.label
  # resolve polytomies randomly 
  tres <- lapply( 1:ntres, function(i) { 
    tr = unroot( multi2di( tr )  )  
    tr$edge.length <- pmax( 1/29e3/5, tr$edge.length  ) #ensures that edge is >0, makes treedater fit a little nicer 
    tr
  })
  tds <- lapply( tres, function(tr){
    dater( unroot(tr), dates[tr$tip.label], s= 906, omega0 = .01, numStartConditions=0, meanRateLimits = c(0.005,0.015), ncpu = ncpu )
  })
  tds
  outtrees = lapply( tds, function(x){
    class(x) <- 'phylo'
    x
  })
  class( outtrees ) <- 'multiPhylo' 
  write.tree( outtrees 
              , file = treeoutfn 
  )
  invisible( tds )
}

# All Sequences -----------------------------------------------------------
seqs <- read.FASTA(paste0("./longlabs_seq_AFGPAK_",latestdate,".fasta"))
write.FASTA(seqs, paste0("./",latestdate,"/longlabs_seq_",latestdate,".fasta"))


# Save starting trees
if(file.exists(paste0("./all_",latestdate))){cat("folder exists")} else {dir.create(paste0("./all_",latestdate))}

tree <- read.tree(paste0(fn, ".treefile"))
make_starting_trees(tree, treeoutfn = paste0("./all_",latestdate,"/startTrees.nwk"), ncpu =8, ntres = 4)


# Stool data only ---------------------------------------------------------

AFP_seqs <- seqs[grep("_AFP_", names(seqs))]

if(file.exists(paste0("./AFP_",latestdate))){cat("folder exists")} else {dir.create(paste0("./AFP_",latestdate))}
write.FASTA(AFP_seqs, paste0("./AFP_",latestdate,"/longlabs_seq_",latestdate,"_AFP.fasta"))

# Starting trees
fn_afp = paste0("./AFP_PkAf_",latestdate,"/longlabs_seq_",latestdate,"_AFP.fasta")
system( paste0( 'iqtree -nt AUTO -redo -m HKY -s ', fn_afp ), intern=FALSE)
tree_afp <- read.tree(paste0(fn, ".treefile"))

make_starting_trees(tree_afp, treeoutfn = paste0("./AFP_",latestdate,"/startTrees.nwk"), ncpu =8, ntres = 4)
