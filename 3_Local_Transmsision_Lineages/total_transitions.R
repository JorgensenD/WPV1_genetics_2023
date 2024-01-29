# Estimating maximum possible number of transitions
pacman::p_load(
  digest,
  ape,
  lubridate,
  treeio,
  plyr,
  dplyr,
  Biostrings,
  tibble
)

# load tree
tree <- read.beast("/Users/dnj13/Documents/PKAF_mid23/all_20230810/MJ_MCC_CA_DS.trees")
# load functions
source("LTL_functions.R")


mrsd=2023.605
mrsd=as.Date(date_decimal(mrsd))

# ancestors
tree@data$ancestor <- mapply(ancestor, tree@data$node)
ancestor_loc <- tree@data[,c("node","state")]
names(ancestor_loc) <- c("ancestor","ancestor_loc")
tree@data <- as_tibble(join(tree@data,ancestor_loc, by="ancestor"))

# midpoint of branch leading to X
tree@data$branching <- as.numeric(tree@data$height) + as.numeric(tree@data$length)/2
tree@data$branching <- date_decimal(decimal_date(mrsd)-tree@data$branching)

# need tip numbers to drop from the dataset which have identical sequences
alignment <- read.fasta("your_fasta_file.fasta")

# Apply function to find all instances of identical sequences and return their names to your sequences
all_identical_sequence_names <- find_all_identical_sequence_names(alignment)

# Find node numbers of these tips
nodenos_drop <- as.numeric(which(tree@phylo$tip.label %in% all_identical_sequence_names))

treedata_2014 <- tree@data %>%
  mutate(node = as.numeric(node)) %>%
  dplyr::filter(!node %in% nodenos_drop,
                branching>=as.POSIXct("2013-01-01") & branching<as.POSIXct("2015-01-01")) 
  
treedata_2020 <- tree@data %>%
  mutate(node = as.numeric(node)) %>%
  dplyr::filter(!node %in% nodenos_drop,
                branching>=as.POSIXct("2019-01-01") & branching<as.POSIXct("2021-01-01"))
