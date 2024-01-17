# load packages -----------------------------------------------------------
pacman::p_load(
  ape
)



# Functions for XML formatting --------------------------------------------
## Format for sequence data
.seq_format <- function( d ){
  y = '<sequence> 
  <taxon idref="NAME"/>
  SEQUENCE
  </sequence>'
  seqd = lapply(  rownames(d) , function(sid){
    y1 = gsub( pattern = 'NAME', replace=sid, y )
    seq = paste( collapse='', as.character( d[sid, ] )[1,]  ) 
    y2 = gsub( pattern ='SEQUENCE', replace= seq, y1 )
    y2 
  })
  paste( seqd, collapse = '\n' )
}

# Format for taxon data
.taxon_format <- function( d ){
  y = '<taxon id="NAME">
  <date value="DATE" direction="forwards" units="years"/>
		</taxon>'
  seqd = lapply(  rownames(d) , function(sid){
    y1 = gsub( pattern = 'NAME', replace=sid, y )
    date = sapply( strsplit( sid, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )})
    y2 = gsub( pattern ='DATE', replace= date, y1 )
    y2 
  })
  paste( seqd, collapse = '\n' )
}



# Add template to working directory ---------------------------------------
directory_interest <- paste0("./AFP_20230810/")
fasta_name <- paste0(directory_interest, "longlabs_seq_20230810_AFP.fasta")
  xml_name <- paste0(directory_interest, "TREEGEN_TEMPLATE.xml")

## Load template - copy to working directory
file.copy("TREEGEN_TEMPLATE.xml", xml_name)

## load trees
trees <- read.tree(paste0(directory_interest, 'startTrees.nwk'))
ntres <- length(trees)

# format xml tags
d = read.dna(paste0(fasta_name), format = 'fasta')
taxon_date <- .taxon_format(d)
seq_data <- .seq_format(d)


## read in skeleton
x = readLines( xml_name ) 
xmlofn = gsub( xml_name, pattern='TEMPLATE', replacement='FINAL' )

## paste data
for ( k in 1:ntres ){
  xk0 = gsub( x , pattern = 'STARTTREE', replacement = write.tree( trees[[k]] )  )
  xk1 = gsub( xk0, pattern = 'SEQUENCES', replacement = seq_data )
  xk2  = gsub( xk1, pattern='TAXONLIST', replacement= taxon_date )

  if ( !grepl( pattern = '\\.xml$', xmlofn )  )
    writeLines( xk2, con =  paste0( xmlofn, '.', k, '.xml' )  )
  else 
    writeLines( xk2, con =  gsub( pattern = '\\.xml$', replacement = paste0('\\.',k,'\\.xml'), xmlofn ) )
}



