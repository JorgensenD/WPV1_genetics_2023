## MJ analysis
pacman::p_load(
  ape
)

# XML format functions ----------------------------------------------------
## These functions take the fasta sequence data as an input and format the taxa, sequence data, discrete trait list and pairwise state data
## how it is expected by BEAST. They are paired with a standard template allowing the same analysis to be carried out on different datasets
## and subsets of the data in the future.

.seq_format <- function( d ){
  y = 
  '<sequence> 
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

.taxon_format_loc <- function( d ){
  y = 
  '<taxon id="NAME">
  <date value="DATE" direction="forwards" units="years"/>
  <attr name="state">
        LOCN
  </attr>
  </taxon>'
  seqd = lapply(  rownames(d) , function(sid){
    y1 = gsub( pattern = 'NAME', replace=sid, y )
    date = sapply( strsplit( sid, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )})
    loc = sapply( strsplit( sid, '\\_' ), function(x){ as.character( tail(x,1)[1] )})
    y2 = gsub( pattern ='DATE', replace= date, y1 )
    y3 = gsub( pattern = "LOCN", replace = loc,y2 )
    y3 
  })
  paste( seqd, collapse = '\n' )
}


.trait_format <- function( d ){
  y = 
  '<state code="LOC_EGS"/>'
  loc = as.character(unique(sapply( strsplit( rownames(d), '\\_' ), function(x){ as.character( tail(x,1)[1] )})))
  seqd = lapply(  loc , function(sid){
    y1 = gsub( pattern = 'LOC_EGS', replace=sid, y )
    y1
  })
  paste( seqd, collapse = '\n' )
}

.state_matrix <- function(d){
  y = 
  '<parameter id="PAIRLOCS" value="MATROW"/>'
  loc = as.character(unique(sapply( strsplit( rownames(d), '\\_' ), function(x){ as.character( tail(x,1)[1] )})))
  grid <- expand.grid(loc,loc)
  sid <- paste0("state.", paste(grid$Var2, grid$Var1, sep = "_"))
  
  locmat <- list()
  for(j in 1:length(sid)){
    locmat[[j]] <- rep(0, length(sid))
    locmat[[j]][j] <- 1
    locmat[[j]] <- paste(locmat[[j]], collapse = " ")
  }
  locmat <- unlist(locmat)
  
  seqd = mapply(function(sid, locmat){
  y1 = gsub( pattern = 'PAIRLOCS', replace=sid, y )
  y2 = gsub( pattern = 'MATROW', replace=locmat, y1 )
  y2
  }, sid, locmat)
  paste( seqd, collapse = '\n' )
}
.state_count <- function(d){
  y=
  '<parameter id= "state.count" value="MATRIX"/>'
  loc = as.character(unique(sapply( strsplit( rownames(d), '\\_' ), function(x){ as.character( tail(x,1)[1] )})))
  ratemat <- matrix(1,length(loc), length(loc))
  diag(ratemat) <- 0
  y1 = gsub( pattern = 'MATRIX', replace=paste(ratemat, collapse = " "), y )
  paste( y1, collapse = '\n' )
  }


# Add data to template ----------------------------------------------------
## Paste data and trees into premade beast template with the other parameters set

directory_interest <- paste0("./all_",latestdate,"/")

fasta_name <- paste0(directory_interest,"longlabs_seq_",latestdate,".fasta")
xml_name <- paste0(directory_interest, "MJ_TEMPLATE.xml")

## name of treeset from treegen beast run after logcombiner used to combine and trim down to 1000 trees
emp_treename <- "combined_newnames.nwk"

## Move template into directory of interest
file.copy("MJ_TEMPLATE.xml", xml_name)

# format xml tags
d = read.dna(fasta_name, format = 'fasta')
taxon_date <- .taxon_format_loc(d)
seq_data <- .seq_format(d)
trait_data <- .trait_format(d)
state_outputs <- .state_matrix(d)
state_count <- .state_count(d)
## read in skeleton
x = readLines( xml_name ) 
xmlofn = gsub( xml_name, pattern='TEMPLATE', replacement='FINAL' )

## paste data

  xk0 = gsub( x , pattern = 'TREEFILE_EMP', replacement = emp_treename  )
  xk1 = gsub( xk0, pattern = 'SEQUENCES', replacement = seq_data )
  xk2  = gsub( xk1, pattern='TAXONLIST', replacement= taxon_date )
  xk3  = gsub( xk2, pattern='LOCATIONS_INCLUDED', replacement= trait_data )
  xk4  = gsub( xk3, pattern='PAIRWISE_MATRICES', replacement= state_outputs )
  xk5  = gsub( xk4, pattern='COUNT_MATRIX', replacement= state_count )
  
  
if ( !grepl( pattern = '\\.xml$', xmlofn )  ){
    writeLines( xk5, con =  paste0( xmlofn, '.', '.xml' )  )
 } else { 
    writeLines( xk5, con =  gsub( pattern = '\\.xml$', replacement = paste0('\\.','\\.xml'), xmlofn ) )
 }

  