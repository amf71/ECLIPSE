###Script to generate eg data to run through ctDNA pipeline ###

home <- "/Volumes/"
project <- "Archer.ctDNA" 
cdt <- paste0(home,"proj-tracerx-lung/tctProjects/frankella/")
Input <- paste0(cdt,project,"/input.data/")
Output <- paste0(cdt,project,"/output.files/")
plots <- paste0(cdt,project,"/plots/")

library(data.table)

load(file=paste0(Input,"/20200717_Archermuts_additional_annotations.RData")) #Archermuts
Archermuts <- data.table::as.data.table(Archermuts)

# choose a sample - highish ctDNA fraction
meanDAOs <- Archermuts[, mean( DAOs ), by = sample_name ]
sample <- meanDAOs[ V1 > 0.5, sample_name ][ 10 ]
sample <- 'GER1800520_LTX109_117027_LP181002027_S83_R1_001' #case with mulitmodal VAF hence CIN

example_data2 <- Archermuts[ sample_name == sample,
                             .(SampleID, sample_name, Chromosome, Position, reference, alternate, 
                               DAOs, Deep_depth, PyCloneCluster, PyCloneClonal, combPyCCF,  Ascat.Cell, 
                               panTotalCN, tnc_error_rate,  DgroupER, `failed filters` )]

names(example_data2) <- c("patient_id", "sample_id", "chr",  "pos", "ref", "alt", 
                          "supporting_reads", "depth", "clone", "is_clonal", "tumour_ccf",  
                          "tumour_cellularity", "total_cpn", "background_error", 
                          "error_group", "error_filter" )

# get an averaege of the cellularity
example_data2[, tumour_cellularity := sapply( strsplit( tumour_cellularity, split = ';' ), function(x){
  mean( sapply( strsplit( x, split = ':'), function(y) as.numeric(y[2])) )
}) ]

# also need the tumour VAFs on there.. get from pipeline output
release <- paste0( home, 'proj-tracerx-lung/tracerx/release_tx421_20190411/')
sample <- gsub( "^.{10}_", "", sample)
sample <- gsub( "_.*$", "", sample)

sample_full <- system( paste( 'ls', release ), intern = TRUE )
sample_full <- sample_full[ grepl( sample, sample_full ) ]

load( file = paste0( release, sample_full, '/exome/mutTable/', sample, '_mutTable.20200616.RData') )

mutTable <- get( ls()[grepl('LTX', ls())])
  
mut_ids <- paste( gsub( 'LTX0', 'LTX', example_data2$patient_id), 
                  example_data2$chr,
                  example_data2$pos,
                  example_data2$ref, sep = ':' )

example_data2[, tumour_vaf := mutTable[ match( mut_ids, rownames(mutTable) ), "RegionSum" ] ]
example_data2[, tumour_vaf := sapply( strsplit( tumour_vaf, split = ';' ), function(x){
  mean( sapply( strsplit( x, split = ':'), function(y) as.numeric( strsplit(y[2], split = '/')[[1]][1] ) / as.numeric( strsplit(y[2], split = '/')[[1]][2]) ))
}) ]

example_data2[, is_clonal := is_clonal == "C"]

# get the tree data for this sample
tree.release <- paste0( home, 'proj-tracerx-lung/tctProjects/lungTx/Tx421/working/trees/trees_release_20200602/' )

load( paste0( tree.release, sample, "/AutoOut.rda") )
tree <- AutoOut$trees[[1]]$auto_tree

tree <- apply( tree, 2, as.numeric )

setwd( paste0(home, "proj-tracerx-lung/tctProjects/frankella/R_packages/eclipse") )
data.table::fwrite(example_data2, 'data/example_data.csv')
data.table::fwrite(tree, 'data/tree.csv')











