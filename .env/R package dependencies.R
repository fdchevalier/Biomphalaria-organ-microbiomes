# System
options("repos" = c(CRAN = "https://cloud.r-project.org"))
library("devtools")
library("BiocManager")

install_github("jalvesaq/colorout", ref="7ea9440", upgrade="never")


install_bioc("3.10/microbiome", upgrade="never")
install_bioc("3.10/phyloseq",   upgrade="never")
install_version("phylogram",  version="2.1.0",  upgrade="never")      # prune function
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", ref="6e09713", upgrade="never")
#install_github("ropensci/taxa",   ref="49969dc", upgrade="never")
install_github("jbisanz/qiime2R", ref="574708a", upgrade="never")


# PCA analysis
install_github("willtownes/glmpca", ref="edc04cc", upgrade="never")

##install_version("doParallel", version="1.0.15", upgrade="never")

##library(rcompanion)
##library(multcompView)


#install("phyloseq", update=FALSE)
#install("picante", update=FALSE)          # PD function
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", ref="6e09713", upgrade="never")
#install_github("ropensci/taxa", ref="49969dc", upgrade="never")
##library("metacoder")
#install_version("phylogram", version="2.1.0", upgrade="never")      # prune function

## Graph
##library("RColorBrewer")
##library("colorspace")
#install_version("gridBase", version="0.4-7", upgrade="never")
#install_version("plotrix", version"3.7-5", upgrade="never")
#install_version("ggpubr", version="0.2.4", upgrade="never")
#install_version("ggdendro", version="0.1-20", upgrade="never")


