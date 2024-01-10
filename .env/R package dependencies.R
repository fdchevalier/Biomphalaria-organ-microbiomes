# System
options("repos" = c(CRAN = "https://cloud.r-project.org", BIOC = "https://bioconductor.org/packages/3.10/bioc"))
library("devtools")
library("BiocManager")

install_github("jalvesaq/colorout", ref="7ea9440", upgrade="never")

# Microbiome
install_version("microbiome", version="1.8.0",  upgrade="never")
install_version("phyloseq",   version="1.30.0", upgrade="never")
install_version("phylogram",  version="2.1.0",  upgrade="never")      # prune function
install_version("picante",    version="1.8",    upgrade="never")      # PD function
install_version("dendextend", version="1.15.2", upgrade="never")
install_version("usedist",    version="0.4.0",  upgrade="never")
install_version("ape",        version="5.6-1",  upgrade="never")      # new version to avoid complex number in matrix
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", ref="6e09713", upgrade="never")
install_github("ropensci/taxa",   ref="49969dc", upgrade="never")
install_github("jbisanz/qiime2R", ref="574708a", upgrade="never")

# Statistics
install_version("multcompView",  version="0.1-8",  upgrade="never")

# Models
install_version("RcppEigen",  version="0.3.3.5.0",  upgrade="never")
install_version("lme4",  version="1.1-21",  upgrade="never")
install_version("lmerTest",  version="3.1-0",  upgrade="never")

## Graph
#install_version("gridBase", version="0.4-7", upgrade="never")
#install_version("plotrix", version"3.7-5", upgrade="never")
#install_version("ggpubr", version="0.2.4", upgrade="never")
#install_version("ggdendro", version="0.1-20", upgrade="never")
install_version("UpSetR", version="1.4.0", upgrade="never")
