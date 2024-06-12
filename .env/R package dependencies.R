# System
options("repos" = c(CRAN = "https://cloud.r-project.org", BIOC = "https://bioconductor.org/packages/3.10/bioc"))
library("devtools")
library("BiocManager")

install_github("jalvesaq/colorout", ref="7ea9440", upgrade="never")

# Microbiome
install_version("microbiome", version="1.8.0",  upgrade="never")
install_version("phylogram",  version="2.1.0",  upgrade="never")      # prune function
install_version("usedist",    version="0.4.0",  upgrade="never")
