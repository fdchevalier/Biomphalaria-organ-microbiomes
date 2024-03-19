#!/usr/bin/env Rscript
# Title: library_stats.R
# Version: 2.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-03-25
# Modified in: 2024-03-19



#==========#
# Comments #
#==========#

# Generate statistic table for 16S libraries processed with Qiime2.



#==========#
# Versions #
#==========#

# v2.0 - 2024-03-19: adapt the code to the organ microbiome project
# v1.0 - 2020-03-25: creation



#==========#
# Packages #
#==========#

library("magrittr")



#===========#
# Variables #
#===========#

# Set working directory
setwd(file.path(getwd(), "scripts"))

# Results and graphs directories
res.d <- "../results/1-qiime/"
graph.d <- "../graphs/"

# Sequencing stat file
myqzv.f <- paste0(res.d, "denoising-stats.qzv")



#=================#
# Data processing #
#=================#

# Check for folder presence
if (! dir.exists(graph.d)) { dir.create(graph.d, recursive=TRUE) }

# Loading data
## Identification of the metadata file within the archive file
myfile <- unzip(myqzv.f, list=TRUE)[1] %>% unlist() %>% grep("metadata.tsv", ., value=T)
## Reading the file from the archive file
mydata <- read.delim(unz(myqzv.f, myfile), stringsAsFactors=FALSE)

# Inditifying numeric columns
num <- which(mydata[1,] == "numeric")

# Clean data frame
mydata <- mydata[-1,]

# Convert values to numeric
mydata[ ,num] <- sapply(mydata[,num], function(x) as.numeric(x))

# Export table
mytb <- mydata[, c("Species", "Type", "id", "input", "filtered", "denoised", "merged", "non.chimeric") ]
colnames(mytb) <- c("Species", "Type", "Sample ID", "Number of input sequences", "Number of sequences after filtering", "Number of sequences after denoising", "Number of sequencing after merging", "Number of non chimeric sequences")
mytb <- mytb[order(mytb[,3]),]
cat("Output: ", paste0(res.d, "lib-seq-stat (supp. table 1).tsv"), "\n")
write.table(mytb, paste0(res.d, "lib-seq-stat (supp. table 1).tsv"), sep="\t", quote=FALSE, row.names=FALSE)

output <- paste0(graph.d, "seq.ct.pdf")
cat("Output: ", output, "\n")
pdf(output, width=10)
boxplot(mydata[,"input"] ~ mydata[,"Type"], ylab="Number of sequences")
invisible(dev.off())

output <- paste0(graph.d,"seq.final.pdf")
cat("Output: ", output, "\n")
pdf(output, width=10)
boxplot(mydata[,"non.chimeric"] ~ mydata[,"Type"], ylab="Number of sequences")
invisible(dev.off())


#--------------#
# Data summary #
#--------------#

# Per library
mylb.mn <- mean(mydata[,"input"])
mylb.sd <- sd(mydata[,"input"])

mycln <- c("input", "filtered", "denoised", "merged", "non.chimeric")
mylb.tb   <- apply(mydata[,mycln], 1, function(x) x/x[1])
mylb.tb.s <- apply(mylb.tb, 1, function(x) c(mean(x), sd(x), max(x), min(x)) )
rownames(mylb.tb.s) <- c("Mean", "SD", "Max", "Min")


# Per population and replicate (mean and proportion)
mytbs     <- split(mydata, mydata[,"Type"])
mytb.s.mn <- mytbs %>% lapply(., function(x) sapply(x[,mycln], mean) ) %>% as.data.frame()
mytb.s.sd <- mytbs %>% lapply(., function(x) sapply(x[,mycln], sd) )   %>% as.data.frame()
mytb.s.p  <- apply(mytb.s.mn, 2, function(x) x/x[1])

mytb.s <- sapply(1:ncol(mytb.s.mn), function(i) paste(round(mytb.s.mn[,i],1), round(mytb.s.sd[,i],1), sep=" ± "))
colnames(mytb.s) <- colnames(mytb.s.mn)
rownames(mytb.s) <- rownames(mytb.s.mn)

cat("Output: ", paste0(res.d, "lib-stat.tsv"), "\n")
write.table(mytb.s, paste0(res.d, "lib-stat.tsv"), sep="\t", quote=FALSE)
