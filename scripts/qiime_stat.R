library("magrittr")

#===========#
# Variables #
#===========#

res.d <- "results/1-qiime/"

# Denoise file
myqzv.f <- paste0(res.d, "denoising-stats.qzv")

## Metadat file
#md.f <- "sample-metadata.tsv"



#=================#
# Data processing #
#=================#

# Load metadata
md <- read.delim(md.f, header=TRUE)

# Loading data
## Identification of the metadata file within the archive file
myfile <- unzip(myqzv.f, list=TRUE)[1] %>% unlist() %>% grep("metadata.tsv", ., value=T)
# Reading the file from the archive file
mydata <- read.delim(unz(myqzv.f, myfile), stringsAsFactors=FALSE)

# Inditifying numeric columns
num <- which(mydata[1,] == "numeric")

# Clean data frame
mydata <- mydata[-1,]

# Convert values to numeric
mydata[ ,num] <- sapply(mydata[,num], function(x) as.numeric(x))

## Merge metadata with read data
#mydata <- merge(mydata, md, by=1)

#for (i in 1:2) {
#    pdf(paste0("seq.ct_",i,".pdf"), width=10)
#    boxplot(mydata[mydata$Replicate == i,"input"] ~ mydata[mydata$Replicate == i,"Population"], ylab="Number of sequence")
#    dev.off()
#    
#    pdf(paste0("seq.final_",i,".pdf"), width=10)
#    boxplot(mydata[mydata$Replicate == i,"non.chimeric"] ~ mydata[mydata$Replicate == i,"Population"], ylab="Number of sequence")
#    dev.off()
#}


#--------------#
# Data summary #
#--------------#

# Per library
mylb.mn <- mean(mydata[,"input"])
mylb.sd <- sd(mydata[,"input"])

mylb.tb <- apply(mydata[,c("input", "filtered", "denoised", "merged", "non.chimeric")], 1, function(x) x/x[1])
mylb.tb.s <- apply(mylb.tb, 1, function(x) c(mean(x), sd(x), range(x)) )
rownames(mylb.tb.s) <- c("Mean", "SD", "Min", "Max")


# Per populations and replicates (mean and proportion)
mytb.s.mn <- split(mydata, paste(mydata[,"Species"], mydata[,"Tissue"])) %>% lapply(., function(x) sapply(x[,c("input", "filtered","denoised","merged","non.chimeric")], mean) ) %>% as.data.frame()
mytb.s.sd <- split(mydata, paste(mydata[,"Species"], mydata[,"Tissue"])) %>% lapply(., function(x) sapply(x[,c("input", "filtered","denoised","merged","non.chimeric")], sd) )   %>% as.data.frame()
mytb.s.p  <- apply(mytb.s.mn, 2, function(x) x/x[1])

mytb.s <- sapply(1:ncol(mytb.s.mn), function(i) paste(round(mytb.s.mn[,i],1), round(mytb.s.sd[,i],1), sep=" ± "))
colnames(mytb.s) <- colnames(mytb.s.mn)
rownames(mytb.s) <- rownames(mytb.s.mn)

write.table(mytb.s, paste0(res.d, "lib-stat.tsv"), sep="\t", quote=FALSE)
