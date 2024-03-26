#!/usr/bin/env Rscript
# Title: microbiome_density.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>, Lauren Carruthers <Lauren.Carruthers@glasgow.ac.uk>
# Created in: 2024-03-15



#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    # library("qiime2R")
    library("magrittr")
    library("ggpubr")

    # library("rcompanion")
    library("plyr")
    library("multcompView")
})



#===========#
# Functions #
#===========#

# Compute geom_boxplot upper whisker
## source: https://stackoverflow.com/a/53795076
q <- function(x) {
    y   <- sort(x)
    iqr <- quantile(y, 0.75) - quantile(y, 0.25)
    myq <- tail(y[which(y <= quantile(y, 0.75) + 1.5*iqr)],1)
    return(myq)
}



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Working directory
setwd(file.path(getwd(), "scripts"))

# # Qiime file
# asv.f  <- "../results/1-qiime/table.qza"

# qPCR file
qpcr.f <- "../data/qPCR/snail_microbiome_qPCR_data.tsv"

# # Picrust file
# picrust.f <- "../results/2-picrust2/marker_predicted_and_nsti.tsv.gz"

# Population order and color
pop.data <- matrix(c("Ba",     "#ffac40",
                     "BgBS90", "#87d687"), byrow=TRUE, ncol=2)

# Tissue order
org.data <- matrix(c(
                "H",  "Hemolymph",
                "S",  "Stomach",
                "G",  "Gut",
                "L",  "Hepatopancreas",
                "O",  "Ovotestis",
                "W",  "Whole snail",
                "TY", "Water tray",
                "TK", "Water tank"), ncol=2, byrow=TRUE)

# Graph folder
graph.d <- "../graphs/"

# GGplot options
theme_set(theme_classic())

# Volumes
elu.vol <- 50
hml.vol <- 40



#=================#
# Data processing #
#=================#

# Load data
mydata  <- read.delim(qpcr.f, header = TRUE, stringsAsFactors = TRUE)
mydata[, "Type.code"] <- factor(mydata[, "Type.code"], org.data[,1])

# Remove qPCR water control
mydata <- mydata[ ! grepl("Control", mydata[, "Species"]), ]

my16S  <- mydata[ mydata[, "Target"] == "16s", ]
mypiwi <- mydata[ mydata[, "Target"] == "piwi4", ]


#--------------------#
# Data normalization #
#--------------------#

# Liquid samples
spl_lqd <- grepl("H|TK|TY", my16S[, "Type.code"])

## Normalize qPCR values with dilution factor
## Dilution factor to reflect sample (hemolymph or water) concentration instead of elution concentration
dil.fct <- elu.vol/hml.vol
my16S[spl_lqd, "Quantity_mean"] <- my16S[spl_lqd, "Quantity_mean"] * dil.fct

# Organs and whole snail
for (i in my16S[! spl_lqd, "Sample_ID"]) {
    my16S[ my16S[, "Sample_ID"] == i, "Quantity_mean" ] <- my16S[ my16S[, "Sample_ID"] == i, "Quantity_mean" ] / (mypiwi[ mypiwi[, "Sample_ID"] == i, "Quantity_mean" ] / 2)
}


# Statistical tests
mytests   <- vector("list", 2)
myletters <- vector("list", 2)
myorder   <- order(match(my16S[,"Type.code"], org.data[,1]))
mygrp     <- apply(my16S[, c("Type.code", "Population")], 1, paste, collapse=" ") %>% factor(., levels = .[myorder] %>% unique())

for (i in 1:2) {
    if (i == 1) { myrows <- ! spl_lqd } else { myrows <- spl_lqd }
    mytest <- pairwise.wilcox.test(my16S[myrows,"Quantity_mean"], mygrp[myrows], p.adjust.method = "none")
    mypv <- as.vector(mytest$p.value)
    names(mypv) <- sapply(colnames(mytest$p.value), function(x) paste0(x, "-", rownames(mytest$p.value))) %>% as.vector()
    mypv <- mypv[!is.na(mypv)]

    mytests[[i]]   <- mytest
    myletters[[i]] <- generate.label.df(mypv)
}

# P-value matrix for each relevant index
p <- vector("list", length(mytests))
mybreaks <- c(10^(seq(log10(mytests %>% unlist() %>% min(., na.rm = TRUE) %>% as.numeric()), log10(0.05), length.out=4)), 1)
for (i in 1:length(mytests)) {
    p[[i]] <- plotPvalues(mytests[[i]]$p.value, border.col="black", pval.col="black", font.size=3) +
        xlab("") + ylab("") +
        scale_fill_gradient2(low = "red", mid = "yellow", high = "white", limits = c(NA, 0.05), guide = "colorbar", na.value = "transparent", trans = "log10", breaks = mybreaks, label = format(mybreaks, digit = 2, scientific = TRUE)) +
        guides(fill = guide_colourbar(barwidth = 9, title.vjust = 0.75)) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        # ggtitle(myidx[i, 2])
}

my16S$Sample.code <- paste(my16S[, "Snail.ID"], my16S[, "Tank.ID"], sep = ".") %>% as.factor()

my16S_pop <- split(my16S, my16S[, "Population"]) %>% lapply(., function(x) x[! grepl("W|T.", x[, "Type.code"]), ] %>% droplevels())

my16S_spl <- lapply(my16S_pop, function(x) split(x, x[, "Sample.code"]))


mytype_spl <- my16S_spl[[1]][[1]][, "Type.code"]
mycomb    <- mytype_spl %>% combn(., 2)


cor_stat <- vector("list", length(pop.data[, 1])) %>% setNames(., pop.data[, 1])

for (i in pop.data[, 1]) {

    mymat <- matrix(NA, ncol = length(mytype_spl), nrow = length(mytype_spl))
    rownames(mymat) <- mytype_spl
    colnames(mymat) <- mytype_spl
    
    myspl_tmp <- my16S_spl[[i]]

    for (j in 1:ncol(mycomb)) {
        type1 <- mycomb[1, j] %>% as.character()
        type2 <- mycomb[2, j] %>% as.character()


        mycomp <- sapply(myspl_tmp, function(x, y = type1, z = type2) cbind(x[ x[, "Type.code"] == y, "Quantity_mean"], x[ x[, "Type.code"] == z, "Quantity_mean"])) %>% t()

        mycor <- cor.test(mycomp[, 1], mycomp[, 2], method = "kendall")

        mymat[type1, type2] <- mycor$p.value
        mymat[type2, type1] <- mycor$estimate
    }

    cor_stat[[i]] <- mymat
}


my16S_pop <- split(my16S, my16S[, "Population"]) %>% lapply(., function(x) x[ grepl("H|T.", x[, "Type.code"]), ] %>% droplevels())

my16S_spl <- lapply(my16S_pop, function(x) split(x, x[, "Tank.ID"]))


mytype_spl <- my16S_spl[[1]][[1]][, "Type.code"] %>% levels()
mycomb    <- mytype_spl %>% combn(., 2)


cor_stat <- vector("list", length(pop.data[, 1])) %>% setNames(., pop.data[, 1])

for (i in pop.data[, 1]) {

    mymat <- matrix(NA, ncol = length(mytype_spl), nrow = length(mytype_spl))
    rownames(mymat) <- mytype_spl
    colnames(mymat) <- mytype_spl

    myspl_tmp <- my16S_spl[[i]]

    for (j in 1:ncol(mycomb)) {
        type1 <- mycomb[1, j] %>% as.character()
        type2 <- mycomb[2, j] %>% as.character()

print(j)
        mycomp <- sapply(myspl_tmp, function(x, y = type1, z = type2) cbind(x[ x[, "Type.code"] == y, "Quantity_mean"], x[ x[, "Type.code"] == z, "Quantity_mean"])) %>% t()

        mycor <- cor.test(mycomp[, 1], mycomp[, 2], method = "kendall")

        mymat[type1, type2] <- mycor$p.value
        mymat[type2, type1] <- mycor$estimate
    }

    cor_stat[[i]] <- mymat
}





# Raw 16S copy number between populations
cat("Raw 16S copy number:\n\n")
my16S_pop <- split(my16S, my16S[, "Population"]) 
my16S_spl <- lapply(my16S_pop, function(x) split(x, x[, "Type.code"]))
lapply(my16S_spl, function(x) {sapply(x, function(y) {z <- y[, "Quantity_mean"]; z <- z[z < 1e5];  c(Mean = mean(z), SE = sd(z)/length(z), Min = min(z), Max = max(z), N = length(z))}) }) %>% print()

lapply(my16S_spl, function(x) { y <- do.call(rbind, x[c("G", "S")]) ; z <- y[, "Quantity_mean"]; z <- z[z < 1e5];  c(Mean = mean(z), SE = sd(z)/length(z), Min = min(z), Max = max(z), N = length(z)) }) %>% print()



# # Note: no data cleaning because contaminant are also quantified with qPCR
# asv <- apply(taxa, 2, function (x) sum(x > 0) ) %>% data.frame(asv=.)

# # Update mrkr.pred with contaminants
# contam    <- ! rownames(taxa) %in% mrkr.pred[,1]
# contam.mt <- cbind(rownames(taxa)[contam], 1, 1)
# colnames(contam.mt) <- colnames(mrkr.pred)
# mrkr.pred <- rbind(mrkr.pred, contam.mt)

# # Reset 16S count prediction for ASV with nsis > 2 because of incertainty of annotation (see PiCRUST documentation)
# mrkr.pred[mrkr.pred[, 3] > 2, 2] <- 1

# # Adjust marker count
# mrkr.count <- apply(mrkr.pred[,1:2], 1, function(x) taxa[ x[1], ] * as.numeric(x[2]))
# mrkr.sum   <- as.data.frame(rowSums(mrkr.count))

# # Bacterial index
# mrkr.idx <- colSums(taxa) / mrkr.sum

# # Merge data
# mydata <- merge(mydata, mrkr.idx, by.x=1, by.y="row.names")


# Correct and normlize data
## Dilution factor to reflect sample (hemolymph or water) concentration (instead of elution concentration)
dil.fct <- elu.vol/hml.vol

## Correct qPCR values with dilution factor
mydata[,2] <- mydata[,2] * dil.fct

## Normalize qPCR count
mydata[,ncol(mydata)+1] <- mydata[,2] * mydata[,3]

## Merge ASV number
mydata <- merge(mydata, asv, by.x=1, by.y="row.names")


# Population labels
mydata[,ncol(mydata)+1] <- lapply(strsplit(mydata[,1], ".", fixed=TRUE), function(x) x[1]) %>% unlist()
mydata[ grepl("Water", mydata[,1]), ncol(mydata)] <- "Water"
colnames(mydata)[ncol(mydata)] <- "Population"
mydata[,"Population"] <- factor(mydata[,"Population"], pop.data[,1])
mydata <- mydata[order(match(mydata$Population,pop.data[,1])),] # Reorder

# Raw 16S copy number between populations
cat("Raw 16S copy number\n\n")
mycln <- 2

## Mean and range
cat("16S copy number:\n\t- mean: ", mean(my16[,mycln]), "\n\t- range (min max): ", range(mydata[,mycln]), "\n\n")

## Test for normal distribution without replicates
not.norm <- any(sapply(pop.data[,1], function(x) shapiro.test(mydata[mydata[,"Population"] == x, mycln])$p.value) < 0.05)
if (not.norm) cat("Data per population is not following normal distribution\n")

## Formating data
qpcr.d <- mydata[, c(7,mycln,6)]

## Kruskal-Wallis test
qpcr.d.k <- kruskal.test(qpcr.d[,2], qpcr.d[,3])
print(qpcr.d.k)

## Statistical comparison
qpcr.d.stat    <- pairwise.wilcox.test(qpcr.d[,2], qpcr.d[,3], p.adjust="none")
qpcr.d.tb      <- fullPTable(qpcr.d.stat$p.value)
qpcr.d.letters <- multcompLetters(qpcr.d.tb, compare="<=", threshold=0.05, Letters=letters)$Letters

## Print results
cat("Pairwise comparisons of 16S copy number between population (mean of the replicates): ", qpcr.d.letters, "\n\n\n")


# Estimated bacteria between populations
cat("Estimateed bacteria\n\n")
mycln <- 4

## Mean and range
cat("Bacteria number:\n\t- mean: ", mean(mydata[,mycln]), "\n\t- range (min max): ", range(mydata[,mycln]), "\n\n")

## Test for normal distribution without replicates
not.norm <- any(sapply(pop.data[,1], function(x) shapiro.test(mydata[mydata[,"Population"] == x, mycln])$p.value) < 0.05)
if (not.norm) cat("Data per population is not following normal distribution\n")

## Formating data
asv.d <- mydata[, c(7,mycln,6)]

## Kruskal-Wallis test
asv.d.k <- kruskal.test(asv.d[,2], asv.d[,3])
print(asv.d.k)

## Statistical comparison
asv.d.stat    <- pairwise.wilcox.test(asv.d[,2], asv.d[,3])
asv.d.tb      <- fullPTable(asv.d.stat$p.value)
asv.d.letters <- multcompLetters(asv.d.tb, compare="<=", threshold=0.05, Letters=letters)$Letters

## Print results
cat("Pairwise comparisons of estimated bacteria between population (mean of the replicates): ", asv.d.letters, "\n\n")

## Mean bacteria in snails
snail.b <- asv.d[asv.d[,"Population"] != "Water",2]
cat("Mean bacteria number in snails:\n\t- mean: ", mean(snail.b), "\n\t- standard error: ", sd(snail.b)/sqrt(length(snail.b)), "\n\t- range (min max): ", range(snail.b), "\n\n")

## Mean bacteria in water
water.b <- asv.d[asv.d[,"Population"] == "Water",2]
cat("Mean bacteria number in waters:\n\t- mean: ", mean(water.b), "\n\t- standard error: ", sd(water.b)/sqrt(length(water.b)), "\n\t- range (min max): ", range(water.b), "\n\n")



#=========#
# Figures #
#=========#

# Check graph folder
if (! dir.exists(graph.d)) { dir.create(graph.d, recursive=TRUE) }

p.ls <- vector("list", 2)

my16S_org <- my16S[! spl_lqd, ]
my16S_lqd <- my16S[  spl_lqd, ]

p.ls[[1]] <- ggplot(my16S_org, aes(x = Type.code, y = Quantity_mean, fill = Population)) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "", y = bquote("16S copy number per host cell (log"[2] * ")")) +
        scale_y_continuous(limits = quantile(my16S_org[, "Quantity_mean"], c(0.1, 0.9)), trans = "log2") +
        scale_fill_manual(values = pop.data[,2]) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p.ls[[2]] <- ggplot(my16S_lqd, aes(x = Type.code, y = Quantity_mean, fill = Population)) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "", y = bquote("16S copy number per µL (log"[2] * ")")) +
        scale_y_continuous(limits = quantile(my16S_lqd[, "Quantity_mean"], c(0.1, 0.9)), trans = "log2") +
        scale_fill_manual(values = pop.data[,2]) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

pdf(NULL)
p <- ggarrange(plotlist = p.ls, labels = LETTERS[1:length(p.ls)], ncol = length(p.ls), nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv")
invisible(dev.off())
ggsave(paste0(graph.d, "Fig. 7 - 16S density.pdf"), p, width = 8, height = 4, useDingbats = FALSE)








#--------------------#
# Microbiome density #
#--------------------#

p.ls <- vector("list", 3)

ylim1 <- boxplot(mydata$V2 ~ mydata$Population, outline=FALSE, plot=FALSE)$stats[1:4,] %>% max() %>% signif(., 1)
ylim2 <- boxplot(mydata$V4 ~ mydata$Population, outline=FALSE, plot=FALSE)$stats[1:4,] %>% max() %>% signif(., 1)
ylim <- max(c(ylim1, ylim2))
# Raw data
p.ls[[1]] <- ggplot(mydata, aes(x = Population, y = V2, fill = Population)) +
        geom_boxplot(outlier.shape=NA) + xlab("") + ylab("16S copy number per µL") +
        scale_y_continuous(limits = c(0, ylim)) +
        scale_fill_manual(values = pop.data[,2]) +
        stat_summary(geom = 'text', label = qpcr.d.letters, fun.y = q, vjust = -1) +
        theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Normalized data
p.ls[[2]] <- ggplot(mydata, aes(x = Population, y = V4, fill = Population)) +
        geom_boxplot(outlier.shape=NA) + xlab("") + ylab("Estimated number of bacteria per µL") + 
        scale_y_continuous(limits = c(0, ylim)) +
        scale_fill_manual(values = pop.data[,2]) +
        stat_summary(geom = 'text', label = asv.d.letters, fun.y = q, vjust = -1) +
        theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Correlation
p.ls[[3]] <- ggplot(mydata, aes(x = V4, y = asv)) +
        geom_point(aes(color = Population)) + 
        stat_smooth(method = lm, se = FALSE) +
        xlab("Bacterial number") + ylab("Number of observed ASVs") +
        scale_color_manual(values = pop.data[,2]) +
        theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
        stat_cor(method = "kendall", cor.coef.name = "tau", label.x.npc = "right", label.y.npc = "top", hjust = 1)

p <- ggarrange(plotlist=p.ls, ncol=1, labels=LETTERS[1:length(p.ls)])
ggsave(paste0(graph.d,"Fig. 5 - qPCR-diversity.pdf"), p, height=15, width=5)


# Clean tmp file
if (file.exists("Rplots.pdf")) { unlink("Rplots.pdf") }

