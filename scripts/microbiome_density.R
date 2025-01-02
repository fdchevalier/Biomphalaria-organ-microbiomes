#!/usr/bin/env Rscript
# Title: microbiome_density.R
# Version: 1.2
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>, Lauren Carruthers <Lauren.Carruthers@glasgow.ac.uk>
# Created in: 2024-03-15
# Modified in: 2024-01-02



#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    library("cli")
    library("magrittr")
    library("plyr")
    library("ggpubr")
    library("multcompView")
})



#===========#
# Functions #
#===========#

cat("Loading functions...\n")

# Working directory
setwd(file.path(getwd(), "scripts"))

# Compute geom_boxplot upper whisker
## source: https://stackoverflow.com/a/53795076
q <- function(x) {
    y   <- sort(x)
    iqr <- quantile(y, 0.75) - quantile(y, 0.25)
    myq <- tail(y[which(y <= quantile(y, 0.75) + 1.5*iqr)],1)
    return(myq)
}

source("microbiome_diversity_func.R")


#===========#
# Variables #
#===========#

cli_h1("Setting variables")

# Suppress warning messages
options(warn=-1)

# Folders, sample type info, and graphic options
source("microbiome_module.R")

# qPCR file
qpcr.f <- "../data/qPCR/snail_microbiome_qPCR_data.tsv"

# Volumes
elu.vol <- 50
hml.vol <- 40



#=================#
# Data processing #
#=================#

cli_h1("Processing data")

# Load data
mydata  <- read.delim(qpcr.f, header = TRUE, stringsAsFactors = TRUE)
mydata[, "Type.code"] <- factor(mydata[, "Type.code"], org.data[,1])

# Remove qPCR water control
mydata <- mydata[ ! grepl("Control", mydata[, "Species"]), ]

# Add new type code
mydata$Code <- sapply(mydata[, "Type.code"], function(x) org.data[org.data[,1] == x, 3]) %>% factor(., levels = org.data[, 3])

my16S  <- mydata[ mydata[, "Target"] == "16s", ]
mypiwi <- mydata[ mydata[, "Target"] == "piwi4", ]


#--------------------#
# Data normalization #
#--------------------#

# Liquid samples
spl_lqd <- grepl("H|TK|TY", my16S[, "Type.code"])

# Normalize qPCR values with dilution factor or piwi
## Dilution factor to reflect sample (hemolymph or water) concentration instead of elution concentration
dil.fct <- elu.vol / hml.vol
my16S[spl_lqd, "Quantity_mean"] <- my16S[spl_lqd, "Quantity_mean"] * dil.fct

## Organs and whole snails
for (i in my16S[! spl_lqd, "Sample_ID"]) {
    my16S[ my16S[, "Sample_ID"] == i, "Quantity_mean" ] <- my16S[ my16S[, "Sample_ID"] == i, "Quantity_mean" ] / (mypiwi[ mypiwi[, "Sample_ID"] == i, "Quantity_mean" ] / 2)
}

# List of hemolymph samples with high normalized quantities
cat("List of hemolymph samples with high normalized quantities (>1e4):\n")
hml_high <- my16S[ my16S[, "Quantity_mean"] > 1e4 & my16S[, "Type.code"] == "H", ]
print(hml_high)

# Remove outlier in hemolymph
cat("Removing outliers among hemolymph samples (>1e5 normalized 16S copies):\n")
outlier <- my16S[, "Quantity_mean"] > 1e5 & my16S[, "Type.code"] == "H"
my16S   <- my16S[ ! outlier,]
spl_lqd <- spl_lqd[ ! outlier]


#------------#
# Statistics #
#------------#

cli_h1("Statistics")

# 16S copy number statistics for each snail populations
cli_h2("16S copy number statistics")
my16S_pop <- split(my16S, my16S[, "Population"])
my16S_spl <- lapply(my16S_pop, function(x) split(x, x[, "Code"]))
lapply(my16S_spl, function(x) {sapply(x, function(y) {z <- y[, "Quantity_mean"]; c(Mean = mean(z), SE = sd(z)/length(z), Min = min(z), Max = max(z), N = length(z))}) }) %>% print()

# Statistical tests
cli_h2("Pairwise comparison on 16S quantities")
mytests   <- vector("list", 2)
myletters <- vector("list", 2)
myorder   <- order(match(my16S[,"Type.code"], org.data[,1]))
mygrp     <- apply(my16S[, c("Code", "Population")], 1, paste, collapse=" ") %>% factor(., levels = .[myorder] %>% unique())

for (i in 1:2) {
    # Organs or hemolymph/water
    if (i == 1) { myrows <- ! spl_lqd } else { myrows <- spl_lqd }

    # Test
    mytest <- pairwise.wilcox.test(my16S[myrows, "Quantity_mean"], mygrp[myrows], p.adjust.method = "none")

    # Extract p-values
    mypv   <- as.vector(mytest$p.value)
    names(mypv) <- sapply(colnames(mytest$p.value), function(x) paste0(x, "-", rownames(mytest$p.value))) %>% as.vector()
    mypv <- mypv[!is.na(mypv)]

    # Store data
    mytests[[i]]   <- mytest
    myletters[[i]] <- generate.label.df(mypv)
}

print(mytests)

# Correlation between hemolymph and organs
cli_h2("Correlation between hemolymph and organs")

## Split data to individuals per population and sample type
my16S$Sample.code <- paste(my16S[, "Snail.ID"], my16S[, "Tank.ID"], sep = ".") %>% as.factor()
my16S_pop <- split(my16S, my16S[, "Population"]) %>% lapply(., function(x) x[! grepl("W|T.", x[, "Type.code"]), ] %>% droplevels())
my16S_spl <- lapply(my16S_pop, function(x) split(x, x[, "Sample.code"]))

mytype_spl <- my16S_spl[[1]][[1]][, "Code"]
mycomb     <- mytype_spl %>% combn(., 2)

## Test for correlation
cor_stat <- vector("list", length(pop.data[, 1])) %>% setNames(., pop.data[, 1])

for (i in pop.data[, 1]) {

    mymat <- matrix(NA, ncol = length(mytype_spl), nrow = length(mytype_spl))
    rownames(mymat) <- mytype_spl
    colnames(mymat) <- mytype_spl

    myspl_tmp <- my16S_spl[[i]]

    mymethod <- "pearson"
    norm_test <- sapply(mytype_spl, function(type) sapply(myspl_tmp, function(x, y = type) { row <- x[, "Code"] == y ; ifelse(all(row), NA, x[ x[, "Code"] == y, "Quantity_mean"]) }) %>% shapiro.test(.)) %>% .["p.value", ] %>% unlist() %>% any(. < 0.05)
    if (norm_test) {
        cli_alert(paste0(i, ": at least one 16S distribution does not follow a normal distribution."))
        mymethod <- "kendall"
    }

    for (j in 1:ncol(mycomb)) {
        type1 <- mycomb[1, j] %>% as.character()
        type2 <- mycomb[2, j] %>% as.character()

        mycomp <- sapply(myspl_tmp, function(x, y = type1, z = type2) x[ match(c(y, z), x[, "Code"]), "Quantity_mean"]) %>% t()

        mycor <- cor.test(mycomp[, 1], mycomp[, 2], method = mymethod)

        mymat[type1, type2] <- mycor$p.value
        mymat[type2, type1] <- mycor$estimate
    }

    cor_stat[[i]] <- mymat
}

print(cor_stat)

## Alternative: correlation between hemolymph and water.
# my16S_pop <- split(my16S, my16S[, "Population"]) %>% lapply(., function(x) x[ grepl("H|T.", x[, "Type.code"]), ] %>% droplevels())
# my16S_spl <- lapply(my16S_pop, function(x) split(x, x[, "Tank.ID"]))



#=========#
# Figures #
#=========#

# Check graph folder
if (! dir.exists(graph.d)) { dir.create(graph.d, recursive=TRUE) }

# 16S density boxplot
p.ls <- vector("list", 2)

my16S_org <- my16S[! spl_lqd, ]
my16S_lqd <- my16S[  spl_lqd, ]

y_lab <- "16S rRNA gene copy number"

p.ls[[1]] <- ggplot(my16S_org, aes(x = Code, y = Quantity_mean, fill = Population)) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "", y = paste(y_lab, "\nper host cell")) +
        scale_y_continuous(trans = "log2") +
        scale_fill_manual(values = pop.data[,2]) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p.ls[[2]] <- ggplot(my16S_lqd, aes(x = Code, y = Quantity_mean, fill = Population)) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "", y = paste(y_lab, "\nper µL")) +
        scale_y_continuous(trans = "log2") +
        scale_fill_manual(values = pop.data[,2]) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

pdf(NULL)
p <- ggarrange(plotlist = p.ls, labels = LETTERS[1:length(p.ls)], ncol = length(p.ls), nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv")
invisible(dev.off())
ggsave(paste0(graph.d, "Fig. 7 - 16S density.pdf"), p, width = 8, height = 4, useDingbats = FALSE)


# P-value matrix for organs and hemolymph comparisons
p.ls <- vector("list", length(mytests))
mybreaks <- sapply(mytests, function(x) min(x$p.value, na.rm = T)) %>% min() %>% { seq(log10(.), log10(0.05), length.out = 4) } %>% { c(10^(.), 1) }
for (i in 1:length(mytests)) {
    x_clr <- rep(pop.data[1, 2], ncol(mytests[[i]]$p.value))
    x_clr[ grepl(pop.data[2, 1], colnames(mytests[[i]]$p.value)) ] <- pop.data[2, 2]
    y_clr <- rev(x_clr)
    p.ls[[i]] <- plotPvalues(mytests[[i]]$p.value, border.col = "black", pval.col = "black", font.size = 3) +
                    xlab("") + ylab("") +
                    scale_fill_gradient2(low = "red", mid = "yellow", high = "white", limits = c(min(mybreaks), 0.05), guide = "colorbar", na.value = "transparent", trans = "log10", breaks = mybreaks, label = format(mybreaks, digit = 2, scientific = TRUE)) +
                    guides(fill = guide_colourbar(barwidth = 9, title.vjust = 0.75)) +
                    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", axis.text.x = element_text(color = x_clr, angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(color = y_clr))
}

pdf(NULL)
p <- ggarrange(plotlist = p.ls, labels = LETTERS[1:length(p.ls)], ncol = length(p.ls), nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv")
invisible(dev.off())
ggsave(paste0(graph.d, "Supp. Fig. 7 - density_p-val.pdf"), p, width = 7 * length(p.ls), height = 7, useDingbats = FALSE)

# Clean tmp file
if (file.exists("Rplots.pdf")) { unlink("Rplots.pdf") }
