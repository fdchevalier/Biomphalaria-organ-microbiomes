#!/usr/bin/env Rscript
# Title: microbiome_diversity.R
# Version: 2.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-03-25
# Modified in: 2025-01-02



#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    # System
    library("cli")
    library("magrittr")
    library("plyr")             # function count
    library("tidyr")
    library("reshape")

    # Microbiome
    library("microbiome")       # function core
    library("qiime2R")
    library("phyloseq")
    library("picante")          # PD function
    library("pairwiseAdonis")   # see github
    library("phylogram")        # prune function
    library("dendextend")
    library("usedist")
    library("lme4")
    library("lmerTest")

    # Graph
    library("colorspace")       # function rainbow_hcl
    library("ggpubr")
    library("ggdendro")
    library("grid")             # function grid.rect
    library("lemon")
    library("UpSetR")
})



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

cat("Loading functions...\n")

# Compute shared ASV proportion
prop_asv <- function(list_df, org.vec, type, private = FALSE) {
    res <- sapply(org.vec[! org.vec %in% type], function(y) sapply(list_df, function(x) {
                                                                    if (private) {
                                                                        z <- rowSums(x[, type, drop = FALSE]  == 1 & x[, y]  == 1) > 0 & rowSums(x[, setdiff(names(x), c(type, y))]) == 0
                                                                    } else {
                                                                        z <- rowSums(x[, type, drop = FALSE]  == 1 & x[, y]  == 1) > 0 }
                                                                    sum(z)/length(z[rowSums(x) > 1]) } ))
    return(res)
}

# External functions
source("microbiome_diversity_func.R")



#===========#
# Variables #
#===========#

cli_h1("Setting variables")

# # Suppress warning messages
# options(warn=-1)

# Qiime files
asv.f  <- "../results/1-qiime/table.qza"
taxa.f <- "../results/1-qiime/rep-seqs_taxa.qza"
tree.f <- "../results/1-qiime/rooted-tree.qza"

# Metadata file
md.f   <- "../data/sample-metadata.tsv"

# Contaminant file
blast.f <- "../results/1-qiime/rep-seqs_unassigned/unassigned.blastn.tsv"

# Folders, sample type info, and graphic options
source("microbiome_module.R")

# Update result folder
result.d <- paste0(result.d, "/2-diversity/")

# Seed for reproducibility
myseed <- 26017



#=================#
# Data processing #
#=================#

cli_h1("Processing data")

if (! dir.exists(result.d)) { dir.create(result.d, recursive = TRUE) }
if (! dir.exists(graph.d)) { dir.create(graph.d, recursive = TRUE) }

#----------------#
# Preparing data #
#----------------#

# Loading ASV data
asv    <- read_qza(asv.f)$data
asv.tb <- otu_table(asv, taxa_are_rows = TRUE)

# Loading contaminants
blast <- read.delim(blast.f, header = FALSE)

# Loading taxa data and reorganizing table
taxa <- read_qza(taxa.f)$data

## Split info
cln.nm <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa   <- separate(taxa, 2, cln.nm, sep = ";")

# Polishing data
rownames(taxa) <- taxa[,1]
taxa <- taxa[,-c(1,ncol(taxa))]

# Assign "Unassigned" to all level of unassigned ASV
taxa[ taxa[,1] == "Unassigned", ] <- "Unassigned"

## Convert table
taxa.tb             <- tax_table(as.matrix(taxa))
taxa_names(taxa.tb) <- rownames(taxa)


# Loading tree
tree <- read_qza(tree.f)$data

# Loading metadata and reorder factors
md <- read.delim(md.f, row.names = 1)
md[,"Population"] <- factor(md[,"Population"], pop.data[,1])
md[,"Type.code"]  <- factor(md[,"Type.code"], org.data[,1])

# Add new type code
md$Code <- sapply(md[, "Type.code"], function(x) org.data[org.data[,1] == x, 3]) %>% factor(., levels = org.data[, 3])

# Building phyloseq object
mybiom <- phyloseq(asv.tb, taxa.tb, tree, sample_data(md))

# Populations
mypop <- md

# Cleaning data from contaminants
nb.asv   <- nrow(tax_table(mybiom))
nb.unass <- (taxa[,1] == "Unassigned") %>% sum()

## Mitochondria and chloroplast
my.mt <- apply(tax_table(mybiom), 1, function(x) any(grepl("mitochondria", x, ignore.case = TRUE)))
my.cp <- apply(tax_table(mybiom), 1, function(x) any(grepl("chloroplast", x, ignore.case = TRUE)))

## Eukaryote contaminants among the unassigned
my.euk <- blast[ blast[,14] == "Eukaryota", 1]
my.euk <- rownames(tax_table(mybiom)) %in% my.euk

## Exporting contaminent list
write.table(rownames(tax_table(mybiom))[(my.mt | my.cp | my.euk)], "../data/contaminants", row.names = FALSE, col.names = FALSE, quote = FALSE)

## Remove contaminants
mybiom <- subset_taxa(mybiom, ! (my.mt | my.cp | my.euk))

## Final numbers
cat("Initial number of ASVs:", nrow(asv.tb), "\n\n")
cat("\nNumber of ASVs removed:",
    "\n\t- mitochondria and chloroplasts:", sum(my.mt | my.cp), "/", signif(sum(my.mt | my.cp) / nb.asv * 100, 2), "% of the total ASVs",
    "\n\t- eukaryotes:", sum(my.euk), "/", signif(sum(my.euk) / nb.unass * 100, 2), "% of the unassigned ASVs", "/", signif(sum(my.euk) / nb.asv * 100, 2), "% of the total ASVs",
    "\n\t- total:", sum(my.mt | my.cp | my.euk), "/", signif(sum(my.mt | my.cp | my.euk)/ nb.asv * 100, 2), "% of the total ASVs", "\n\n")

# Number per snail species
otu_nb_sp <- vector("list", nrow(pop.data)) %>% setNames(., pop.data[,1])

for (i in pop.data[, 1]) {
    mysp  <- paste0(i, ".")
    mycln <- otu_table(mybiom) %>% colnames() %>% grep(org.tp[,1] %>% paste0(mysp, ., collapse = "|"), ., value = TRUE)
    otu_nm <- ((otu_table(mybiom)[, mycln] %>% rowSums()) > 0) %>% which() %>% names()
    otu_unass <- otu_nm[otu_nm %in% rownames(tax_table(mybiom)[ tax_table(mybiom)[,1] == "Unassigned", ])]
    otu_nb_sp[[i]] <- list(otu_nm, otu_unass)
}

cat("\nNumber of ASVs per snail species:\n")
for (i in pop.data[, 1]) {
    cat("\t- ", i, ": ", length(otu_nb_sp[[i]][[1]]), " (", length(otu_nb_sp[[i]][[2]]), " unassigned)\n", sep = "")
}
cat("\t- Common ASVs: ", intersect(otu_nb_sp[[1]][[1]], otu_nb_sp[[2]][[1]]) %>% length(), " (",  intersect(otu_nb_sp[[1]][[2]], otu_nb_sp[[2]][[2]]) %>% length(), " unassigned)\n", sep = "")


#--------------------#
# Rarefaction curves #
#--------------------#

cli_h1("Statistics")
cli_h2("Rarefaction curves")

p <- ggrare(mybiom, step = 1000, color = "Population", se = FALSE, plot = FALSE) +
        scale_color_manual(values = pop.data[,2]) +
        facet_wrap(~Code)

# Adjust legend
pdf(NULL)
p <- shift_legend(p)
invisible(dev.off())

ggsave(paste0(graph.d, "Supp. Fig. 1 - rarefaction.pdf"), p, height = 7, width = 7, useDingbats = FALSE)


#-----------------#
# Alpha-diversity #
#-----------------#

cli_h2("Alpha-diversity")

myalpha <- microbiome::alpha(mybiom)
mypd    <- pd(t(asv), tree)
myalpha <- merge(myalpha, mypd, by = "row.names")
myalpha <- merge(myalpha, md, by.x = 1, by.y = "row.names")

myidx <- matrix(c("observed",           "Observed richness",
                  "PD",                 "Faith's Phylogenetic diversity",
                  "evenness_simpson",   "Simpson evenness"
                  ), ncol = 2, byrow = TRUE)

# Statistical tests
mya.test <- vector("list", length(myidx[,1]))
myorder  <- order(match(myalpha[,"Type.code"], org.data[,1]))
mygrp    <- apply(myalpha[, c("Type", "Population")], 1, paste, collapse=" ") %>% factor(., levels = .[myorder] %>% unique())
for (i in myidx[,1]) {
    mya.test[[ match(i,myidx[,1]) ]] <- pairwise.wilcox.test(myalpha[,i], mygrp, p.adjust.method = "fdr")
}

mya.letters <- vector("list", length(myidx[,1]))
for (i in 1:length(mya.test)) {
    d <- as.vector(mya.test[[i]]$p.value)
    names(d) <- sapply(colnames(mya.test[[i]]$p.value), function(x) paste0(x, "-", rownames(mya.test[[i]]$p.value))) %>% as.vector()
    d <- d[!is.na(d)]
    mya.letters[[i]] <- generate.label.df(d)
}

# Linear mixed-effect models
myalpha2 <- myalpha[ myalpha[, "Species"] != "Water", ]
myalpha2[, "Snail.ID"] <- paste0(myalpha2[, "Snail.ID"], myalpha2[, "Tank.ID"])
myalpha2[, "Tank.ID"] %<>% as.factor()

## If test without whole snail needed
#myalpha2 <- myalpha2[ myalpha2[, "Type"] != "Whole snail ", ]

## Indices to be tested
indices <- myidx[,1]

## Species and sample types
mylms  <- vector("list", length(indices)) %>% setNames(., indices)
for (i in indices) {
    # mylm   <- lmer(paste(a, "~ Species * Type + (1|Snail.ID) + (1|Tank.ID)"), data = myalpha2) ## Note: This model gives the same results has the one below.
    mylm   <- lmer(paste(i, "~ Species * Type + (1|Snail.ID)"), data = myalpha2)
    # print(BIC(mylm))
    myanova <- anova(mylm)
    mycoef <- mylm %>% summary() %>% coef()

    # Output ANOVA results
    append  <- ifelse(match(i, indices) == 1, FALSE, TRUE)
    outfile <- paste0(result.d, "a-div_anova_results.txt")
    write(i, outfile, append = append)
    write.table(myanova, outfile, append = TRUE, sep = "\t")

    # Output LMER results
    outfile <- paste0(result.d, "a-div_lmer_results.txt")
    write(i, outfile, append = append)
    write.table(mycoef, outfile, append = TRUE, sep = "\t")

    # Store results
    mylms[[i]] <- list(lm = mylm, anova = myanova, coef = mycoef)
}

## Conditions (period and date)
myconds <- c("Collection.time", "Collection.date")
mylms2 <- vector("list", length(pop.data[,1])) %>% setNames(., pop.data[,1])
for (p in pop.data[,1]) {
    myrows  <- myalpha2[, "Population"] == p
    mylms2[[p]] <- vector("list", length(indices)) %>% setNames(., indices)
    for (i in indices) {
        mylms2[[p]][[i]] <- vector("list", length(myconds)) %>% setNames(., myconds)
        for (j in myconds) {
            mylm   <- lmer(paste(i, "~ Type *", j, "+ (1|Snail.ID)"), data = myalpha2[myrows, ])
            # print(BIC(mylm))
            myanova <- anova(mylm)
            mycoef <- mylm %>% summary() %>% coef()

            # Output ANOVA results
            append <- ifelse(match(i, indices) == 1 & match(j, myconds) == 1, FALSE, TRUE)
            outfile <- paste0(result.d, "a-div_cond_anova_results.txt")
            write(paste(i, j), outfile, append = append)
            write.table(myanova, outfile, append = TRUE, sep = "\t")

            # Output LMER results
            outfile <- paste0(result.d, "a-div_cond_lmer_results.txt")
            write(paste(i, j), outfile, append = append)
            write.table(mycoef, outfile, append = TRUE, sep = "\t")

            # Store results
            mylms2[[p]][[i]][[j]] <- list(lm = mylm, anova = myanova, coef = mycoef)
        }
    }
}

cli_alert(paste0("Linear mixed-model effect results output in ", result.d, "."))


#~~~~~~~~~#
# Figures #
#~~~~~~~~~#

# Check model residuals
pdf(paste0(graph.d, "lm residuals - sp-spl.pdf"))
for (i in indices) {
    myresid <- resid(mylms[[i]]$lm) %>% data.frame(x = .)
    p <- ggqqplot(myresid, "x") + ggtitle(i)
    print(p)
}
invisible(dev.off())

pdf(paste0(graph.d, "lm residuals - conditions.pdf"))
for (p in pop.data[,1]) {
    for(i in indices) {
        for(j in myconds) {
            myresid <- resid(mylms2[[p]][[i]][[j]]$lm) %>% data.frame(x = .)
            ggp <- ggqqplot(myresid, "x") + ggtitle(paste(p, "-", i, "-", j))
            print(ggp)
        }
    }
}
invisible(dev.off())


# Boxplot for each relevant index
## No stat added on graph because unreadable
p.ls <- vector("list", nrow(myidx))
for (i in myidx[,1]) {
    df.tmp <- data.frame(i = myalpha[,i], Tissue = myalpha[,"Code"], Species = myalpha[,"Population"])
    p.ls[[match(i,myidx)]] <- ggplot(df.tmp, aes(x = Tissue, y = i, fill = Species)) +
                                geom_boxplot() +
                                labs(x = "", y = "", title = myidx[myidx[,1] == i,2]) +
                                scale_fill_manual(values = pop.data[,2]) +
                                theme(plot.title = element_text(hjust =0.5), legend.position = "none")
}
pdf(NULL)
p <- ggarrange(plotlist = p.ls, labels = LETTERS[1:length(p.ls)], ncol = length(p.ls), nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv")
invisible(dev.off())
ggsave(paste0(graph.d, "Fig. 2 - alpha-div.pdf"), p, width = 10, height = 3, useDingbats = FALSE)


# P-value matrix for each relevant index
p.ls <- vector("list", nrow(myidx))
mybreaks <- sapply(mya.test, function(x) min(x$p.value, na.rm = T)) %>% min() %>% { seq(log10(.), log10(0.05), length.out = 4) } %>% { c(10^(.), 1) }
for (i in 1:nrow(myidx)) {
    # Colors
    x_clr <- rep(pop.data[1, 2], ncol(mya.test[[i]]$p.value))
    x_clr[ grepl(pop.data[2, 1], colnames(mya.test[[i]]$p.value)) ] <- pop.data[2, 2]
    y_clr <- rep(pop.data[1, 2], nrow(mya.test[[i]]$p.value))
    y_clr[ grepl(pop.data[2, 1], rownames(mya.test[[i]]$p.value)) ] <- pop.data[2, 2]

    p.ls[[i]] <- plotPvalues(mya.test[[i]]$p.value, border.col = "black", pval.col = "black", font.size = 3) +
                    labs(x = "", y = "", title = myidx[i, 2]) +
                    scale_fill_gradient2(low = "red", mid = "yellow", high = "white", limits = c(NA, 0.05), guide = "colorbar", na.value = "transparent", trans = "log10", breaks = mybreaks, label = format(mybreaks, digit = 2, scientific = TRUE)) +
                    guides(fill = guide_colourbar(barwidth = 9, title.vjust = 0.75)) +
                    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", axis.text.x = element_text(color = x_clr, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(color = y_clr))
}

pdf(NULL)
p <- ggarrange(plotlist = p.ls, labels = LETTERS[1:length(p.ls)], ncol = length(p.ls), nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv")
invisible(dev.off())
ggsave(paste0(graph.d, "Supp. Fig. 2 - alpha-div_p-val.pdf"), p, width = 10 * length(p.ls), height = 10, useDingbats = FALSE)


#----------------#
# Beta-diversity #
#----------------#

cli_h2("Beta-diversity")

# mydist <- c("jaccard", "Unweighted UniFrac", "Weighted UniFrac", "bray")
mydist <- c("Unweighted UniFrac", "Weighted UniFrac", "bray")

# This links could be use for non-parametric test instead of PERMANOVA but does not work properly.
# source("https://raw.githubusercontent.com/alekseyenko/Tw2/27520182e704ec28324e74ddeab0fb584c7bbf10/code/Tw2.R")
# source("https://raw.githubusercontent.com/alekseyenko/WdStar/32053d4e8ee2cacb7a8d6e4a0ba223e0b27527ca/Wd.R")

pop.ordr <- pop.data[, 1]

# Lists for saving statistics results
mybiom.r.cdis   <- vector("list", length(pop.ordr))
mybiom.r.pcoa   <- vector("list", length(pop.ordr))
mybiom.r.stat   <- vector("list", length(pop.ordr))
mybiom.r.stat.p <- vector("list", length(pop.ordr))
mybiom.r.stat.b <- vector("list", length(pop.ordr))
mybiom.r.stat.T <- vector("list", length(pop.ordr))
mybiom.r.stat.m <- vector("list", length(pop.ordr))

for (r in 1:length(pop.ordr)) {
    # Data
    mybiom.tmp <- subset_samples(mybiom, Population == pop.ordr[[r]])
    mybiom.r   <- rarefy_even_depth(mybiom.tmp, rngseed=myseed, verbose = FALSE)

    # Correcting tree
    ## see: https://github.com/joey711/phyloseq/issues/936
    phy_tree(mybiom.r) %<>% ape::multi2di(.)

    # Lists for saving results
    mybiom.r.cdis[[r]]   <- vector("list", length(mydist))
    mybiom.r.pcoa[[r]]   <- vector("list", length(mydist))
    mybiom.r.stat[[r]]   <- vector("list", length(mydist))
    mybiom.r.stat.p[[r]] <- vector("list", length(mydist))
    mybiom.r.stat.b[[r]] <- vector("list", length(mydist))
    mybiom.r.stat.T[[r]] <- vector("list", length(mydist))
    mybiom.r.stat.m[[r]] <- vector("list", length(mydist))

    # Test for each distance index
    for (i in mydist) {
        myls.idx <- match(i, mydist)

        if (i == "jaccard") {
            mydist.tmp <- distance(mybiom.r, i, binary=TRUE)
        } else {
            mydist.tmp <- distance(mybiom.r, i)
        }

        # Centroid distance between organs/tissue
        mybiom.r.cdis[[r]][[myls.idx]] <- dist_multi_centroids(mydist.tmp, sample_data(mybiom.r)$Type)

        # PCoA
        mybiom.r.pcoa[[r]][[myls.idx]] <- ordinate(mybiom.r, "PCoA", mydist.tmp)

        # Permanova stat
        ## Reorder table and select only the current samples
        mycln     <- "Code"
        mypop.tmp <- mypop[ order(mypop[, mycln]), ]
        mypop.tmp <- mypop.tmp[ rownames(mypop.tmp) %in% labels(mydist.tmp),]

        ## Reorder the distance matrix
        mydist.tmp2 <- dist_subset(mydist.tmp, rownames(mypop.tmp))

        ## Run test
        mybiom.r.stat[[r]][[myls.idx]]   <- adonis(mydist.tmp ~ Type, data = get_variable(mybiom.r), permutations = 1000)
        set.seed(myseed)
        mybiom.r.stat.p[[r]][[myls.idx]] <- pairwise.adonis(mydist.tmp2, mypop.tmp[, mycln],  p.adjust.m = "fdr", perm = 1000)

        ## Test impact of sampling time
        # adonis(mydist.tmp ~ Collection.date + Type, data = mypop.tmp[,c("Type", "Snail.ID", "Collection.time", "Collection.date")], permutations = 1000)

        # Dispersion stat
        mybiom.r.stat.b[[r]][[myls.idx]] <- betadisper(mydist.tmp, get_variable(mybiom.r)[, mycln], type = "median", bias.adjust = TRUE)
        mybiom.r.stat.T[[r]][[myls.idx]] <- TukeyHSD(mybiom.r.stat.b[[r]][[myls.idx]])

        # P-values matrix
        mynames <- levels(mypop[, mycln])
        mymat <- matrix(NA, ncol=length(mynames), nrow=length(mynames))
        colnames(mymat) <- rownames(mymat) <- mynames
        mymat[lower.tri(mymat)] <- mybiom.r.stat.p[[r]][[myls.idx]][,7]
        mymat <- t(mymat)
        mymat[lower.tri(mymat)] <- mybiom.r.stat.T[[r]][[myls.idx]]$group[,4]
        mybiom.r.stat.m[[r]][[myls.idx]] <- mymat
    }

    # Naming slots
    names(mybiom.r.cdis[[r]])   <- mydist
    names(mybiom.r.pcoa[[r]])   <- mydist
    names(mybiom.r.stat[[r]])   <- mydist
    names(mybiom.r.stat.p[[r]]) <- mydist
    names(mybiom.r.stat.b[[r]]) <- mydist
    names(mybiom.r.stat.T[[r]]) <- mydist
    names(mybiom.r.stat.m[[r]]) <- mydist
}


# UPGMA Tree
myhc <- vector("list", length(pop.ordr)) %>% setNames(., pop.ordr)
for (r in 1:length(pop.ordr)) {
    myhc[[r]] <- vector("list", length(mydist)) %>% setNames(., mydist)
    for (i in 1:length(mydist)) {
        myhc[[r]][[i]] <- hclust(mybiom.r.cdis[[r]][[i]], method = "average")
    }
}


#~~~~~~~~~#
# Figures #
#~~~~~~~~~#

# Axes to plot
myax <- list(c(1,2), c(1,3), c(3,2))

# List to save plots
p.ls   <- vector("list", length(pop.ordr))
p.ls.r <- vector("list", length(pop.ordr))
l      <- 0

for (r in 1:length(pop.ordr)) {

    mybiom.tmp <- subset_samples(mybiom, Population == pop.ordr[r])
    mybiom.r   <- rarefy_even_depth(mybiom.tmp, rngseed = myseed, verbose = FALSE)

    p.ls[[r]] <- vector("list", length(mydist))

    k <- 0
    for (i in 1:length(mydist)) {
        for (j in 1:length(myax)) {
            mytitle <- ""
            if (j == 1) {
                l %<>% + 1
                mytitle <- paste0(LETTERS[l], ". ", pop.ordr[[r]], " - ", tools::toTitleCase(mydist[i]))
            }

            p.ls.r[[j]] <- plot_ordination(mybiom.r, mybiom.r.pcoa[[r]][[i]], color = "Code", axes=myax[[j]]) +
                                    scale_color_manual(values = org.data[,4], labels = org.data[,2]) +
                                    stat_ellipse(type = "norm") +
                                    ggtitle(mytitle) +
                                    theme(axis.text = element_text(size = rel(0.45)), axis.title = element_text(size = rel(0.6)), legend.position="none")
        }

    lgd.com <- FALSE
    lgd.pos <- NULL
    if (r == length(pop.ordr) & i == length(mydist)) {
        lgd.com <- TRUE
        lgd.pos <- "bottom"
    }

    pdf(NULL)
    p.ls[[r]][[i]] <- ggarrange(plotlist = p.ls.r, ncol = length(myax), nrow = 1, hjust = 0, common.legend = lgd.com, legend = lgd.pos)
    invisible(dev.off())

    }
}

# Flatten lists for plotting
p.ls %<>% purrr::flatten()

pdf(NULL)
p <- ggarrange(plotlist = p.ls, ncol = 1, nrow = length(p.ls), heights = c(1, 1, 1, 1.15))
invisible(dev.off())

ggsave(paste0(graph.d, "Supp. Fig. 3 - b-div_PCoA.pdf"), p, width = 3.5 * 3, height = 3.5 * length(p.ls), useDingbats = FALSE)

# List to save plots
p.ls <- vector("list", length(pop.ordr))
l    <- 0

# Breaks
mybreaks <- mybiom.r.stat.m %>% purrr::flatten() %>% sapply(., min, na.rm = TRUE) %>% min() %>% { seq(log10(.), log10(0.05), length.out=4) } %>% { c(10^(.), 1) }

for (r in 1:length(pop.ordr)) {

    # Data
    mybiom.tmp <- subset_samples(mybiom, Population == pop.ordr[[r]])
    mybiom.r   <- rarefy_even_depth(mybiom.tmp, rngseed=myseed, verbose=FALSE)

    p.ls[[r]] <- vector("list", length(mydist)) %>% setNames(., mydist)

    for (i in 1:length(mydist)) {
        myletter <- ""
        if (grepl("UniFrac", mydist[i])) {
            l %<>% +1
            myletter <- paste0(LETTERS[l], ". ")
        }
        mytitle <- paste0(myletter, pop.ordr[[r]], " - ", tools::toTitleCase(mydist[i]))

        # Plot tree
        myhc_tmp <- myhc[[r]][[i]]
        myordr <- c(grep("Tank|Tray", labels(myhc_tmp), value = TRUE), grep("Tank|Tray", labels(myhc_tmp), invert = TRUE, value = TRUE))
        myhc_tmp <- as.dendrogram(myhc_tmp) %>% dendextend::rotate(., myordr) %>% set("labels_cex", 0.75) %>% set("branches_lwd", 0.5)

        if (! grepl("Tank|Tray", labels(myhc_tmp))[1]) {
            mygrp <- partition_leaves(myhc_tmp)
            mygrp_wt    <- sapply(mygrp, function(x) grepl("Tank|Tray", x) %>% any())
            mygrp_lg    <- lapply(mygrp, length)
            mygrp_lg_wt <- (mygrp_lg[mygrp_wt] %>% unlist() %>% order(., decreasing = TRUE))[2]
            mygrp5      <- mygrp[mygrp_wt][mygrp_lg_wt] %>% unlist() %>% paste0(., collapse = "|")
            myordr      <- c(grep(mygrp5, labels(myhc_tmp), value = TRUE), grep(mygrp5, labels(myhc_tmp), value = TRUE, invert = TRUE))
            myhc_tmp    <- myhc_tmp %>% dendextend::rotate(., myordr)
        }

        # # Reroot tray with environment if needed
        # tree_tb <- cutree(myhc_tmp, k=1:attributes(myhc_tmp)$members)
        # if (! (grepl("Tank|Tray", rownames(tree_tb)[tree_tb[,2] == 2]) %>% any())) {
        #     root     <- phylogram::prune(myhc_tmp, pattern = "Tank|Tray", keep = TRUE)
        #     subtree  <- phylogram::prune(myhc_tmp, pattern = "Tank|Tray", keep = FALSE)
        #     myhc_tmp <- list(root, subtree) %>% remidpoint() %>% as.cladogram()
        # }

        # Update labels
        labels(myhc_tmp) %<>% sapply(., function(x) paste0(" ", x, " (", org.data[grep(trimws(x), org.data[,2], ignore.case = TRUE), 3], ")"))


        # mylim <- range(pretty(c(0, attributes(myhc_tmp)$height)))
        p.ls.d <- as.ggdend(myhc_tmp) %>% ggplot(., horiz = TRUE) +
                    coord_flip(clip = "off") +
                    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(), plot.margin = unit(c(1.2, 5.5, 1, 1), "lines"))

        mymat    <- mybiom.r.stat.m[[r]][[i]]
        # mybreaks <- c(10^(seq(log10(min(mymat, na.rm=T)), log10(0.05), length.out=4)), 1)

        p.ls.t <- plotPvalues(mymat, border.col = "black", pval.col = "black", font.size = 2, scientific = FALSE, digits = 5) +
                    xlab(expression(paste(beta,"-dispersion test"))) + ylab("Permanova test") +
                    scale_fill_gradient2(low = "red", mid = "yellow", high = "white", limits = c(min(mybreaks), 0.05), guide = "colorbar", na.value = "transparent", trans = "log10", breaks = mybreaks, label = format(mybreaks, digit = 2, scientific = TRUE)) +
                    guides(fill = guide_colourbar(barwidth = 9, title.vjust = 0.75)) +
                    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", plot.margin = margin(20, 1, 5.5, 5.5))

        lgd.com <- FALSE
        lgd.pos <- NULL
        if (r == length(pop.ordr) & i == length(mydist)) {
            lgd.com <- TRUE
            lgd.pos <- "bottom"
        }

        pdf(NULL)
        p.ls[[r]][[i]] <- ggarrange(p.ls.d, p.ls.t, labels = mytitle, ncol = 2, hjust = 0, common.legend = lgd.com, legend = lgd.pos)
        invisible(dev.off())

    }
}

# Flatten lists for plotting
p.ls.1 <- lapply(p.ls, function(x) x[mydist[1:2]])
p.ls.1 %<>% purrr::flatten()

pdf(NULL)
p <- ggarrange(plotlist = p.ls.1, ncol = 1, nrow = length(p.ls.1), heights = c(1, 1, 1, 1.15))
invisible(dev.off())

ggsave(paste0(graph.d, "Fig. 3 - b-div_FDR.pdf"), p, width = 5 * 2, height = 3.5 * length(p.ls.1), useDingbats = FALSE)

# Flatten lists for plotting
p.ls.2 <- lapply(p.ls, function(x) x[mydist[3]])
p.ls.2 %<>% purrr::flatten()

pdf(NULL)
p <- ggarrange(plotlist = p.ls.2, ncol = 1, nrow = length(p.ls.2), heights = c(1, 1.15))
invisible(dev.off())

ggsave(paste0(graph.d, "b-div_bray_FDR.pdf"), p, width = 5 * 2, height = 3.5 * length(p.ls.2), useDingbats = FALSE)


#---------------------#
# Taxonomic diversity #
#---------------------#

cli_h2("Taxonomic diversity")

# Data storage for graphs
p.ls   <- vector("list", length(pop.ordr))
pop.wd <- vector("list", length(pop.ordr))

#~~~~~~~~#
# Phylum #
#~~~~~~~~#

# Color
mylevel <- "Phylum"
set.seed(myseed)
mybiom.P <- mytax_glom(mybiom, mylevel) %>% core(., detection=0, prevalence=0)
myclr.P <- tax_table(mybiom.P) %>% nrow() %>% distinctColorPalette(., seed=myseed - 70)
names(myclr.P) <- tax_table(mybiom.P)[, mylevel] %>% sort()

p.ls   <- vector("list", length(pop.ordr))

for (r in 1:length(pop.ordr)) {

    mybiom.tmp <- subset_samples(mybiom, Population == pop.ordr[[r]])

    # Distance matrix
    mybiom.d <- apply(otu_table(mybiom.tmp), 1, function(x) { x <- x + 1 ; log(x) - mean(log(x)) })
    mydm     <- dist(mybiom.d, method="euclidian")

    # Cluster the data
    myhc <- hclust(mydm, method="ward.D2")

    # Generate the tree
    mydd <- as.dendrogram(myhc)
    mydd_tmp <-  mydd

    # Remove non present taxa
    mybiom.P <- mytax_glom(mybiom.tmp, mylevel)
    mybiom.P <- core(mybiom.P, detection=0, prevalence=0)

    # Sample biome
    mybiom.P.s <- transform(mybiom.P, transform="compositional")

    # Organ biome
    mybiom.P.p <- mysample_glom(mybiom.P, "Type")
    mybiom.P.p <- transform(mybiom.P.p, transform="compositional")

    # Vector color
    myclr <- myclr.P[ tax_table(mybiom.P.p)[,mylevel] ]

    # List to store graphs
    p <- vector("list", 2)

    # Tree
    # mylabels <- as.vector(dendro_data(mydd.tmp)$labels$label)
    mylabels <- labels(mydd_tmp)
    p[[1]] <- ggdendrogram(mydd_tmp) +
                scale_x_continuous(expand = c(1/(42*2),0), labels=mylabels, breaks=1:length(mylabels)) +
                theme(axis.text.y = element_blank(), plot.margin=margin(b=-10))

    # Barplot at sample level
    mybiom.P.m <- psmelt(mybiom.P.s)
    mybiom.P.m <- mybiom.P.m[ order(mybiom.P.m$Abundance, decreasing=TRUE), ]
    p[[2]] <- ggplot(mybiom.P.m, aes(x = Sample, y = Abundance, fill = !! ensym(mylevel))) +
                  geom_bar(stat = "identity") +
                  scale_fill_manual(values = myclr) +
                  scale_x_discrete(limits=mylabels) +
                  # Remove x axis title and legend
                  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="none") +
                  ylab("Relative Abundance\n")

    pdf(NULL)
    p.ls[[r]] <- ggarrange(plotlist = p, labels = paste0(LETTERS[r], ". ", pop.ordr[r]), ncol = 1, common.legend = TRUE, legend = "bottom", align = "v", hjust = 0)
    invisible(dev.off())
}

pdf(NULL)
p <- ggarrange(plotlist = p.ls, ncol = 1, align = "v")
invisible(dev.off())

ggsave(paste0(graph.d, "taxo-div.pdf"), p, width = 15, height = 9 * length(p.ls), useDingbats = FALSE)

# Per Organ
for (r in 1:length(pop.ordr)) {
    mybiom.tmp <- subset_samples(mybiom, Population == pop.ordr[[r]])

    # Remove non present taxa
    mybiom.P <- mytax_glom(mybiom.tmp, mylevel)
    mybiom.P <- core(mybiom.P, detection=1, prevalence=0)

    # Organ biome
    mybiom.P.p <- mysample_glom(mybiom.P, "Type")

    # Transform data to compositional
    mybiom.P.p <- transform(mybiom.P.p, transform="compositional")

    # List to store graphs
    p <- vector("list", 3)

    # Barplot at population level
    mybiom.P.m <- psmelt(mybiom.P.p)

    # Vector color
    myclr <- myclr.P[ tax_table(mybiom.P.p)[,mylevel] ]

    p.ls[[r]] <- ggplot(mybiom.P.m, aes(x = Code, y = Abundance, fill = !! ensym(mylevel))) +   # Use of quasiquotation !!r otherwise lastest slot of the list taken for plotting
                    geom_bar(stat = "identity", color = "black", size = 0.05) +
                    scale_fill_manual(values = myclr) +
                    ## Update legend look (because of the bar borders)
                    guides(fill = guide_legend(override.aes = list(colour = NA, size = 0.5))) +
                    # Remove x axis title and legend
                    labs(title = pop.ordr[[r]], y = "Relative Abundance", x = "") +
                    theme(legend.position = "none")
}

pdf(NULL)
p <- ggarrange(plotlist = p.ls, labels = LETTERS[1:length(p.ls)], ncol = 1, common.legend = TRUE, legend = "bottom", align = "v")
invisible(dev.off())

ggsave(paste0(graph.d, "Fig. 4 - Organ taxo-div.pdf"), p, width = 9, height = 7, useDingbats = FALSE)

#~~~~~~~#
# Class #
#~~~~~~~#

# Color Phylum
mylevel <- "Class"
set.seed(myseed)
mybiom.C <- mytax_glom(mybiom, mylevel) %>% core(., detection=0, prevalence=0)
myclr.C  <- tax_table(mybiom.C) %>% nrow() %>% distinctColorPalette(., seed=myseed - 70)
names(myclr.C) <- tax_table(mybiom.C)[, mylevel] %>% sort()

p.ls   <- vector("list", length(pop.ordr))

# Per Organ
for (r in 1:length(pop.ordr)) {
    mybiom.tmp <- subset_samples(mybiom, Population == pop.ordr[[r]])

    # Remove non present taxa
    mybiom.C <- mytax_glom(mybiom.tmp, mylevel)
    mybiom.C <- core(mybiom.C, detection=1, prevalence=0)

    # Organ biome
    mybiom.C.p <- mysample_glom(mybiom.C, "Type")

    # Transform data to compositional
    mybiom.C.p <- transform(mybiom.C.p, transform="compositional")

    # List to store graphs
    p <- vector("list", 3)

    # Barplot at population level
    mybiom.C.m <- psmelt(mybiom.C.p)

    # Vector color
    myclr <- myclr.C[ tax_table(mybiom.C.p)[,mylevel] ]

    p.ls[[r]] <- ggplot(mybiom.C.m, aes(x = Code, y = Abundance, fill = !! ensym(mylevel))) +   # Use of quasiquotation !!r otherwise lastest slot of the list taken for plotting
                    geom_bar(stat = "identity", color = "black", size = 0.05) +
                    scale_fill_manual(values = myclr) +
                    ## Update legend look (because of the bar borders)
                    guides(fill = guide_legend(override.aes = list(colour = NA, size = 0.5))) +
                    # Remove x axis title and legend
                    labs(title = pop.ordr[[r]], y = "Relative Abundance", x = "") +
                    theme(legend.position = "none")
}

pdf(NULL)
p <- ggarrange(plotlist = p.ls, labels = LETTERS[1:length(p.ls)], ncol = 1, common.legend = TRUE, legend = "bottom", align = "v")
invisible(dev.off())

ggsave(paste0(graph.d, "Supp. Fig. 4 - Organ taxo-div_class.pdf"), p, width = 9, height = 8, useDingbats = FALSE)

#-----------------------------#
# Shared taxa between samples #
#-----------------------------#

cli_h2("Shared taxa")

# pop.ordr <- md[,"Species"] %>% levels
rm.water <- FALSE

# Thresholds of ASV per intersection
tsh.ls <- c(0, 4)

ct.ls  <- vector("list", length(tsh.ls))
low.ls <- vector("list", length(tsh.ls))

# Table of ASV occurence for each sample type and species
for (i in 1:length(tsh.ls)) {

    ct.ls[[i]]  <- vector("list", length(pop.ordr)) %>% setNames(., pop.ordr)
    low.ls[[i]] <- vector("list", length(pop.ordr)) %>% setNames(., pop.ordr)

    # Minimal count in sample type
    mytsh <- 0

    # Invert prevalence of the select minimal count [0-1] (0 = all, 1 = none)
    mypv <- 0

    # Minimum number of intersection
    is.tsh <- tsh.ls[i]

    for (s in pop.ordr) {
        md.tmp <- md[ md[,"Population"] == s, ]
        myotu.tmp <- otu_table(mybiom)[,rownames(md.tmp)]
        org.data.tmp <- org.data

        # Remove irrelevant ASV
        myotu.tmp <- myotu.tmp[ rowSums(myotu.tmp) > 0, ]

        # Remove water
        if (rm.water) { org.data.tmp <- org.data.tmp[ ! grepl("T.*", org.data.tmp[,1]) ,] }

        # Filter data per sample type
        otu_type <- sapply(org.data.tmp[,1], function(x) myotu.tmp[,rownames(md.tmp[ md.tmp[,"Type.code"] == x,])] )
        otu_type <- lapply(otu_type, function(x) rowSums(x > mytsh) > ncol(x) * mypv )
        otu_type <- sapply(otu_type, c) %>% apply(., 2, function(x) as.numeric(x) %>% setNames(., names(x)))

        # Remove irrelevant ASV (after sample type filtering)
        otu_type <- otu_type[ rowSums(otu_type) > 0, ]

        # Filter on minimal intersection
        # otu_type <- otu_type[ rowSums(otu_type) > is.tsh, ]
        mycomb     <- apply(otu_type, 1, paste0, collapse = "")
        mycomb_sel <- which(table(mycomb) > is.tsh) %>% names()
        otu_sel    <- mycomb[mycomb %in% mycomb_sel] %>% names()
        otu_type   <- otu_type[otu_sel, ]

        # Store final data
        ct.ls[[i]][[s]]  <- otu_type %>% as.data.frame()
        low.ls[[i]][[s]] <- table(mycomb) > is.tsh
    }
}

# Ubiquitous ASVs (snail + water) shared between species
ubq_asv <- lapply(ct.ls[[1]], function(x) which((x %>% rowSums(.)) == 8) %>% names())
com_asv <- ubq_asv[[1]] %in% ubq_asv[[2]]

cat("\nNumber of ubiquitous ASVs:",
    "\n\t- ", names(ubq_asv)[[1]], ": ", length(ubq_asv[[1]]), " ASVs",
    "\n\t- ", names(ubq_asv)[[2]], ": ", length(ubq_asv[[2]]), " ASVs",
    "\n\t- Common:", sum(com_asv), " ASVs.\n\n", sep = "")

# Poportion of ASV shared per species
res <- sapply(ct.ls[[1]], function(x) nrow(x[rowSums(x) > 1, ]) / nrow(x))
cat("\nProportion of shared ASVs within each species:\n")
print(res)

# Poportion of ASV shared per sample type
res <- sapply(org.data[,1], function(x) lapply(ct.ls[[1]], function(y) { y1 <- y[ y[,x] == 1, ] ; nrow(y1[ rowSums(y1) > 1, ]) / nrow(y1)}))
cat("\nProportion of shared ASVs per sample type:\n")
print(res)

# Poportion of ASV per sample type
res <- sapply(org.data[,1], function(x) lapply(ct.ls[[1]], function(y) { y1 <- y[ y[,x] == 1, ] ; nrow(y1) / nrow(y)}))
cat("\nProportion of ASVs per sample type:\n")
print(res)

# Proportion of shared ASVs between H and another sample type
# lapply(ct.ls[[1]], function(x) { a <- (x[, "H"]  == 1 & (x[, "W"]  == 1 | x[, "H"]  == 1 & x[, "O"]  == 1)) ; sum(a)/length(a[rowSums(x) > 1]) })
cat("\nProportion of shared ASVs between hemolymph and another sample types (but not limited to the intersection):\n")
prop_asv(ct.ls[[1]], org.data[,1], "H") %>% print()

# Proportion of shared ASVs between whole snail and another sample type
cat("\nProportion of private shared ASVs between whole snail and another sample type (but not found outside of the intersections):\n")
prop_asv(ct.ls[[1]], org.data[,1], "W", private = TRUE) %>% print()

# Proportion of shared ASVs between TK|TY and another sample type
# lapply(ct.ls[[1]], function(x) { a <- (x[, "TK"]  == 1 | x[, "TY"]  == 1) & (x[, "H"]  == 1 | x[, "O"]  == 1 | x[, "W"]  == 1) ; sum(a)/length(a[rowSums(x) > 1]) })
cat("\nProportion of shared ASVs between water (TK|TY) and another sample type (but not limited to the intersection):\n")
prop_asv(ct.ls[[1]], org.data[,1], c("TK", "TY")) %>% print()

# Proportion of shared ASVs between TK|TY and G|S
cat("\nProportion of shared ASVs between water (TK|TY) and digestive tract (G|S) (but not limited to the intersection):\n")
sapply(ct.ls[[1]], function(x) { a <- (x[, "TK"]  == 1 | x[, "TY"]  == 1) & (x[, "G"]  == 1 | x[, "S"] == 1) ; sum(a)/length(a[rowSums(x) > 1]) }) %>% print()

# Proportion of ASVs not found in water (including fully private ASVs)
cat("\nProportion of ASVs not found in  water (TK|TY) (including fully private ASVs to sample type):\n")
sapply(org.data[, 1], function(y) lapply(ct.ls[[1]], function(x) { a <- (x[, "TK"]  == 1 | x[, "TY"]  == 1) ; b <- rowSums(x[! a & x[, y]  == 1,]) >= 1 ; sum(b) / sum(!a) }))

# Proportion of shared ASVs not found in water (excluing fully private ASVs)
cat("\nProportion of shared ASVs not found in  water (TK|TY) (excluding fully private ASVs to sample type):\n")
sapply(org.data[, 1], function(y) lapply(ct.ls[[1]], function(x) { a <- (x[, "TK"]  == 1 | x[, "TY"]  == 1) ; b <- rowSums(x[! a & x[, y]  == 1,]) > 1 ; sum(b) / sum(!a) }))

# Proportion of ASV found only in snails (amoung all ASVs)
cat("\nProportion of ASVs only found in snails:\n")
sapply(ct.ls[[1]], function(x) sum(! rowSums(x[, grepl("^T", colnames(x))]) > 0) / nrow(x))

# Proportion of shared ASV found only in snails (amoung all shared ASVs)
#lapply(ct.ls[[1]], function(x) {x1 <- x[rowSums(x) > 1, ] ; sum(! rowSums(x1[, grepl("^T", colnames(x1))]) > 0) / nrow(x1) })

# Proportion of ubiquitous shared ASVs across snail sample types
cat("\nProportion of ubiquitous shared ASVs only found in all snail sample types:\n")
sapply(ct.ls[[1]], function(x) {x <- x[rowSums(x) > 1, ]
                                xT <- x[! rowSums(x[, grepl("^T", colnames(x))]) > 0,  ! grepl("^T", colnames(x))]
                                sum(rowSums(xT) == ncol(xT)) / nrow(xT) }) %>% print()

# Proportion of ubiquitous shared ASVs across all sample types
cat("\nProportion of ubiquitous shared ASVs found in all sample types (snail and water):\n")
sapply(ct.ls[[1]], function(x) {sum(rowSums(x) == ncol(x)) / nrow(x[rowSums(x) > 1,]) }) %>% print()

# Proportion of intersection with low sharing
cat("\nProportion of intersections with low ASV sharing (<= ", tsh.ls[-1], "):\n", sep = "")
sapply(2:length(tsh.ls), function(y) lapply(low.ls[[y]], function(x) 1 - sum(x) / length(x))) %>% print()

# Proportion of sample types associated with intersection with low ASVs
cat("\nProportion of sample types associated with low ASVs intersections:\n")
sapply(low.ls[[2]], function(x) {names(x)[x] %>% sapply(., strsplit, "") %>% do.call("rbind", .) %>% apply(., 2, as.numeric) %>% colSums() %>% setNames(., org.data[,1])} ) %>% t()

# Pairwise sharing
org.ordr <- org.data[, 1]
pw.ls <- vector("list", length(pop.ordr)) %>% setNames(., pop.ordr)
for (s in pop.ordr) {
    x <- ct.ls[[1]][[s]]
    cb <- combn(org.ordr, 2)
    # mydf <- data.frame(a = cb[1,], b = cb[2,], c = rep(NA, ncol(cb)))
    mydf <- matrix(NA, nrow = length(org.ordr), ncol = length(org.ordr))
    colnames(mydf) <- org.ordr
    # rownames(mydf) <- rev(org.ordr)
    rownames(mydf) <- org.ordr
    for (n in 1:ncol(cb)) {
        a <- (x[, cb[1,n]]  == 1 & x[, cb[2,n]]  == 1)
        mydf[cb[1,n], cb[2,n]] <- (sum(a) / length(a[rowSums(x) > 1])) %>% round(., digit = 2)
    }
    pw.ls[[s]] <- mydf
}

pw.shr <- pw.ls[[1]]
pw.shr[lower.tri(pw.shr)] <- (pw.ls[[2]] %>% t())[lower.tri(pw.shr)]
colnames(pw.shr) %<>% match(., org.data[,1]) %>% org.data[., 3]
rownames(pw.shr) %<>% match(., org.data[,1]) %>% org.data[., 3]


# ASV shared by species for each sample type
rm.water <- FALSE
org.ordr <- org.data[, 1]

ct.sp.ls <- vector("list", length(org.ordr)) %>% setNames(., org.ordr)

for (s in org.ordr) {
    md.tmp <- md[ md[,"Type.code"] == s, ]
    myotu.tmp <- otu_table(mybiom)[,rownames(md.tmp)]
    org.data.tmp <- org.data

    # Remove irrelevant ASV
    myotu.tmp <- myotu.tmp[ rowSums(myotu.tmp) > 0, ]

    # Remove water
    if (rm.water) { org.data.tmp <- org.data.tmp[ ! grepl("T.*", org.data.tmp[,1]) ,] }

    # myclr <- pop.clr[ ! grepl("T.*", org.data.tmp[,1]) ]

    # Filter data per sample type
    mytsh <- 0
    mypv  <- 0
    otu_type <- sapply(pop.ordr, function(x) myotu.tmp[,rownames(md.tmp[ md.tmp[,"Population"] == x,])], simplify = FALSE )
    otu_type <- lapply(otu_type, function(x) rowSums(x > mytsh) > ncol(x) * mypv )
    otu_type <- sapply(otu_type, c) %>% apply(., 2, function(x) as.numeric(x) %>% setNames(., names(x)))

    # Remove irrelevant ASV (after sample type filtering)
    otu_type <- otu_type[ rowSums(otu_type) > 0, ]

    # Filter on minimal intersection
    mycomb     <- apply(otu_type, 1, paste0, collapse = "")
    mycomb_sel <- which(table(mycomb) > is.tsh) %>% names()
    otu_sel    <- mycomb[mycomb %in% mycomb_sel] %>% names()
    otu_type   <- otu_type[otu_sel, ]

    # Store final data
    ct.sp.ls[[s]] <- otu_type %>% as.data.frame()
}


#~~~~~~~~~#
# Figures #
#~~~~~~~~~#

# Shared ASVs by sample types for each snail species
for (i in 1:length(tsh.ls)) {

    # Threshold values
    # mytsh <- tsh.ls[[i]][[1]]
    # mypv  <- tsh.ls[[i]][[2]]
    mytsh <- 0
    mypv  <- 0

    fig_pfx <- "Fig. 5"
    myh <- 5
    myw <- 10
    if (i == 1) {
        fig_pfx <- "Supp. Fig. 5"
        myw <- 20
    }

    for (s in pop.ordr) {
        pdf(paste0(graph.d, fig_pfx, " - ", s," ASV upset - tsh ", mytsh, " pv ", mypv, ".pdf"), height = myh, width = myw, useDingbats = FALSE, onefile = FALSE)
        otu_type_fn <- ct.ls[[i]][[s]]
        clr <- pop.data[ pop.data[, 1] == s, 2]
        ups <- upset(otu_type_fn, sets = colnames(otu_type_fn), nintersects = NA, main.bar.color = clr, mainbar.y.label = "ASVs per intersection", sets.x.label = "ASVs per sample type", keep.order = TRUE)
        print(ups)
        grid.text(s, x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
        invisible(dev.off())
    }
}


# Shared ASVs by species for each sample type
p.ls <- vector("list", length(org.ordr)) %>% setNames(., org.ordr)

myw <- 2.75 * length(org.ordr) / 2
myh <- 2.75 * length(org.ordr) / 4

myclr <- pop.data[,2] %>% c(., midcol(.[1], .[2]))

for (s in org.ordr) {
    otu_type_fn <- ct.sp.ls[[s]]

    z <- apply(otu_type_fn, 1, function(x) {paste(names(x)[as.logical(x)], collapse="-")}) %>% table() %>% as.data.frame()
    z[, 1] %<>% levels() %>% {list(grep("-", ., value = TRUE, invert = TRUE), grep("-", ., value = TRUE))} %>% unlist() %>% factor(z[, 1], levels = .)

    p.ls[[s]] <- ggplot(data = z, aes(x = ., y = Freq)) +
                    geom_bar(stat = "identity", fill = myclr) +
                    geom_text(aes(label = Freq), vjust = -0.3, color = myclr[match(z[,1], levels(z[,1]))]) +
                    labs(x = "", y = "ASVs per intersection", title = org.data[org.data[,1] == s, 2]) +
                    coord_cartesian(clip = 'off') +
                    theme(plot.title = element_text(hjust = 0.5))
}


# Panel with replicate
pdf(paste0(graph.d, "Fig. 6 - ASV upset by spl type.pdf"), width = myw, height = myh)
p <- ggarrange(plotlist = p.ls, nrow = length(p.ls) / 4, ncol = length(p.ls) / 2)
print(p)
invisible(dev.off())


# Pairwise ASV sharing
pdf(NULL)
myclr <- pw.shr
myclr[lower.tri(myclr)] <- pop.data[names(pw.ls)[1] == pop.data[,1], 2]
myclr[upper.tri(myclr)] <- pop.data[names(pw.ls)[2] == pop.data[,1], 2]

p <- plotPvalues(pw.shr, border.col = "black", pval.col = as.vector(myclr), font.size = 3, scientific = FALSE) +
        xlab("") + ylab("") +
        scale_fill_gradient2(low = "white", mid = "skyblue", high = "steelblue4", limits = c(0.2, 0.6), guide = "colorbar", na.value = "transparent") +
        guides(fill = guide_colourbar(barwidth = 9, title.vjust = 0.75)) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", axis.text.x = element_text(vjust = 0.5, hjust = 1), plot.margin = margin(2, 2, 2, 2))
invisible(dev.off())
ggsave(paste0(graph.d, "Supp. Fig. 6 - pairwise shared ASV.pdf"), p, width = 5, height = 5, useDingbats = FALSE)
