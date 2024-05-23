#!/usr/bin/env Rscript
# Title: microbiome_diversity.R
# !! Version: 1.0  !!
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# !! Created in: 2020-03-25 !!

#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    # System
    library("magrittr")
    library("plyr")             # function count
    library("tidyr")
    library("reshape")

    # Microbiome
    # library("DESeq2")         # !! TO INSTALL !!
    library("microbiome")       # function core
    library("qiime2R")
    library("phyloseq")
    library("picante")          # PD function
    library("pairwiseAdonis")   # see github
    library("phylogram")        # prune function
    # # library("pvclust")        # !! TO INSTALL !!
    library("dendextend")     # pvclust_edges function    !! TO INSTALL !!
    library("usedist")        # dist_multi_centroids function !! TO INSTALL !!
    library("lme4")
    library("lmerTest")

    # Graph
    library("colorspace")       # function rainbow_hcl
    # library("cowplot")
    library("ggpubr")
    library("ggdendro")
    library("grid")             # function grid.rect
    library("lemon")
    library("UpSetR")

# Maybe add vegan needed by shiftlegend
})

#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

cat("Loading functions...\n")
source("microbiome_diversity_func.R")

#===========#
# Variables #
#===========#

cat("Setting variables...\n")

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

# Seed for reproducibility
myseed <- 26017
set.seed(myseed)

#=================#
# Data processing #
#=================#

if (! dir.exists(graph.d)) { dir.create(graph.d, recursive=TRUE) }

#----------------#
# Preparing data #
#----------------#

# Loading ASV data
asv    <- read_qza(asv.f)$data
asv.tb <- otu_table(asv, taxa_are_rows=TRUE)

# Loading contaminants
blast <- read.delim(blast.f, header=FALSE)

# Loading taxa data and reorganizing table
taxa <- read_qza(taxa.f)$data

## Split info
cln.nm <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa   <- separate(taxa, 2, cln.nm, sep=";")

## Cleaning names
taxa[, (1:length(cln.nm))+1] <- apply(taxa[,(1:length(cln.nm))+1], 2, function(x) gsub(".*__","", x))

## Polishing data
rownames(taxa) <- taxa[,1]
taxa <- taxa[,-c(1,ncol(taxa))]

## Assign "Unassigned" to all level of unassigned ASV
taxa[ taxa[,1] == "Unassigned", ] <- "Unassigned"

## Convert table
taxa.tb             <- tax_table(as.matrix(taxa))
taxa_names(taxa.tb) <- rownames(taxa)

# Loading tree
tree <- read_qza(tree.f)$data

# Loading metadata and reorder factors
md <- read.delim(md.f, row.names=1)
md[,"Population"] <- factor(md[,"Population"], pop.data[,1])
md[,"Type.code"] <- factor(md[,"Type.code"], org.data[,1])

# Building phyloseq object
mybiom <- phyloseq(asv.tb, taxa.tb, tree, sample_data(md))

# Populations
mypop <- md

# Cleaning data from contaminants
nb.asv   <- nrow(tax_table(mybiom))
nb.unass <- (taxa[,1] == "Unassigned") %>% sum()

## Mitochondria and chloroplast
my.mt <- apply(tax_table(mybiom), 1, function(x) any(grepl("mitochondria", x, ignore.case=TRUE)))
my.cp <- apply(tax_table(mybiom), 1, function(x) any(grepl("chloroplast", x, ignore.case=TRUE)))

## Eukaryote contaminants among the unassigned
my.euk <- blast[ blast[,14] == "Eukaryota", 1]
my.euk <- rownames(tax_table(mybiom)) %in% my.euk

## Exporting contaminent list
write.table(rownames(tax_table(mybiom))[(my.mt | my.cp | my.euk)], "../data/contaminants", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Remove contaminants
mybiom <- subset_taxa(mybiom, ! (my.mt | my.cp | my.euk))

## Final numbers
cat("\nNumber of ASVs removed:",
    "\n\t- mitochondria and chloroplasts:", sum(my.mt | my.cp), "/", signif(sum(my.mt | my.cp) / nb.asv * 100, 2), "% of the total ASVs",
    "\n\t- eukaryotes:", sum(my.euk), "/", signif(sum(my.euk) / nb.unass * 100, 2), "% of the unassigned ASVs", "/", signif(sum(my.euk) / nb.asv * 100, 2), "% of the total ASVs",
    "\n\t- total:", sum(my.mt | my.cp | my.euk), "/", signif(sum(my.mt | my.cp | my.euk)/ nb.asv * 100, 2), "% of the total ASVs", "\n\n")

# Number per snail species
otu_nb_sp <- vector("list", nrow(pop.data)) %>% setNames(., pop.data[,1])

for (i in pop.data[, 1]) {
    mysp  <- paste0(i, ".")
    mycln <- otu_table(mybiom) %>% colnames() %>% grep(org.tp[,1] %>% paste0(mysp, ., collapse="|"), ., value = TRUE)
    otu_nm <- ((otu_table(mybiom)[, mycln] %>% rowSums()) > 0) %>% which() %>% names()
    otu_unass <- otu_nm[otu_nm %in% rownames(tax_table(mybiom)[ tax_table(mybiom)[,1] == "Unassigned", ])]
    otu_nb_sp[[i]] <- list(otu_nm, otu_unass)
}

cat("\nNumber of ASVs per snail species:\n")
for (i in pop.data[, 1]) {
    cat("\t- ", i, ": ", length(otu_nb_sp[[i]][[1]]), " (", length(otu_nb_sp[[i]][[2]]), " unassigned)\n", sep = "")
}
cat("\t- Common ASVs: ", intersect(otu_nb_sp[[1]][[1]], otu_nb_sp[[2]][[1]]) %>% length(), " (",  intersect(otu_nb_sp[[1]][[2]], otu_nb_sp[[2]][[2]]) %>% length(), ")\n", sep = "")

#--------------------#
# Rarefaction curves #
#--------------------#

p <- ggrare(mybiom, step = 1000, color = "Population", se = FALSE, plot=FALSE) +
        scale_color_manual(values = pop.data[,2]) +
        facet_wrap(~Type.code)

# Adjust legend
pdf(NULL)
p <- shift_legend(p)
invisible(dev.off())

ggsave(paste0(graph.d, "Supp. Fig. 1 - rarefaction.pdf"), p, height=7, width=7, useDingbats=FALSE)

#-----------#
# Diversity #
#-----------#

myalpha <- microbiome::alpha(mybiom)
mypd    <- pd(t(asv), tree)
myalpha <- merge(myalpha, mypd, by="row.names")
myalpha <- merge(myalpha, md, by.x=1, by.y="row.names")

myidx <- matrix(c("observed",           "Observed richness",
                  "PD",                 "Faith's Phylogenetic diversity",
                  "evenness_simpson",   "Simpson evenness"
                  ), ncol=2, byrow=TRUE)

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

# GLM
myalpha2 <- myalpha[ myalpha[, "Species"] != "Water", ]
myalpha2[, "Snail.ID"] <- paste0(myalpha2[, "Snail.ID"], myalpha2[, "Tank.ID"])
myalpha2[, "Tank.ID"] %<>% as.factor()

#myalpha3 <- myalpha2[ myalpha2[, "Type"] != "Whole snail ", ]

##indices <- c("observed", "chao1", "diversity_shannon", "evenness_simpson", "dominance_simpson", "PD")
#indices <- myidx[,1]
## myconds <- c("Collection.time", "Collection.date", "Tank.ID")
#myconds <- c("Collection.time", "Collection.date")
#myglms  <- vector("list", length(indices)) %>% setNames(., indices)
## for (a in 2:24) {
#for (a in indices) {
#    myglms[[a]] <- vector("list", length(myconds)) %>% setNames(., myconds)
#    for (i in myconds) {
#        # myglm <- lmer(paste(colnames(myalpha2)[a], "~ Species * Type *", i, "+ (1|Snail.ID)"), data = myalpha2)
#        myglm   <- lmer(paste(a, "~ Species * Type *", i, "+ (1|Snail.ID)"), data = myalpha2)
#        print(BIC(myglm))
#        myanova <- anova(myglm)
#        myglms[[a]][[i]] <- list(glm = myglm, anova = myanova)
#    }
#}

indices <- myidx[,1]
myglms  <- vector("list", length(indices)) %>% setNames(., indices)
for (a in indices) {
    ## myglms[[a]] <- vector("list", length(myconds)) %>% setNames(., myconds)
    # myglm   <- lmer(paste(a, "~ Species * Type + (1|Snail.ID) + (1|Tank.ID)"), data = myalpha2) ## Note: This model gives the same results has the one below.
    myglm   <- lmer(paste(a, "~ Species * Type + (1|Snail.ID)"), data = myalpha2)
    print(BIC(myglm))
    myanova <- anova(myglm)
    mycoef <- myglm %>% summary() %>% coef()

    # Output results
    append <- ifelse(match(a, indices) == 1, FALSE, TRUE)
    outfile <- "a-div_lmer_results.txt"
    write(a, outfile, append = append)
    write.table(mycoef, outfile, append = TRUE, sep = "\t")

    # Store results
    myglms[[a]] <- list(glm = myglm, anova = myanova, coef = mycoef)
}

myconds <- c("Collection.time", "Collection.date")
myglms2 <- vector("list", length(pop.data[,1])) %>% setNames(., pop.data[,1])
for (p in pop.data[,1]) {
    myrows  <- myalpha2[, "Population"] == p
    myglms2[[p]] <- vector("list", length(indices)) %>% setNames(., indices)
    for (a in indices) {
        myglms2[[p]][[a]] <- vector("list", length(myconds)) %>% setNames(., myconds)
        for (i in myconds) {
            myglm   <- lmer(paste(a, "~ Type *", i, "+ (1|Snail.ID)"), data = myalpha2[myrows, ])
            print(BIC(myglm))
            myanova <- anova(myglm)
            mycoef <- myglm %>% summary() %>% coef()

            # Output results
            append <- ifelse(match(a, indices) == 1 & match(i, myconds) == 1, FALSE, TRUE)
            outfile <- "a-div_cond_lmer_results.txt"
            write(paste(a, i), outfile, append = append)
            write.table(mycoef, outfile, append = TRUE, sep = "\t")

            # Store results
            myglms2[[p]][[a]][[i]] <- list(glm = myglm, anova = myanova, coef = mycoef)
        }
    }
}

#~~~~~~~~~#
# Figures #
#~~~~~~~~~#

# Boxplot for each relevant index
## No stat added on graph because unreadable
p <- vector("list", nrow(myidx))
for (i in myidx[,1]) {
    df.tmp <- data.frame(i=myalpha[,i], Tissue=myalpha[,"Type.code"], Species=myalpha[,"Population"])
    p[[match(i,myidx)]] <- ggplot(df.tmp, aes(x=Tissue, y=i, fill=Species)) +
#            geom_blank(data=df.tmp, aes(x=Tissue, y=i*1.1)) +
            geom_boxplot() + xlab("") + ylab("") + ggtitle(myidx[myidx[,1] == i,2]) +
#                stat_summary(geom = 'text', label = mya.letters[[r]][[match(i,myidx[,1])]][,1], fun.y = max, vjust = -1) +
#            stat_summary(aes(group=Species), geom="text", fun.y = max, vjust = -1) +
#            stat_summary(aes(group=Species)) +
            scale_fill_manual(values=pop.data[,2]) +
            theme(plot.title = element_text(hjust =0.5), legend.position="none")
}

pdf(NULL)
p <- ggarrange(plotlist = p, labels = LETTERS[1:length(p)], ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv")
invisible(dev.off())
ggsave(paste0(graph.d, "Fig. 2 - alpha-div.pdf"), p, width = 10, height = 3, useDingbats = FALSE)

# P-value matrix for each relevant index
p <- vector("list", nrow(myidx))
mybreaks <- sapply(mya.test, function(x) min(x$p.value, na.rm=T)) %>% min() %>% c(10^(seq(log10(.), log10(0.05), length.out=4)), 1)
for (i in 1:nrow(myidx)) {
    # mybreaks <- c(10^(seq(log10(min(mya.test[[i]]$p.value, na.rm = T)), log10(0.05), length.out = 4)), 1)
    x_clr <- rep(pop.data[1, 2], ncol(mya.test[[i]]$p.value))
    x_clr[ grepl(pop.data[2, 1], colnames(mya.test[[i]]$p.value)) ] <- pop.data[2, 2]
    y_clr <- rev(x_clr)
    p[[i]] <- plotPvalues(mya.test[[i]]$p.value, border.col="black", pval.col="black", font.size=3) +
        xlab("") + ylab("") +
        scale_fill_gradient2(low = "red", mid = "yellow", high = "white", limits = c(NA, 0.05), guide = "colorbar", na.value = "transparent", trans = "log10", breaks = mybreaks, label = format(mybreaks, digit = 2, scientific = TRUE)) +
        guides(fill = guide_colourbar(barwidth = 9, title.vjust = 0.75)) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", axis.text.x = element_text(color = x_clr, angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(color = y_clr)) +
        ggtitle(myidx[i, 2])
}

# pdf(paste0(graphs.d, "Supp. Fig. 1 - alpha-div_p-val.pdf"), height = 10, width = 10)
# lapply(p, print) %>% invisible()
# invisible(dev.off())
pdf(NULL)
p <- ggarrange(plotlist = p, labels = LETTERS[1:length(p)], ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv")
invisible(dev.off())
ggsave(paste0(graph.d, "Supp. Fig. 2 - alpha-div_p-val_v3.pdf"), p, width = 30, height = 10, useDingbats = FALSE)

#------#
# PCoA #
#------#

# For PCoA
# mydist <- c("jaccard", "Unweighted UniFrac", "Weighted UniFrac", "bray")
mydist <- c("Unweighted UniFrac", "Weighted UniFrac")

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
        mypop.tmp <- mypop[ order(mypop[, "Type"]), ]
        mypop.tmp <- mypop.tmp[ rownames(mypop.tmp) %in% labels(mydist.tmp),]

        ## Reorder the distance matrix
        mydist.tmp2 <- dist_subset(mydist.tmp, rownames(mypop.tmp))

        ## Run test
        mybiom.r.stat[[r]][[myls.idx]]   <- adonis(mydist.tmp ~ Type, data = get_variable(mybiom.r), permutations = 1000)
        set.seed(myseed)
        # mybiom.r.stat.p[[r]][[myls.idx]] <- pairwise.adonis(mydist.tmp2, mypop.tmp[, "Type"], perm = 1000)
        mybiom.r.stat.p[[r]][[myls.idx]] <- pairwise.adonis(mydist.tmp2, mypop.tmp[, "Type"],  p.adjust.m = "fdr", perm = 1000)

        ## Test impact of sampling time TO BE DONE AND INCLUDE IN PAPER
        adonis(mydist.tmp ~ Collection.date + Type, data = mypop.tmp[,c("Type", "Snail.ID", "Collection.time", "Collection.date")], permutations = 1000)

        # Dispersion stat
        mybiom.r.stat.b[[r]][[myls.idx]] <- betadisper(mydist.tmp, get_variable(mybiom.r)[, "Type"], type = "centroid")
        mybiom.r.stat.T[[r]][[myls.idx]] <- TukeyHSD(mybiom.r.stat.b[[r]][[myls.idx]])

        # P-values matrix
        mynames <- sort(mypop[, "Type"]) %>% unique
        mymat <- matrix(NA, ncol=length(mynames), nrow=length(mynames))
        colnames(mymat) <- rownames(mymat) <- mynames
        # mymat[lower.tri(mymat)] <- mybiom.r.stat.p[[r]][[myls.idx]][,6]
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

# # Centroid distances between species and populations

# ## List species and remove water
# sp.vec <- sample_data(mybiom)$Species %>% as.character() %>% unique()
# sp.vec <- sp.vec[! sp.vec %in% "Water"]

# ## List samples per species
# spl.ls        <- sapply(sp.vec, function(x) sample_data(mybiom)[ sample_data(mybiom)$Species == x, "Population"] %>% unique())
# names(spl.ls) <- sp.vec

# for (r in 1:length(pop.ordr)) {
#     cat("\nCentroid distances for replicate ", r, ":\n", sep="")

#     for(i in mydist) {
#         mymt <- mybiom.r.cdis[[r]][[i]] %>% as.matrix()

#         # Add NA
#         mymt[upper.tri(mymt, diag = TRUE)] <- NA

#         # Remove water
#         myspl <- ! rownames(mymt) %in% "Water"
#         mymt  <- mymt[myspl, myspl]

#         # Test selection
#         norm <- lapply(spl.ls, function(x) shapiro.test(as.vector(mymt[, x]))$p.value > 0.05) %>% unlist() %>% all()

#         # Test
#         if (norm) {
#             p.val <- t.test(mymt[, spl.ls[[1]]], mymt[, spl.ls[[2]]])$p.value %>% signif(., 2)
#             test  <- "T test"
#         } else {
#             p.val <- wilcox.test(mymt[, spl.ls[[1]]], mymt[, spl.ls[[2]]])$p.value %>% signif(., 2)
#             test  <- "Wilcoxon test"
#         }

#         # Average distance
#         average.d <- lapply(spl.ls, function(x) mean(mymt[, x], na.rm=TRUE) %>% signif(., 2))

#         cat("\t- ", tools::toTitleCase(i), " distance: the average distance between Ba and Bg is ", average.d[[1]], " and the average distance within Bg is ", average.d[[2]], " (", test, "; p-value :", p.val, ").\n", sep="")

#     }
# }

# UPGMA Tree
myhc <- vector("list", length(pop.ordr)) %>% setNames(., pop.ordr)
for (r in 1:length(pop.ordr)) {
    myhc[[r]] <- vector("list", length(mydist)) %>% setNames(., mydist)
    for (i in 1:length(mydist)) {
        myhc[[r]][[i]] <- hclust(mybiom.r.cdis[[r]][[i]], method = "average")
    }
}

# GLM
myalpha2 <- myalpha[ myalpha[, "Species"] != "Water", ]
myalpha2[, "Snail.ID"] <- paste0(myalpha2[, "Snail.ID"], myalpha2[, "Tank.ID"])
myalpha2[, "Tank.ID"] %<>% as.factor()

##indices <- c("observed", "chao1", "diversity_shannon", "evenness_simpson", "dominance_simpson", "PD")
#indices <- myidx[,1]
## myconds <- c("Collection.time", "Collection.date", "Tank.ID")
#myconds <- c("Collection.time", "Collection.date")
#myglms  <- vector("list", length(indices)) %>% setNames(., indices)
## for (a in 2:24) {
#for (a in indices) {
#    myglms[[a]] <- vector("list", length(myconds)) %>% setNames(., myconds)
#    for (i in myconds) {
#        # myglm <- lmer(paste(colnames(myalpha2)[a], "~ Species * Type *", i, "+ (1|Snail.ID)"), data = myalpha2)
#        myglm   <- lmer(paste(a, "~ Species * Type *", i, "+ (1|Snail.ID)"), data = myalpha2)
#        print(BIC(myglm))
#        myanova <- anova(myglm)
#        myglms[[a]][[i]] <- list(glm = myglm, anova = myanova)
#    }
#}

indices <- myidx[,1]
myglms  <- vector("list", length(indices)) %>% setNames(., indices)
for (a in indices) {
    myglms[[a]] <- vector("list", length(myconds)) %>% setNames(., myconds)
    # myglm   <- lmer(paste(a, "~ Species * Type + (1|Snail.ID) + (1|Tank.ID)"), data = myalpha2) ## Note: This model gives the same results has the one below.
    myglm   <- lmer(paste(a, "~ Species * Type + (1|Snail.ID)"), data = myalpha2)
    print(BIC(myglm))
    myanova <- anova(myglm)
    myglms[[a]] <- list(glm = myglm, anova = myanova)
}

#~~~~~~~~~#
# Figures #
#~~~~~~~~~#

# Axes to plot
myax <- list(c(1,2), c(1,3), c(3,2))

# # List to save plots
# p.ls.r <- vector("list", length(pop.ordr))

# for (r in 1:length(pop.ordr)) {

#     # Data
#     mybiom.tmp <- subset_samples(mybiom, Population == pop.ordr[[r]])
#     mybiom.r <- rarefy_even_depth(mybiom.tmp, rngseed=myseed, verbose=FALSE)

#     p.ls.r[[r]] <- vector("list", length(mydist) * length(myax))

#     k <- 0
#     for (i in 1:length(mydist)) {
#         for (j in 1:length(myax)) {
#             mytitle <- ""
#             if (j == 1) { mytitle <- tools::toTitleCase(mydist[i]) }

#             k <- k + 1
#             p.ls.r[[r]][[k]] <- plot_ordination(mybiom.r, mybiom.r.pcoa[[r]][[i]], color = "Type", axes=myax[[j]]) +
#                            # scale_color_manual(values = org.data[,3]) +
#                            stat_ellipse(type = "norm") +
#                            # stat_ellipse(type = "euclid") +
#                            # stat_ellipse(type = "t", geom = "polygon") +
#                            ggtitle(mytitle) +
#                            theme(axis.text = element_text(size = rel(0.45))) +
#                            theme(axis.title = element_text(size = rel(0.6))) +
#                            theme(legend.position="none")
#         }
#     }
# }

# myletters <- lapply(1:length(mydist), function(x) c(LETTERS[x], rep("",length(myax)-1))) %>% unlist()

# mylg <- length(p.ls.r[[1]])
# pdf(NULL)
# p <- ggarrange(plotlist=p.ls.r[[1]], nrow=length(p.ls.r[[1]])/3, ncol=3, common.legend=TRUE, legend="bottom", labels=myletters)
# invisible(dev.off())

# ggsave(paste0(graph.d,"b-div.pdf"), p, width=3.5*3, height=3.5*length(p.ls.r[[1]])/3, useDingbats=FALSE)

# # Panel with replicates
# pdf(NULL)
# p <- ggarrange(ggarrange(plotlist=p.ls.r[[1]], nrow=length(p.ls.r[[1]])/3, ncol=3, common.legend=TRUE, legend="bottom"),
#                # as_ggplot(grid.lines(c(0.5,0.5), c(0,1), gp=gpar(lwd=1.5, col="black"))),
#                ggarrange(plotlist=p.ls.r[[2]], nrow=length(p.ls.r[[2]])/3, ncol=3, common.legend=TRUE, legend="bottom"),
#                ncol=length(pop.ordr)+1, nrow=1, widths=c(1,0.01,1), labels=c(LETTERS[1], "", LETTERS[2]))
# invisible(dev.off())

# ggsave(paste0(graph.d,"Supp. Fig. 4 - b-div_all.pdf"), p, width=3.5*3*2, height=3.5*length(p.ls.r[[1]])/3, useDingbats=FALSE)

# List to save plots
p.ls.t <- vector("list", length(pop.ordr))
p.ls.d <- vector("list", length(pop.ordr))
mytitle <- vector("list", length(pop.ordr))

for (r in 1:length(pop.ordr)) {

    # Data
    mybiom.tmp <- subset_samples(mybiom, Population == pop.ordr[[r]])
    mybiom.r   <- rarefy_even_depth(mybiom.tmp, rngseed=myseed, verbose=FALSE)

    p.ls.t[[r]] <- vector("list", length(mydist))
    p.ls.d[[r]] <- vector("list", length(mydist))
    mytitle[[r]] <- vector("list", length(mydist))

    for (i in 1:length(mydist)) {
        # mytitle <- paste(" ", pop.ordr[[r]], "-", tools::toTitleCase(mydist[i]))
        mytitle[[r]][[i]] <- paste(pop.ordr[[r]], "-", tools::toTitleCase(mydist[i]))

        # Plot tree
        myhc_tmp <- myhc[[r]][[i]]
        myordr <- c(grep("Tank|Tray", labels(myhc_tmp), value = TRUE), grep("Tank|Tray", labels(myhc_tmp), invert = TRUE, value = TRUE))
        myhc_tmp <- as.dendrogram(myhc_tmp) %>% dendextend::rotate(., myordr) %>% set("labels_cex", 0.75) %>% set("labels", paste(" ", labels(.))) %>% set("branches_lwd", 0.5)

        if (! grepl("Tank|Tray", labels(myhc_tmp))[1]) {
            myordr <- c(grep("Whole snail|Hemolymph", labels(myhc_tmp), invert = TRUE, value = TRUE), grep("Whole|Hemolymph", labels(myhc_tmp), value = TRUE))
            myhc_tmp <- myhc_tmp %>% dendextend::rotate(., myordr)
        }

        # # Reroot tray with environment if needed
        # tree_tb <- cutree(myhc_tmp, k=1:attributes(myhc_tmp)$members)
        # if (! (grepl("Tank|Tray", rownames(tree_tb)[tree_tb[,2] == 2]) %>% any())) {
        #     root     <- phylogram::prune(myhc_tmp, pattern = "Tank|Tray", keep = TRUE)
        #     subtree  <- phylogram::prune(myhc_tmp, pattern = "Tank|Tray", keep = FALSE)
        #     myhc_tmp <- list(root, subtree) %>% remidpoint() %>% as.cladogram()
        # }

        # mylim <- range(pretty(c(0, attributes(myhc_tmp)$height)))
        p.ls.d[[r]][[i]] <- as.ggdend(myhc_tmp) %>% ggplot(., horiz = TRUE) +
                                # ylab("Heights") +
                                coord_flip(clip = "off") +
                                # coord_flex_flip(bottom = capped_horizontal("both")) +
                                # labs(title = mytitle) +
                                theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(), plot.margin = unit(c(1.2, 5, 1, 1), "lines"))
                                # theme(axis.line.x = element_line(size = 0.5), axis.ticks.x = element_line(size = 0.5), axis.text.x = element_text(), plot.margin = unit(c(0, 4, 6, 0), "lines"))

        mymat     <- mybiom.r.stat.m[[r]][[i]]
        mybreaks <- c(10^(seq(log10(min(mymat, na.rm=T)), log10(0.05), length.out=4)), 1)

        p.ls.t[[r]][[i]] <- plotPvalues(mymat, border.col="black", pval.col="black", font.size=2) +
                    xlab(expression(paste(beta,"-dispersion test"))) + ylab("Permanova test") +
                    # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
                    scale_fill_gradient2(low = "red", mid = "yellow", high = "white", limits = c(NA, 0.05), guide = "colorbar", na.value = "transparent", trans = "log10", breaks = mybreaks, label = format(mybreaks, digit = 2, scientific = TRUE)) +
                    guides(fill = guide_colourbar(barwidth = 9, title.vjust = 0.75)) +
                    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

    }
}

# Flatten lists for plotting
p.ls.d %<>% purrr::flatten()
p.ls.t %<>% purrr::flatten()

# mylg <- length(p.ls.d)
mylg <- paste0(LETTERS[1:length(p.ls.d)], ". ",  unlist(mytitle))
pdf(NULL)
# p <- ggarrange(ggarrange(plotlist=p.ls.d, nrow=length(p.ls.d), ncol=1, common.legend=TRUE, legend="bottom", labels=LETTERS[1:mylg]),
p <- ggarrange(ggarrange(plotlist=p.ls.d, nrow=length(p.ls.d), ncol=1, common.legend=TRUE, legend="bottom", labels=mylg, hjust = 0),
               ggarrange(plotlist=p.ls.t, nrow=length(p.ls.t), ncol=1, common.legend=TRUE, legend="bottom"),
               ncol=2, nrow=1)
invisible(dev.off())

ggsave(paste0(graph.d, "Fig. 3 - b-div_FDR.pdf"), p, width = 5 * 2, height = 3.5 * length(p.ls.d), useDingbats = FALSE)

## vvv to delete
myax <- list(c(1,2), c(1,3), c(2,3))
p.ls.r <- vector("list", length(pop.ordr))

for (r in 1:length(pop.ordr)) {

    mybiom.tmp <- subset_samples(mybiom, Species == pop.ordr[[r]])
    mybiom.r <- rarefy_even_depth(mybiom.tmp, rngseed=myseed)

    p.ls.r[[r]] <- vector("list", length(mydist) * length(myax))

    k <- 0
    for (i in 1:length(mydist)) {
        for (j in 1:length(myax)) {
            mytitle <- ""
            if (j == 1) { mytitle <- tools::toTitleCase(mydist[i]) }

            k <- k + 1
            p.ls.r[[r]][[k]] <- plot_ordination(mybiom.r, mybiom.r.pcoa[[r]][[i]], color = "Tissue", axes=myax[[j]]) +
                           scale_color_manual(values=pop.clr) +
                           stat_ellipse(type = "norm") +
                           ggtitle(mytitle) +
                           theme(axis.text = element_text(size = rel(0.45))) +
                           theme(axis.title = element_text(size = rel(0.6))) +
                           theme(legend.position="none")

        }
    }
}

myletters <- lapply(1:length(mydist), function(x) c(LETTERS[x], rep("",length(myax)-1))) %>% unlist()

mylg <- length(p.ls.r[[1]])
pdf(NULL)
p <- ggarrange(plotlist=p.ls.r[[1]], nrow=length(p.ls.r[[1]])/3, ncol=3, common.legend=TRUE, legend="bottom", labels=myletters)
dev.off()

pdf(paste0(graph.d,"b-div.pdf"), width=3.5*3, height=3.5*length(p.ls.r[[1]])/3)
print(p)
dev.off()

# Panel with replicate
pdf(NULL)
p <- ggarrange(ggarrange(plotlist=p.ls.r[[1]], nrow=length(p.ls.r[[1]])/3, ncol=3, common.legend=TRUE, legend="bottom"),
               ggarrange(plotlist=p.ls.r[[2]], nrow=length(p.ls.r[[2]])/3, ncol=3, common.legend=TRUE, legend="bottom"),
               ncol=2, nrow=1, labels=LETTERS[1:2])
dev.off()

pdf(paste0(graph.d,"b-div_all.pdf"), width=3.5*3*2, height=3.5*length(p.ls.r[[1]])/3)
print(p)
dev.off()

#---------------------#
# Taxonomic diversity #
#---------------------#

# Data storage for graphs
p.ls   <- vector("list", length(pop.ordr))
pop.wd <- vector("list", length(pop.ordr))

# Color
mylevel <- "Phylum"
set.seed(myseed)

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
    myclr <- myclr.P[ tax_table(mybiom.P)[,2] ]

    # List to store graphs
    p <- vector("list", 3)

    # Tree
    # mylabels <- as.vector(dendro_data(mydd.tmp)$labels$label)
    mylabels <- labels(mydd_tmp)
    p[[1]] <- ggdendrogram(mydd_tmp) +
    # p[[1]] <-  as.ggdend(mydd_tmp) %>% ggplot(.) +
        scale_x_continuous(expand = c(1/(42*2),0), labels=mylabels, breaks=1:length(mylabels)) +
        theme(axis.text.y = element_blank(), plot.margin=margin(b=-10))

    # Barplot at sample level
    mybiom.P.m <- psmelt(mybiom.P.s)
    mybiom.P.m <- mybiom.P.m[ order(mybiom.P.m$Abundance, decreasing=TRUE), ]
    p[[2]] <- ggplot(mybiom.P.m, aes(x = Sample, y = Abundance, fill = Phylum)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = myclr) +
      scale_x_discrete(limits=mylabels) +
      # Remove x axis title and legend
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="none") +
      ylab("Relative Abundance\n")

    # Barplot at population level
    mybiom.P.m <- psmelt(mybiom.P.p)
    mybiom.P.m <- mybiom.P.m[ order(mybiom.P.m$Abundance, decreasing=TRUE), ]

    ## Number of samples per population
    spl.count   <- count(sample_data(mybiom.P)[,"Population"])
    pop.ls      <- as.vector(mybiom.P.m$Type) %>% sapply(., function(x) spl.count[ spl.count[,1] == x, 2])
    pop.wd[[r]] <- pop.ls

    ## Correction factor to make sure to have odd number under bars
    if (r == 1) myftr <- 2
    if (r == 2) myftr <- 2

    ## Label breaks
    mylabels.b <- rep("", nsamples(mybiom.tmp) / myftr)

    #
    mylabels.u  <- sapply(as.vector(mylabels), function(x) sample_data(mybiom.tmp)[x, "Type"]) %>% unlist() %>% as.vector() %>% unique()
    j <- 0
    for (l in mylabels.u) {
        i <- spl.count[ spl.count[,1] == l,2]/2 / myftr
        mylabels.b[ j+ceiling(i) ] <- l
        j <- j+i*2
    }

    pop.wd[[r]] <- pop.wd[[r]]/myftr

    p[[3]] <- ggplot(mybiom.P.m, aes(x = Type, y = Abundance, fill = Phylum, width=pop.wd[[!!r]]*0.98)) +   # Use of quasiquotation !!r otherwise lastest slot of the list taken for plotting
      geom_bar(stat = "identity") +
      scale_x_discrete(limits=mylabels.b, labels=mylabels.b) +
      scale_fill_manual(values = myclr) +
      # Remove x axis title and legend
      theme(axis.title.x = element_blank(), legend.position="none", axis.ticks = element_line(linetype=c("blank", rep("solid", length(mylabels.u))))) +     # axis.ticks needed to remove the first tick
      ylab("Relative Abundance\n")

    pdf(NULL)
    p.ls[[r]] <- ggarrange(plotlist=p[1:3], labels=LETTERS[1:length(p[[1:3]])], ncol = 1, common.legend = TRUE, legend = "bottom", align = "v")
    dev.off()
}

pdf(paste0(graph.d,pop.ordr[[1]]," taxo-div.pdf"), width=15, height=15, useDingbats=FALSE)
p.ls[[1]]
dev.off()

pdf(paste0(graph.d,pop.ordr[[2]]," taxo-div.pdf"), width=15, height=9, useDingbats=FALSE)
p.ls[[2]]
dev.off()

# Per Organ
# for (i in seq(151, 250)) {
i <- 150
p.ls   <- vector("list", length(pop.ordr))
set.seed(myseed)
mybiom.P <- mytax_glom(mybiom, mylevel) %>% core(., detection=0, prevalence=0)
myclr.P <- tax_table(mybiom.P) %>% nrow() %>% distinctColorPalette(., seed=myseed - i)
names(myclr.P) <- tax_table(mybiom.P)[, mylevel] %>% sort()

for (r in 1:length(pop.ordr)) {
    mybiom.tmp <- subset_samples(mybiom, Population == pop.ordr[[r]])

    # Remove non present taxa
    mybiom.P <- mytax_glom(mybiom.tmp, "Phylum")
    mybiom.P <- core(mybiom.P, detection=1, prevalence=0)
    # mybiom.P <- core(mybiom.P, detection=200, prevalence=0.1)

    # # Sample biome
    # mybiom.P.s <- transform(mybiom.P, transform="compositional")

    # Organ biome
    mybiom.P.p <- mysample_glom(mybiom.P, "Type")

    # Transform data to compositional
    mybiom.P.p <- transform(mybiom.P.p, transform="compositional")

    # Vector color
    myclr <- myclr.P[ tax_table(mybiom.P.p)[,2] ]

    # List to store graphs
    p <- vector("list", 3)

    # Barplot at population level
    mybiom.P.m <- psmelt(mybiom.P.p)
    mybiom.P.m <- mybiom.P.m[ order(mybiom.P.m$Abundance, decreasing=TRUE), ]

    # Vector color
    myclr <- myclr.P[ tax_table(mybiom.P)[,2] ]

    p.ls[[r]] <- ggplot(mybiom.P.m, aes(x = Type.code, y = Abundance, fill = Phylum,)) +   # Use of quasiquotation !!r otherwise lastest slot of the list taken for plotting
        geom_bar(stat = "identity", color = "black", size = 0.05) +
        # scale_x_discrete(limits=mylabels.b, labels=mylabels.b) +
        scale_fill_manual(values = myclr) +
        ## Update legend look (because of the bar borders)
        guides(fill = guide_legend(override.aes = list(colour = NA, size = 0.5))) +
        # Remove x axis title and legend
        # theme(axis.title.x = element_blank(), legend.position="none", axis.ticks = element_line(linetype=c("blank", rep("solid", length(mylabels.u))))) +     # axis.ticks needed to remove the first tick
        labs(title = pop.ordr[[r]], y = "Relative Abundance", x = "") +
        theme(legend.position = "none")
}

# pdf(paste0(graph.d, pop.ordr[[1]], " organ taxo-div.pdf"), width = 8, height = 7, useDingbats = FALSE)
# ggarrange(plotlist=p.ls, labels=LETTERS[1:length(p.ls)], ncol = 1, common.legend = TRUE, legend = "bottom", align = "v")
# dev.off()

pdf(NULL)
p <- ggarrange(plotlist=p.ls, labels=LETTERS[1:length(p.ls)], ncol = 1, common.legend = TRUE, legend = "bottom", align = "v")
invisible(dev.off())

ggsave(paste0(graph.d, "Fig. 4 - Organ taxo-div ", i, ".pdf"), p, width = 8, height = 7, useDingbats = FALSE)
# }

# pdf(paste0(graph.d,pop.ordr[[2]]," organ taxo-div.pdf"), width=8, height=5, useDingbats=FALSE)
# p.ls[[2]][[3]]
# dev.off()

#-----------------------------#
# Shared taxa between samples #
#-----------------------------#

# pop.ordr <- md[,"Species"] %>% levels
rm.water <- FALSE

# # Test several thresholding (minimal count in sample type, minimal prevalance))
# tsh.ls <- list(c(0, 0), c(100, 0.2))
tsh.ls <- c(0, 4)

ct.ls <- vector("list", length(tsh.ls))
low.ls <- vector("list", length(tsh.ls))

for (i in 1:length(tsh.ls)) {

    ct.ls[[i]]  <- vector("list", length(pop.ordr)) %>% setNames(., pop.ordr)
    low.ls[[i]] <- vector("list", length(pop.ordr)) %>% setNames(., pop.ordr)

    # Minimal count in sample type
    # mytsh <- tsh.ls[[i]][[1]]
    mytsh <- 0

    # Invert prevalence of the select minimal count [0-1] (0 = all, 1 = none)
    # mypv <- tsh.ls[[i]][[2]]
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

        # myclr <- pop.clr[ ! grepl("T.*", org.data.tmp[,1]) ]

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

for (i in 1:length(tsh.ls)) {

    # Threshold values
    # mytsh <- tsh.ls[[i]][[1]]
    # mypv  <- tsh.ls[[i]][[2]]
    mytsh <- 0
    mypv  <- 0

    myh <- 5
    myw <- 10
    if (i == 1) myw <- 20

    for (s in pop.ordr) {
        pdf(paste0(graph.d, "Fig.5 - ", i, s," ASV upset test - tsh ", mytsh, " pv ", mypv, ".pdf"), height = myh, width = myw, useDingbats = FALSE)
        otu_type_fn <- ct.ls[[i]][[s]]
        clr <- pop.data[ pop.data[, 1] == s, 2]
        ups <- upset(otu_type_fn, sets = colnames(otu_type_fn), nintersects = NA, main.bar.color = clr, mainbar.y.label = "ASVs per intersection", sets.x.label = "ASVs per sample type", keep.order = TRUE)
        print(ups)
        grid.text(s, x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
        dev.off()
    }
}

# Ubiquitous ASVs (snail + water) shared between species
a <- which((ct.ls[[1]][[1]] %>% rowSums(.)) == 8) %>% names()
b <- which((ct.ls[[1]][[2]] %>% rowSums(.)) == 8) %>% names()

a %in% b

# Poportion of ASV shared per species
lapply(ct.ls[[1]], function(x) nrow(x[rowSums(x) > 1, ]) / nrow(x))

# Poportion of ASV shared per sample type
lapply(ct.ls[[1]], function(a) sapply(1:8, function(x) { a1 <- a[ a[,x] == 1, ] ; nrow(a1[ rowSums(a1) > 1, ]) / nrow(a1)}) %>% setNames(., colnames(a))) %>% melt()

# Proportion of shared ASVs between H and O|W
# lapply(ct.ls[[1]], function(x) { a <- (x[, "H"]  == 1 & (x[, "W"]  == 1 | x[, "H"]  == 1 & x[, "O"]  == 1)) ; sum(a)/length(a[rowSums(x) > 1]) })
spl <- "H"
lapply(org.data[org.data[,1] != spl, 1], function(y) sapply(ct.ls[[1]], function(x) { a <- (x[, spl]  == 1 & x[, y]  == 1) ; sum(a)/length(a[rowSums(x) > 1]) })) %>% setNames(., paste(spl, "and", org.data[org.data[,1] != spl, 1]))

# Proportion of shared ASVs between whole snail and another sample type
spl <- "W"
lapply(org.data[org.data[,1] != spl, 1], function(y) sapply(ct.ls[[1]], function(x) { a <- (x[, spl]  == 1 & x[, y]  == 1) & rowSums(x[, setdiff(names(x), c(spl, y))])  == 0 ; sum(a)/length(a[rowSums(x) > 1]) })) %>% setNames(., paste(spl, "and", org.data[org.data[,1] != spl, 1]))

# Proportion of shared ASVs between TK|TY and H
# lapply(ct.ls[[1]], function(x) { a <- (x[, "TK"]  == 1 | x[, "TY"]  == 1) & (x[, "H"]  == 1 | x[, "O"]  == 1 | x[, "W"]  == 1) ; sum(a)/length(a[rowSums(x) > 1]) })
lapply(ct.ls[[1]], function(x) { a <- (x[, "TK"]  == 1 | x[, "TY"]  == 1) & (x[, "H"]  == 1) ; sum(a)/length(a[rowSums(x) > 1]) })

# Proportion of shared ASVs between TK|TY and G|S
# lapply(ct.ls[[1]], function(x) { a <- (x[, "TK"]  == 1 | x[, "TY"]  == 1) & (x[, "G"]  == 1 | x[, "S"] == 1) ; sum(a)/length(a[rowSums(x) > 1]) })

lapply(org.data[, 1], function(y) lapply(ct.ls[[1]], function(x) { a <- (x[, "TK"]  == 1 | x[, "TY"]  == 1) ; b <- rowSums(x[! a & x[, y]  == 1,]) >= 1 ; sum(b) / sum(!a) })) %>% setNames(., org.data[, 1])

# # Proportion of ASV found only in snails (amoung all ASVs)
# a <- ct.ls[[1]][[2]]
# sum(! rowSums(a[, grepl("^T", colnames(a))]) > 0) / nrow(a)
lapply(ct.ls[[1]], function(x) sum(! rowSums(x[, grepl("^T", colnames(x))]) > 0) / nrow(x))

# Proportion of shared ASV found only in snails (amoung all shared ASVs)
lapply(ct.ls[[1]], function(x) {x1 <- x[rowSums(x) > 1, ] ; sum(! rowSums(x1[, grepl("^T", colnames(x1))]) > 0) / nrow(x1) })

## Proportion of ubiquitous shared ASVs across snail sample types
# aT <- a[, ! grepl("^T", colnames(a))]
lapply(ct.ls[[1]], function(x) {x <- x[rowSums(x) > 1, ]
                                xT <- x[! rowSums(x[, grepl("^T", colnames(x))]) > 0,  ! grepl("^T", colnames(x))]
                                sum(rowSums(xT) == ncol(xT)) / nrow(xT) })

## Proportion of ubiquitous shared ASVs across all sample types
lapply(ct.ls[[1]], function(x) {sum(rowSums(x) == ncol(x)) / nrow(x[rowSums(x) > 1,]) })

# Proportion of shared ASV found only in snails (amoung all shared ASVs)
lapply(ct.ls[[1]], function(x) {x1 <- x[rowSums(x) > 1, ] ; sum(rowSums(x1) > 4) / nrow(x1) })

# ## Proportion of shared ASVs low sharing
# for (i in 1:length(ct.ls)) {
#     a <- ct.ls[[1]][[i]]
#     b <- ct.ls[[2]][[i]]

#     a <- a[rowSums(a) > 1,]
#     b <- b[rowSums(b) > 1,]

#     print(1 - nrow(b) / nrow(a))
# }

## Proportion of intersection with low sharing
for (i in 2:length(tsh.ls)) {
    for (s in pop.ordr) {
        mylow <- low.ls[[i]][[s]]
        print(1 - sum(mylow) / length(mylow))
    }
}

## Proportion of sample types associated with intersection with low ASVs
lapply(low.ls[[2]], function(x) {names(x)[x] %>% sapply(., strsplit, "") %>% do.call("rbind", .) %>% apply(., 2, as.numeric) %>% colSums() %>% setNames(., org.data[,1])} )

# pairwise sharing
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

a <- pw.ls[[1]]
a[lower.tri(a)] <- (pw.ls[[2]] %>% t())[lower.tri(a)]

# P-value matrix for each relevant index
# p <- vector("list", nrow(myidx))
# mybreaks <- c(10^(seq(log10(min(mymat, na.rm=T)), log10(0.05), length.out=4)), 1)
pdf(NULL)
p <- plotPvalues(a, border.col="black", pval.col="black", font.size=3, scientific = FALSE) +
    xlab("") + ylab("") +
    # scale_fill_gradient2(low = "red", mid = "yellow", high = "white", limits = c(NA, 0.05), guide = "colorbar", na.value = "transparent", trans = "log10", breaks = mybreaks, label = format(mybreaks, digit = 2, scientific = TRUE)) +
    scale_fill_gradient2(low = "white", mid = "yellow", high = "red", limits = c(0.2, 0.6), guide = "colorbar", na.value = "transparent") +
    guides(fill = guide_colourbar(barwidth = 9, title.vjust = 0.75)) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
    # ggtitle(myidx[i, 2])
invisible(dev.off())
ggsave(paste0(graph.d, "Supp. Fig. 4 - pairwise shared ASV.pdf"), p, width = 10, height = 10, useDingbats = FALSE)

#~~~~~~~~~~
# By organ

rm.water <- FALSE

org.ordr <- org.data[, 1]

ct.ls <- vector("list", length(org.ordr)) %>% setNames(., org.ordr)

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
    ct.ls[[s]] <- otu_type %>% as.data.frame()
}

p.ls <- vector("list", length(org.ordr)) %>% setNames(., org.ordr)

myw <- 2.75 * length(org.ordr) / 2
myh <- 2.75 * length(org.ordr) / 4

myclr <- pop.data[,2] %>% c(., midcol(.[1], .[2]))

for (s in org.ordr) {
    otu_type_fn <- ct.ls[[s]]

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
pdf(paste0(graph.d, "1 ASV upset by spl_tp - test.pdf"), width = myw, height = myh)
p <- ggarrange(plotlist = p.ls, nrow = length(p.ls) / 4, ncol = length(p.ls) / 2) #, common.legend=TRUE, legend="bottom", labels=LETTERS[1:2])
print(p)
dev.off()

#=========#
# Figures #
#=========#

if (! dir.exists(graph.dir) ) { dir.create(graph.dir, recursive=TRUE) }

