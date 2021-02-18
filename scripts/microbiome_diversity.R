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
#    library("magrittr")
    library("plyr")             # function count
    library("tidyr")

    # Microbiome
    library("microbiome")       # function core
    library("qiime2R")
    library("phyloseq")
    library("picante")          # PD function
    library("pairwiseAdonis")   # see github
    library("phylogram")        # prune function
    library("pvclust")
    library("dendextend")        # pvclust_edges function

    # Graph
    library("colorspace")       # function rainbow_hcl
    library("cowplot")
    library("ggpubr")
    library("ggdendro")
    library("grid")             # function grid.rect
})



#===========#
# Functions #
#===========#

# Working directory
#setwd(file.path(getwd(), "scripts"))

cat("Loading functions...\n")
source("microbiome_diversity_func.R")



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

asv.f  <- "../results/1-qiime/table.qza"
taxa.f <- "../results/1-qiime/rep-seqs_taxa.qza"
tree.f <- "../results/1-qiime/rooted-tree.qza"

md.f   <- "../data/sample-metadata.tsv"

# Population order and color
pop.data <- matrix(c("Ba",   "#ffac40",
                     "Bg90", "#87d687"), byrow=TRUE, ncol=2) # Change to BgBS90

# Tissue order
ts.data <- matrix(c(
		"H", "Hemolymph",
		"S", "Stomach",
		"G", "Gut",
		"L", "Hepatopancreas",
		"O", "Ovotestis",
        "W", "Whole snail",
		"TY", "Water tray",
		"TK", "Water tank"), ncol=2, byrow=TRUE)


# Graph output directory
graph.d <- "../graphs/"

# GGplot options
theme_set(theme_classic())

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
# vvvv To be removed when md file updated
colnames(md)[,1] <- "Population"
# ^^^^
md[,"Population"] <- factor(md[,"Population"], pop.data[,1])
md[,"Tissue"] <- factor(md[,"Tissue"], ts.data[,1])

# Building phyloseq object
mybiom <- phyloseq(asv.tb, taxa.tb, tree, sample_data(md))

# Cleaning data from contaminants
my.mt  <- apply(tax_table(mybiom), 1, function(x) any(grepl("mitochondria", x, ignore.case=TRUE)))
my.cp  <- apply(tax_table(mybiom), 1, function(x) any(grepl("chloroplast", x, ignore.case=TRUE)))
mybiom <- subset_taxa(mybiom, ! (my.mt | my.cp))


# Determine populations and their order
##mypop     <- cbind(rownames(md), md[,c(2,4)])
#mypop     <- cbind(rownames(md), md[,2:4])
#mypop[,2] <- factor(mypop[,2], pop.ordr) 
mypop <- md



#mybiom.rep <- mybiom
rep    <- unclass(sample_data(mybiom)[,"Replicate"])[[1]]
nb.rep <- unique(rep)
min.seq <- rep(NA, length(nb.rep))
#for (i in 1:length(nb.rep)) {
#    min.seq[i] <- min(sample_sums(mybiom)[rep == nb.rep[i]])
#}
#myrows <- rep == nb.rep[which.max(min.seq)]
#mybiom <- subset_samples(mybiom, myrows)


#--------------------#
# Rarefaction curves #
#--------------------#
 
p <- ggrare(mybiom, step = 1000, color = "Population", se = FALSE, plot=FALSE) +
        scale_color_manual(values = pop.data[,2]) +
        theme(legend.position="none") +
        facet_wrap(~Population)

ggsave(paste0(graph.d,"Fig. 2 - rarefaction.pdf"), p, height=7, width=7, useDingbats=FALSE)
#p <- ggrare(mybiom, step = 1000, color = "Population", label = "Sample", se = FALSE, plot=FALSE)
p <- ggrare(mybiom, step = 1000, color = "Species", se = FALSE, plot=FALSE)
p <- p + facet_wrap(~Tissue)

pdf(NULL)
p <- shift_legend(p)
dev.off()

pdf(paste0(graph.d,"rarefaction.pdf"), height=7, width=14)
plot(p)
dev.off()


#-----------#
# Diversity #
#-----------#

## TO DO ##
# Remove duplicated samples

myalpha <- microbiome::alpha(mybiom)
mypd    <- pd(t(asv), tree)
myalpha <- merge(myalpha, mypd, by="row.names")
myalpha <- merge(myalpha, mypop, by.x=1, by.y="row.names")

#myidx <- c("observed", "chao1", "diversity_shannon", "evenness_simpson", "dominance_simpson", "PD")
myidx <- matrix(c("observed",           "Observed richness",
                  "chao1",              "Chao1 richness",
                  "diversity_shannon",  "Shannon diversity",
                  "PD",                 "Faith's Phylogenetic diversity",
                  "evenness_simpson",   "Simpson evenness",
                  "dominance_simpson",  "Simpson dominance"), ncol=2, byrow=TRUE)





# Statistical tests
#mya.test    <- vector("list", length(nb.rep))
#mya.letters <- vector("list", length(nb.rep))

#for (r in nb.rep) { 
#    mya.test[[r]]   <- vector("list", length(myidx[,1]))
#    for (i in myidx[,1]) {
#        mya.test[[r]][[ match(i,myidx[,1]) ]] <- pairwise.wilcox.test(myalpha[myalpha["Replicate"]==r,i], myalpha[myalpha["Replicate"]==r,"Population"], p.adjust.method="none")
#    }
#
#    mya.letters[[r]] <- vector("list", length(myidx[,1]))
#    for (i in 1:length(mya.test[[r]])) {
#        d <- as.vector(mya.test[[r]][[i]]$p.value)
#        names(d) <- sapply(colnames(mya.test[[r]][[i]]$p.value), function(x) paste0(x, "-", rownames(mya.test[[r]][[i]]$p.value))) %>% as.vector()
#        d <- d[!is.na(d)]
#        mya.letters[[r]][[i]] <- generate.label.df(d)
#    }
#}
#


## Graph
p <- vector("list", nrow(myidx))
for (i in myidx[,1]) {
    df.tmp <- data.frame(i=myalpha[,i], Tissue=myalpha[,"Tissue"], Species=myalpha[,"Species"])
    p[[match(i,myidx)]] <- ggplot(df.tmp, aes(x=Tissue, y=i, fill=Species)) +
#            geom_blank(data=df.tmp, aes(x=Tissue, y=i*1.1)) +
            geom_boxplot() + xlab("") + ylab("") + ggtitle(myidx[myidx[,1] == i,2]) +
#                stat_summary(geom = 'text', label = mya.letters[[r]][[match(i,myidx[,1])]][,1], fun.y = max, vjust = -1) +
#            stat_summary(aes(group=Species), geom="text", fun.y = max, vjust = -1) +
#            stat_summary(aes(group=Species)) +
            scale_fill_manual(values=pop.clr) +
            theme(plot.title = element_text(hjust =0.5), legend.position="none")

}
pdf(paste0(graph.d,"alpha-div.pdf"), width=10, height=10)
ggarrange(plotlist=p, labels=LETTERS[1:length(p)], ncol=2, nrow=ceiling(length(p)/2), common.legend=TRUE, legend="bottom", align="hv")
dev.off()



### Rarity
#myidx <- matrix(c("rarity_low_abundance", "Low abundance"), ncol=2, byrow=TRUE)
#
#myr.test    <- vector("list", length(nb.rep))
#myr.letters <- vector("list", length(nb.rep))
#
#for (r in nb.rep) { 
#    myr.test[[r]]   <- vector("list", length(myidx[,1]))
#    for (i in myidx[,1]) {
#        myr.test[[r]][[ match(i,myidx[,1]) ]] <- pairwise.wilcox.test(myalpha[myalpha["Replicate"]==r,i], myalpha[myalpha["Replicate"]==r,"Population"], p.adjust.method="none")
#    }
#
#    myr.letters[[r]] <- vector("list", length(myidx[,1]))
#    for (i in 1:length(myr.test[[r]])) {
#        d <- as.vector(myr.test[[r]][[i]]$p.value)
#        names(d) <- sapply(colnames(myr.test[[r]][[i]]$p.value), function(x) paste0(x, "-", rownames(myr.test[[r]][[i]]$p.value))) %>% as.vector()
#        d <- d[!is.na(d)]
#        myr.letters[[r]][[i]] <- generate.label.df(d)
#    }
#}
#
#
#
#p <- vector("list", length(nb.rep))
#for (r in nb.rep) {
#    p[[r]] <- vector("list", length(myidx[,1]))
#    for (i in myidx[,1]) {
#        df.tmp <- data.frame(i=myalpha[myalpha["Replicate"]==r,i], Population=myalpha[myalpha["Replicate"]==r,"Population"])
#
#        p[[r]][[match(i,myidx)]] <- ggplot(df.tmp, aes(x=Population, y=i, fill=Population)) +
#                geom_blank(data=df.tmp, aes(x=Population, y=i*1.1)) +
#                geom_boxplot() + xlab("") + ylab("") + #gtitle(myidx[myidx[,1] == i,2]) +
#                stat_summary(geom = 'text', label = myr.letters[[r]][[match(i,myidx[,1])]][,1], fun.y = max, vjust = -1) +
#                scale_fill_manual(values=pop.clr) +
#                theme(plot.title = element_text(hjust =0.5), legend.position="none")
#
#    }
#}
#
#myletters <- NULL
#if (nrow(myidx) > 1) { myletters <- LETTERS[1:nrow(myidx)] }
#
#pdf(paste0(graph.d,"rarity.pdf"), width=10, height=10)
#ggarrange(plotlist=p[[1]], labels=myletters, ncol=length(p[[1]]), nrow=ceiling(length(p[[r]])/2))
#dev.off()
#
#pdf(paste0(graph.d,"rarity_all.pdf"), width=10*length(nb.rep), height=10)
#ggarrange( ggarrange(plotlist=p[[1]], labels=LETTERS[1], ncol=length(p[[r]]), nrow=ceiling(length(p[[r]])/2)),
#           ggarrange(plotlist=p[[2]], labels=LETTERS[2], ncol=length(p[[r]]), nrow=ceiling(length(p[[r]])/2)),
#           ncol=length(nb.rep), nrow=1)
#dev.off()



#------#
# PCoA #
#------#

# For PCoA
mydist <- c("jaccard", "Unweighted UniFrac", "Weighted UniFrac", "bray")

#mybiom.r.pcoa   <- vector("list", length(mydist))
#mybiom.r.stat   <- vector("list", length(mydist))
#mybiom.r.stat.p <- vector("list", length(mydist))
#mybiom.r.stat.b <- vector("list", length(mydist))
#mybiom.r.stat.T <- vector("list", length(mydist))
#
#mybiom.tmp <- mybiom
#mybiom.r <- rarefy_even_depth(mybiom.tmp, rngseed=myseed)
#
#
#for (i in mydist) {
#    myls.idx <- match(i, mydist)
#
#    if (i == "jaccard") {
#        mydist.tmp <- distance(mybiom.r, i, binary=TRUE)
#    } else {
#        mydist.tmp <- distance(mybiom.r, i)
#    }
#    
#    mybiom.r.pcoa[[myls.idx]]   <- ordinate(mybiom.r, "PCoA", mydist.tmp)
#
##    # Permanova stat
##    mybiom.r.stat[[myls.idx]]   <- adonis(mydist.tmp ~ Population, data = get_variable(mybiom.r), permutations = 1000)
##    mybiom.r.stat.p[[myls.idx]] <- pairwise.adonis(mydist.tmp, get_variable(mybiom.r)[,"Population"], perm = 1000)
##
##    # Dispersion stat
##    mybiom.r.stat.b[[myls.idx]] <- betadisper(mydist.tmp, get_variable(mybiom.r)[,"Population"])
##    mybiom.r.stat.T[[myls.idx]] <- TukeyHSD(mybiom.r.stat.b[[myls.idx]])
#}
#
#names(mybiom.r.pcoa)   <- mydist
#names(mybiom.r.stat)   <- mydist
#names(mybiom.r.stat.p) <- mydist
#names(mybiom.r.stat.b) <- mydist
#names(mybiom.r.stat.T) <- mydist
#
#
#
#myax <- list(c(1,2), c(1,3), c(2,3))
#
#
#p.ls.r <- vector("list", length(mydist) * length(myax))
#
#k <- 0
#for (i in 1:length(mydist)) {
#    for (j in 1:length(myax)) {
#        mytitle <- ""
#        if (j == 1) { mytitle <- tools::toTitleCase(mydist[i]) }
#
#        k <- k + 1
#        p.ls.r[[k]] <- plot_ordination(mybiom.r, mybiom.r.pcoa[[i]], color = "Population", axes=myax[[j]]) +
#                       scale_color_manual(values=pop.clr) +
#                       stat_ellipse(type = "norm") +
#                       ggtitle(mytitle) +
#                       theme(axis.text = element_text(size = rel(0.45))) +
#                       theme(axis.title = element_text(size = rel(0.6))) +
#                       theme(legend.position="none")
#
#    }
#}
#
#myletters <- lapply(1:length(mydist), function(x) c(LETTERS[x], rep("",length(myax)-1))) %>% unlist()
#
#mylg <- length(p.ls.r[[1]])
#pdf(NULL)
#p <- ggarrange(plotlist=p.ls.r[[1]], nrow=length(p.ls.r[[1]])/3, ncol=3, common.legend=TRUE, legend="bottom", labels=myletters)
#dev.off()
#
#
#pdf(paste0(graph.d,"b-div.pdf"), width=3.5*3, height=3.5*length(p.ls.r[[1]])/3)
#print(p)
#dev.off()
#
#
## Panel with replicate
#pdf(NULL)
#p <- ggarrange(ggarrange(plotlist=p.ls.r[[1]], nrow=length(p.ls.r[[1]])/3, ncol=3, common.legend=TRUE, legend="bottom"),
#               ggarrange(plotlist=p.ls.r[[2]], nrow=length(p.ls.r[[2]])/3, ncol=3, common.legend=TRUE, legend="bottom"),
#               ncol=2, nrow=1, labels=LETTERS[1:2])
#dev.off()
#
#
#pdf(paste0(graph.d,"b-div_all.pdf"), width=3.5*3*2, height=3.5*length(p.ls.r[[1]])/3)
#print(p)
#dev.off()

mybiom.r.pcoa   <- vector("list", length(pop.ordr))
mybiom.r.stat   <- vector("list", length(pop.ordr))
mybiom.r.stat.p <- vector("list", length(pop.ordr))
mybiom.r.stat.b <- vector("list", length(pop.ordr))
mybiom.r.stat.T <- vector("list", length(pop.ordr))

for (r in 1:length(pop.ordr)) {
    mybiom.tmp <- subset_samples(mybiom, Species == pop.ordr[[r]])
    mybiom.r <- rarefy_even_depth(mybiom.tmp, rngseed=myseed)

    mybiom.r.pcoa[[r]]   <- vector("list", length(mydist))
    mybiom.r.stat[[r]]   <- vector("list", length(mydist))
    mybiom.r.stat.p[[r]] <- vector("list", length(mydist))
    mybiom.r.stat.b[[r]] <- vector("list", length(mydist))
    mybiom.r.stat.T[[r]] <- vector("list", length(mydist))

    for (i in mydist) {
        myls.idx <- match(i, mydist)

        if (i == "jaccard") {
            mydist.tmp <- distance(mybiom.r, i, binary=TRUE)
        } else {
            mydist.tmp <- distance(mybiom.r, i)
        }
        
        mybiom.r.pcoa[[r]][[myls.idx]]   <- ordinate(mybiom.r, "PCoA", mydist.tmp)

#        # Permanova stat
#        mybiom.r.stat[[r]][[myls.idx]]   <- adonis(mydist.tmp ~ Population, data = get_variable(mybiom.r), permutations = 1000)
#        mybiom.r.stat.p[[r]][[myls.idx]] <- pairwise.adonis(mydist.tmp, get_variable(mybiom.r)[,"Tissue"], perm = 1000)
#
#        # Dispersion stat
#        mybiom.r.stat.b[[r]][[myls.idx]] <- betadisper(mydist.tmp, get_variable(mybiom.r)[,"Tissue"])
#        mybiom.r.stat.T[[r]][[myls.idx]] <- TukeyHSD(mybiom.r.stat.b[[r]][[myls.idx]])
    }

    names(mybiom.r.pcoa[[r]])   <- mydist
    names(mybiom.r.stat[[r]])   <- mydist
    names(mybiom.r.stat.p[[r]]) <- mydist
    names(mybiom.r.stat.b[[r]]) <- mydist
    names(mybiom.r.stat.T[[r]]) <- mydist

}


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
set.seed(myseed)
mybiom.P <- mytax_glom(mybiom, "Phylum") %>% core(., detection=0, prevalence=0)
myclr.P  <- tax_table(mybiom.P) %>% nrow() %>% rainbow_hcl() %>% sample()
names(myclr.P) <- tax_table(mybiom.P)[,2]

for (r in 1:length(pop.ordr)) {

#    if (r == 1 ) {
#        mybiom.tmp <- subset_samples(mybiom, Replicate == r)
#        mydd.tmp <-  prune(mydd, pattern = ".2$")
#    } else {
#        mybiom.tmp <- mybiom
#        mydd.tmp   <- mydd
#    }

    mybiom.tmp <- subset_samples(mybiom, Species == pop.ordr[[r]])

    # Distance matrix
    mybiom.d <- apply(otu_table(mybiom.tmp), 1, function(x){ x <- x+1; log(x) - mean(log(x))})
    mydm     <- dist(mybiom.d, method="euclidian")

    # Cluster the data
    myhc <- hclust(mydm, method="ward.D2")

    # Generate the tree
    mydd <- as.dendrogram(myhc)

    ## Reroot the tree with water samples
    #myroot  <- prune(mydd, pattern = "Water", keep = TRUE)
    #subtree <- prune(mydd, pattern = "Water")
    #mydd    <- list(subtree, myroot) %>% remidpoint() %>% as.cladogram()

#    mydd.tmp <-  prune(mydd, pattern = pop.ordr[[r]], keep=TRUE)
    mydd.tmp <-  mydd

    # Remove non present taxa
    mybiom.P <- mytax_glom(mybiom.tmp, "Phylum")
    mybiom.P <- core(mybiom.P, detection=0, prevalence=0)

    # Sample biome
    mybiom.P.s <- transform(mybiom.P, transform="compositional")

    # Population biome
    mybiom.P.p <- mysample_glom(mybiom.P, "Tissue")
    mybiom.P.p <- transform(mybiom.P.p, transform="compositional")
    
    # Vector color
    myclr <- myclr.P[ tax_table(mybiom.P)[,2] ]

    # List to store graphs
    p <- vector("list", 3)

    # Tree
    mylabels <- as.vector(dendro_data(mydd.tmp)$labels$label)
    p[[1]] <- ggdendrogram(mydd.tmp) +
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
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position="none") + 
      ylab("Relative Abundance\n")


#    # Barplot at population level
#    mybiom.P.m <- psmelt(mybiom.P.p)
#    mybiom.P.m <- mybiom.P.m[ order(mybiom.P.m$Abundance, decreasing=TRUE), ]
#
#    ## Number of samples per population
#    spl.count   <- count(sample_data(mybiom.P)[,2])
#    pop.ls      <- as.vector(mybiom.P.m$Tissue) %>% sapply(., function(x) spl.count[ spl.count[,1] == x, 2])
#    pop.wd[[r]] <- pop.ls
#
#    ## Correction factor to make sure to have odd number under bars
#    if (r == 1) myftr <- 2
#    if (r == 2) myftr <- 2
#
#    ## Label breaks
#    mylabels.b <- rep("", nsamples(mybiom.tmp) / myftr)
#
#    #
#    mylabels.u  <- sapply(as.vector(mylabels), function(x) sample_data(mybiom.tmp)[x, "Tissue"]) %>% unlist() %>% as.vector() %>% unique()
#    j <- 0
#    for (l in mylabels.u) {
#        i <- spl.count[ spl.count[,1] == l,2]/2 / myftr
#        mylabels.b[ j+ceiling(i) ] <- l
#        j <- j+i*2
#    }
#
#    pop.wd[[r]] <- pop.wd[[r]]/myftr
#
#
#    p[[3]] <- ggplot(mybiom.P.m, aes(x = Population, y = Abundance, fill = Phylum, width=pop.wd[[!!r]]*0.98)) +   # Use of quasiquotation !!r otherwise lastest slot of the list taken for plotting
#      geom_bar(stat = "identity") +
#      scale_x_discrete(limits=mylabels.b, labels=mylabels.b) +
#      scale_fill_manual(values = myclr) +
#      # Remove x axis title and legend
#      theme(axis.title.x = element_blank(), legend.position="none", axis.ticks = element_line(linetype=c("blank", rep("solid", length(mylabels.u))))) +     # axis.ticks needed to remove the first tick
#      ylab("Relative Abundance\n")

    pdf(NULL)
    p.ls[[r]] <- ggarrange(plotlist=p[1:2], labels=LETTERS[1:length(p[[1:2]])], ncol=1, common.legend=T, legend="bottom", align="v")
    dev.off()
}


pdf(paste0(graph.d,pop.ordr[[1]]," taxo-div.pdf"), width=15, height=9, useDingbats=FALSE)
p.ls[[1]]
dev.off()


pdf(paste0(graph.d,pop.ordr[[2]]," taxo-div.pdf"), width=15, height=9, useDingbats=FALSE)
p.ls[[2]]
dev.off()


#-----------------------------#
# Common and specific species #
#-----------------------------#

library("venn")

ct.ls <- vector("list", length(pop.ordr))
p.ls  <- vector("list", length(pop.ordr))


for (s in pop.ordr) {
    myidx <- match(s, pop.ordr)

    md.tmp <- md[ md[,"Species"] == s, ]
    myotu.tmp <- otu_table(mybiom)[,rownames(md.tmp)]
    # Remove irrelevant ASV
    myotu.tmp <- myotu.tmp[ rowSums(myotu.tmp) > 0, ]

    myts.tmp <- ts.ordr[ ! grepl("T.*", ts.ordr[,1]) ,]

    myclr <- pop.clr[ ! grepl("T.*", ts.ordr[,1]) ]

    # Threshold to select minimal count in tissues
    mytsh <- 0

    # Invert prevalence of the select minimal count [0-1] (0 = all, 1 = none)
    mypv <- 0

    a <- sapply(myts.tmp[,1], function(x) myotu.tmp[,rownames(md.tmp[ md.tmp[,"Tissue"] == x,])] )
#    b <- lapply(a, function(x) rowSums(x) > 0)
    b <- lapply(a, function(x) rowSums(x > mytsh) > ncol(x)*mypv )
    c <- data.frame(sapply(b,c))

    # Remove irrelevant ASV (after tissue filtration)
    c <- c[ rowSums(c) > 0, ]

    # Rename colums (Not a good idea if image too small)
#    for (i in 1:ncol(c)) { colnames(c)[i] <- ts.ordr[ ts.ordr[,1] == colnames(c)[i],2]}

    ct.ls[[myidx]] <- nrow(c)


#    # Common species between tissues and whole snail
#    c.com <- c[ c[,"W"], ]
#
#    # Specific species to tissues
#    c.spec <- c[ ! c[,"W"], ]
#
##    c.com.ct  <- signif( colSums(c.com) / nrow(c.com) , 2)
##    c.spec.ct <- signif( colSums(c.spec) / nrow(c.spec) , 2)
##    c.com.ct  <- signif( c(colSums(c.com), Total=nrow(c.com)) , 2)
##    c.spec.ct <- signif( c(colSums(c.spec), Total=nrow(c.spec)) , 2)
##    ct.ls[[myidx]] <- list(c.com.ct, c.spec.ct)
#    ct.ls[[myidx]] <- data.frame( Common=signif( c(colSums(c.com), Total=nrow(c.com)) , 2),
#                                  Specific=signif( c(colSums(c.spec), Total=nrow(c.spec)) , 2),
#                                  Tissue=c(colnames(c.com), "Total"))
#
#    ct.ls[[myidx]][,3] <- factor(ct.ls[[myidx]][,3], c(myts.tmp[,1], "Total"))
#
#    p[[1]] <- ggplot(ct.ls[[myidx]], aes(x=Tissue, y=Specific)) +
#            geom_bar(stat="identity") + ylab("Number of ASV observed")
#    
#    p[[2]] <- ggplot(ct.ls[[myidx]], aes(x=Tissue, y=Common)) +
#            geom_bar(stat="identity") + ylab("Number of ASV observed")
#
#    pdf(paste0(graph.d,s," ASV specific whole snail.pdf"), useDingbats=FALSE)
#    print(p[[1]])
#    dev.off()
#    
#    pdf(paste0(graph.d,s," ASV common whole snail.pdf"), useDingbats=FALSE)
#    print(p[[2]])
#    dev.off()


    pdf(paste0(graph.d,s," ASV venn - tsh ",mytsh," pv ",mypv,".pdf"), height=2.5, width=2.5, useDingbats=FALSE)
    venn(c, zcolor=myclr, box=FALSE)
    dev.off()
}




#---------#
# Heatmap #
#---------#

p.ls <- vector("list", length(nb.rep))
for (r in nb.rep) {
    
    if (r ==1 ) {
        mybiom.tmp <- subset_samples(mybiom, Replicate == r)
    } else {
        mybiom.tmp <- mybiom
    }
    
    mygrp <- split(get_variable(mybiom.tmp), get_variable(mybiom.tmp, "Population")) %>% lapply(., rownames)

    mygrp.otu <- lapply(mygrp, function(x) otu_table(mybiom.tmp)[,x] %>% rowSums())

    # Top 50 OTU from each group
    mygrp.otu.s <- lapply(mygrp.otu, function(x) sort(x, decreasing=TRUE) %>% names() %>% head(., 50))

    mygrp.otu.uniq <- mygrp.otu.s %>% unlist() %>% unique()

    # Prune table
    mybiom.o  <- prune_taxa(mygrp.otu.uniq, mybiom.tmp)

    # Identify interval to add blank bars
    spl.count <- count(sample_data(mybiom.o)[,2])
    myint     <- cumsum(spl.count[-nrow(spl.count),2])

    # Plot
    p.ls[[r]] <- plot_heatmap(mybiom.o, method=NULL, distance="bray", sample.order=unlist(mygrp), taxa.order=mygrp.otu.uniq) +
                    geom_vline(xintercept = myint+0.5, size = 2, color="white") +
                    labs(y="ASV") +
                    theme(axis.ticks.y=element_blank(), axis.line=element_blank(), axis.text.y=element_blank())
}

pdf(paste0(graph.d,"taxo-div_hm.pdf"), width=10, height=20, useDingbats=FALSE)
p.ls[[1]]
dev.off()


pdf(paste0(graph.d,"taxo-div_hm_all.pdf"), width=20, height=20, useDingbats=FALSE)
p.ls[[2]]
dev.off()


## Top 50 OTU from each group
#mygrp.otu.s <- lapply(mygrp.otu, function(x) sort(x[x > 100]) %>% names() %>% head(., 50))
#
#mygrp.otu.uniq <- mygrp.otu.s %>% unlist() %>% unique()
#
#d <- prune_taxa(mygrp.otu.uniq, mybiom)
#
#pdf("heatmap_rare.pdf", width=10, height=20)
#plot_heatmap(d, method=NULL, distance="bray", sample.order=unlist(mygrp), taxa.order=mygrp.otu.uniq)
#dev.off()




# Close the registered backend cluster
stopCluster(cl)
rm(cl)
registerDoSEQ()



#=========#
# Figures #
#=========#

if (! dir.exists(graph.dir) ) { dir.create(graph.dir, recursive=TRUE) }

#clr <- rainbow(unique(mypop[,2]) %>% length())
clr <- colorspace::rainbow_hcl(unique(mypop[,2]) %>% length())
names(clr) <- unique(mypop[,2]) 

#pdf(paste0(graph.dir,"Nb_OTU.pdf"), width=4, height=3)
#    par(mar=c(3,2,0,0)+0.1)
#    boxplot(mypop.otu[,3] ~ mypop.otu[,2], ylab="Number of OTUs", xlab="Popuations")
#dev.off()
#
#
#pdf(paste0(graph.dir,"Phyla_by_pop.pdf"), width=5)
#    plot_composition(mybiom.core.ss.P.t, average_by = "Population", sample.sort=sort(sample_names(mybiom.core.ss.P.t))) +
#    scale_fill_discrete(name="Phylum", breaks=taxa_names(mybiom.core.ss.P.t), labels=( tax_table(mybiom.core.ss.P.t)[,2] %>% strsplit(., "_") %>% lapply(., function(x) rev(x)[1]) %>% unlist() ) )
#dev.off()
#
##color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
##color <- grDevices::terrain.colors(length(taxa_names(mybiom.core.ss.P.t)))
##color <- cm.colors(length(taxa_names(mybiom.core.ss.P.t)))
#
#pdf(paste0(graph.dir,"Phyla_by_samp.pdf"), width=15)
#    plot_composition(mybiom.core.ss.P.t, sample.sort=sort(sample_names(mybiom.core.ss.P.t))) +
#    scale_fill_discrete(name="Phylum", breaks=taxa_names(mybiom.core.ss.P.t), labels=( tax_table(mybiom.core.ss.P.t)[,2] %>% strsplit(., "_") %>% lapply(., function(x) rev(x)[1]) %>% unlist() ))
##scale_fill_manual(values=color, name="Phylum", breaks=taxa_names(mybiom.core.ss.P.t), labels=( tax_table(mybiom.core.ss.P.t)[,2] %>% strsplit(., "_") %>% lapply(., function(x) rev(x)[1]) %>% unlist() ) ) 
#dev.off()
#
#pdf(paste0(graph.dir,"Phyla_rep.pdf"), width=5)
#    plot_composition(mybiom.core.ss.P.t.sp, sample.sort=sort(sample_names(mybiom.core.ss.P.t.sp))) +
#    scale_fill_discrete(name="Phylum", breaks=taxa_names(mybiom.core.ss.P.t.sp), labels=( tax_table(mybiom.core.ss.P.t.sp)[,2] %>% strsplit(., "_") %>% lapply(., function(x) rev(x)[1]) %>% unlist() ))
##scale_fill_manual(values=color, name="Phylum", breaks=taxa_names(mybiom.core.ss.P.t), labels=( tax_table(mybiom.core.ss.P.t)[,2] %>% strsplit(., "_") %>% lapply(., function(x) rev(x)[1]) %>% unlist() ) ) 
#dev.off()


#pdf(paste0(graph.dir,"Fig 2.pdf"), width=7.5, height=3.5)
pdf(paste0(graph.dir,"Fig 3.pdf"), width=6.5, height=3.5, useDingbats=FALSE)

mycex <- 0.6

#par(cex.axis=0.55, cex.lab=0.75)
par(cex.axis=mycex, cex.lab=mycex*1.2)
par(mar=c(3,3,1,0)+0.1)

layout(matrix(1:4, ncol=2, byrow=TRUE))

#------#
# OTUs #
#------#

# Graph
b <- boxplot(mypop.otu[,3] ~ mypop.otu[,2], ylab="", outline=FALSE, frame=FALSE, lwd=0.5, axes=FALSE, col=clr[order(names(clr))])

# x axis
axis(1, at=1:length(b$n), labels=FALSE, lwd=0.5)
text(1:length(b$n), par("usr")[3] - diff(par("yaxp")[1:2])*0.1, labels=b$names, srt=45, pos=1, cex=mycex, xpd=TRUE)

# y axis
title(ylab="Observed OTUs", line=2)
axis(2, lwd=0.5)

# Annotation
myshift <- par("pin")[2] / diff(par("usr")[3:4]) * (diff(par("yaxp")[1:2]) / par("yaxp")[3]) * 1000
text(1:length(b$n), b$stats[5,]+myshift, otu.letters, xpd=TRUE, cex=mycex)

#------#
# qPCR #
#------#

# Graph
b <- boxplot(qpcr[,2] ~ qpcr[,3], ylab="", outline=FALSE, frame=FALSE, lwd=0.5, axes=FALSE, col=clr[order(names(clr))])

# x axis
axis(1, at=1:length(b$n), labels=FALSE, lwd=0.5)
text(1:length(b$n), par("usr")[3] - diff(par("yaxp")[1:2])*0.1, labels=b$names, srt=45, pos=1, cex=mycex, xpd=TRUE)

# y axis
axis(2, lwd=0.5)
title(ylab=expression("Number of 16S per "* mu * "L"), line=2)
#title(ylab=("Number of 16S per µL"), line=2)

#-------------------#
# Rarefaction curve #
#-------------------#

table.rc <- table.rc[ table.rc[,"Depth"] <= 80000, ]
plot(0, ylim=c(0, max(table.rc[,3])), xlim=c(0, max(table.rc[,1])), axes=FALSE, bty="n", type="n", xlab="", ylab="")
#plot(0, ylim=c(0, max(table.rc[table.rc[,"Depth"] <= 80000,3])), xlim=c(0, 8e4), axes=FALSE, bty="n", type="n", xlab="", ylab="")
title(xlab="Number of sequences sampled", ylab="Observed OTUs", line=2)
axis(1, lwd=0.5)
axis(2, lwd=0.5)
for (p in unique(mypop[,2])) {
    b.tmp <- table.rc[ table.rc[,2] == p, ]
    plotCI(b.tmp[,1], b.tmp[,3], ui=b.tmp[,3]+b.tmp[,4], li=b.tmp[,3]-b.tmp[,4], pch=19, col=clr[match(p,unique(mypop[,2]))], add=TRUE)
#   fit <- lm(V1 ~ poly(Depth,2,raw=TRUE), b.tmp)
#   fit <- loess(V1 ~ Depth, b.tmp)
   fit <- lm(V1 ~ log10(Depth), b.tmp)
#    fit <- gam(V1 ~ s(Depth,k=nrow(b.tmp)), data=b.tmp)
    xx  <- seq(1, max(table.rc[,1]), (max(table.rc[,1]) - 1)/1000)#length=nrow(b.tmp))
    pred <- predict(fit, data.frame(Depth=xx), interval="confidence")
#    pred <- predict(fit, data.frame(Depth=xx))
    lines(xx, pred[,1], col=clr[match(p,unique(mypop[,2]))])
#    polygon(c(xx, rev(xx)), c(pred[, 2], rev(pred[, 3])), col=gsub("#","#80",clr[match(p,unique(mypop[,2]))]), border=NA)
}
#legend("topleft", legend=unique(mypop[,2]) %>% sort, col=clr,lty=1, lwd=0.5, bty="n", cex=0.5)


#------#
# PCoA #
#------#

# How to integrate ggplot in layout: https://stackoverflow.com/a/14125565

plot.new()
vps <- baseViewports()
pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
vp1 <-plotViewport(c(0,0,0,0)) ## create new vp with margins, you play with this values 

#pdf(paste0(graph.dir,"PCoA_2.pdf"), width=3.5, height=3.5)
#    p <- plot_ordination(mb, e, color = "Population")$data
#    p <- p[!grepl("Ba|water",p[,"Species"]),]
#    mycol <- sapply(as.vector(p[,6]), function(x) clr[match(x, unique(mypop[,2]))])
#    plot(p[,2] ~ p[,1], xlab=paste0("Axis 1 (",round(e$value[1,2]*100, 2), "%)"), ylab=paste0("Axis 2 (",round(e$value[2,2]*100, 2), "%)"), col=mycol, pch=19)
p <- plot_ordination(mb, e, color = "Population") +
        scale_color_manual(values=clr[order(names(clr))]) +
        stat_ellipse(type = "norm") + theme_bw() +
        theme(axis.text = element_text(size = rel(0.45))) +
        theme(axis.title = element_text(size = rel(0.6))) +
        theme(legend.position="none")
p$data <- p$data[!grepl("Ba|water",p$data[,"Species"]),]
#p <- p + stat_ellipse(type = "norm") + theme_bw()
#print(p + theme(legend.position="none"), vp=vp1)
print(p, vp=vp1)
    
dev.off()