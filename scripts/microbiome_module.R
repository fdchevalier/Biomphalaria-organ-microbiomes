#!/usr/bin/env Rscript
# Title: microbiome_module.R
# Version: 0.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2024-05-23
# Modified in: 2024-12-22



#===========#
# Variables #
#===========#

# Population order and color
pop.data <- matrix(c("Ba",     "#ffac40",
                     "BgBS90", "#87d687"), byrow=TRUE, ncol=2)

# Sample type order
org.data <- matrix(c(
                "H",  "Hemolymph",      "Hm", "#da2b2b",
                "S",  "Stomach",        "S",  "#7cae00",
                "G",  "Gut",            "G",  "#00be67",
                "L",  "Hepatopancreas", "Hp", "#000000",
                "O",  "Ovotestis",      "O",  "#9b9b9b",
                "W",  "Whole snail",    "W",  "#683131",
                "TY", "Water tray",     "Ty", "#22dadf",
                "TK", "Water tank",     "Tk", "#00a9ff"
                ), ncol = 4, byrow = TRUE)

org.tp <- org.data[ grep("TY|TK", org.data[,1], invert = TRUE), ]

# Folders
graph.d  <- "../graphs/"
result.d <- "../results/"

# GGplot options
theme_set(theme_classic())

