#!/usr/bin/env Rscript
# Title: microbiome_module.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2024-05-23



#===========#
# Variables #
#===========#

# Population order and color
pop.data <- matrix(c("Ba",     "#ffac40",
                     "BgBS90", "#87d687"), byrow=TRUE, ncol=2)

# Sample type order
org.data <- matrix(c(
                "H",  "Hemolymph",      "Hm",
                "S",  "Stomach",        "S",
                "G",  "Gut",            "G",
                "L",  "Hepatopancreas", "Hp",
                "O",  "Ovotestis",      "O",
                "W",  "Whole snail",    "W",
                "TY", "Water tray",     "Ty",
                "TK", "Water tank",     "Tk"
                ), ncol = 3, byrow = TRUE)

org.tp <- org.data[ grep("TY|TK", org.data[,1], invert = TRUE), ]

# Folders
graph.d  <- "../graphs/"
result.d <- "../results/"

# GGplot options
theme_set(theme_classic())

