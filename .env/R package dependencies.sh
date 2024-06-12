#!/bin/bash

function R_install {
    wget -q -O /tmp/R_pckg.arx $1

    mkdir -p /tmp/R_pckg
    if [[ ${1##*.} == "zip" ]]
    then
        unzip -qq -d /tmp/R_pckg /tmp/R_pckg.arx
    else
        tar -xzf /tmp/R_pckg.arx -C /tmp/R_pckg
    fi

    folder=$(dirname $(find /tmp/R_pckg/ -name DESCRIPTION))
    R CMD INSTALL "$folder"
    rm -R /tmp/R_pckg*
}

# Colorout
R_install https://github.com/jalvesaq/colorout/archive/7ea9440.zip

# Pairwise Adonis package
R_install https://github.com/pmartinezarbizu/pairwiseAdonis/archive/6e09713.zip

# Taxa package
R_install https://github.com/ropensci/taxa/archive/49969dc.zip

# qiime2R package
R_install https://github.com/jbisanz/qiime2R/archive/574708a.zip

# Updated version of ape
R_install https://cran.r-project.org/src/contrib/Archive/ape/ape_5.6-1.tar.gz
