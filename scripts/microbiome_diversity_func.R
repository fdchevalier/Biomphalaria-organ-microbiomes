# Title: microbiome_diversity_function.R
# Version: 1.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-03-25
# Modified in: 2020-08-31



#=======================#
# Computation functions #
#=======================#

# Standard error
se <- function(x) { sd(x) / sqrt(length(x)) }


# Edited tax_glom function from phyloseq
## This removed NA columns from the taxa table
mytax_glom <- function(x, taxrank=NULL, ...) {
    if (is.null(taxrank)) { stop("taxrank required") }
    if (! any(grepl(taxrank, rank_names(x))) ) { stop("taxrank unknown") }

    idx <- grep(taxrank, rank_names(x))
    x   <- tax_glom(x, taxrank, ...)
    tax_table(x) <- tax_table(x)[,1:idx]

    return(x)
}


# Agglomerate taxa from the same sample category
mysample_glom <- function(x, cln) {

    myls <- split(get_variable(x), get_variable(x)[,cln])

    myotu <- myls %>% lapply(., rownames) %>% sapply(., function(y) rowSums(otu_table(x)[,y]))

    myspl <- myls %>% lapply(., function(y) y[1,]) %>% unsplit(., levels(get_variable(x)[,cln]))
    rownames(myspl) <- myspl[,cln]

    x <- phyloseq(otu_table(myotu, taxa_are_rows=TRUE), tax_table(x), phy_tree(x), sample_data(myspl))

    return(x)
}


# Generate letter labels from post-hoc test
generate.label.df <- function(x){

    # x         Tukey test
    # variable  variable used for the Anova

    # Dependency
    library(multcompView)

    # Extract labels and factor levels from Tukey post-hoc
    Tukey.levels <- x
    Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])

    # Put the labels in the same order as in the boxplot
    Tukey.labels$treatment <- rownames(Tukey.labels)
#    Tukey.labels           <- Tukey.labels[ order(Tukey.labels$treatment) , ]
    return(Tukey.labels)
}



#=====================#
# Graphical functions #
#=====================#

#-------------------#
# Rarefaction curve #
#-------------------#

## source: https://github.com/mahendra-mariadassou/phyloseq-extended/blob/9efcd88995d9e60104592d1edcee784cf59c3381/R/graphical_methods.R
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
    ## Args:
    ## - physeq: phyloseq class object, from which abundance data are extracted
    ## - step: Step size for sample size in rarefaction curves
    ## - label: Default `NULL`. Character string. The name of the variable
    ##          to map to text labels on the plot. Similar to color option
    ##          but for plotting text.
    ## - color: (Optional). Default ‘NULL’. Character string. The name of the
    ##          variable to map to colors in the plot. This can be a sample
    ##          variable (among the set returned by
    ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
    ##          returned by ‘rank_names(physeq)’).
    ##
    ##          Finally, The color scheme is chosen automatically by
    ##          ‘link{ggplot}’, but it can be modified afterward with an
    ##          additional layer using ‘scale_color_manual’.
    ## - color: Default `NULL`. Character string. The name of the variable
    ##          to map to text labels on the plot. Similar to color option
    ##          but for plotting text.
    ## - plot:  Logical, should the graphic be plotted.
    ## - parallel: should rarefaction be parallelized (using parallel framework)
    ## - se:    Default TRUE. Logical. Should standard errors be computed.

    ## Require vegan and ggplot2
    library("vegan")
    library("ggplot2")

    x <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(physeq)) { x <- t(x) }

    ## This script is adapted from vegan `rarecurve` function
    tot <- rowSums(x)
    S <- rowSums(x > 0)
    nr <- nrow(x)

    rarefun <- function(i) {
        cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i]) {
            n <- c(n, tot[i])
        }
        y <- rarefy(x[i, ,drop = FALSE], n, se = se)
        if (nrow(y) != 1) {
            rownames(y) <- c(".S", ".se")
            return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
        } else {
            return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
        }
    }

    if (parallel) {
        out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
    } else {
        out <- lapply(seq_len(nr), rarefun)
    }
    df <- do.call(rbind, out)

    ## Get sample data
    if (!is.null(sample_data(physeq, FALSE))) {
        sdf <- as(sample_data(physeq), "data.frame")
        sdf$Sample <- rownames(sdf)
        data <- merge(df, sdf, by = "Sample")
        labels <- data.frame(x = tot, y = S, Sample = rownames(x))
        labels <- merge(labels, sdf, by = "Sample")
    }

    ## Add, any custom-supplied plot-mapped variables
    if( length(color) > 1 ){
        data$color <- color
        names(data)[names(data)=="color"] <- deparse(substitute(color))
        color <- deparse(substitute(color))
    }
    if( length(label) > 1 ){
        labels$label <- label
        names(labels)[names(labels)=="label"] <- deparse(substitute(label))
        label <- deparse(substitute(label))
    }

    p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
    p <- p + labs(x = "Sample Size", y = "Species Richness")

    if (!is.null(label)) {
        p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                           size = 4, hjust = 0)
    }
    p <- p + geom_line()

    if (se) { ## add standard error if available
        p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
    }

    if (plot) {
        plot(p)
    }
    invisible(p)
}


#--------------#
# Shift legend #
#--------------#

## source: https://stackoverflow.com/a/58734961
shift_legend <- function(p, position = "center") {
    # Dependencies
    suppressMessages({
        library("cowplot")
        library("lemon")
    })

    pnls <- cowplot::plot_to_gtable(p) %>% gtable::gtable_filter("panel") %>%
        with(setNames(grobs, layout$name)) %>% purrr::keep(~identical(.x,zeroGrob()))

    if( length(pnls) == 0 ) stop( "No empty facets in the plot" )

    lemon::reposition_legend( p, position, panel=names(pnls) )
}


#---------------------#
# Plot p-value matrix #
#---------------------#

## source: https://github.com/b0rxa/scmamp/blob/d74d0c085f82ed4bfef45ed537f773b20d9eaf85/R/plotting.R
plotPvalues <- function(pvalue.matrix, alg.order=NULL, show.pvalue=TRUE, font.size=NULL, border.col="white", pval.col="white", scientific = TRUE) {

    # Dependencies
    if (!requireNamespace("ggplot2", quietly=TRUE)) {
        stop("This function requires the ggplot2 package. Please install it.", call.=FALSE)
    }
    suppressMessages(library("reshape2"))

    # Convert the matrix into a data frame and order the algorithms according to
    # the desired order.
    df <- melt(pvalue.matrix)
    colnames(df) <- c("X", "Y", "p.value")
    if (!is.null(alg.order)) {
        l <- colnames(pvalue.matrix)[alg.order]
        df$X <- factor(df$X, levels=l)
        df$Y <- factor(df$Y, levels=l)
    }

    gplot <- ggplot2::ggplot(df, ggplot2::aes(x=X, y=Y, fill=p.value)) + ggplot2::geom_tile(col=border.col) +
    ggplot2::scale_fill_continuous("p-value") + ggplot2::labs(x="Algorithm" , y="Algorithm")

    if (show.pvalue) {
        ## geom_text size are not font size
        ## source :https://community.rstudio.com/t/why-does-ggplot-size-parameter-not-behave-consistently/21619/4
        if (is.null(font.size)) font.size <- GeomLabel$default_aes$size / .pt
        p.value.f <- df$p.value
        if (scientific) {
            p.value.f[ ! is.na(p.value.f) ] <- format(p.value.f[ ! is.na(p.value.f) ], digits=2, scientific = TRUE, na.encode = FALSE)
        } else {
            p.value.f <- round(p.value.f, 2)
        }
        gplot <- gplot + ggplot2::geom_text(ggplot2::aes(label = p.value.f),
                                        size=font.size, col=pval.col)
    }
    return(gplot)
}


#------------------------#
# Distinct color palette #
#------------------------#

## source: https://github.com/ronammar/randomcoloR/blob/e8bf2c21c4e9e6b14a95e20e0e19853840f41b7e/R/randomcolor.R
distinctColorPalette <-function(k=1, altCol=FALSE, runTsne=FALSE, seed=42) {
    # k         number of colors (>= 1). May be ineffective for k > 40.
    # altCol    Use an alternate color space
    # runTsne   Preprocess color space with t-SNE to obtain distinct colors. Reduces performance.
    # seed      number to be repeatable.

    # Dependencies
    library("cluster")
    library("colorspace")
    library("Rtsne")

    # Compute a 2000 color spectrum and convert to LAB
    runif(1)
    old_seed <- .Random.seed ## Keep this seed work only locally https://stackoverflow.com/a/14324316/1608734
    on.exit({.Random.seed <<- old_seed})
    set.seed(seed)
    n <- 2e3
    currentColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
    currentColorSpace <- as(currentColorSpace, "LAB")
    currentColorSpace <- currentColorSpace@coords
    if (altCol) {
        currentColorSpace <- t(unique(grDevices::col2rgb(scales::hue_pal(l=60:100)(n)))) ## Note: hue_pal no longer accepts multiple l
    }

    if (runTsne) {
        # Run 2D t-SNE before clustering
        tsne <- Rtsne(currentColorSpace, perplexity=50, check_duplicates=FALSE, pca=FALSE, max_iter=500)
        pamx <- pam(tsne$Y, k)  # k-medoids
        if (altCol) {
            colors <- rgb(currentColorSpace[pamx$id.med, ], maxColorValue=255)
        } else {
            colors <- hex(LAB(currentColorSpace[pamx$id.med, ]))
        }
    } else {
        # Set iter.max to 20 to avoid convergence warnings.
        km <- kmeans(currentColorSpace, k, iter.max=20)
        if (altCol) {
            colors <- rgb(round(km$centers), maxColorValue=255)
        } else {
            colors <- unname(hex(LAB(km$centers)))
        }
    }

    return(colors)
}

