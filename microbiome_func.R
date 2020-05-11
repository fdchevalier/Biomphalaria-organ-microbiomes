se <- function(x) { sd(x) / sqrt(length(x)) }

tsv2csv <- function(x) {
    require(magrittr)
    require(rio)
    
    myfn <- basename(x) %>% tools::file_path_sans_ext()
    out.path <- paste0(tempdir(),"/",myfn,".csv")
    convert(x, out.path, in_opts = list(format = "tsv"))
    return(out.path)
}


mytax_glom <- function(x, taxrank=NULL, ...) {
    if (is.null(taxrank)) { stop("taxrank required") }
    if (! any(grepl(taxrank, rank_names(x))) ) { stop("taxrank unknown") }

    idx <- grep(taxrank, rank_names(x))
    x   <- tax_glom(x, taxrank, ...)
    tax_table(x) <- tax_table(x)[,1:idx]

    return(x)
}

mysample_glom <- function(x, cln) {

    myls <- split(get_variable(x), get_variable(x)[,cln])

    myotu <- myls %>% lapply(., rownames) %>% sapply(., function(y) rowSums(otu_table(x)[,y]))

#    otu_table(x) <- otu_table(myotu, taxa_are_rows=TRUE)
    myspl <- myls %>% lapply(., function(y) y[1,]) %>% unsplit(., levels(get_variable(x)[,cln]))
    rownames(myspl) <- myspl[,cln]

    x <- phyloseq(otu_table(myotu, taxa_are_rows=TRUE), tax_table(x), phy_tree(x), sample_data(myspl))

    return(x)
}

calculate_rarefaction_curves <- function(psdata, measures, depths, parallel=FALSE) {
    require('plyr') # ldply
    require('reshape2') # melt
    require('doParallel')

    # Set parallel options if required
    if (parallel) {
        paropts  <- list(.packages=c("phyloseq", "reshape2"))
    } else {
        paropts  <- NULL
    }
    
    estimate_rarified_richness <- function(psdata, measures, depth) {
        if(max(sample_sums(psdata)) < depth) return()
        psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)

        rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)

        alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)

        # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
        molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')

        molten_alpha_diversity
    } 

    names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
    rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive() && ! parallel, 'text', 'none'), .parallel=parallel, .paropts=paropts)

    # convert Depth from factor to numeric
    rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]

    rarefaction_curve_data
}



#generate.label.df <- function(x, variable){
#
#    # x         Tukey test
#    # variable  variable used for the Anova
#
#    # Dependency
#    library(multcompView)
#
#     # Extract labels and factor levels from Tukey post-hoc
#     Tukey.levels <- x[[variable]][,4]
#     Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
#
#     # Put the labels in the same order as in the boxplot
#     Tukey.labels$treatment <- rownames(Tukey.labels)
#     Tukey.labels           <- Tukey.labels[ order(Tukey.labels$treatment) , ]
#     return(Tukey.labels)
#}

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
#     Tukey.labels           <- Tukey.labels[ order(Tukey.labels$treatment) , ]
     return(Tukey.labels)
}











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


## source: https://stackoverflow.com/a/58734961
shift_legend <- function(p) {
    # Dependencies
    library("cowplot")
    library("lemon")
    
    pnls <- cowplot::plot_to_gtable(p) %>% gtable::gtable_filter("panel") %>%
        with(setNames(grobs, layout$name)) %>% purrr::keep(~identical(.x,zeroGrob()))

    if( length(pnls) == 0 ) stop( "No empty facets in the plot" )

    lemon::reposition_legend( p, "center", panel=names(pnls) )
}
