
#' Plot clustered heatmaps
#'
#' Plot clustered heatmaps from SummarizedExperiment with pheatmap and return object as ggplot
#'
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param scaledAssay an scaled assay name from SummarizedExperiment::assayNames(se)
#' @param ntop number of most-variable genes to select. Igored if "features" is specified.
#' @param features character vector features/genes to be used to measure distance between the samples
#' @param annotation_col a character vector of colnames of colData(se)
#' @param annotation_row a list of character vector of colnames of rowData(se)
#' @param show_geneAs a character vector of colnames of rowData(se)
#' @param ... other arguments like color or shape whose values should be similar
#' to colData columns names passed to \code{\link{pheatmap}}
#'
#' @return ggplot object
#' @export
#' @import pheatmap
#' @importFrom ggplotify as.ggplot
#' @importFrom dplyr select
#'
#' @examples
#' se <- readRDS(system.file("extdata","rat_vole_mouseSE_salmon.rds", package = "broadSeq"))
#'
#' se <- broadSeq::normalizeEdgerCPM(se ,method = "none",cpm.log = TRUE )
#'
#' broadSeq::plotHeatmapCluster(
#'     se,
#'     scaledAssay = "logCPM",
#'     annotation_col = c("species", "stage"),
#'     annotation_row = c("Class","gene_biotype"),
#'     ntop = 30, show_geneAs = "symbol",
#'     cluster_cols = TRUE, cluster_rows = FALSE,
#'     show_rownames = TRUE, show_colnames = FALSE,
#'     main = "Top 30 variable gene vst"
#' )
plotHeatmapCluster <- function(se, scaledAssay="vst", ntop = 500L,
                               features = NULL, show_geneAs = NULL,
                               annotation_col=NA, annotation_row=NA,...){
    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot("scaledAssay must be an assay of se"= (scaledAssay %in% assayNames(se)))

    selected_genes <- getSelectedGene(features, se, ntop, scaledAssay)

    dottedArg <- list( ...)
    if(is.null(dottedArg$main) ) main = paste("Top ",length(selected_genes) ," variable genes")

    data_mat <- assays(se)[[scaledAssay]][selected_genes,]
    if(is.na2(annotation_col) == FALSE){
        stopifnot(
            "annotation_col must be a subset of colnames of colData(se)"= ((annotation_col %in% colnames(colData(se))))
        )
        df <- as.data.frame(colData(se)[,annotation_col])
        colnames(df) <- annotation_col
        rownames(df) <- colnames(se)
        annotation_col <- df
    }

    if(is.na2(annotation_row) == FALSE){
        stopifnot(
            "annotation_row must be a subset of colnames of rowData(se)"= ((annotation_row %in% colnames(rowData(se))))
        )
        df <- rowData(se)[selected_genes,] %>% as.data.frame()  %>%
            dplyr::select(all_of(annotation_row))
        colnames(df) <- annotation_row
        annotation_row <- df
    }

    if(!is.null(show_geneAs)){
        stopifnot(
            "show_geneAs must be a subset of colnames of rowData(se)"= ( (show_geneAs %in% colnames(rowData(se))))
        )
        gene_name <- as.data.frame(rowData(se)[selected_genes,])  %>% dplyr::select(all_of(show_geneAs))
        rownames(data_mat) <- gene_name[rownames(data_mat),]
        if(is.data.frame(annotation_row)){
             rownames(annotation_row) <- gene_name[rownames(annotation_row),]
        }
    }

    p<-pheatmap(data_mat, annotation_col=annotation_col,
                annotation_row=annotation_row,
                ...)
    return(as.ggplot(p))
}

#' Classical multidimensional scaling
#'
#' Classical multidimensional scaling is based on measuring the distance between the samples.
#'
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param scaledAssay an scaled assay name from SummarizedExperiment::assayNames(se)
#' @param ntop number of most-variable genes to select. Igored if "features" is specified.
#' @param features character vector features/genes to be used to measure distance between the samples
#' @param ... other arguments like color or shape whose values should be similar
#' to colData columns names passed to ggpubr::\code{\link{ggscatter}}
#'
#' @return ggplot object
#' @export
#' @importFrom ggpubr ggscatter ggpar
#' @examples
#' se <- readRDS(system.file("extdata","rat_vole_mouseSE_salmon.rds", package = "broadSeq"))
#'
#' se <- broadSeq::transformDESeq2(se,method = "vst"  )
#'
#' broadSeq::plot_MDS(se, scaledAssay = "vst", ntop=500,
#'                     color = "species", shape = "stage",
#'                     ellipse=TRUE, legend = "bottom")
plot_MDS <- function(se, scaledAssay="vst", ntop = 500L, features = NULL,...){
    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot("scaledAssay must be an assay of se"= (scaledAssay %in% assayNames(se)))

    selected_genes <- getSelectedGene(features, se, ntop, scaledAssay)

    # Get the data for those genes
    selected_expr <- assays(se)[[scaledAssay]][selected_genes, ]
    sampleDists <- dist( t( selected_expr ) )
    sampleDistMatrix <- as.matrix( sampleDists )
    mds <- data.frame(cmdscale(sampleDistMatrix))
    mds <- cbind(mds, colData(se))

    return(mds %>%
               ggscatter(x = "X1", y = "X2", ...) %>%
               ggpar(title = paste("Top ",length(selected_genes) ," variable genes"))
        )
}


#' Perform Principal Components Analysis
#'
#'
#' This function returns the results of stats::\code{\link{prcomp}} in a tidy list format.
#' This is more flexible for further custom PCA , biplot and exploring gene(factor) loading of the PCA.
#'
#' [Reused code](https://gist.github.com/tavareshugo/5ca8a5e18fedc3f23f5ec4b09a9fc906)
#'
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param scaledAssay an scaled assay name from SummarizedExperiment::assayNames(se)
#' @param ntop number of most-variable genes to select. Igored if "features" is specified.
#' @param features character vector features/genes to be used for PCA
#' @param ... other arguments to be passed to stats::\code{\link{prcomp}}
#'
#' @return a list with four `data.frame` objects: pc_scores, eigen_values,
#' loadings (eigen vectors) and the original data.
#' @export
#' @importFrom dplyr %>%
#' @examples
#' se <- readRDS(system.file("extdata","rat_vole_mouseSE_salmon.rds", package = "broadSeq"))
#'
#' se <- broadSeq::normalizeEdgerCPM(se ,method = "none",cpm.log = TRUE )
#' computedPCA_logCPM <- broadSeq::prcompTidy(se, scaledAssay = "logCPM", ntop = 500)
#'
#' plotAnyPC(computedPCA = computedPCA_logCPM, x = 1, y = 2, color = "species",
#'          shape = "stage",legend = "bottom")
#' plotAnyPC(computedPCA = computedPCA_logCPM, x = 2, y = 3, color = "species",
#'          shape = "stage",legend = "bottom")
#'
#' computedPCA_logCPM$eigen_values %>%
#'  dplyr::filter(var_exp >= 0.5) %>% # Selecting PC explaining more than 1% variance
#'     ggbarplot(x="PC",y="var_exp", label = TRUE, label.pos = "out")
#'
prcompTidy <- function(se, scaledAssay="vst", ntop = 500L, features = NULL, ...){
    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot("scaledAssay must be an assay of se"= (scaledAssay %in% assayNames(se)))
    stopifnot("features should be either NULL or subset of rownames"=
                  (is.null(features) | features %in% SummarizedExperiment::rownames(se)))

    selected_genes <- getSelectedGene(features, se, ntop, scaledAssay)

    # Get the data for those genes
    selected_expr <- assays(se)[[scaledAssay]][selected_genes, ]

    # perform a PCA on the data in assay(x) for the selected genes
    ## Need to transpose the matrix as prcomp clusters by rows
    pca <- prcomp(t(selected_expr), ...)

    # Get sample info
    sample_info <- as.data.frame(SummarizedExperiment::colData(se))

    #### Eigen scores table ####
    pc_scores <- sample_info %>%
        dplyr::bind_cols(as.data.frame(pca$x))

    #### Eigen values table ####
    eigen_values <- data.frame(PC = colnames(pca$x), stdev = pca$sdev) %>%
        dplyr::mutate(var = stdev^2,
                      var_pct = var/sum(var),
                      cum_var = cumsum(var_pct),
                      var_exp =round(var_pct,digits = 4) *100 , ## Variance Explained
                      PC = forcats::fct_inorder(PC))

    #### Factor loadings table ####
    factor_loadings <- pca$rotation %>%
        as.data.frame() %>%
        dplyr::mutate(gene = row.names(.)) %>%
        dplyr::select(gene, dplyr::everything())
    rownames(factor_loadings) <- factor_loadings$gene

    gene_info <- as.data.frame(SummarizedExperiment::rowData(se))

    factor_loadings <- purrr::reduce(lapply(list(factor_loadings,gene_info),
                              function(x) data.frame(x, rn = row.names(x))),
                       merge, all.x=TRUE)
    rownames(factor_loadings) <- factor_loadings$rn
    factor_loadings$rn <- NULL

    #### Convert the original data to a data.frame ####
    selected_expr <- selected_expr %>%
        as.data.frame() %>%
        dplyr::rename_all(dplyr::funs(paste0("sample", .))) %>%
        dplyr::mutate(gene = rownames(.)) %>%
        dplyr::select(gene, dplyr::everything())

    # Return a list with each of these xs
    return(list(pc_scores = pc_scores,
                eigen_values = eigen_values,
                loadings = factor_loadings,
                original = selected_expr))

}

#' @importFrom genefilter rowVars
getSelectedGene <- function(genes, se, ntop, scaledAssay) {
    if(!is.null(genes)){
        message("Only using ", length(genes), " genes as requested.")

        if(!all(genes %in% rownames(se))) stop("Not all provided genes are in the gene count matrix.")

        selected_genes <- genes # which(genes %in% rownames(se))
    } else if(is.numeric(ntop) & ntop < nrow(se)){
        ntop <- round(ntop)
        message("Only using ", ntop, " most variable genes.")

        # calculate the variance for each gene
        rv <- rowVars(assays(se)[[scaledAssay]])

        # select the ntop genes by variance
        selected_genes <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    } else {
        message("Using all ", nrow(se), " genes.")
        selected_genes <- 1:nrow(se)
    }
    return(selected_genes)
}


#' Plots principal components
#'
#' @param computedPCA a list of data.frame returned by \code{\link{prcompTidy}}
#' @param x PC number for x-axis default 1
#' @param y PC number for y-axis default 2
#' @param ... other arguments like color or shape whose values should be similar
#' to colData columns names passed to ggpubr::\code{\link{ggscatter}}
#'
#' @return ggplot object
#' @export
#'
#' @importFrom ggpubr ggscatter ggpar
#' @importFrom ggplot2 scale_shape_manual
#' @rdname prcompTidy
plotAnyPC <- function(computedPCA,x =1 ,y = 2, ...){
    pc_x = paste("PC",x,sep = "")
    pc_y = paste("PC",y,sep = "")
    pct <- computedPCA$eigen_values %>%
        dplyr::filter(PC %in% c(pc_x,pc_y)) %>%
        dplyr::pull(var_pct) %>%
        round(digits = 4) *100
    dottedArg <- list( ...)

    p <- computedPCA$pc_scores %>% ggscatter(x = pc_x, y = pc_y,...)
    if(!is.null(dottedArg$shape) ) {
        computedPCA$pc_scores[,dottedArg$shape] <- factor(computedPCA$pc_scores[,dottedArg$shape])
        p <- p + scale_shape_manual(values=1:nlevels(computedPCA$pc_scores[,dottedArg$shape]))
    }
    return(
        p %>%
            ggpar(title = paste(nrow(computedPCA$original) ," variable genes"),
                  xlab = paste(pc_x,pct[1]), ylab = paste(pc_y,pct[2]))
        )
}


#' biplotAnyPC
#'
#' @param computedPCA a list of data.frame returned by \code{\link{prcompTidy}}
#' @param genes if genes is NULL then top max and min loaded genes of each PCs are plotted
#' @param genesLabel one of rowData column names
#' @param ... other arguments like color or shape whose values should be similar
#' to colData columns names passed to ggpubr::\code{\link{ggscatter}}
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' @importFrom ggpubr ggscatter ggpar
#' @importFrom ggplot2 geom_segment geom_text scale_shape_manual sym
#' @rdname prcompTidy
biplotAnyPC <- function(computedPCA,x = 1,y = 2,
                        genes=NULL, genesLabel = NULL,...){
    pc_x = paste("PC",x,sep = "")
    pc_y = paste("PC",y,sep = "")
    pct <- computedPCA$eigen_values %>%
        dplyr::filter(PC %in% c(pc_x,pc_y)) %>%
        dplyr::pull(var_pct) %>%
        round(digits = 4) *100
    dottedArg <- list( ...)

    p <- computedPCA$pc_scores %>% ggscatter(x = pc_x, y = pc_y,...)
    if(!is.null(dottedArg$shape) ) {
        computedPCA$pc_scores[,dottedArg$shape] <- factor(computedPCA$pc_scores[,dottedArg$shape])
        p <- p + scale_shape_manual(values=1:nlevels(computedPCA$pc_scores[,dottedArg$shape]))
    }

    if(is.null(genes)){
        pc1_genes <- computedPCA$loadings %>%
            filter(!!sym(pc_x) == max(!!sym(pc_x)) | !!sym(pc_x) == min(!!sym(pc_x))) %>% pull(gene)
        pc2_genes <- computedPCA$loadings %>%
            filter(!!sym(pc_y) == max(!!sym(pc_y)) | !!sym(pc_y) == min(!!sym(pc_y))) %>% pull(gene)
        genes <- c(pc1_genes,pc2_genes)
    }
    dA <- computedPCA$loadings[genes,] %>%
        dplyr::select(!!sym(pc_x),!!sym(pc_y), !!sym(genesLabel)) %>%
        dplyr::mutate_if( is.numeric, ~ . * 100)
    p <- p +
        geom_segment(data=dA, aes(x=0, y=0,xend=!!sym(pc_x), yend=!!sym(pc_y)), size = 1,
                     arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")+
        geom_text(data=dA, aes(x=!!sym(pc_x), y=!!sym(pc_y),label= !!sym(genesLabel)),
                  size = 3, vjust=1, color="blue")
    return(
        p %>%
            ggpar(title = paste(nrow(computedPCA$original) ," variable genes"),
                  xlab = paste(pc_x,pct[1]), ylab = paste(pc_y,pct[2]))
        )
}

#' Useful for PCA
#'
#' @param computedPCA a list of data.frame returned by \code{\link{prcompTidy}}
#' @param pcs The numbers of PCs
#' @param topN Number of features per PC
#' @param keep the column names of rowData to keep the corresponding information
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' @rdname prcompTidy
getFeatureLoadRanking <- function(computedPCA, pcs=1:5, topN = 10, keep ){
    x <- plyr::ldply(paste("PC",pcs,sep = ""), .fun = extract_topGeneLoadings,
                loadings = computedPCA$loadings,
                topN = topN, keep = keep)
    x$PC <- x$PC %>% factor(levels =  paste("PC",1:10,sep = ""))
    return(x)
}

extract_topGeneLoadings <- function(loadings,whichpc,topN,keep){
    geneloadings_sorted <- dplyr::arrange(loadings, desc(abs( !!sym(whichpc)))) %>%
        dplyr::select(all_of(keep), dplyr::all_of(whichpc))
    geneloadings_extreme <- geneloadings_sorted %>% dplyr::slice(1:topN)#, with_ties = FALSE)
    colnames(geneloadings_extreme) <- c(dplyr::all_of(keep), "loading")
    geneloadings_extreme$PC <- whichpc
    geneloadings_extreme$Rank <- 1:topN
    return(geneloadings_extreme)
}
