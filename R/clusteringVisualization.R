
#' Title
#'
#' @param se
#' @param scaledAssay
#' @param ntop
#' @param genes
#' @param annotation_col a character vector of colnames of colData(se)
#'
#' @return
#' @export
#' @import pheatmap
#'
#' @examples
plotHeatmapCluster <- function(se, scaledAssay="vst", ntop = 500L, genes = NULL,
                               annotation_col=c(""),...){
    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot("scaledAssay must be an assay of se"= (scaledAssay %in% assayNames(se)))
    stopifnot(
        "annotation_col must be a subset of colnames of colData(se)"= (annotation_col %in% colnames(colData(se)))
        )
    selected_genes <- getSelectedGene(genes, se, ntop, scaledAssay)

    dottedArg <- list( ...)
    if(is.null(dottedArg$main) ) main = paste("Top ",length(selected_genes) ," variable genes")
    df <- as.data.frame(colData(se)[,annotation_col])
    p<-pheatmap(assays(se)[[scaledAssay]][selected_genes,],
             annotation_col=df,...)
    return(p)
}
#' Title
#'
#' @param se
#' @param scaledAssay
#' @param ntop
#' @param genes
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_MDS <- function(se, scaledAssay="vst", ntop = 500L, genes = NULL,
                     shapeBY=NULL, colorBY=NULL, sizeBY=NULL, labelBY=NULL){
    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot("scaledAssay must be an assay of se"= (scaledAssay %in% assayNames(se)))

    selected_genes <- getSelectedGene(genes, se, ntop, scaledAssay)

    # Get the data for those genes
    selected_expr <- assays(se)[[scaledAssay]][selected_genes, ]
    sampleDists <- dist( t( selected_expr ) )
    sampleDistMatrix <- as.matrix( sampleDists )
    mds <- data.frame(cmdscale(sampleDistMatrix))
    mds <- cbind(mds, colData(se))
    # mds %>% ggplot(aes(X1,X2,color=species,shape=stage))+geom_point()
    if(!is.null(shapeBY)){
        mds[,shapeBY] <- factor(mds[,shapeBY])
    }
    return(mds %>%
    ggplot(aes_string("X1","X2",
                      text=labelBY,
                      shape=shapeBY,
                      color=colorBY,
                      size=sizeBY))+
        geom_point()+
        guides(shape=guide_legend(title = shapeBY,title.position = "top"),
               color=guide_legend(title.position = "top",title = colorBY),
               size=guide_legend(title.position = "top",title = sizeBY))+
        scale_shape_manual(values=1:nlevels(mds[,shapeBY]))+
        ggtitle(paste("Top ",length(selected_genes) ," variable genes"))+
        theme_bw())
}


#' Perform Principal Components Analysis on a DESeqTransform object
#' [reused code from https://gist.github.com/tavareshugo/5ca8a5e18fedc3f23f5ec4b09a9fc906]
#'
#' This function is based on the `DESeq2::plotPCA()` function, but returns the
#' results of `prcomp` in a tidy list format. This is more flexible for further
#' custom plotting and exploring factor loadings of the PCA.
#'
#' @param se
#' @param scaledAssay
#' @param ntop
#' @param genes
#' @param ...
#'
#' @ntop number of most-variable genes to select. Igored if "genes" is specified.
#' @genes character vector of specific genes to use
#'
#' @return a list with four `data.frame` objects: pc_scores, eigen_values,
#' loadings (eigen vectors) and the original data.
#' @export
#' @importFrom dplyr %>%
prcompTidy <- function(se, scaledAssay="vst", ntop = 500L, genes = NULL, ...){
    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot("scaledAssay must be an assay of se"= (scaledAssay %in% assayNames(se)))

    selected_genes <- getSelectedGene(genes, se, ntop, scaledAssay)

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

#' Title
#'
#' @param genes
#' @param se
#' @param ntop
#' @param scaledAssay
#'
#' @return
#' @importFrom genefilter rowVars
#'
#' @examples
getSelectedGene <- function(genes, se, ntop, scaledAssay) {
    if(!is.null(genes)){
        message("Only using ", genes, " genes as requested.")

        if(!all(genes %in% rownames(se))) stop("Not all provided genes are in the gene count matrix.")

        selected_genes <- which(genes %in% rownames(se))
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


#' Title
#'
#' @param computedPCA
#' @param x
#' @param y
#' @param shapeBY
#' @param colorBY
#' @param sizeBY
#' @param labelBY
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggplot2 ggplot aes_string geom_point guide_legend guides guide_legend scale_shape_manual ggtitle labs theme_bw
plotAnyPC <- function(computedPCA,x,y, shapeBY=NULL, colorBY=NULL, sizeBY=NULL, labelBY=NULL){
    pc_x = paste("PC",x,sep = "")
    pc_y = paste("PC",y,sep = "")
    pct <- computedPCA$eigen_values %>%
        dplyr::filter(PC %in% c(pc_x,pc_y)) %>%
        dplyr::pull(var_pct) %>%
        round(digits = 4) *100
    if(!is.null(shapeBY)){
        computedPCA$pc_scores[,shapeBY] <- factor(computedPCA$pc_scores[,shapeBY])
    }
    return(
        computedPCA$pc_scores %>%
            ggplot(aes_string(pc_x,pc_y,
                              text=labelBY,
                              shape=shapeBY,
                              color=colorBY,
                              size=sizeBY))+
            geom_point()+
            guides(shape=guide_legend(title = shapeBY,title.position = "top"),
                   color=guide_legend(title.position = "top",title = colorBY),
                   size=guide_legend(title.position = "top",title = sizeBY))+
            scale_shape_manual(values=1:nlevels(computedPCA$pc_scores[,shapeBY]))+
            ggtitle(paste("Top ",nrow(computedPCA$original) ," variable genes"))+
            labs(x=paste(pc_x,pct[1]),y=paste(pc_y,pct[2]))+
            theme_bw())
}
