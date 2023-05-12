#' Perform Principal Components Analysis on a DESeqTransform object
#' [reused code from https://gist.github.com/tavareshugo/5ca8a5e18fedc3f23f5ec4b09a9fc906]
#'
#' This function is based on the `DESeq2::plotPCA()` function, but returns the
#' results of `prcomp` in a tidy list format. This is more flexible for further
#' custom plotting and exploring factor loadings of the PCA.
#'
#' @param x an object of class DESeqTransform
#' @ntop number of most-variable genes to select. Igored if "genes" is specified.
#' @genes character vector of specific genes to use
#'
#' @return a list with four `data.frame` objects: pc_scores, eigen_values,
#' loadings (eigen vectors) and the original data.
#' @export
#' @importFrom dplyr %>%
prcompTidy <- function(x, ntop = 500L, genes = NULL, ...){
    # suppressPackageStartupMessages(library(magrittr))

    # Get sample info
    sample_info <- as.data.frame(SummarizedExperiment::colData(x))

    # Get counts
    x <- SummarizedExperiment::assay(x)

    if(!is.null(genes)){
        message("Only using ", genes, " genes as requested.")

        if(!all(genes %in% rownames(x))) stop("Not all provided genes are in the gene count matrix.")

        selected_genes <- which(genes %in% rownames(x))

    } else if(is.numeric(ntop) & ntop < nrow(x)){
        ntop <- round(ntop)
        message("Only using ", ntop, " most variable genes.")

        # calculate the variance for each gene
        rv <- genefilter::rowVars(x)

        # select the ntop genes by variance
        selected_genes <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

    } else {
        message("Using all ", nrow(x), " genes.")
        selected_genes <- 1:nrow(x)
    }

    # Get the data for those genes
    selected_expr <- x[selected_genes, ]

    # perform a PCA on the data in assay(x) for the selected genes
    ## Need to transpose the matrix as prcomp clusters by rows
    pca <- prcomp(t(selected_expr), ...)

    #### Eigen scores table ####
    # Get sample information from DESeq x
    # and bind the PC scores
    pc_scores <- sample_info %>%
        dplyr::bind_cols(as.data.frame(pca$x))

    #### Eigen values table ####
    eigen_values <- data.frame(PC = colnames(pca$x), stdev = pca$sdev) %>%
        dplyr::mutate(var = stdev^2,
                      var_pct = var/sum(var),
                      cum_var = cumsum(var_pct),
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


# ggPCA <- function(vsd) {
#    combo_pcaData <- plotPCA(vsd,
#                              intgroup = c("Temperature","Population","Sex"),
#                              returnData=TRUE)
#     percentVar <- round(100 * attr(combo_pcaData, "percentVar"))
#     p<-ggplot(combo_pcaData, aes(PC1, PC2,  text=name,color=source)) +
#                  geom_point(aes(shape=factor_organism_part,size=factor_dev_stage)) +# geom_line()+
#                  scale_shape_manual(values=1:nlevels(combo_pcaData$factor_organism_part)) +
#                  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#                  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#                  coord_fixed()
#     return(p)
# }

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
    return(
        computedPCA$pc_scores %>%
            ggplot(aes_string(pc_x,pc_y,
                              text=labelBY,
                              shape=shapeBY,#"Sex",
                              color=colorBY, #"Population",
                              size=sizeBY))+ #"Temperature"
            geom_point()+
            guides(shape=guide_legend(title = shapeBY,title.position = "top"),
                   color=guide_legend(title.position = "top",title = colorBY),
                   size=guide_legend(title.position = "top",title = sizeBY))+
            scale_shape_manual(values=1:nlevels(computedPCA$pc_scores[,shapeBY]))+
            ggtitle("Top 500 variable genes")+
            labs(x=paste(pc_x,pct[1]),y=paste(pc_y,pct[2]))+
            theme_bw())
}
