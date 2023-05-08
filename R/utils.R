#' Title
#'
#' @param packageName
#'
#' @return
#'
#' @examples
checkNameSpace <- function(packageName) {
    if (!requireNamespace(packageName, quietly = TRUE)) {
        stop(
            paste("Package \"",  packageName , "\" must be installed to use this function.", sep = ""),
            call. = FALSE
        )
    }else{
        return(TRUE)
    }
}

#' Title
#'
#' @param smrExpt
#' @param features
#' @param x_factor
#' @param y_value
#'
#' @return
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter facet_wrap
#' @examples
genes_plot <- function(smrExpt, features,x_factor,y_value){
    checkNameSpace("sechm")
    checkNameSpace("ggplot2")

    d <- sechm::meltSE(smrExpt,features )
    ggplot2::ggplot(d, aes(.data[[x_factor]], .data[[y_value]], fill=.data[[x_factor]])) +
        geom_boxplot() +
        geom_jitter()+
        facet_wrap(~feature, scale="free")
}

rnaSeq_rank <- function(Methods_list){
    Methods_list$limma <- Methods_list$limma[ order(Methods_list$limma$adj.P.Val),]
    Methods_list$DEseq <- Methods_list$DEseq[ order(Methods_list$DEseq$padj),]
    Methods_list$edgeR <- Methods_list$edgeR[ order(Methods_list$edgeR$FDR),]
    Methods_list$DELocal <- Methods_list$DELocal[ order(-abs(Methods_list$DELocal$logFC)),]
    Methods_list$DELocal_TAD <- Methods_list$DELocal_TAD[ order(-abs(Methods_list$DELocal_TAD$logFC)),]

    makeRank <- function(result){
        result %>%
            dplyr::mutate(rank = 1:n()) %>%
            dplyr::filter(tooth_genes==TRUE) %>%
            dplyr::mutate(sequence_rank = 1:dplyr::n()) %>%
            dplyr::select(ensembl_gene_id,rank,sequence_rank)
    }
    tables <- lapply(Methods_list, makeRank)

    combined_r <- do.call(rbind , tables)

    combined_r$method <- rownames(combined_r) %>% stringr::str_extract(pattern = "[a-zA-Z_]*")
    return(combined_r)
}


