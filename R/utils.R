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
#' @importFrom ggpubr ggboxplot facet
#' @examples
genes_plot <- function(smrExpt, features,...){
    checkNameSpace("sechm")
    checkNameSpace("ggpubr")

    d <- sechm::meltSE(smrExpt,features )
    d %>% ggboxplot( ... ) %>%
        facet(facet.by = "feature", scale="free")
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

#' Applies round function on numeric columns of a dataframe.
#'
#' @param df
#' @param digits
#'
#' @return
#' @export
#'
#' @examples
round_df <- function(df, digits) {
    nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
    df[,nums] <- round(df[,nums], digits = digits)
    (df)
}

