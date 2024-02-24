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

#' Useful to visualize transformed/normalized assay with count assay. Plots 'boxplot' of any assay for each sample. Aesthetic can be added from colData.
#'
#' @param smrExpt
#' @param assayName
#' @param ...
#'
#' @return
#' @export
#' @importFrom ggpubr ggboxplot
#' @importFrom sechm meltSE
#' @importFrom SummarizedExperiment assayNames
#' @examples
sampleAssay_plot <- function(smrExpt, assayName = "counts",...){
    checkNameSpace("sechm")
    stopifnot("assayName not found"=(assayName %in% SummarizedExperiment::assayNames(smrExpt)))

    d <- sechm::meltSE(smrExpt,features=rownames(smrExpt),assayName = assayName )
    d %>% ggpubr::ggboxplot(y = assayName, x= "sample", ...)+rotate_x_text()
}

#' Title
#'
#' @param smrExpt
#' @param features
#' @param x
#' @param y
#' facet.by must be one of the column names of rowData(smrExpt). default "feature" which is equivalent to rownames of rowData
#'
#' @return
#' @export
#' @importFrom ggpubr ggboxplot facet
#' @importFrom dplyr filter
#' @examples
genes_plot <- function(smrExpt, features,facet.by = "feature",...){
    checkNameSpace("sechm")
    stopifnot("features is not character or list"=(is.character(features) | is.list(features)))
    if(is.list(features)){
        stopifnot("list length must be 1"=(length(features)==1))
        stopifnot("name of list must be one of colnames of rowData(smrExpt)"=( all(names(features) %in% colnames(rowData(smrExpt)) ) ) )
        features <- as.data.frame( rowData(smrExpt)) %>% filter(get(names(features)) %in%  features[[1]]) %>% rownames()
    }

    d <- sechm::meltSE(smrExpt,features )
    d %>% ggboxplot( ... ) %>%
        facet(facet.by , scale="free")
}

#' Title
#'
#' @param smrExpt
#' @param feature
#' @param assays
#' @param x
#' @param y
#'
#' @return
#' @export
#' @importFrom ggpubr ggboxplot facet
#' @examples
assay_plot <- function(smrExpt, feature,assays,...){
    checkNameSpace("sechm")

    d <- sechm::meltSE(smrExpt,features = feature,assayName = assays )

    listPlot <- list()
    for(i in 1:length(assays)){
        listPlot[[i]] <- d %>% ggboxplot( y=assays[i],...)
    }

    ggarrange(plotlist = listPlot,common.legend = TRUE, legend =  "bottom") %>%
        annotate_figure(top = text_grob(feature, color = "red", face = "bold", size = 14))
}


#' Provides GO gene set enrichment and over-representation analysis
#'
#' This wrapper function combines clusterProfiler::gseGO and
#' clusterProfiler::enrichGO. The input type of thes two methods are different;
#' order ranked geneList and a vector of entrez gene id. Here combinedEnrichment
#' function internally generates these two data types from user defined DEG_table
#' (differentially expresssed genes).
#'
#'
#' @param DEG_table A data.frame atleast with two columns.
#' @param geneCol The column name of DEG_table which provides gene ids and should
#' be compatible with keytype parameter.
#' @param logCol The column name of DEG_table which provides logfold(numeric) values to
#' create a order ranked geneList for gseGO funtion.
#' @param logfoldCut to filter genes based on parameter logCol
#' @param OrgDB OrgDb; passed to clusterProfiler functions
#' @param keyType keytype of input gene(geneCol). One of the keytypes(OrgDB); passed to clusterProfiler functions
#' @param universe background genes; passed to clusterProfiler::enrichGO.
#' @param ont one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.; passed to clusterProfiler functions
#' @param pvalueCutoff ; passed to clusterProfiler functions
#' @param qvalueCutoff ; passed to clusterProfiler functions
#'
#' @return a named list of three data.frames which are output of gseGO("gseResult")
#' and enrichGO ("oraUP" and "oraDOWN").
#' @export
#'
#' @importFrom clusterProfiler enrichGO gseGO
#' @examples
combinedEnrichment <- function(DEG_table, geneCol = "ID", logCol = "logFoldChange",
                               OrgDB = "org.Hs.eg.db", keyType, universe, ont  = "BP",
                               logfoldCut = 1, pvalueCutoff  = 0.05, qvalueCutoff  = 0.05 ){

    library(OrgDB,quietly = TRUE,character.only = T)
    stopifnot(keyType %in% keytypes(get(OrgDB)))
    ## feature 1: numeric vector
    DEGlist = DEG_table[[logCol]]

    ## feature 2: named vector
    names(DEGlist) = as.character(DEG_table[,geneCol])

    DEGlist = sort(DEGlist, decreasing = TRUE)
    print("Performing: Gene Set Enrichment Analysis")
    gseResult <- gseGO(DEGlist,
                       OrgDb         = OrgDB,
                       ont           = ont,
                       # minGSSize    = 100,
                       # maxGSSize    = 500,
                       pvalueCutoff = pvalueCutoff,
                       verbose      = FALSE,
                       keyType = keyType)

    print("Performing: Over-representation analysis of UP regulated genes")
    oraUP <- enrichGO(subset(DEG_table, get(logCol) >= logfoldCut)[,geneCol],
                      universe      = universe,
                      OrgDb         = OrgDB,
                      ont           = ont,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = pvalueCutoff,
                      qvalueCutoff  = qvalueCutoff,
                      readable      = TRUE,
                      keyType = keyType)

    print("Performing: Over-representation analysis of DOWN regulated genes")
    oraDOWN <- enrichGO(subset(DEG_table, get(logCol) <= -1*logfoldCut)[,geneCol],
                        universe      = universe,
                        OrgDb         = OrgDB,
                        ont           = ont,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = pvalueCutoff,
                        qvalueCutoff  = qvalueCutoff,
                        readable      = TRUE,
                        keyType = keyType)

    return(list(gseResult=gseResult, oraUP= oraUP, oraDOWN = oraDOWN))
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

is.na2 = function(x){
    if(is.list(x) | length(x) > 1){
        return(FALSE)
    }
    if(length(x) == 0){
        return(TRUE)
    }

    return(is.na(x))
}
