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

#' Useful to visualize distribution of assay values for each sample. Plots 'boxplot'
#' of any assay for each sample. Aesthetic can be added from colData.
#'
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param assayName One of the values from SummarizedExperiment::assayNames(se)
#' @param ... other arguments to be passed to ggpubr::\code{\link{ggboxplot}}
#'
#' @return
#' @export
#' @importFrom ggpubr ggboxplot rotate_x_text
#' @importFrom sechm meltSE
#' @importFrom SummarizedExperiment assayNames
#' @examples
#' se <- readRDS(system.file("extdata","rat_vole_mouseSE_salmon.rds", package = "broadSeq"))
#'
#' sampleAssay_plot(se, assayName = "counts",
#' fill="stage", # stage is a column name of colData(se)
#' yscale="log2")
#'
#' se <- broadSeq::normalizeEdgerCPM(se ,method = "none",cpm.log = TRUE )
#'
#' sampleAssay_plot(se, assayName = "logCPM", fill="stage")
sampleAssay_plot <- function(se, assayName = "counts",...){
    checkNameSpace("sechm")
    stopifnot("assayName not found"=(assayName %in% SummarizedExperiment::assayNames(se)))

    d <- sechm::meltSE(se,features=rownames(se), assayName = assayName )
    d %>% ggpubr::ggboxplot(y = assayName, x= "sample", ...)+rotate_x_text()
}

#' Expression of multiple genes/features from a single assay as boxplot (or added dotplot)
#'
#'
#'
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param features a character vector of rownames or named list of character vectors
#'  where name is one of the colnames of rowdata.
#' @param assayName One of the values from SummarizedExperiment::assayNames(se);
#' default is "counts" assay
#' @param facet.by must be one of the column names of rowData(se). default
#' "feature" which is equivalent to rownames of rowData
#' @param x a column name of colData which will be used in x-axis
#' @param ... other arguments to be passed to ggpubr::\code{\link{ggboxplot}}
#'
#' @return
#' @export
#' @importFrom ggpubr ggboxplot facet
#' @importFrom dplyr filter
#' @examples
#' se <- readRDS(system.file("extdata", "mouseSE_dev_tooth_count_length_geneData.rds",
#'     package = "broadSeq"))
#' # The normalized values are added with the assay name "logCPM"
#' se <- broadSeq::normalizeEdgerCPM(se ,method = "none",cpm.log = TRUE )
#'
#' broadSeq::genes_plot(se,
#'                      features = c("ENSMUSG00000000003", "ENSMUSG00000000103"),
#'                      facet.by = "mgi_symbol", # column of rowData
#'                      x = "stage",  fill="stage")
#'
#' broadSeq::genes_plot(se,
#'                      features = list(mgi_symbol=c("Shh","Edar") ),
#'                      facet.by = "mgi_symbol", # column of rowData
#'                      x = "stage",  fill="stage")
#'
#'broadSeq::assay_plot(se, feature = c("ENSMUSG00000002633"),
#'                     assays =  c("counts","logCPM"),
#'                     x = "stage", fill="stage", add="dotplot", palette = "npg")
genes_plot <- function(se, features, assayName = "counts", facet.by = "feature",
                       x, ...){
    checkNameSpace("sechm")
    stopifnot("assayName not found"=(assayName %in% SummarizedExperiment::assayNames(se)))
    stopifnot("features is not character or list"=(is.character(features) | is.list(features)))
    stopifnot("facet.by is not a column name of rowData"=(facet.by %in% colnames(rowData(se))))

    if(is.list(features)){
        stopifnot("list length must be 1"=(length(features)==1))
        stopifnot("name of list must be one of colnames of rowData(se)"=(
            all(names(features) %in% colnames(rowData(se)) ) ) )
        features <- as.data.frame( rowData(se)) %>%
            filter(get(names(features)) %in%  features[[1]]) %>% rownames()
    }

    d <- sechm::meltSE(se,features )
    d %>% ggboxplot(y=assayName, ... ) %>%
        facet(facet.by , scale="free")
}

#' Boxplot of a single gene/feature from multiple assays
#'
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param feature a character vector of rownames or named list of character vectors
#'  where name is one of the column of rowdata.
#' @param assayNames  names from SummarizedExperiment::assayNames(se); default
#' value is "counts"
#' @param x a column name of colData which will be used in x-axis
#' @param ... other arguments to be passed to ggpubr::\code{\link{ggboxplot}}
#'
#' @return return an object of class ggarrange, which is a ggplot or a list of ggplot.
#' @export
#' @importFrom ggpubr ggboxplot facet annotate_figure
#' @rdname genes_plot
assay_plot <- function(se, feature, assayNames = c("counts"), x, ...){
    checkNameSpace("sechm")
    stopifnot("assayName not found"=(all(assayNames %in% SummarizedExperiment::assayNames(se))))
    stopifnot("features is not character or list"=(feature %in% SummarizedExperiment::rownames(se)))

    d <- sechm::meltSE(se,features = feature,assayName = assayNames )

    listPlot <- list()
    for(i in 1:length(assayNames)){
        listPlot[[i]] <- d %>% ggboxplot( y=assayNames[i],...)
    }

    ggarrange(plotlist = listPlot, common.legend = TRUE, legend =  "bottom") %>%
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

    # library(OrgDB,quietly = TRUE,character.only = T)
    requireNamespace(OrgDB,quietly = TRUE)
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

#' Applies round function on numeric columns of a data.frame.
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
