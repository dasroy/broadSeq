#' Title
#'
#' @param se
#' @param cpm.log
#' @param normLibSizes.method "none" ,"TMM"
#' @param ... passed to normLibSizes function
#'
#' @return
#' @export
#' @import edgeR
#'
#' @examples
normalizeEdgerCPM <- function(se, normLibSizes.method="TMM",cpm.log = TRUE, ...){
    stopifnot(is(se, "SummarizedExperiment"))
    x<- normLibSizes(se,method=normLibSizes.method,...)
    x<- cpm(x,log=cpm.log)
    assay_name <- ifelse(cpm.log,"log","")
    assay_name <- paste(assay_name,
                        ifelse(normLibSizes.method != "none",
                               normLibSizes.method,
                               "CPM"),sep = "")
    SummarizedExperiment::assays(se)[[assay_name]] <- x
    return(se)
}

#' Title
#'
#' @param se
#' @param method  "vst" or "rlog"
#' @param ...
#'
#' @return
#' @export
#' @import DESeq2
#'
#' @examples
transformDESeq2 <- function(se, method="vst", ...){
    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot("method must be either vst or rlog"= (method %in% c("vst", "rlog")))
    # Need to convert to DESeqDataSet first from SummarizedExperiment object where 'counts' assay should be the first in assays list
    if (SummarizedExperiment::assayNames(se)[1] != "counts") {
        nuOrder <-
            c("counts",
              setdiff(SummarizedExperiment::assayNames(se), "counts"))
        SummarizedExperiment::assays(se) <-
            SummarizedExperiment::assays(se)[nuOrder]
    }

    if (mode(assays(se)[["counts"]]) != "integer") {
        counts <- round(as.matrix(assays(se)[["counts"]]))
        mode(counts) <- "integer"
        assays(se)[["counts"]] <- counts
    }
    x <- DESeqDataSet(se, design = ~ 1)
    if (method == "rlog") {
        rd <- rlog(x)
    } else{
        rd <- varianceStabilizingTransformation(x, ...)
    }
    assays(se)[[method]] <- assay(rd)
    return(se)
}
