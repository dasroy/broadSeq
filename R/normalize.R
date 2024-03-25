#' Use of edgeR package to normalize count data
#'
#'
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param cpm.log value for edgeR::\code{\link{cpm}} function. default TRUE
#' @param method value for edgeR::\code{\link{normLibSizes}} function. default "TMM"
#' @param ... passed to normLibSizes function
#'
#' @return Object of \code{\link{SummarizedExperiment}} class where a new assay
#' is added to the input object.
#' @export
#' @import edgeR
#'
#' @examples
normalizeEdgerCPM <- function(se, method="TMM", cpm.log = TRUE, ...){
    stopifnot(is(se, "SummarizedExperiment"))
    x<- normLibSizes(se,method=method,...)
    x<- cpm(x,log=cpm.log)
    assay_name <- ifelse(cpm.log,"log","") #if_else is better to handle na values
    assay_name <- paste(assay_name,
                        ifelse(method != "none",
                               method, "CPM"), sep = "")
    message("New assay name ",assay_name)
    SummarizedExperiment::assays(se)[[assay_name]] <- x
    return(se)
}

#' Transform SummarizedExperiment with DESeq2 package
#'
#' To use SummarizedExperiment with DESeq2, this function makes sure that 'counts'
#' assay should be the first in assays list and the mode is integer.
#'
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param method "vst", "normTransform" or "rlog" to choose either of DESeq2::\code{\link{varianceStabilizingTransformation}}
#' DESeq2::\code{\link{normTransform}} and DESeq2::\code{\link{rlog}} function.
#' default is "vst"
#' @param ... arguments passed to \code{\link{varianceStabilizingTransformation}}
#' \code{\link{normTransform}} and \code{\link{rlog}}
#'
#' @return Object of \code{\link{SummarizedExperiment}} class where a new assay
#' is added to the input object.
#'
#' @export
#' @import DESeq2
#' @importFrom SummarizedExperiment assays
#'
#' @examples
transformDESeq2 <- function(se, method="vst", ...){
    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot("method must be either vst or rlog"= (method %in% c("vst", "rlog","normTransform")))
    warning("For length correction assayname must match with avgTxLength\n ")
    ## Should be checked
    # Need to convert to DESeqDataSet first from SummarizedExperiment object where
    # 'counts' assay should be the first in assays list
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
        rd <- rlog(x,...)
    } else if(method == "vst"){
        rd <- varianceStabilizingTransformation(x, ...)
    }else if(method == "normTransform"){
        rd <- normTransform(x, ...)
    }
    assays(se)[[method]] <- assay(rd)
    return(se)
}
