## Based on https://montilab.github.io/BS831/articles/docs/DiffanalysisRNAseqComparison.html
## https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-020-76881-x/MediaObjects/41598_2020_76881_MOESM1_ESM.pdf

#' @export
#' @rdname use_limma
use_limma_trend <- function(se, colData_id, control, treatment, rank=FALSE,...){
    showPlot <- FALSE
    dottedArg <- list( ...)
    if(!is.null(dottedArg$showPlot) ) showPlot <- dottedArg$showPlot

    return(use_limma(se = se, colData_id = colData_id,
                     control = control, treatment = treatment, rank = rank,
                     useVoom = FALSE, showPlot = showPlot,...))
}

#' @export
#' @rdname use_limma
use_limma_voom <- function(se, colData_id, control, treatment,rank=FALSE,...){
    showPlot <- FALSE
    dottedArg <- list( ...)
    if(!is.null(dottedArg$showPlot) ) showPlot <- dottedArg$showPlot

    return(use_limma(se = se, colData_id = colData_id,
                     control = control, treatment = treatment, rank = rank,
                     useVoom = TRUE, showPlot = showPlot,...))
}

#' To use SummarizedExperiment with package limma
#'
#' A wrapper function of limma where input is an object of \code{\link{SummarizedExperiment}}
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param colData_id One of the columns of colData(se). It should be factors of more than one value.
#' @param control Base level and one of the factor values of `colData(se)[[colData_id]]`
#' @param treatment one of the factor values of `colData(se)[[colData_id]]`
#' @param rank Logical value default FALSE. If true the result will have an
#' additional column named "rank"
#' @param useVoom whether to use limma::voom or edgeR::cpm
#' @param showPlot whether to use limma::plotSA ; default FALSE
#' @param ... other arguments to be passed to main function edgeR::calcNormFactors .
#' @param limma.adjust argument for limma::topTable
#' @param limma.sort.by argument for limma::topTable
#' @param limma.number argument for limma::topTable
#'
#' @return a data.frame of output from limma::topTable
#'
#' @export
#' @importFrom dplyr %>% left_join
#' @examples
#'
#' se <- readRDS(system.file("extdata",
#'         "rat_vole_mouseSE_salmon.rds",
#'         package = "broadSeq"))
#'
#' # To reduce runtime
#' se <- se[rowData(se)$chromosome_name == 2,colData(se)$species == "Mouse"]
#' result <-
#'     use_limma(se = se,
#'            colData_id = "stage", control = "Bud", treatment = "Cap",
#'            rank = TRUE)
use_limma <- function(se, colData_id, control, treatment,
                      rank=FALSE, useVoom=TRUE, showPlot=FALSE,
                      limma.adjust="BH", limma.sort.by = "p", limma.number=Inf,
                      ...) {
    checkNameSpace("limma")
    checkNameSpace("edgeR")
    se <- se[,se[[colData_id]] %in% c(control,treatment)]
    # dim(se)
    condition <- factor(as.character(SummarizedExperiment::colData(se)[,colData_id]),
                        levels = c(control,treatment))

    design <- stats::model.matrix(~ 0 + condition)
    colnames(design) <- c(control,treatment)

    command_str <- paste("limma::makeContrasts(",treatment , "-", control,",levels = design)"
                         , sep = "")
    contrast.matrix <- eval(parse(text=command_str))

    ##Normalizing
    dge <- edgeR::calcNormFactors(se, ...)

    if(showPlot & useVoom){
        par(mfrow=c(1,2))
    } else {
        par(mfrow=c(1,1))
    }
    if(useVoom){
        v <- limma::voom(dge, design, plot=showPlot)
        fit <- limma::lmFit(v, design)
    }else{
        ## Differential expression: Limma trend
        ## convert counts to logCPM
        logCPM <- edgeR::cpm(dge,log=TRUE,priot.count=3)
        fit <- limma::lmFit(logCPM, design )
    }

    fit <- limma::contrasts.fit(fit, contrast.matrix)
    fit <- limma::eBayes(fit,trend = !useVoom)

    res <- limma::topTable(fit, coef=1, adjust=limma.adjust,
                           sort.by = limma.sort.by, number=limma.number,
                           ...)
    if(isTRUE(showPlot)){
        limma::plotSA(fit, main="Final model: Mean-variance trend")
    }
    if(rank){
        return(res %>% dplyr::select(logFC:B) %>%
                   dplyr::mutate(rank := 1:dplyr::n()) )
    }else{
        return(res)
    }
}

#' @export
#' @rdname use_edgeR
use_edgeR_GLM <- function(se, colData_id, control, treatment, rank=FALSE,...){
    return(use_edgeR(se = se, colData_id = colData_id,
                     control = control, treatment = treatment, rank = rank,
                     option="GLM",...))
}

#' @export
#' @rdname use_edgeR
use_edgeR_exact <- function(se, colData_id, control, treatment, rank=FALSE,...){
    return(use_edgeR(se = se, colData_id = colData_id,
                     control = control, treatment = treatment, rank = rank,
                     option = "exact",...))
}


#' To use SummarizedExperiment with package edgeR
#'
#' A wrapper function of DESeq2 where input is an object of \code{\link{SummarizedExperiment}}
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param colData_id One of the columns of colData(se). It should be factors of more than one value.
#' @param control Base level and one of the factor values of `colData(se)[[colData_id]]`
#' @param treatment one of the factor values of `colData(se)[[colData_id]]`
#' @param rank Logical value default FALSE. If true the result will have an
#' additional column named "rank"
#' @param edgeR.n argument for edgeR::\code{\link{topTags}}
#' @param edgeR.adjust.method argument for edgeR::\code{\link{topTags}}
#' @param edgeR.sort.by argument for edgeR::\code{\link{topTags}}
#' @param option "GLM" or "exact" to indicate to use either edgeR::\code{\link{glmLRT}}
#' or edgeR::\code{\link{exactTest}}
#' @param ... other arguments to be passed to edgeR::\code{\link{glmLRT}}
#' or edgeR::\code{\link{exactTest}}
#'
#' @return a data.frame of output from edgeR::\code{\link{topTags}}
#' @export
#' @examples
#'
#' se <- readRDS(system.file("extdata",
#'         "rat_vole_mouseSE_salmon.rds",
#'         package = "broadSeq"))
#'
#' # To reduce runtime
#' se <- se[rowData(se)$chromosome_name == 2,colData(se)$species == "Mouse"]
#'
#' result <-
#'     use_edgeR(se = se,
#'            colData_id = "stage", control = "Bud", treatment = "Cap",
#'            rank = TRUE)
use_edgeR <- function(se, colData_id, control, treatment, rank=FALSE,
                      edgeR.n = Inf, edgeR.adjust.method = "BH", edgeR.sort.by = "PValue",
                      option="GLM",...){
    checkNameSpace("edgeR")
    stopifnot("option must be either 'GLM' or 'exact'"=(option=="GLM" | option=="exact"))
    se <- se[,se[[colData_id]] %in% c(control,treatment)]
    condition <- factor(as.character(SummarizedExperiment::colData(se)[,colData_id]),
                        levels = c(control,treatment))

    ##Normalizing
    dt <- edgeR::calcNormFactors(se,method="TMM")
    ##Normalization

    if(option=="GLM"){
        design <- stats::model.matrix(~ 0 + condition)
        colnames(design) <- c(control,treatment)

        d2 <- edgeR::estimateGLMCommonDisp(dt,design)
        d2 <- edgeR::estimateGLMTrendedDisp(d2,design, method="auto")
        # Possible method to "auto", "bin.spline", "power", "spline", "bin.loess".
        # "auto" chooses "bin.spline" when > 200 tags and "power" otherwise
        d2 <- edgeR::estimateGLMTagwiseDisp(d2,design)
        fit <- edgeR::glmFit(d2,design)

        command_str <- paste("limma::makeContrasts(",treatment , "-", control,
                             ",levels = design)" , sep = "")
        contrast.matrix <- eval(parse(text=command_str))
        t12 <- edgeR::glmLRT(fit,contrast=contrast.matrix,...)
    }else if(option=="exact"){
        ##Estimate dispersion
        dt$samples$group <- condition
        d1 <- edgeR::estimateCommonDisp(dt, verbose=TRUE)
        d1 <- edgeR::estimateTagwiseDisp(d1)
        ##Compare groups (exact test)
        t12 <- edgeR::exactTest(d1, pair=c(1,2),...)
    }
    res <- edgeR::topTags(t12, n=edgeR.n, sort.by = edgeR.sort.by ,
                          adjust.method = edgeR.adjust.method,
                          ...) %>% as.data.frame()
    if(rank){
        return(res %>% dplyr::select(logFC:FDR)%>%
                   dplyr::mutate(rank := 1:dplyr::n()) )
    }else{
        return(res)
    }
}

#' To use SummarizedExperiment with package DESeq2
#'
#' A wrapper function of DESeq2 where input is an object of \code{\link{SummarizedExperiment}}
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param colData_id One of the columns of colData(se). It should be factors of more than one value.
#' @param control Base level and one of the factor values of `colData(se)[[colData_id]]`
#' @param treatment one of the factor values of `colData(se)[[colData_id]]`
#' @param rank Logical value default FALSE. If true the result will have an
#' additional column named "rank" and the results are ranked on "padj"
#' @param ...  other arguments to be passed to main function DESeq2::\code{\link{results}}.
#'
#' @return a data.frame converted from DESeq2::\code{\link{DESeqResults}}
#' @export
#' @importFrom dplyr %>% left_join
#' @examples
#'
#' se <- readRDS(system.file("extdata",
#'         "rat_vole_mouseSE_salmon.rds",
#'         package = "broadSeq"))
#'
#' # To reduce runtime
#' se <- se[rowData(se)$chromosome_name == 2,colData(se)$species == "Mouse"]
#'
#' result <-
#'     use_deseq2(se = se,
#'            colData_id = "stage", control = "Bud", treatment = "Cap",
#'            rank = TRUE)
use_deseq2 <- function(se, colData_id, control, treatment,rank=FALSE,...){
    checkNameSpace("DESeq2")
     if (SummarizedExperiment::assayNames(se)[1] != "counts") {
         nuOrder <-
             c("counts",
               setdiff(SummarizedExperiment::assayNames(se), "counts"))
         SummarizedExperiment::assays(se) <-
             SummarizedExperiment::assays(se)[nuOrder]
     }
     se <- se[,se[[colData_id]] %in% c(control,treatment)]

    formula_str <- paste("formula( ~ ",colData_id, ")", sep = "")
    dds <- DESeq2::DESeqDataSet(se, design = eval(parse(text=formula_str)))

    ## set reference to control, otherwise default is alphabetical order
    dds[[colData_id]] <- factor(dds[[colData_id]], levels=c(control,treatment))

    ##   1. estimate size factors
    ##   2. estimate dispersion
    ##   3. negative binomial GLM fitting and wald test
    dds_res <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds_res,...)
    # res$dispersion <- DESeq2::dispersions(dds_res)
    if(rank){
        return(res %>% as.data.frame() %>%
               dplyr::arrange(padj) %>%
               dplyr::mutate(rank = 1:dplyr::n()))
    }else{
        return(res)
    }
}

#' To use SummarizedExperiment with package DELocal
#'
#' A wrapper function of DELocal where input is an object of \code{\link{SummarizedExperiment}}
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param colData_id One of the columns of colData(se). It should be factors of
#' more than one value.
#' @param control Base level and one of the factor values of `colData(se)[[colData_id]]`
#' @param treatment one of the factor values of `colData(se)[[colData_id]]`
#' @param rank Logical value default FALSE. If true the result will have an
#' additional column named "rank" and the results are ranked on "relative.logFC"
#'
#' @param ... other arguments to be passed to main function DELocal::\code{\link{DELocal}}.
#'
#' @return a data.frame from DELocal
#' @export
#' @importFrom dplyr %>% left_join
#' @examples
#'
#' se <- readRDS(system.file("extdata",
#'         "rat_vole_mouseSE_salmon.rds",
#'         package = "broadSeq"))
#'
#' # To reduce runtime
#' se <- se[rowData(se)$chromosome_name == 2,colData(se)$species == "Mouse"]
#'
#' result <-
#'     use_DELocal(se = se,
#'            colData_id = "stage", control = "Bud", treatment = "Cap",
#'            rank = TRUE)
use_DELocal <- function(se, colData_id, control, treatment,rank=FALSE,...){
    checkNameSpace("DELocal")
    formula_str <- paste("formula( ~ ",colData_id, ")", sep = "")
    if (SummarizedExperiment::assayNames(se)[1] != "counts") {
        nuOrder <-
            c("counts",
              setdiff(SummarizedExperiment::assayNames(se), "counts"))
        SummarizedExperiment::assays(se) <-
            SummarizedExperiment::assays(se)[nuOrder]
    }
    se <- se[,se[[colData_id]] %in% c(control,treatment)]

    DELocal_result <- DELocal::DELocal(pSmrExpt = se, # Genes without neighbours are missing
                              nearest_neighbours = 5, pDesign =eval(parse(text=formula_str)),
                              pValue_cut = 1, pLogFold_cut = 0,...)
    if(rank){
        return(DELocal_result %>%
                   dplyr::arrange(desc(abs(relative.logFC))) %>%
                   dplyr::mutate(rank = 1:dplyr::n()))
    }else{
        return(DELocal_result)
    }
}



#' To use SummarizedExperiment with package EBSeq
#'
#' A wrapper function of EBSeq where input is an object of \code{\link{SummarizedExperiment}}
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param colData_id One of the columns of colData(se). It should be factors of
#' more than one value.
#' @param control Base level and one of the factor values of `colData(se)[[colData_id]]`
#' @param treatment one of the factor values of `colData(se)[[colData_id]]`
#' @param rank Logical value default FALSE. If true the result will have an
#' additional column named "rank" and the results are ranked on "PPDE"
#' @param ... other arguments to be passed to main function EBSeq::\code{\link{GetDEResults}}.
#'
#' @return a data.frame object converted from the output of EBSeq::\code{\link{GetDEResults}}.
#' @export
#' @importFrom dplyr %>% left_join
#' @examples
#'
#' se <- readRDS(system.file("extdata",
#'         "rat_vole_mouseSE_salmon.rds",
#'         package = "broadSeq"))
#'
#' # To reduce runtime
#' se <- se[rowData(se)$chromosome_name == 2,colData(se)$species == "Mouse"]
#'
#' result <-
#'     use_EBSeq(se = se,
#'            colData_id = "stage", control = "Bud", treatment = "Cap",
#'            rank = TRUE)
use_EBSeq <- function(se, colData_id, control, treatment, rank=FALSE,...){
    checkNameSpace("EBSeq")
    control_names <- se[,se[[colData_id]]==control] %>% colnames()
    treatment_names <- se[,se[[colData_id]]==treatment] %>% colnames()
    tabla <- SummarizedExperiment::assays(se)[["counts"]][,c(control_names,treatment_names)]
    Sizes <- EBSeq::MedianNorm(tabla)

    EBOut <- EBSeq::EBTest(Data=tabla,
                 Conditions=factor(c(rep(control,length(control_names)),
                                     rep(treatment,length(treatment_names))),
                                      levels = c(control,treatment)),
                 sizeFactors=Sizes, maxround=5)
    EBDERes <- EBSeq::GetDEResults(EBOut, ...)
    ##Calculate FC
    FC <- EBSeq::PostFC(EBOut, SmallNum = 0.01)

    ##Obtain probabilities
    res <- data.frame(PPEE=EBDERes$PPMat[,"PPEE"] ,PPDE=EBDERes$PPMat[,"PPDE"],Status=EBDERes$Status,
                      Direction = FC$Direction) #,PostFC=FC$PostFC,RealFC=FC$RealFC
    if(rank){
        res %>% dplyr::arrange(desc(PPDE)) %>%
            dplyr::mutate(rank = 1:dplyr::n())
    }else{
        res
    }
}


#' Differential expression method for NOISeq
#'
#' This is a wrapper function of NOISeq::\code{\link{noiseqbio}} whose input class is ´eSet´
#' and output class is `Output` which are not widely used. We can use as(se, "ExpressionSet")
#' to get an eSet easily but then it will be hard to refer the treatment and control.
#' The order of factors influence the log fold change sign. To keep it comparable
#' to other methods the `NOISeq::readData()` is used internally.
#'
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param colData_id One of the columns of colData(se). It should be factors of more than one value.
#' @param control Base level and one of the factor values of `colData(se)[[colData_id]]`
#' @param treatment one of the factor values of `colData(se)[[colData_id]]`
#' @param rank Logical value default FALSE. If true the result will have an
#' additional column named "rank" which is ordered by ´prob´ values returned by
#' function NOISeq::\code{\link{noiseqbio}}.
#' @param ...  other arguments to be passed to main function NOISeq::\code{\link{noiseqbio}}.
#' The 'input' and 'factor' argument should not be used.
#'
#' @return A data.frame object from the results of NOISeq::noiseqbio(). For details check
#' the documentation of ´NOISeq´
#' @export
#' @importFrom dplyr %>% left_join
#' @examples
#'
#' se <- readRDS(system.file("extdata",
#'         "rat_vole_mouseSE_salmon.rds",
#'         package = "broadSeq"))
#'
#' # To reduce runtime
#' se <- se[rowData(se)$chromosome_name == 2,colData(se)$species == "Mouse"]
#'
#' result_Noiseq <-
#'     use_NOIseq(se = se,
#'            colData_id = "stage", control = "Bud", treatment = "Cap",
#'            rank = TRUE,
#'            r = 10) # r is an argument of NOISeq::noiseqbio
use_NOIseq <- function(se, colData_id, control, treatment, rank=FALSE, ...){
    control_names <- se[,se[[colData_id]]==control] %>% colnames()
    treatment_names <- se[,se[[colData_id]]==treatment] %>% colnames()
    tabla <- SummarizedExperiment::assays(se)[["counts"]][,c(control_names,treatment_names)] %>% as.data.frame()

    expt_factors <- data.frame(Condition = factor(c(
        rep(treatment, length(treatment_names)),
        rep(control, length(control_names))
    )))
    colnames(expt_factors) <- colData_id
    expt_data <- NOISeq::readData(data = tabla, factors = expt_factors)

    ##Available normalization methods: norm = "rpkm", "uqua", "tmm2, "none"
    expt_noiseqbio <- NOISeq::noiseqbio(input= expt_data, factor=colData_id, ...)

    if(rank){
        expt_noiseqbio@results[[1]] %>%
            dplyr::arrange(desc(prob)) %>%
            dplyr::mutate(rank = 1:dplyr::n())
    }else{
        expt_noiseqbio@results[[1]]
    }
}

#' To use SummarizedExperiment with package samr
#'
#'
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param colData_id One of the columns of colData(se). It should be factors of more than one value.
#' @param control Base level and one of the factor values of `colData(se)[[colData_id]]`
#' @param treatment one of the factor values of `colData(se)[[colData_id]]`
#' @param rank Logical value default FALSE. If true the result will have an
#' @param ... other arguments to be passed to samr::SAMseq
#'
#' @return a data.frame object as a result
#' @importFrom dplyr %>% left_join
#' @importFrom stringr str_remove
#' @examples
#'
use_SAMseq <- function(se, colData_id, control, treatment,rank = FALSE,...){
    checkNameSpace("samr")
    control_names <- se[,se[[colData_id]]==control] %>% colnames()
    treatment_names <- se[,se[[colData_id]]==treatment] %>% colnames()
    tabla <- SummarizedExperiment::assays(se)[["counts"]][,c(control_names,treatment_names)] %>% as.data.frame()
    y <- factor(c(rep(control,length(control_names)),
                 rep(treatment,length(treatment_names))),
               levels = c(control,treatment))

    # it uses cat() for lots of text, annoying
    samfit <- samr::SAMseq(tabla, y, resp.type = "Two class unpaired",
                           geneid = rownames(tabla), genenames = rownames(tabla),
                           random.seed=1, fdr.output = 1,...)
    if(is.null(samfit$siggenes.table$genes.up) & is.null(samfit$siggenes.table$genes.lo)){
        return(NULL)
    }else if(is.null(samfit$siggenes.table$genes.up)) {
        res <- as.data.frame(samfit$siggenes.table$genes.lo)
        res[,"Fold Change"] <- as.numeric(res[,"Fold Change"]) * (-1)
    }else if(is.null(samfit$siggenes.table$genes.lo)) {
        res <- as.data.frame(samfit$siggenes.table$genes.up)
    }else{
        down <- as.data.frame(samfit$siggenes.table$genes.lo)
        down[,"Fold Change"] <- as.numeric(down[,"Fold Change"]) * (-1)
        res <- rbind(as.data.frame(samfit$siggenes.table$genes.up),down)
    }
    colnames(res) <- colnames(res) %>% stringr::str_remove(" ")
    missing_id <- setdiff(rownames(se),res$GeneID)
    missing_data <- data.frame(GeneID=missing_id, GeneName=missing_id,
                               `Score(d)` = NA, 'FoldChange' = NA, `q-value(%)` = NA)
    colnames(missing_data) <- colnames(res)

    res <- rbind(res,missing_data)
    row.names(res) <- res$GeneID

    res <- res %>% as.data.frame() %>% dplyr::mutate_at(c('Score(d)', 'FoldChange', 'q-value(%)'), as.numeric)

    res <- res %>% dplyr::select('Score(d)': 'q-value(%)')

    if(rank){
        res <- res %>%
         dplyr::arrange(as.numeric(`q-value(%)`))  %>%
         dplyr::mutate(rank = 1:dplyr::n())
    }

    return(res)
}

#' To identify differentially expressed genes by multiple methods
#'
#' @param deFun_list a list of function which can perform differential expression analysis
#' @param return.df whether to return all results aggregated form of data.frame or
#' a list of results. Default is FALSE
#' @param se Object of \code{\link{SummarizedExperiment}} class
#' @param colData_id One of the columns of colData(se). It should be factors of more than one value.
#' @param control Base level and one of the factor values of `colData(se)[[colData_id]]`
#' @param treatment one of the factor values of `colData(se)[[colData_id]]`
#' @param ... other arguments to be passed to functions listed in deFun_list
#'
#' @importFrom purrr reduce
#' @return a list or data.frame
#' @export
#'
#' @examples
#'
#' se <- readRDS(system.file("extdata",
#'         "rat_vole_mouseSE_salmon.rds",
#'         package = "broadSeq"))
#'
#' # To reduce runtime
#' se <- se[rowData(se)$chromosome_name == 2,colData(se)$species == "Mouse"]
#'
#' # First define a named list of functions
#' funs <- list(limma_trend = use_limma_trend, limma_voom = use_limma_voom,
#'              edgeR_exact = use_edgeR_exact, edgeR_glm = use_edgeR_GLM,
#'              deseq2 = use_deseq2,
#'              DELocal = use_DELocal, noiseq = use_NOIseq,
#'              EBSeq = use_EBSeq)
#'
#'
#' multi_result <- broadSeq::use_multDE(
#'     se <- se[rowData(se)$chromosome_name == 2,colData(se)$species == "Mouse"],
#'     deFun_list = funs, return.df = TRUE,
#'     colData_id = "stage", control = "Bud", treatment = "Cap",
#'     rank = TRUE)
use_multDE <- function(deFun_list, return.df= FALSE ,se ,
                       colData_id , control , treatment , ... ) {
    wrap_it <- function( f, ...) {
        message(paste("Now executing >> ",f))
        message("####")
        deFun_list[[f]](... )
    }
    SummarizedExperiment::colData(se)[,colData_id] %>% summary()

    y <- lapply(names(deFun_list), wrap_it, se = se, colData_id = colData_id,
                control = control, treatment = treatment, ...)

    named_result <- list()
    for (i in seq_len(length(y))) {
        named_result[[names(deFun_list)[i]]] <- y[[i]]
    }
    y <- named_result

    if(return.df){
        for (prefix in names(y)) {
            colnames(y[[prefix]]) <- paste(prefix,colnames(y[[prefix]]), sep = "_")
        }

        y <- purrr::reduce(lapply(y, function(x) data.frame(x, rn = row.names(x))),
                           merge, all=TRUE)
        gene_info <- as.data.frame(SummarizedExperiment::rowData(se))
        rownames(y) <- y$rn
        y$rn <- NULL
        y <- purrr::reduce(lapply(list(y,gene_info),
                                function(x) data.frame(x, rn = row.names(x))),
                         merge, all.x=TRUE)
        rownames(y) <- y$rn
        y$rn <- NULL
    }
    return(y)
}

#' Volcano plot with formatted x and y axis label.
#'
#' @param df a data.frame object
#' @param pValName column name of df which provides p-values
#' @param lFCName column name of df which provides log fold change values
#' @param sigThreshold Threshold for p-values
#' @param logFCThreshold Threshold for log fold change values
#' @param labelName column name of df to label the dots
#' @param selectedLabel which dots to highlight
#' @param palette one of "npg" ,"aaas", "lancet", "jco", "ucscgb", "uchicago",
#' "simpsons" and "nejm" or similar to viridis::cividis(3)
#'
#' @return ggplot object
#' @export
#' @importFrom dplyr %>% left_join if_else
#' @importFrom ggpubr ggscatter
#' @examples
volcanoPlot <- function(df,pValName,lFCName, sigThreshold=0.05, logFCThreshold = 1,
                        labelName=NULL, selectedLabel = NULL, palette = "nejm"){
    df <- df %>% dplyr::mutate(padj=-log10(!!sym(pValName)),
                              Significant = if_else((!!sym(pValName) < sigThreshold & !!sym(lFCName) > logFCThreshold),"UP",
                                                    if_else((!!sym(pValName) < sigThreshold & !!sym(lFCName) < -logFCThreshold),"DOWN",
                                                          "Not", missing="Not"), missing="Not"
                              ))
    df$Significant <- factor(df$Significant,levels = c("DOWN","UP","Not"))
    plot <- df %>%  ggscatter(
            x = lFCName, y = "padj",
            color = "Significant", palette =palette,
            title = "Volcano plot",
            label = labelName, repel = TRUE,label.rectangle=TRUE, show.legend=FALSE,
            label.select = selectedLabel)+
        labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))+
        geom_abline(intercept = -log10(sigThreshold), slope = 0, linetype="dashed")
    return(plot)
}
