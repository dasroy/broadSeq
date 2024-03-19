## Based on https://montilab.github.io/BS831/articles/docs/DiffanalysisRNAseqComparison.html
## https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-020-76881-x/MediaObjects/41598_2020_76881_MOESM1_ESM.pdf


#'use_limma_trend
#'
#' @param smrExpt a `SummarizedExperiment` object
#' @param class_id one of the colnames of colData(smrExpt) and must be a factor.
#' @param control one of the factor values of class_id
#' @param treatment another factor values of class_id
#' @param rank (TRUE/FALSE) whether to rank genes
#' @param ... additional arguments for limma::topTable
#'
#' @return results of limma::topTable
#' @export
#'
#' @examples
use_limma_trend <- function(smrExpt, class_id, control, treatment, rank=FALSE,...){
    showPlot = FALSE
    dottedArg <- list( ...)
    if(!is.null(dottedArg$showPlot) ) showPlot = dottedArg$showPlot

    return(use_limma(smrExpt = smrExpt, class_id = class_id,
                     control = control, treatment = treatment, rank = rank,
                     useVoom = FALSE, showPlot = showPlot,...))
}

#' Title
#'
#' @param smrExpt
#' @param class_id
#' @param control
#' @param treatment
#' @param rank
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
use_limma_voom <- function(smrExpt, class_id, control, treatment,rank=FALSE,...){
    showPlot = FALSE
    dottedArg <- list( ...)
    if(!is.null(dottedArg$showPlot) ) showPlot = dottedArg$showPlot

    return(use_limma(smrExpt = smrExpt, class_id = class_id,
                     control = control, treatment = treatment, rank = rank,
                     useVoom = TRUE, showPlot = showPlot,...))
}

#' Title
#'
#' @param smrExpt
#' @param class_id
#' @param control
#' @param treatment
#' @param useVoom
#' @param rank
#' @param showPlot
#' @param ...
#' @param limma.adjust
#' @param limma.sort.by
#' @param limma.number
#'
#' @return
#'
#' @importFrom dplyr %>% left_join
#' @examples
use_limma <- function(smrExpt, class_id, control, treatment,
                      rank=FALSE, useVoom=TRUE, showPlot=FALSE,
                      limma.adjust="BH", limma.sort.by = "p", limma.number=Inf,
                      ...) {
    checkNameSpace("limma")
    checkNameSpace("edgeR")
    smrExpt <- smrExpt[,smrExpt[[class_id]] %in% c(control,treatment)]
    # dim(smrExpt)
    condition <- factor(as.character(SummarizedExperiment::colData(smrExpt)[,class_id]),
                        levels = c(control,treatment))

    design <- stats::model.matrix(~ 0 + condition)
    colnames(design) <- c(control,treatment)

    command_str <- paste("limma::makeContrasts(",treatment , "-", control,",levels = design)"
                         , sep = "")
    contrast.matrix <- eval(parse(text=command_str))

    ##Normalizing
    dge <- edgeR::calcNormFactors(smrExpt, ...)

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

    res <- limma::topTable(fit, coef=1, adjust=limma.adjust, sort.by = limma.sort.by, number=limma.number,
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

#' Title
#'
#' @param smrExpt
#' @param class_id
#' @param control
#' @param treatment
#' @param rank
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
use_edgeR_GLM <- function(smrExpt, class_id, control, treatment, rank=FALSE,...){
    return(use_edgeR(smrExpt = smrExpt, class_id = class_id,
                     control = control, treatment = treatment, rank = rank,
                     option="GLM",...))
}


#' Title
#'
#' @param smrExpt
#' @param class_id
#' @param control
#' @param treatment
#' @param rank
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
use_edgeR_exact <- function(smrExpt, class_id, control, treatment, rank=FALSE,...){
    return(use_edgeR(smrExpt = smrExpt, class_id = class_id,
                     control = control, treatment = treatment, rank = rank,
                     option = "exact",...))
}

#' Title
#'
#' @param smrExpt
#' @param class_id
#' @param control
#' @param treatment
#' @param rank
#' @param option
#' @param ...
#' @param edgeR.n
#' @param edgeR.adjust.method
#' @param edgeR.sort.by
#'
#' @return
#' @importFrom dplyr %>% left_join
#' @examples
use_edgeR <- function(smrExpt, class_id, control, treatment, rank=FALSE,
                      edgeR.n = Inf, edgeR.adjust.method = "BH", edgeR.sort.by = "PValue",
                      option="GLM",...){ # "exact"
    checkNameSpace("edgeR")
    stopifnot("option must be either 'GLM' or 'exact'"=(option=="GLM" | option=="exact"))
    smrExpt <- smrExpt[,smrExpt[[class_id]] %in% c(control,treatment)]
    condition <- factor(as.character(SummarizedExperiment::colData(smrExpt)[,class_id]),
                        levels = c(control,treatment))

    ##Normalizing
    dt <- edgeR::calcNormFactors(smrExpt,method="TMM")
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

        command_str <- paste("limma::makeContrasts(",treatment , "-", control,",levels = design)"
                             , sep = "")
        contrast.matrix <- eval(parse(text=command_str))
        t12 <- edgeR::glmLRT(fit,contrast=contrast.matrix)
    }else if(option=="exact"){
        ##Estimate dispersion
        dt$samples$group <- condition
        d1 <- edgeR::estimateCommonDisp(dt,  verbose=T)
        d1 <- edgeR::estimateTagwiseDisp(d1)
        ##Compare groups (exact test)
        t12 <- edgeR::exactTest(d1, pair=c(1,2))
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

#' Title
#'
#' @param smrExpt
#' @param class_id
#' @param control
#' @param treatment
#' @param rank
#' @param ...
#'
#' @return
#' @export
#' @importFrom dplyr %>% left_join
#' @examples
 use_deseq2 <- function(smrExpt, class_id, control, treatment,rank=FALSE,...){
    checkNameSpace("DESeq2")
     if (SummarizedExperiment::assayNames(smrExpt)[1] != "counts") {
         nuOrder <-
             c("counts",
               setdiff(SummarizedExperiment::assayNames(smrExpt), "counts"))
         SummarizedExperiment::assays(smrExpt) <-
             SummarizedExperiment::assays(smrExpt)[nuOrder]
     }
     smrExpt <- smrExpt[,smrExpt[[class_id]] %in% c(control,treatment)]

    formula_str <- paste("formula( ~ ",class_id, ")", sep = "")
    dds <- DESeq2::DESeqDataSet(smrExpt, design = eval(parse(text=formula_str)))

    ## set reference to control, otherwise default is alphabetical order
    dds[[class_id]] <- factor(dds[[class_id]], levels=c(control,treatment))

    ##   1. estimate size factors
    ##   2. estimate dispersion
    ##   3. negative binomial GLM fitting and wald test
    dds_res <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds_res)
    # res$dispersion <- DESeq2::dispersions(dds_res)
    if(rank){
        return(res %>% as.data.frame() %>%
               dplyr::arrange(padj) %>%
               dplyr::mutate(rank = 1:dplyr::n()))
    }else{
        return(res)
    }
}

#' Title
#'
#' @param smrExpt
#' @param class_id
#' @param control
#' @param treatment
#' @param rank
#' @param ...
#'
#' @return
#' @export
#' @importFrom dplyr %>% left_join
#' @examples
use_DELocal <- function(smrExpt, class_id, control, treatment,rank=FALSE,...){
    checkNameSpace("DELocal")
    formula_str <- paste("formula( ~ ",class_id, ")", sep = "")
    if (SummarizedExperiment::assayNames(smrExpt)[1] != "counts") {
        nuOrder <-
            c("counts",
              setdiff(SummarizedExperiment::assayNames(smrExpt), "counts"))
        SummarizedExperiment::assays(smrExpt) <-
            SummarizedExperiment::assays(smrExpt)[nuOrder]
    }
    smrExpt <- smrExpt[,smrExpt[[class_id]] %in% c(control,treatment)]

    DELocal_result <- DELocal::DELocal(pSmrExpt = smrExpt, # Genes without neighbours are missing
                              nearest_neighbours = 5,pDesign =eval(parse(text=formula_str)),
                              pValue_cut = 1, pLogFold_cut = 0)
    if(rank){
        return(DELocal_result %>%
                   dplyr::arrange(desc(abs(relative.logFC))) %>%
                   dplyr::mutate(rank = 1:dplyr::n()))
    }else{
        return(DELocal_result)
    }
}



#' Title
#'
#' @param smrExpt
#' @param class_id
#' @param control
#' @param treatment
#' @param rank
#' @param ...
#'
#' @return
#' @export
#' @importFrom dplyr %>% left_join
#' @examples
use_EBSeq <- function(smrExpt, class_id, control, treatment, rank=FALSE,...){
    checkNameSpace("EBSeq")
    control_names <- smrExpt[,smrExpt[[class_id]]==control] %>% colnames()
    treatment_names <- smrExpt[,smrExpt[[class_id]]==treatment] %>% colnames()
    tabla <- SummarizedExperiment::assays(smrExpt)[["counts"]][,c(control_names,treatment_names)]
    Sizes=EBSeq::MedianNorm(tabla)

    EBOut=EBSeq::EBTest(Data=tabla,
                 Conditions=factor(c(rep(control,length(control_names)),
                                     rep(treatment,length(treatment_names))),
                                      levels = c(control,treatment)),
                 sizeFactors=Sizes, maxround=5)
    EBDERes=EBSeq::GetDEResults(EBOut, ...)
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
#' and output class is `Output` which are not widely used. We can use as(smrExpt, "ExpressionSet")
#' to get an eSet easily but then it will be hard to refer the treatment and control.
#' The order of factors influence the log fold change sign. To keep it comparable
#' to other methods the `NOISeq::readData()` is used internally.
#'
#' @param smrExpt Object of \code{\link{SummarizedExperiment}} class
#' @param class_id One of the columns of colData(smrExpt). It should be factors of more than one value.
#' @param control Base level and one of the factor values of `colData(smrExpt)[["condition"]]`
#' @param treatment one of the factor values of `colData(smrExpt)[["condition"]]`
#' @param rank Logical value default FALSE. If true the result will have an
#' additional column named "rank" which is ordered by ´prob´ values returned by
#' function NOISeq::\code{\link{noiseqbio}}.
#' @param ...  other arguments to be passed to main function NOISeq::\code{\link{noiseqbio}}.
#' The 'input' and 'factor' argument should not be used.
#'
#' @return A data.frame object from the "NOISeq::noiseqbio()@results[[1]]". For details check
#' the documentation of ´NOISeq´
#' @export
#' @importFrom dplyr %>% left_join
#' @examples
use_NOIseq <- function(smrExpt, class_id, control, treatment, rank=FALSE, ...){
    control_names <- smrExpt[,smrExpt[[class_id]]==control] %>% colnames()
    treatment_names <- smrExpt[,smrExpt[[class_id]]==treatment] %>% colnames()
    tabla <- SummarizedExperiment::assays(smrExpt)[["counts"]][,c(control_names,treatment_names)] %>% as.data.frame()

    expt_factors = data.frame(Condition = factor(c(
        rep(treatment, length(treatment_names)),
        rep(control, length(control_names))
    )))
    colnames(expt_factors) <- class_id
    expt_data <- NOISeq::readData(data = tabla, factors = expt_factors)

    ##Available normalization methods: norm = "rpkm", "uqua", "tmm2, "none"
    expt_noiseqbio = NOISeq::noiseqbio(input= expt_data, factor=class_id, ...)

    if(rank){
        expt_noiseqbio@results[[1]] %>%
            dplyr::arrange(desc(prob)) %>%
            dplyr::mutate(rank = 1:dplyr::n())
    }else{
        expt_noiseqbio@results[[1]]
    }
}

#' Title
#'
#' @param smrExpt
#' @param class_id
#' @param control
#' @param treatment
#' @param rank
#' @param ...
#'
#' @return
#' @export
#' @importFrom dplyr %>% left_join
#' @examples
use_SAMseq <- function(smrExpt, class_id, control, treatment,rank = FALSE,...){
    checkNameSpace("samr")
    control_names <- smrExpt[,smrExpt[[class_id]]==control] %>% colnames()
    treatment_names <- smrExpt[,smrExpt[[class_id]]==treatment] %>% colnames()
    tabla <- SummarizedExperiment::assays(smrExpt)[["counts"]][,c(control_names,treatment_names)] %>% as.data.frame()
    y <- factor(c(rep(control,length(control_names)),
                 rep(treatment,length(treatment_names))),
               levels = c(control,treatment))

    # it uses cat() for lots of text, annoying
    samfit <- samr::SAMseq(tabla, y, resp.type = "Two class unpaired",
                           geneid = rownames(tabla), genenames = rownames(tabla),
                           random.seed=1, fdr.output = 1)
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
    missing_id <- setdiff(rownames(smrExpt),res$GeneID)
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

#' Title
#'
#' @param deFun_list
#' @param return.df
#' @param smrExpt
#' @param class_id
#' @param control
#' @param treatment
#' @param ...
#'
#' @importFrom purrr reduce
#' @return
#' @export
#'
#' @examples
use_multDE <- function(deFun_list, return.df= FALSE ,smrExpt , class_id , control , treatment , ... ) {
    wrap_it <- function( f, ...) {
        message(paste("Now executing >> ",f))
        message("####")
        deFun_list[[f]](... )
    }
    SummarizedExperiment::colData(smrExpt)[,class_id] %>% summary() #%>% message()

    y <- lapply(names(deFun_list), wrap_it, smrExpt = smrExpt, class_id = class_id,
                control = control, treatment = treatment, ...)

    named_result <- list()
    for (i in 1:length(y)) {
        named_result[[names(deFun_list)[i]]] <- y[[i]]
    }
    y <- named_result

    if(return.df){
        for (prefix in names(y)) {
            colnames(y[[prefix]]) <- paste(prefix,colnames(y[[prefix]]), sep = "_")
        }

        y <- purrr::reduce(lapply(y, function(x) data.frame(x, rn = row.names(x))),
                           merge, all=TRUE)
        gene_info <- as.data.frame(SummarizedExperiment::rowData(smrExpt))
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

#' Title
#'
#' @param df
#' @param pValName
#' @param lFCName
#' @param sigThreshold
#' @param logFCThreshold
#' @param labelName
#' @param selectedLabel
#' @param palette one of "npg" ,"aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "nejm" or similar to viridis::cividis(3)
#'
#' @return
#' @export
#' @importFrom dplyr %>% left_join if_else
#' @importFrom ggpubr ggscatter
#' @examples
volcanoPlot <- function(df,pValName,lFCName, sigThreshold=0.05, logFCThreshold = 1,
                        labelName=NULL, selectedLabel = NULL, palette = "nejm"){
    df <- df %>% dplyr::mutate(padj=-log10(!!sym(pValName)),
                              Significant = if_else((!!sym(pValName) < sigThreshold & !!sym(lFCName) > logFCThreshold),"UP",
                                                    if_else((!!sym(pValName) < sigThreshold & !!sym(lFCName) < -logFCThreshold),"DOWN",
                                                          "Not",missing="Not"),missing="Not"
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
