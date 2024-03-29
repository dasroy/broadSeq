% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{use_limma_trend}
\alias{use_limma_trend}
\alias{use_limma_voom}
\alias{use_limma}
\title{To use SummarizedExperiment with package limma}
\usage{
use_limma_trend(se, colData_id, control, treatment, rank = FALSE, ...)

use_limma_voom(se, colData_id, control, treatment, rank = FALSE, ...)

use_limma(
  se,
  colData_id,
  control,
  treatment,
  rank = FALSE,
  useVoom = TRUE,
  showPlot = FALSE,
  limma.adjust = "BH",
  limma.sort.by = "p",
  limma.number = Inf,
  ...
)
}
\arguments{
\item{se}{Object of \code{\link{SummarizedExperiment}} class}

\item{colData_id}{One of the columns of colData(se). It should be factors of more than one value.}

\item{control}{Base level and one of the factor values of \code{colData(se)[[colData_id]]}}

\item{treatment}{one of the factor values of \code{colData(se)[[colData_id]]}}

\item{rank}{Logical value default FALSE. If true the result will have an
additional column named "rank"}

\item{...}{other arguments to be passed to main function edgeR::calcNormFactors .}

\item{useVoom}{whether to use limma::voom or edgeR::cpm}

\item{showPlot}{whether to use limma::plotSA ; default FALSE}

\item{limma.adjust}{argument for limma::topTable}

\item{limma.sort.by}{argument for limma::topTable}

\item{limma.number}{argument for limma::topTable}
}
\value{
a data.frame of output from limma::topTable
}
\description{
A wrapper function of limma where input is an object of \code{\link{SummarizedExperiment}}
}
\examples{

se <- readRDS(system.file("extdata",
        "rat_vole_mouseSE_salmon.rds",
        package = "broadSeq"))

# To reduce runtime
se <- se[rowData(se)$chromosome_name == 2,colData(se)$species == "Mouse"]
result <-
    use_limma(se = se,
           colData_id = "stage", control = "Bud", treatment = "Cap",
           rank = TRUE)
}
