% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{use_edgeR_GLM}
\alias{use_edgeR_GLM}
\alias{use_edgeR_exact}
\alias{use_edgeR}
\title{To use SummarizedExperiment with package edgeR}
\usage{
use_edgeR_GLM(se, colData_id, control, treatment, rank = FALSE, ...)

use_edgeR_exact(se, colData_id, control, treatment, rank = FALSE, ...)

use_edgeR(
  se,
  colData_id,
  control,
  treatment,
  rank = FALSE,
  edgeR.n = Inf,
  edgeR.adjust.method = "BH",
  edgeR.sort.by = "PValue",
  option = "GLM",
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

\item{...}{other arguments to be passed to edgeR::\code{\link{glmLRT}}
or edgeR::\code{\link{exactTest}}}

\item{edgeR.n}{argument for edgeR::\code{\link{topTags}}}

\item{edgeR.adjust.method}{argument for edgeR::\code{\link{topTags}}}

\item{edgeR.sort.by}{argument for edgeR::\code{\link{topTags}}}

\item{option}{"GLM" or "exact" to indicate to use either edgeR::\code{\link{glmLRT}}
or edgeR::\code{\link{exactTest}}}
}
\value{
a data.frame of output from edgeR::\code{\link{topTags}}
}
\description{
A wrapper function of DESeq2 where input is an object of \code{\link{SummarizedExperiment}}
}
\examples{

se <- readRDS(system.file("extdata",
        "rat_vole_mouseSE_salmon.rds",
        package = "broadSeq"))

# To reduce runtime
se <- se[rowData(se)$chromosome_name == 2,colData(se)$species == "Mouse"]

result <-
    use_edgeR(se = se,
           colData_id = "stage", control = "Bud", treatment = "Cap",
           rank = TRUE)
}
