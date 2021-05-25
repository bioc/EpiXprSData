\name{EpiXprSData-package}
\alias{EpiXprSData-package}
\alias{EpiXprSData}
\title{Tissue Specific Gene Expression Prediction Models Constructed
    From Paired DNA Methylation and RNA Counts}
\description{
    Models constructed on paired DNA methylation and gene expression
    counts for the purpose of predicting expression in samples with only DNA
    methylation. This data package is designed to be used in conjunction with
    the \code{EpiXprS} Bioconductor package. The data are available in
    \code{ExperimentHub}.
}
\details{
    See the vignette for examples of using this data in the \code{EpiXprS}
    Bioconductor package.

    browseVignettes("EpiXprSData")

    A brief description for how the models were constructed are presented in the
    scripts/ directory of the source package. For a full description, see the
    \code{EpiXprS} package.
}
\examples{
library(ExperimentHub)
hub <- ExperimentHub()
x <- query(hub, c("EpiXprSData","PRAD"))

}
