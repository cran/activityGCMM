% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelFitFunctions.R
\name{compareGCMMfit}
\alias{compareGCMMfit}
\title{Compare fit of GCMM models based on circular residuals}
\usage{
compareGCMMfit(model1, model2, sample = 10000)
}
\arguments{
\item{model1}{Object of class \code{GCMM} with output from \code{\link{GCMM}} function to compare with model2; residuals must be saved when running the \code{\link{GCMM}} function using \code{saveResids=TRUE}}

\item{model2}{Object of class \code{GCMM} with output from \code{\link{GCMM}} function to compare with model1; residuals must be saved when running the \code{\link{GCMM}} function using \code{saveResids=TRUE}}

\item{sample}{Number of posterior samples; default=10000}
}
\value{
Returns object of class  \code{GCMMestimate} with list of output
}
\description{
Compare fit of GCMM models by comparing posterior distributions of the summed circular residuals
}
\examples{
\donttest{ 
FoxVMGCMM<-GCMM(data=redfoxsample$Radians, RE1=redfoxsample$CameraTrapID, family="vonmises", 
    saveResids=TRUE, scale=c("2pi"), autorun=FALSE, adapt=0, sample=1000, burnin=500, thin=1) 
FoxWCGCMM<-GCMM(data=redfoxsample$Radians, RE1=redfoxsample$CameraTrapID, family="wrappedcauchy",
    saveResids=TRUE, scale=c("2pi"), autorun=FALSE, adapt=0, sample=1000, burnin=500, thin=1) 
FoxModelCompare<-compareGCMMfit(FoxVMGCMM, FoxWCGCMM)  }

}
