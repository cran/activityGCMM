% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/APDfunctions.R
\name{APDatPeak}
\alias{APDatPeak}
\title{Activity Probability Density at Peak Time of Another}
\usage{
APDatPeak(model1, model2, sample = 1000, HDIprob = 0.95)
}
\arguments{
\item{model1}{Object of class \code{GCMM} with output from \code{\link{GCMM}} function}

\item{model2}{Object of class \code{GCMM} with output from \code{\link{GCMM}} function}

\item{sample}{Number of posterior samples from which to build the HDIs}

\item{HDIprob}{Value for probability mass to use for HDI; default=95\%}
}
\value{
Returns matrix with the mean and HDI of activity probability density estimated from both GCMM models at the peak activity time of the other. Posterior distributions of activity probability density and peak activity times for both GCMM models are also saved.
}
\description{
Calculates activity probability density from one GCMM model at the peak activity time of a second GCMM model
}
