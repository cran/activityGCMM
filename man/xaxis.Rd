% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supportfunctions.R
\name{xaxis}
\alias{xaxis}
\title{Axis labels for temporal activity plots}
\usage{
xaxis(
  axisunits = c("radians", "solar", "sun", "time", "none"),
  lines = TRUE,
  cex.axis = 0.8
)
}
\arguments{
\item{axisunits}{Scale to use for the xaxis, either "radians", "time", "solar", "sun", or "none"; default="radians"}

\item{lines}{Whether to include lines on the graph for the x axis labels; default=TRUE}

\item{cex.axis}{Font size for axis labels}
}
\value{
Prints axis
}
\description{
Support function for xaxis labels for graphing temporal activity curves
}
