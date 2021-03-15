#' Global variables
#' @title Global variables
#' @description List of variables
#' @details Variables List of variables
#' @name Variables
#' @return Summary of model output, depending on printvars & mixture plot of components; also saves to global environment the model intercepts, concentration parameters, mixture component weights, model as a runjags object ("GCMM") which can be extended with extend.jags or monitored parameters added/dropped
#' @return out Model output
#' @return Mu Component means
#' @return P Component weights
#' @return K Component concentration
#' @return GCMMcomponents Separate mixture components
#' @return GCMMmixture Mixture
	utils::globalVariables(names(c("out","Mu","P","K","GCMMcomponents","GCMMmixture")))#,"mixtureplot")))
