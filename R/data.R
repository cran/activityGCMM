#' Sample data of camera trap observations of humans
#' 
#'
#' @description Example dataset for fitting circular mixed effect mixture models with activityGCMM package
#' 
#' @name humanssample
#' @docType data
#' @title Sample data of camera trap observations of humans
#' 
#'
#'
#' @format Dataframes with 3 variables
#'    Radians Time of observations, in radians (0 to 2pi)
#'    CameraTrapID Variable identifying camera traps
#'    SamplingPeriod Variable identifying sampling period during which camera traps were recording
#'
#' 
#' @source \ Campbell L.A.D. 2017
#' 
#' @keywords datasets
#' 
#' @examples
#' data(humanssample)
#' \dontrun{ GCMM(data=humanssample$Radians, RE1=humanssample$SamplingPeriod, 
#'     scale=c("2pi"), family="vonmises", autojags=TRUE, thin=3) }
#'
"humanssample"

#' Sample data of camera trap observations of humans
#' 
#'
#' @description Example dataset for fitting circular mixed effect mixture models with activityGCMM package
#' 
#' @name redfoxsample
#' @docType data
#' @title Sample data of camera trap observations of red fox
#' 
#'
#'
#' @format Dataframes with 3 variables
#'    Radians Time of observations, in radians (0 to 2pi)
#'    CameraTrapID Variable identifying camera traps
#'    SamplingPeriod Variable identifying sampling period during which camera traps were recording
#'
#' 
#' @source \ Campbell L.A.D. 2017
#' 
#' @keywords datasets
#'
#' @examples
#' data(redfoxsample)
#' \dontrun{ GCMM(data=redfoxsample$Radians, RE1=redfoxsample$SamplingPeriod, 
#'     scale=c("2pi"), family="vonmises", autojags=FALSE,
#'     adapt=0, sample=300, burnin=300, thin=1, n.chains=2  ) }
#'
"redfoxsample"
