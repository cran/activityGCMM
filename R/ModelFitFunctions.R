#####################################################
#' Calculate sum of absolute circular residuals
#' @description Calculate posterior probability distribution of summed absolute circular residuals for assessing GCMM model fit
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function; residuals must be saved when running the \code{\link{GCMM}} function using \code{saveResids=TRUE}
#' @return Returns list with output summary with mean and 95% HDI and posterior distribution of summed absolute circular residuals
#' @examples
#' \donttest{ FoxGCMMresids<-GCMM(data=redfoxsample$Radians, RE1=redfoxsample$SamplingPeriod, saveResids=TRUE,
#'     scale=c("2pi"), family="vonmises", autorun=FALSE, adapt=0, sample=300, burnin=300, thin=1) 
#'     FoxResids<-sumCircResids(FoxGCMMresids) }
#' @export

	sumCircResids<-function(model) {
		MCMC<-combineMCMC(model)
		names<-c()
		for (i in 1:length(model$data)) { n<-paste("Resid[",i,"]",sep=""); names<-c(names,n) }
		k<-length(model$Mu)+length(model$c)+length(model$P)+model$Nclust+as.numeric(length(model$RE2)>0)*model$Nclust
			RawRes<-MCMC[,((k+1):(ncol(MCMC)))]
				for (i in 1:nrow(RawRes)) { 
					for (j in 1:ncol(RawRes)) {
						if (RawRes[i,j]>pi) { RawRes[i,j]<-RawRes[i,j]-(2*pi) }
						if (RawRes[i,j]<pi*-1) { RawRes[i,j]<-RawRes[i,j]+(2*pi) }
				}	}

		SumAbs<-rowSums(abs(RawRes))
		SumAbsHDI<-HDI(SumAbs)
				print(SumAbsHDI)

		out<-list(summary=SumAbsHDI,PD=list(SumCircResids=SumAbs))
		class(out)<-"GCMMestimate"
		return(out)
	}



############################################
#' Compare fit of GCMM models based on circular residuals
#' @description Compare fit of GCMM models by comparing posterior distributions of the summed circular residuals
#' @param model1 Object of class \code{GCMM} with output from \code{\link{GCMM}} function to compare with model2; residuals must be saved when running the \code{\link{GCMM}} function using \code{saveResids=TRUE}
#' @param model2 Object of class \code{GCMM} with output from \code{\link{GCMM}} function to compare with model1; residuals must be saved when running the \code{\link{GCMM}} function using \code{saveResids=TRUE}
#' @param sample Number of posterior samples; default=10000
#' @return Returns object of class  \code{GCMMestimate} with list of output
#' 
#' @examples
#' \donttest{ 
#' FoxVMGCMM<-GCMM(data=redfoxsample$Radians, RE1=redfoxsample$CameraTrapID, family="vonmises", 
#'     saveResids=TRUE, scale=c("2pi"), autorun=FALSE, adapt=0, sample=1000, burnin=500, thin=1) 
#' FoxWCGCMM<-GCMM(data=redfoxsample$Radians, RE1=redfoxsample$CameraTrapID, family="wrappedcauchy",
#'     saveResids=TRUE, scale=c("2pi"), autorun=FALSE, adapt=0, sample=1000, burnin=500, thin=1) 
#' FoxModelCompare<-compareGCMMfit(FoxVMGCMM, FoxWCGCMM)  }
#' 
#' @export

	compareGCMMfit<-function(model1, model2,sample=10000) {

	message("");message("Sum of circular residuals for model 1:")
		E1<-sumCircResids(model1)
	message("");message("Sum of circular residuals for model 2:")
		E2<-sumCircResids(model2)
	message("")
		comp<-compareGCMM(model1=E1,p1="SumCircResids",model2=E2,p2="SumCircResids",sample=sample,plot=FALSE)

		graphics::par(mfrow=c(3,1),mar=c(4.1,4.1,2,2))
		graphics::hist(comp$PD$PD1,breaks=100, prob=TRUE, main="",xlab="Model1 SumCircResids Estimate", col="aquamarine3",border="lightseagreen")
			HDI1<-HDI(comp$PD$PD1)
			graphics::segments(x0=HDI1[2], y0=-.001, x1=HDI1[3], y1=-.001, lwd=5,col="grey15"); graphics::box(lwd=2)
		graphics::hist(comp$PD$PD2,breaks=100,prob=TRUE,main="", xlab="Model2 SumCircResids Estimate", col="aquamarine3",border=col)
			HDI2<-HDI(comp$PD$PD2)
			graphics::segments(x0=HDI2[2], y0=-.001, x1=HDI2[3], y1=-.001, lwd=5,col="grey15"); graphics::box(lwd=2)
		graphics::hist(comp$PD$Difference,breaks=100,prob=TRUE,main="",xlab="Difference Estimate",col="lightseagreen",border="cyan3")
			HDIdiff<-HDI(comp$PD$Difference)
			graphics::segments(x0=HDIdiff[2], y0=.0005, x1=HDIdiff[3], y1=.0005, lwd=5,col="red2");graphics::box(lwd=2)
			graphics::abline(v=0,lwd=5,col="red3")			

	return(comp)
	}


############################################
#' Posterior predictive check of GCMM model
#' @description Conduct posterior predictive check (PPC) by simulating data from fitted GCMM model and plotting against observed data
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function; if YExp and clustIDs are not provided as vectors, the GCMM model
#'     must contain this information using the arguments \code{saveYExp=TRUE} and \code{saveclustIDs=TRUE}
#' @param YExp Vector of YExp values from GCMM function; see also \code{\link{GCMM}}
#' @param clustID Vector of clustID values from GCMM function; see also \code{\link{GCMM}}
#' @return Returns vector of simulated values and prints plot of simulated and raw values
#' @examples \donttest{ 
#'   FoxGCMMPPC<-GCMM(data=redfoxsample$Radians, RE1=redfoxsample$CameraTrapID, family="vonmises", 
#'     saveclustIDs=TRUE, saveYExp=TRUE,
#'     scale=c("2pi"), autorun=FALSE, adapt=0, sample=300, burnin=300, thin=1) 
#'   FoxPPC<-GCMMppc(FoxGCMMPPC)  }
#' @export

	GCMMppc<-function(model,YExp=NULL,clustID=NULL) {

	if (length(YExp)==0) { YExp<-as.numeric(model$YExp[,4]) }
	if (length(clustID)==0) { clustID<-model$clustIDs }
	df<-data.frame(YExp,clustID); Yexp<-list()
		for (c in 1:model$Nclust) { Yexp[[c]]<-df[which(df$clustID==c),1] }

	YNew<-c()
	for (c in 1:model$Nclust) {
		for (i in 1:length(Yexp[[c]]) ) {
			if (model$family=="wrappedcauchy") {
				s<-circular::rwrappedcauchy(1, circular::circular(Yexp[[c]][i]), as.numeric(model$c[c]))  
				YNew<-c(YNew,s)
				}
		if (model$family=="vonmises") { 
				s<-circular::rvonmises(1, circular::circular(Yexp[[c]][i]), as.numeric(model$c[c])) 
				YNew<-c(YNew,s)
				}
			}
		}

	YNew<-convertRad(YNew,to=model$scale)
	graphics::plot(model$data)
		graphics::points(YNew, col="red")
		
	return(YNew)

	}



############################################
