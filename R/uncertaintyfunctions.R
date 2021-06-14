###############################################
#' Plot GCMM Activity Curve Posterior Samples
#' 
#' @description Plot GCMM activity curve posterior samples for visualizing estimate uncertatinty
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param sample Number of posterior samples to plot; default=100 
#' @param scale Scale for which to plot the activity curve, either "pi" for -pi, pi or "2pi" for 0, 2pi; default is that which is recommended by the GCMM function 
#' @param ymax Value to use as upper limit for y-axis
#' @param plotmean Logical argument for whether to plot activity curve from posterior distribution mean; default=TRUE
#' @param RGB Vector of RBG values for line colour
#' @param alpha Value for line transparency, between 0 (completely transparent) to 1 (completely opaque); default=0.05
#' @param axisunits Scale to use for the xaxis, either "radians", "time", "solar", "sun", or "non"; default="radians"
#' @param lines Whether to include lines on the graph for the x axis labels; default=TRUE
#' @param cex.axis Font size for axis labels 
#' @return No return value; prints plot of activity curve posterior samples
#' @export


plotGCMMsamples<-function(model, sample=100, scale="NULL", ymax="NULL", plotmean=TRUE, RGB=c(200,200,200), alpha=0.05,
		axisunits="radians",lines=TRUE,cex.axis=0.8) {

	if (scale=="NULL") { scale<-model$scale }

	MCMC<-combineMCMC(model)
	MCMCsample<-samplerows(df=MCMC,n=sample)

	if (model$Nclust==1) {
		Mu1<-c(MCMCsample[,1])
		K1<-c(MCMCsample[,2])
		P1<-c(MCMCsample[,4])
		}
	if (model$Nclust==2) {
		Mu1<-c(MCMCsample[,1])
		Mu2<-c(MCMCsample[,2])
		K1<-c(MCMCsample[,3])
		K2<-c(MCMCsample[,4])
		P1<-c(MCMCsample[,7])
		P2<-c(MCMCsample[,8])
		}
	if (model$Nclust==3) {
		Mu1<-c(MCMCsample[,1])
		Mu2<-c(MCMCsample[,2])
		Mu3<-c(MCMCsample[,3])
		K1<-c(MCMCsample[,4])
		K2<-c(MCMCsample[,5])
		K3<-c(MCMCsample[,6])
		P1<-c(MCMCsample[,10])
		P2<-c(MCMCsample[,11])
		P3<-c(MCMCsample[,12])
		}
	if (model$Nclust==4) {
		Mu1<-c(MCMCsample[,1])
		Mu2<-c(MCMCsample[,2])
		Mu3<-c(MCMCsample[,3])
		Mu4<-c(MCMCsample[,4])
		K1<-c(MCMCsample[,5])
		K2<-c(MCMCsample[,6])
		K3<-c(MCMCsample[,7])
		K4<-c(MCMCsample[,8])
		P1<-c(MCMCsample[,13])
		P2<-c(MCMCsample[,14])
		P3<-c(MCMCsample[,15])
		P4<-c(MCMCsample[,16])
		}

	if (ymax=="NULL") { ymax<-yMax(list(model),type="mixture")*1.5 }


	graphics::par(mar=c(4.2,4.2,1,3))
		if (scale=="2pi") { xlim1<-0; xlim2<-2*pi; xscale="noon" }
		if (scale=="pi") { xlim1<- -pi; xlim2<-pi; xscale="midnight" }
		col <- grDevices::rgb(RGB[1],RGB[2],RGB[3], max = 255, alpha = alpha*255)

	 overlap::densityPlot(as.numeric(model$GCMMmixture),lwd=3,xaxt="n",axt="n",xlim=c(xlim1,xlim2),cex=.8,xcenter=xscale,
			main="",ylim=c(-.02,ymax),adjust=5,col="white",extend=NA,xscale=NA,xaxs="i",yaxs="i")
			xaxis(axisunits, lines, cex.axis)			
			graphics::par(new=TRUE) 


	for (s in 1:sample) {

	    if (model$Nclust==1) {
			if (model$family=="vonmises") { 
			gcmmcomponentsC1<-circular::rvonmises(as.numeric(P1[s])*10000, circular::circular(Mu1[s]), K1[s])   
				gcmmmixture<-c(gcmmcomponentsC1)
					}
			if (model$family=="wrappedcauchy") { 
			gcmmcomponentsC1<-circular::rwrappedcauchy(as.numeric(P1[s])*10000, circular::circular(Mu1[s]), K1[s])  
				gcmmmixture<-c(gcmmcomponentsC1)
					}
				}
	    if (model$Nclust==2) {
			if (model$family=="vonmises") { 
			gcmmcomponentsC1<-circular::rvonmises(as.numeric(P1[s])*10000, circular::circular(Mu1[s]), K1[s])  
			gcmmcomponentsC2<-circular::rvonmises(as.numeric(P2[s])*10000, circular::circular(Mu2[s]), K2[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2)
					}
			if (model$family=="wrappedcauchy") { 
			gcmmcomponentsC1<-circular::rwrappedcauchy(as.numeric(P1[s])*10000, circular::circular(Mu1[s]), K1[s])  
			gcmmcomponentsC2<-circular::rwrappedcauchy(as.numeric(P2[s])*10000, circular::circular(Mu2[s]), K2[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2)
					}
				}
	    if (model$Nclust==3) {
			if (model$family=="vonmises") { 
			gcmmcomponentsC1<-circular::rvonmises(as.numeric(P1[s])*10000, circular::circular(Mu1[s]), K1[s])  
			gcmmcomponentsC2<-circular::rvonmises(as.numeric(P2[s])*10000, circular::circular(Mu2[s]), K2[s])  
			gcmmcomponentsC3<-circular::rvonmises(as.numeric(P3[s])*10000, circular::circular(Mu3[s]), K3[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2,gcmmcomponentsC3)
					}
			if (model$family=="wrappedcauchy") { 
			gcmmcomponentsC1<-circular::rwrappedcauchy(as.numeric(P1[s])*10000, circular::circular(Mu1[s]), K1[s])  
			gcmmcomponentsC2<-circular::rwrappedcauchy(as.numeric(P2[s])*10000, circular::circular(Mu2[s]), K2[s])  
			gcmmcomponentsC3<-circular::rwrappedcauchy(as.numeric(P3[s])*10000, circular::circular(Mu3[s]), K3[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2,gcmmcomponentsC3)
					}
				}
	    if (model$Nclust==4) {
			if (model$family=="vonmises") { 
			gcmmcomponentsC1<-circular::rvonmises(as.numeric(P1[s])*10000, circular::circular(Mu1[s]), K1[s])  
			gcmmcomponentsC2<-circular::rvonmises(as.numeric(P2[s])*10000, circular::circular(Mu2[s]), K2[s])  
			gcmmcomponentsC3<-circular::rvonmises(as.numeric(P3[s])*10000, circular::circular(Mu3[s]), K3[s])  
			gcmmcomponentsC4<-circular::rvonmises(as.numeric(P4[s])*10000, circular::circular(Mu4[s]), K4[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2,gcmmcomponentsC3,gcmmcomponentsC4)
					}
			if (model$family=="wrappedcauchy") { 
			gcmmcomponentsC1<-circular::rwrappedcauchy(as.numeric(P1[s])*10000, circular::circular(Mu1[s]), K1[s])  
			gcmmcomponentsC2<-circular::rwrappedcauchy(as.numeric(P2[s])*10000, circular::circular(Mu2[s]), K2[s])  
			gcmmcomponentsC3<-circular::rwrappedcauchy(as.numeric(P3[s])*10000, circular::circular(Mu3[s]), K3[s])  
			gcmmcomponentsC4<-circular::rwrappedcauchy(as.numeric(P4[s])*10000, circular::circular(Mu4[s]), K4[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2,gcmmcomponentsC3,gcmmcomponentsC4)
					}
				}

	graphics::par(new=TRUE)
	overlap::densityPlot(as.numeric(gcmmmixture),lwd=2,xaxt="n",axt="n",xlim=c(xlim1,xlim2),xcenter=xscale,
			main="",ylim=c(-.02,ymax),adjust=5,col=col,extend=NA,xscale=NA,xaxs="i",yaxs="i") 

	progress(s=s, sample=sample)
	  }

	if (plotmean==TRUE) {
		graphics::par(new=TRUE)
		overlap::densityPlot(as.numeric(model$GCMMmixture),lwd=3,xaxt="n",axt="n",xlim=c(xlim1,xlim2),xcenter=xscale,
			main="",ylim=c(-.02,ymax),adjust=5,col="grey15",extend=NA,xscale=NA,xaxs="i",yaxs="i") 
			}

		datrug<-convertRad(model$data,to=scale)
		graphics::abline(h=0,col="white"); graphics::rug(datrug,lwd=2)
		graphics::box(lwd=2)

	}

###########################################################
#' Plot histogram of posterior distribution
#' 
#' @description Plot histogram of posterior samples 
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param param Parameter for which to plot the posterior samples
#' @param col Colour of histogram
#' @return Returns matrix with the mean and HDI of activity probability density estimated from both GCMM models at the peak activity time of the other. Posterior distributions of activity probability density and peak activity times for both GCMM models are also saved.
#' @export

	posteriorhistplot<-function(model, param, col="cyan4") {

	PD<-extractparam(model,x=param)
	hdi<-HDI(PD)
	graphics::hist(PD,breaks=100,prob=TRUE, main="",xlab="Estimate",col=col,border=col)
		graphics::segments(x0=hdi[2], y0=-.001, x1=hdi[3], y1=-.001, lwd=7,col="grey15")
		graphics::abline(v=hdi[1],lwd=3,col="grey15")
		graphics::box(lwd=2)
	} 
