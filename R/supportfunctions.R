#' Convert Radians Scale
#' 
#' @description Converts scale of observations in radians to (0,2pi) or (-pi,pi) or hours
#' @param x Vector of observations in radians
#' @param to Character vector to specify conversion, either "2pi" for (0, 2pi), "pi" for (-pi, pi), or "hours" for hours
#' @return Returns vector of observations on desired scale
#' @examples Rad2<-convertRad(redfoxsample$Radians,to="2pi")
#' @export
	convertRad<-function(x, to) {
	y<-x
		if(to=="2pi") {
			if (min(y)<0) { repeat { 
				if (min(y)>=0) { break }
				for(i in 1:length(y)) { if(y[i]<0) { y[i]<-y[i]+2*pi } } 
			}	}
			if (max(y)>2*pi) { repeat { 
				if (max(y)<=2*pi) { break }
				for(i in 1:length(y)) { if(y[i]>2*pi) { y[i]<-y[i]-2*pi } } 
			}	}
		}

		if(to=="pi") {
			if (min(y)<(-1*pi)) { repeat { 
				if (min(y)>=(-1*pi)) { break }
				for(i in 1:length(y)) { if(y[i]<(-1*pi)) { y[i]<-y[i]+2*pi } } 
			}	}
			if (max(y)>pi) { repeat { 
				if (max(y)<=pi) { break }
				for(i in 1:length(y)) { if(y[i]>pi) { y[i]<-y[i]-2*pi } } 
			}	}
		}

		if(to=="hours") {
		# first converting to 2pi:
			if (min(y)<0) { repeat { 
				if (min(y)>=0) { break }
				for(i in 1:length(y)) { if(y[i]<0) { y[i]<-y[i]+2*pi } } 
			}	}
			if (max(y)>2*pi) { repeat { 
				if (max(y)<=2*pi) { break }
				for(i in 1:length(y)) { if(y[i]>2*pi) { y[i]<-y[i]-2*pi } } 
			}	}
		# then converting to hours:
		y<-y/(2*pi)*24
		}

	return(y)
	}


###################################################################
#' Calculate highest density interval
#' 
#' @description Calculates the highest density interval
#' @param x Vector of data
#' @param prob Value for probability mass of HDI; default=95%
#' @return Returns matrix of the mean and upper and lower bounds of the HDI
#' @export


	HDI<-function(x, prob=.95) {

	func<-stats::approxfun(density(x))
	density<-data.frame(x=x,y=func(x))	
	density<-density[order(-density$y),]
	int<-density[1:(prob*length(x)),]

	low<-paste("Lower",prob*100,sep="")
	high<-paste("Upper",prob*100,sep="")

	out<-round(c(mean(x),min(int$x),max(int$x)),3)
		names(out)<-c("Mean",low,high)

	return(out)
	}


###################################################################
#' Mode
#' 
#' @description Returns the mode of a vector
#' @param x Vector of data
#' @return Returns the mode of x
#' @export

	mode<-function(x) { 
		t<-table(x)
		n<-names(t)
			for(i in 1:length(t)) { if (max(t)==t[i]) { m<-as.numeric(n[i]) } } 
		return(m)
		}


###################################################################

#' Posterior distribution summaries and support
#' 
#' @description Calculates the proportion > 0, < 0, and posterior distribution support (PDS)
#' @param x Vector of data
#' @return Returns matrix with data summary
#' @export

	PDS<-function(x) {

	output<-c(sum(as.numeric(x>0))/length(x), sum(as.numeric(x<0))/length(x))
	max<-max(output)
	output<-c(output, max)
		output<-round(output,3)
		names(output)<-c("PD>0","PD<0","PDS")
		return(output)
	}

###################################################################

#' Calculate y-axis limit for plotting multiple activity curves
#' 
#' @description Identifies maximum probability density for multiple activity curves to select y-axis limit when plotting multiple curves
#' @param models List of objects of class GCMM containing output from GCMM function
#' @param type Identifier for whether to use maximum probability density from GCMM mixture, using "mixture", or components, using "components"; default is to use the GCMM mixture density
#' @return Returns a value of the maximum probability density plus buffer space to be used as the y-axis limit in activity curve plots
#' @export

	yMax<-function(models,type="mixture") {

	max<-c()
	if (type=="mixture") {
		for (i in 1:length(models) ) { 
			APD<-calcAPD(x=seq(min(models[[i]]$GCMMmixture),max(models[[i]]$GCMMmixture),length=100),
				curve=models[[i]]$GCMMmixture)
			max<-c(max,max(APD))
		}	}

	if (type=="components") {
		for (i in 1:length(models) ) {
			for (c in 1:models[[i]]$Nclust)  { 

				APD<-calcAPD(x=seq(min(models[[i]]$GCMMcomponents[[c]]),max(models[[i]]$GCMMcomponents[[c]]),length=100),
					curve=models[[i]]$GCMMcomponents[[c]])
				max<-c(max,max(APD))
		}	}	}

	return(max(max)*1.4)
	}



###################################################################
#' Axis labels for temporal activity plots
#' 
#' @description Support function for xaxis labels for graphing temporal activity curves
#' @param axisunits Scale to use for the xaxis, either "radians", "time", "solar", "sun", or "none"; default="radians"
#' @param lines Whether to include lines on the graph for the x axis labels; default=TRUE
#' @param cex.axis Font size for axis labels 
#' @return Prints axis
#' @export

	xaxis<-function(axisunits=c("radians","solar","sun","time","none"), lines=TRUE, cex.axis=0.8) {

		if(axisunits=="time") { labels<-c("12:00","18:00","24:00","6:00","12:00","18:00","24:00") }
		if(axisunits=="solar") { labels<-c("solar noon","sunset","solar midnight","sunrise","solar noon","sunset","solar midnight") }
		if(axisunits=="sun") { labels<-c("noon","sunset","midnight","sunrise","noon","sunset","midnight") }
		if(axisunits=="radians") { labels<-c("-pi","-3pi/2","0","pi/2","pi","3pi/2","2pi") }
		if(axisunits=="none") { labels<-c("","","","","","","") }

		if (lines==TRUE) { graphics::abline(v=c(-pi, -pi/2, 0, pi/2, pi, 3*pi/2, 2*pi),col="grey55", lty=2) }
			graphics::axis(side=1, at=c(-pi, -pi/2, 0, pi/2, pi, 3*pi/2, 2*pi), 
				labels=c(labels), tck=-0.02,col="grey15",cex.axis=cex.axis)
	}

###################################################################
#' Axis labels for circular temporal data plots
#' 
#' @description Support function for axis labels for circular plots
#' @param axisunits Scale to use for the xaxis, either "radians", "time", "sun", or "none"; default="radians"
#' @return Prints axis 
#' @export

	circaxis<-function(axisunits=c("radians","time","sun","none")) {

	if(axisunits[1]=="radians") { labels<-c(0,"pi/2","pi","3pi/2",0.6) }
	if(axisunits[1]=="time") { labels<-c("24:00","6:00","12:00","18:00",0.5) }
	if(axisunits[1]=="sun") { labels<-c("midnight","sunrise","noon","sunset",0.4) }
	if(axisunits[1]=="none") { labels<-c("","","","",0.4) }

		graphics::text(x=0, y=.7, labels[1], cex=as.numeric(labels[5]))
		graphics::text(x=.7, y=0, labels[2], cex=as.numeric(labels[5]))
		graphics::text(x=0, y=-.7, labels[3], cex=as.numeric(labels[5]))
		graphics::text(x=-.7, y=.0, labels[4], cex=as.numeric(labels[5]))

	}


###################################################################

#' Combine MCMC chains for posterior simulations
#' 
#' @description Support function that extracts MCMC chains for creating posterior simulations of activity curves
#' @param model Object of class GCMM containing output from GCMM function
#' @return Returns a list of MCMC chains
#' @export

		combineMCMC<-function(model) { 

	return( rbind(model$runjags$mcmc[[1]],model$runjags$mcmc[[2]],model$runjags$mcmc[[3]] ) ) }


###################################################################

#' Extract parameters for posterior simulations
#' 
#' @description Support function that extracts parameter estimates for creating posterior simulations of activity curves
#' @param model Object of class GCMM containing output from GCMM function
#' @param x Name of parameter to be extracted
#' @return Returns posterior samples of the parameter
#' @export

	extractparam<-function(model,x) {

	if(class(model)[1]=="GCMM") {
		MCMC<-combineMCMC(model)
		posterior<-c()
			if(sum(as.numeric(x==colnames(MCMC)))>0) {
				for (i in 1:ncol(MCMC) ) {
					if (x==colnames(MCMC)[i]) {
						posterior<-MCMC[,i] } } 
			}
			if(sum(as.numeric(x==names(model)))>0) {
				for (i in 1:length(model) ) {
					if (x==names(model)[i]) {
						posterior<-model[i] } }
			}
	} else {
		if(class(model)[1]=="GCMMestimate") {
			posterior<-c()
			if(sum(as.numeric(x==names(model$PD)))>0) {
				for (i in 1:length(model$PD)) {
					if ( x==names(model$PD[i])) {
						posterior<-model$PD[i] } } } 
	} else {
		posterior<-c()
			for (i in 1:ncol(model) ) {
				if (x==colnames(model)[i]) {
					posterior<-model[,i] } }
	} } 

	if(class(posterior)=="list") { posterior<-unlist(posterior) }

	return(posterior) 
	} 


###################################################################

#' Calculate activity probability density
#' 
#' @description Support function for calculating activity probability density at a specified time from an activity curve
#' @param x Value in radians for which to predict probability density
#' @param curve Temporal data to predict density from
#' @return Returns activity probability density value
#' @export

	calcAPD<-function(x, curve) {

    bw<-overlap::getBandWidth(curve, kmax=3)/5
    dens<-overlap::densityFit(curve, seq(0,2*pi,length=128), bw)
    plot<-cbind(x=seq(0,2*pi,length=128), y=dens)
    apd<-stats::approxfun(plot)
    APDvals<-round(apd(x),3)
	return(APDvals)
}


###################################################################

#' Calculate proportions of circular variable within an interval
#' 
#' @description Support function that calculates proportion of a vector of circular data within an interval
#' @param x Vector of data
#' @param p1 Number identifying start of interval
#' @param p2 Number identifying end of interval
#' @return Returns proportion of the vector within the interval
#' @export

		calcprop<-function(x, p1, p2) {
			count<-0
			for (c in 1:length(x)) { 
				if (as.numeric(x[c]) > p1) {
				if (as.numeric(x[c]) < p2) {
					count<-count+1
			} } }
		return(count/length(x))
		}


###################################################################

#' Sample rows of dataframe
#' 
#' @description Support function that samples rows of data from a dataframe
#' @param df Dataframe
#' @param n Number of samples
#' @return Returns sample of dataframe with the number of specified rows
#' @export


	samplerows<-function(df,n){
		   return(df[sample(nrow(df),n),])
		}



###################################################################

#' Extract GCMM parameters for running GCMM simulations
#' 
#' @description Support function that extracts posterior samples of GCMM parameters for posterior simulations of GCMM activity curves
#' @param model Object of class GCMM containing output from GCMM function
#' @param sample Number of posterior samples
#' @return Returns a vector posterior draws
#' @export


	GCMMsimsparams<-function(model, sample) {

		MCMC<-combineMCMC(model)
		MCMCs<-samplerows(df=MCMC,n=sample)

	Mu1<-extractparam(MCMCs,"CircularIntercept[1]")
	P1<-extractparam(MCMCs,"ClustProb[1]")
	if (model$family=="wrappedcauchy") { 
		K1<-extractparam(MCMCs,"rhoC[1]") } else { 
		K1<-extractparam(MCMCs,"K[1]") } 
	if(model$Nclust>1) {
		Mu2<-extractparam(MCMCs,"CircularIntercept[2]")
		P2<-extractparam(MCMCs,"ClustProb[2]")
		if (model$family=="wrappedcauchy") { 
				K2<-extractparam(MCMCs,"rhoC[2]") } else { 
				K2<-extractparam(MCMCs,"K[2]") } 
		} else { Mu2<-NULL; P2<-NULL; K2<-NULL }
	if(model$Nclust>2) {
		Mu3<-extractparam(MCMCs,"CircularIntercept[3]")
		P3<-extractparam(MCMCs,"ClustProb[3]")
		if (model$family=="wrappedcauchy") { 
				K3<-extractparam(MCMCs,"rhoC[3]") } else { 
				K3<-extractparam(MCMCs,"K[3]") } 
		} else { Mu3<-NULL; P3<-NULL; K3<-NULL }
	if(model$Nclust>3) {
		Mu4<-extractparam(MCMCs,"CircularIntercept[4]")
		P4<-extractparam(MCMCs,"ClustProb[4]")
		if (model$family=="wrappedcauchy") { 
				K4<-extractparam(MCMCs,"rhoC[4]") } else { 
				K4<-extractparam(MCMCs,"K[4]") } 
			} else { Mu4<-NULL; P4<-NULL; K4<-NULL }

	simspd<-list(Mu1=Mu1,Mu2=Mu2,Mu3=Mu3,Mu4=Mu4, 
			P1=P1,P2=P2,P3=P3,P4=P4, 
			K1=K1,K2=K2,K3=K3,K4=K4,
			Nclust=model$Nclust, family=model$family	)

		return(simspd)
	}


###################################################################

#' Create GCMM simulations
#' 
#' @description Support function that creates posterior simulations of GCMM activity curves
#' @param PD Posterior samples of GCMM parameters; output from GCMMsimsparams function
#' @param s Index value for running simulations
#' @return Returns a vector of data simulated from the GCMM mixture
#' @export

	GCMMsims<-function(PD, s) { 

	    if (PD$Nclust==1) {
			if (PD$family=="vonmises") { 
				gcmmcomponentsC1<-circular::rvonmises(as.numeric(PD$P1[s])*10000, circular::circular(PD$Mu1[s]), PD$K1[s])  
				gcmmmixture<-c(gcmmcomponentsC1)
				}
			if (PD$model$family=="wrappedcauchy") { 
				gcmmcomponentsC1<-circular::rwrappedcauchy(as.numeric(PD$P1[s])*10000, circular::circular(PD$Mu1[s]), PD$K1[s])  
				gcmmmixture<-c(gcmmcomponentsC1)
				}
			}
	    if (PD$Nclust==2) {
			if (PD$family=="vonmises") { 
				gcmmcomponentsC1<-circular::rvonmises(as.numeric(PD$P1[s])*10000, circular::circular(PD$Mu1[s]), PD$K1[s])  
				gcmmcomponentsC2<-circular::rvonmises(as.numeric(PD$P2[s])*10000, circular::circular(PD$Mu2[s]), PD$K2[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2)
				}
			if (PD$family=="wrappedcauchy") { 
				gcmmcomponentsC1<-circular::rwrappedcauchy(as.numeric(PD$P1[s])*10000, circular::circular(PD$Mu1[s]), PD$K1[s])  
				gcmmcomponentsC2<-circular::rwrappedcauchy(as.numeric(PD$P2[s])*10000, circular::circular(PD$Mu2[s]), PD$K2[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2)
				}
			}
	    if (PD$Nclust==3) {
			if (PD$family=="vonmises") { 
				gcmmcomponentsC1<-circular::rvonmises(as.numeric(PD$P1[s])*10000, circular::circular(PD$Mu1[s]), PD$K1[s])  
				gcmmcomponentsC2<-circular::rvonmises(as.numeric(PD$P2[s])*10000, circular::circular(PD$Mu2[s]), PD$K2[s])  
				gcmmcomponentsC3<-circular::rvonmises(as.numeric(PD$P3[s])*10000, circular::circular(PD$Mu3[s]), PD$K3[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2,gcmmcomponentsC3)
				}
			if (PD$family=="wrappedcauchy") { 
				gcmmcomponentsC1<-circular::rwrappedcauchy(as.numeric(PD$P1[s])*10000, circular::circular(PD$Mu1[s]), PD$K1[s])  
				gcmmcomponentsC2<-circular::rwrappedcauchy(as.numeric(PD$P2[s])*10000, circular::circular(PD$Mu2[s]), PD$K2[s])  
				gcmmcomponentsC3<-circular::rwrappedcauchy(as.numeric(PD$P3[s])*10000, circular::circular(PD$Mu3[s]), PD$K3[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2,gcmmcomponentsC3)
				}
			}
	    if (PD$Nclust==4) {
			if (PD$family=="vonmises") { 
				gcmmcomponentsC1<-circular::rvonmises(as.numeric(PD$P1[s])*10000, circular::circular(PD$Mu1[s]), PD$K1[s])  
				gcmmcomponentsC2<-circular::rvonmises(as.numeric(PD$P2[s])*10000, circular::circular(PD$Mu2[s]), PD$K2[s])  
				gcmmcomponentsC3<-circular::rvonmises(as.numeric(PD$P3[s])*10000, circular::circular(PD$Mu3[s]), PD$K3[s])  
				gcmmcomponentsC4<-circular::rvonmises(as.numeric(PD$P4[s])*10000, circular::circular(PD$Mu4[s]), PD$K4[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2,gcmmcomponentsC3,gcmmcomponentsC4)
				}
			if (PD$family=="wrappedcauchy") { 
				gcmmcomponentsC1<-circular::rwrappedcauchy(as.numeric(PD$P1[s])*10000, circular::circular(PD$Mu1[s]), PD$K1[s])  
				gcmmcomponentsC2<-circular::rwrappedcauchy(as.numeric(PD$P2[s])*10000, circular::circular(PD$Mu2[s]), PD$K2[s])  
				gcmmcomponentsC3<-circular::rwrappedcauchy(as.numeric(PD$P3[s])*10000, circular::circular(PD$Mu3[s]), PD$K3[s])  
				gcmmcomponentsC4<-circular::rwrappedcauchy(as.numeric(PD$P4[s])*10000, circular::circular(PD$Mu4[s]), PD$K4[s])  
				gcmmmixture<-c(gcmmcomponentsC1,gcmmcomponentsC2,gcmmcomponentsC3,gcmmcomponentsC4)
				}
			}
	
		return(gcmmmixture)

	}



###################################################################

#' Report function progress
#' 
#' @description Support function that prints progress of functions with long computation times at 10% intervals
#' @param s Number of iteration in the loop
#' @param sample Total number of iterations in the loop
#' @return No return value; prints progress of function at 10% intervals
#' @export

	progress<-function(s, sample) {

	if (s==1) { message("  Starting analysis...") }
	if (s==floor(sample)*.1 ) { message("    10% completed...") }
	if (s==floor(sample)*.2 ) { message("    20% completed...") }
	if (s==floor(sample)*.3 ) { message("    30% completed...") }
	if (s==floor(sample)*.4 ) { message("    40% completed...") }
	if (s==floor(sample)*.5 ) { message("    50% completed...") }
	if (s==floor(sample)*.6 ) { message("    60% completed...") }
	if (s==floor(sample)*.7 ) { message("    70% completed...") }
	if (s==floor(sample)*.8 ) { message("    80% completed...") }
	if (s==floor(sample)*.9 ) { message("    90% completed...") }
	if (s==floor(sample))  { message("    100% completed") }

	}


