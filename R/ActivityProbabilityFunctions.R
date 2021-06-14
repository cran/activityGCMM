#########################################################

#' Probability of Activity during Time Period
#' 
#' @description Calculate activity probability estimates from a GCMM model during a specific period of time
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param timestart Start of time period (in radians)
#' @param timeend End of time period (in radians)
#' @param sample Number of posterior samples from which to calculate the HDI; default=1000
#' @return Returns object of class \code{GCMMestimate} containing list with summary of output with mean and 95% HDI and the posterior distribution of predicted activity probability
#' @export

		GCMMprob<-function(model, timestart, timeend, sample=1000) {

	t1<-timestart; t2<-timeend
	if(timeend<timestart) { t1<-timeend; t2<-timestart }

	MCMC<-combineMCMC(model)
	MCMCsample<-samplerows(df=MCMC,n=sample)

	proportion<-c()
		PD<-GCMMsimsparams(model, sample)
		for (s in 1:sample) {
			gcmmmixture<-GCMMsims(PD, s)

			count<-0
			for (c in 1:length(gcmmmixture)) { 
				if (as.numeric(gcmmmixture[c]) > t1) {
					if (as.numeric(gcmmmixture[c]) < t2) {
						count<-count+1
					} } }
		proportion[s]<-count/length(gcmmmixture)

		progress(s=s,sample=sample)

		}

	if(t2<t1) { proportion<- 1-proportion }
	out<-HDI(proportion)
	print(out)
	
	output<-list(summary=out,PD=list(Probability=proportion))
	class(output)<-"GCMMestimate"
	return(output)

	}



#########################################################

#' Activity highest posterior density interval from mean activity curve
#' 
#' @description Estimates activity highest posterior density interval (HPD), HPD duration, number of activity peaks, 
#'         peak activity times and maximum activity probability density for a given probability density mass from the 
#'         GCMM activity curve predicted by the GCMM parameter posterior distribution means
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param prob Value to use for probability density mass; default=0.50
#' @param silent Logical vector for whether to print output and plot; default=FALSE
#' @param plot Logical argument for whether to plot activity curve with HPD; default=FALSE
#' @param scale Scale of the data, either "2pi" for 0,2pi or "pi" for -pi,pi; default="2pi"
#' @param ymax Value to use as upper limit of y axis on activity curve plot
#' @param col Colour for HPD on plot
#' @return Returns object of class \code{GCMMestimate} with list with estimated peak activity times, maximum activity probability density,  
#'         HPD interval, HPD duration, and number of activity peaks. 
#' @export

	activityHPDmean<-function(model, prob=0.5, scale="2pi", silent=FALSE, plot=FALSE, col="cyan4", ymax="NULL") {

		if (class(model)[1]=="GCMM") { X<-model$GCMMmixture } else { X<-as.numeric(model) }
		X<-convertRad(X,to="2pi")

		apd<-calcAPD(x=X,curve=X)
			df<-data.frame(apd,X) 						# dataframe of density
			df<-df[order(df$X),]; df$index<-1:nrow(df) 	  	# sorted by time; added index
			df<-df[order(-df$apd),] 					# sorted by highest apd
			HD<-df[1:floor(prob*nrow(df)),] 				# taking proportion according to prob
			HD<-HD[order(HD$X),]						# sorting by time
			int<-matrix(ncol=3);colnames(int)<-colnames(HD)
			for(i in 2:length(HD$X)) { 					# finding breaks in continuity, i.e. new interval	
				if ((HD$index[i]-HD$index[i-1])>1) {
					int<-rbind(int,HD[i-1,],HD[i,])
				} }
				int[1,]<-HD[1,]
				int<-rbind(int,HD[length(HD$X),]) 
					mat<-matrix(int$X,ncol=2,byrow=TRUE)
					colnames(mat)<-c("Start","End")
 
		maxAPD<-c();maxAPDtime<-c(); Start<-c(); End<-c(); intlength<-c()
		for (i in 1:nrow(mat) ) { 
			Start<-c(Start,mat[i,1])
			End<-c(End,mat[i,2])
			tmp<-HD[which(HD$X>=mat[i,1] & HD$X<=mat[i,2]),]	
			maxAPD<-c(maxAPD,max(tmp$apd))
				for(z in 1:nrow(tmp)) { #for (l in 1:length(maxAPD)) {
					if (max(tmp$apd)==tmp$apd[z]) { maxAPDtime[i]<-tmp$X[z] } }
			intlength<-c(intlength,nrow(tmp))
			}
		Duration<-End-Start
		maximums<-data.frame(Time=round(maxAPDtime,3),APD=round(maxAPD,3),
			TimeStart=round(Start,3),TimeEnd=round(End,3),Duration=round(Duration,3),Prop=round(intlength/length(X),3))

			if (	min(df$X) == min(HD$X) ) {
				newend<-maximums$TimeEnd[1]
				newstart<-maximums$TimeStart[nrow(maximums)]
				newAPD<-max(maximums$APD[1],maximums$APD[nrow(maximums)])
				if(maximums$APD[1]==newAPD){newTime<-maximums$Time[1]} else {
					newTime<-maximums$Time[nrow(maximums)] }
				newprop<-sum(maximums$Prop[1],maximums$Prop[nrow(maximums)])
				newDur<-sum(maximums$Duration[1],maximums$Duration[nrow(maximums)])
			maximums$Time[1]<-newTime
			maximums$APD[1]<-newAPD
			maximums$TimeStart[1]<-newstart
			maximums$TimeEnd[1]<-newend
			maximums$Duration[1]<-newDur
			maximums$Prop[1]<-newprop
			maximums<-maximums[-nrow(maximums),]
			} 
		names(maximums)<-c("PeakTime","MaxAPD","StartTimeHDE","EndTimeHDE","Duration","Prop")
		output<-maximums[order(-maximums$MaxAPD),]
		if (min(output$Prop)<.01) { 
			output<-output[order(output$Prop),]
			output<-output[-1,]
			}
		output<-output[order(-output$MaxAPD),]
		Npeaks<-nrow(output)


		hpd<-c(HD$X,convertRad(HD$X,to="pi"))
		if(silent==FALSE) {
			if(plot==TRUE) {
			if (class(model)[1]=="GCMM") { 
				mixtureplot(model,scale=scale,ymax=ymax)
					graphics::abline(v=c(hpd),col=col)
					graphics::par(new=TRUE); mixtureplot(model,scale=scale,ymax=ymax)
		}	}
	} 
		if(silent==FALSE) { print(output) 

		PeakTimeHrs<-round(convertRad(output$PeakTime,to="hours"),2)
		StartTimeHrs<-round(convertRad(output$StartTimeHDE,to="hours"),2)
		EndTimeHrs<-round(convertRad(output$EndTimeHDE,to="hours"),2)
		DurationHrs<-round(convertRad(output$Duration,to="hours"),2)
		Hours<-data.frame(PeakTime=PeakTimeHrs,StartTimeHDE=StartTimeHrs,
			EndTimeHDE=EndTimeHrs,Duration=DurationHrs)
		message(""); message("Results Converted to Hours:")
		print(Hours)
		 }

		HDE<-list(summary=output,Npeaks=Npeaks,DurationTotal=sum(output$Duration),
			HD=HD$X,lines=hpd,intervals=mat)
		return(HDE)
	}


#########################################################

#' Activity highest posterior density interval estimates (activityHPD)
#' 
#' @description Calculate activity highest posterior density interval (activityHPD), HPD duration, number of activity peaks, 
#'         peak activity times and maximum activity probability density for a given probability density mass
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param prob Value to use for probability density mass; default=0.50
#' @param scale Scale of the data for plotting, either "2pi" for 0,2pi or "pi" for -pi,pi; default="2pi"
#' @param sample Number of posterior samples for which to calculate 95% HDIs
#' @return Returns object of class \code{GCMMestimate} with list including estimated peak activity times, maximum activity  
#'         probability density, HPD interval, HPD duration, and number of activity peaks
#' @export

	activityHPD<-function(model, sample=1000, prob=0.5, scale="2pi") {
		PD<-GCMMsimsparams(model, sample=sample)
		Peaks<-c(); MaxAPD<-c(); Duration<-c(); PeakTime<-c()
		for (s in 1:sample) {
			gcmmmixture<-GCMMsims(PD, s=s) 
				gcmmmixture<-convertRad(as.numeric(gcmmmixture),to="2pi")
			HDE<-activityHPDmean(gcmmmixture,prob=prob,scale="2pi",silent=TRUE)
				Peaks[s]<-HDE$Npeaks
				MaxAPD[s]<-max(HDE$summary$MaxAPD)
				Duration[s]<-HDE$DurationTotal
				PeakTime[s]<-HDE$summary[1,1]
		progress(s=s,sample=sample)
		}

	p<-c(mode(Peaks),HDI(Peaks)[2],HDI(Peaks)[3])
	apd<-HDI(MaxAPD)
	d<-HDI(Duration)
		t<-mean(circular::as.circular(PeakTime,type="angles",units="radians",template="none", zero=0,rotation="clock",modulo="asis"))
			t<-round(c(as.numeric(t),HDI(PeakTime)[2],HDI(PeakTime)[3]),3)
			PeakTimepi<-convertRad(PeakTime,to="pi")
			t2<-mean(circular::as.circular(PeakTimepi,type="angles",units="radians",template="none", zero=0,rotation="clock",modulo="asis"))
			t2<-c(as.numeric(t2),HDI(PeakTimepi)[2],HDI(PeakTimepi)[3])
				if((t[3]-t[2])>(t2[3]-t2[2])) { t[2]<-t2[2]; t[3]<-t2[3] }
			t<-round(t,3)

	output<-rbind(p,apd,t,d);colnames(output)<-c("Mode/Mean","Lower95","Upper95")
		rownames(output)<-c("NPeaks","MaxAPD","PeakTime","Duration")
		message(""); print(output)

	message(""); message("Mean Estimate:")
	meanHPD<-activityHPDmean(model,prob=prob,scale=scale,silent=FALSE)

	out<-list(summary=output,mean=meanHPD$summary,HPDoutput=meanHPD,PD=list(Peaks=Peaks,MaxAPD=MaxAPD,PeakTime=PeakTime,Duration=Duration),prob=prob,sample=sample, model=model)
	class(out)<-"GCMMestimate"

	return(out)
	
	}

	
############################

#' Plot estimated number of activity peaks
#' @description Plot of posterior samples of estimated number of activity peaks from \code{\link{activityHPD}} function
#' @param x Object of class \code{GCMMestimate} with output from \code{\link{activityHPD}} function
#' @param col Colour for plot
#' @return No return value; prints histogram plot of posterior estimates
#' @export

		peaksPDplot<-function(x,col="cyan4") {
			t<-table(x$PD$Peaks)
			p<-t/sum(t)
			graphics::barplot(p,col=col,border=col,xlab="Estimated Number of Activity Peaks",ylab="Proportion")
		print(p)
		}


############################

#' Plot estimated time of activity peaks
#' @description Plot mean GCMM activity curve with peak activity time from \code{\link{activityHPD}} function
#' @param x Object of class \code{GCMMestimate} with output from \code{\link{activityHPD}} function
#' @param scale Scale for plotting the activity curve, either "2pi" for 0,2pi or "pi" for -pi,pi; default="2pi"
#' @param ymax Value for upper limit of y axis
#' @return No return value; prints plot activity curve and peak activity time
#' @export

	GCMMpeaksplot<-function(x,scale="2pi",ymax="NULL") {
		mixtureplot(x$model,scale=scale,ymax=ymax)
		times<-convertRad(x$mean[,1],to="scale")
		graphics::abline(v=c(times),lty=5,lwd=3,col="red2")
		}


###########################################################

#' Plot activity curve with activityHPDs
#' @description Plot GCMM activity curve with activityHPDs
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param prob Vector of values of probability density mass for activityHPD; default=c(0.50,0.25)
#' @param scale Scale for plotting the activity curve, either "2pi" for 0,2pi or "pi" for -pi,pi; default="2pi"
#' @param col Vector of colours for plot
#' @param ymax Value for upper limit of y axis
#' @param axisunits Units for x axis, either "radians", "time", "solar", "sun", or "none"; default="radians"
#' @return No return value; prints plot of activity curve with activityHPDs
#' @export

	plotactivityHPD<-function(model,prob=c(0.75,0.50,0.25),col=c("lightseagreen","aquamarine3","aquamarine"),scale="2pi",ymax="NULL",axisunits="radians") {
		mixtureplot(model,scale=scale,ymax=ymax,axisunits=axisunits)
		prob<-sort(prob,decreasing=TRUE)
		for(i in 1:length(prob)) {
			graphics::par(new=TRUE);x<-activityHPDmean(model,prob=prob[i],col=col[i],plot=TRUE)
			}
	}


##################################################################

#' @title Activity HPD Overlap
#' @description Calculate whether there is overlap between two HPDs, the amount of overlap, and the probability of activity
#'    during the HPD of the other
#' @param model1 Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param model2 Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param prob Value of probability density mass for activityHPD; default=0.50
#' @param sample Number of posterior samples for calculating 95% HDI; default=1000
#' @param scale Scale for plotting the activity curve, either "2pi" for 0,2pi or "pi" for -pi,pi; default="2pi"
#' @return Returns object of class  \code{GCMMestimate} containing list of posterior distributions and summary information of mean and 95% HDI
#' @export

	HPDoverlap<-function(model1, model2, prob=0.5, sample=1000, scale="2pi") {

	m1PD<-GCMMsimsparams(model1, sample=sample) 
	m2PD<-GCMMsimsparams(model2, sample=sample) 
		probm1<-c(); probm2<-c(); Diff<-c(); Overlap<-c()
		m1Dur<-c(); m2Dur<-c(); m1probZ<-c(); m2probZ<-c()
		HDEDiff1<-c(); HDEOverlap1<-c(); HDEDiff2<-c()
		HDEOverlap2<-c()

	for (s in 1:sample) {
		m1gcmmmixture<-GCMMsims(m1PD, s=s) 
			m1gcmmmixture<-convertRad(as.numeric(m1gcmmmixture),to="2pi")
		m2gcmmmixture<-GCMMsims(m2PD, s=s) 
			m2gcmmmixture<-convertRad(as.numeric(m2gcmmmixture),to="2pi")

		m1HDE<-activityHPDmean(m1gcmmmixture,prob=prob,scale="2pi",silent=TRUE)
		m2HDE<-activityHPDmean(m2gcmmmixture,prob=prob,scale="2pi",silent=TRUE)
			m1Dur[s]<-m1HDE$DurationTotal
			m2Dur[s]<-m2HDE$DurationTotal

			m1prop<-c(); hdeDiff1<-0; hdeOverlap1<-c()
				for(e in 1:nrow(m2HDE$intervals)) { 
					prop1<-calcprop(x=m1gcmmmixture, p1=m2HDE$intervals[e,1], p2=m2HDE$intervals[e,2])
					m1prop<-sum(m1prop,prop1)		
					hded1<-calcprop(x=m1HDE$HD, p1=m2HDE$intervals[e,1], p2=m2HDE$intervals[e,2])
					hdeOverlap1<-sum(hdeOverlap1,hded1)
					if (hdeOverlap1>0) { hdeDiff1<-1}
					}
			m2prop<-c(); hdeDiff2<-0; hdeOverlap2<-c()
				for(e in 1:nrow(m1HDE$intervals)) { 
					prop2<-calcprop(x=m2gcmmmixture, p1=m1HDE$intervals[e,1], p2=m1HDE$intervals[e,2])
					m2prop<-sum(m2prop,prop2)
					hded2<-calcprop(x=m2HDE$HD, p1=m1HDE$intervals[e,1], p2=m1HDE$intervals[e,2])
					hdeOverlap2<-sum(hdeOverlap2,hded2)
					if (hdeOverlap2>0) { hdeDiff2<-1 }
					}
			probm1[s]<-m1prop
			probm2[s]<-m2prop
			HDEDiff1[s]<-hdeDiff1
			HDEOverlap1[s]<-hdeOverlap1
			HDEDiff2[s]<-hdeDiff2
			HDEOverlap2[s]<-hdeOverlap2
		Diff[s]<-mean(hdeDiff1,hdeDiff2)
		Overlap[s]<-mean(hdeOverlap1,hdeOverlap2)
				m1probZ[s]<-m1prop/(m2HDE$DurationTotal/(2*pi))
				m2probZ[s]<-m2prop/(m1HDE$DurationTotal/(2*pi))

	progress(s=s,sample=sample)

	  }


	M1<-HDI(probm1)
	M2<-HDI(probm2)
		M1z<-HDI(m1probZ)
		M2z<-HDI(m2probZ)

	Probability<-rbind(M1,M2,M1z,M2z);colnames(Probability)<-c("Mean","Lower95","Upper95")
		rownames(Probability)<-c("Model 1 Raw","Model 2 Raw","Model 1 Standardized","Model 2 Standardized")
		message("")
		message(paste("Probability of activity during ",prob*100,"% HDE of the other:",sep=""))
		print(Probability)

	HDEDiff<-HDI(Diff)
		message("")
		message(paste("Probability of overlap between the ",prob*100,"% HDEs:",sep=""))
		print(HDEDiff)

	HDEOverlap<-HDI(Overlap)
		message("")
		message(paste("Estimated amount of overlap between the two ",prob*100,"% HDEs:",sep=""))
		print(HDEOverlap)


	out<-list(summary=Probability, PD=list(AP1=probm1, AP2=probm2, AP1z=m1probZ, AP2z=m2probZ, Duration1=m1Dur,Duration2=m2Dur,
			HPDOverlap1=HDEOverlap1, HPDOverlap2=HDEOverlap2, MeanOverlap=Overlap), 
			HPDDiff=HDEDiff, HPDDiff1=HDEDiff1, HPDDiff2=HDEDiff2, MeanDiff=Diff,
			HPDOverlap=HDEOverlap)
	class(out)<-"GCMMestimate"
	return(out)

	}



##########################################

#' @title Circular plot of activity HPD intervals
#' @description Circular plot of activity HPD intervals from GCMM activity curves
#' @param models List of one or more objects of class \code{GCMM} containing output from the \code{\link{GCMM}} function
#' @param prob Value for activityHPD probability density mass; default=0.5 (i.e. 50% HPD)
#' @param col Vector of colours to use in the plot
#' @param axisunits Units to be used for the axis, either "radians", "sun", "time", or "none"
#' @return No return value; prints plot
#' @examples
#' \donttest{ FoxActivityGCMM<-GCMM(data=redfoxsample$Radians, 
#'         RE1=redfoxsample$SamplingPeriod, family="vonmises", autorun=FALSE)
#'      HumanActivityGCMM<-GCMM(data=humanssample$Radians, RE1=humanssample$SamplingPeriod, 
#'         family="vonmises", autorun=FALSE)
#'     circplotHPD(models=list(FoxActivityGCMM,HumanActivityGCMM)) }
#' @export

	circplotHPD<-function(models, prob=0.5, col=c("cyan3","orchid","deeppink","dodgerblue"), axisunits=c("radians","sun","time") )	{

	graphics::par(mfrow=c(1,1))
	for (m in 1:length(models) ) { 
		HPD<-activityHPDmean(models[[m]],prob=prob,silent=TRUE)
			graphics::plot(circular::as.circular(HPD$lines, type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),
					cex=1.4,shrink=1.7,start.sep=.15*m,axes=FALSE,tol=0.501,main="")
				graphics::points(circular::as.circular(HPD$lines,type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),
					shrink=1.7,start.sep=.15*m,cex=.8,col=col[m])
				graphics::par(new=TRUE)
		}
			circaxis(axisunits)
			graphics::text(x=0, y=0, paste(prob*100,"% HPD",sep=""), font=2, cex=0.6)

	}

