#' Predict activity probability density at a time point
#' 
#' @description Calculates predicted activity probability density from GCMM model at a specific time point
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param timepoint Time point for which to predict activity probability density
#' @param HDI Logical argument for whether to calculate 95% HDI; default=TRUE
#' @param sample Number of posterior samples from which to build HDI; default=1000
#' @param scale Scale of the data, either "2pi" for 0,2pi or "pi" for -pi,pi; default="2pi"
#' @return Returns numeric vector of estimated probability density if HDI=FALSE or posterior probability distribution if HDI=TRUE
#' @export


	GCMMpdens<-function(model,timepoint,HDI=TRUE,sample=1000,scale="2pi") {

	if (HDI==FALSE) {
		if (class(model)[1]=="GCMM") {
			mix<-convertRad(as.numeric(model$GCMMmixture),to=scale)
			pdens<-calcAPD(x=timepoint,curve=mix)
			} 
		print(pdens)
	}

	if (HDI==TRUE) { 
		pdens<-c()
		PD<-GCMMsimsparams(model, sample=sample) 
			for (s in 1:sample) {
				gcmmmixture<-GCMMsims(PD, s=s) 
				mix<-convertRad(as.numeric(gcmmmixture),to="scale")
				p<-calcAPD(x=timepoint,curve=mix)#AnimalAPD::APDraw(focal=x,contingent=mix,adjust=5)
				pdens<-c(pdens,p)					
					progress(s=s,sample=sample)
				}
		out<-HDI(pdens)
			print(out)
		} else { HDI<-"NULL" }

		PD<-list(pdens=pdens)
		output<-list(summary=out,PD=PD)
		class(output)<-"GCMMestimate"
		return(output)

	}


#####################################################################

#' Activity Probability Density at Peak Time of Another
#' 
#' @description Calculates activity probability density from one GCMM model at the peak activity time of a second GCMM model
#' @param model1 Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param model2 Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param sample Number of posterior samples from which to build the HDIs
#' @param HDIprob Value for probability mass to use for HDI; default=95%
#' @return Returns matrix with the mean and HDI of activity probability density estimated from both GCMM models at the peak activity time of the other. Posterior distributions of activity probability density and peak activity times for both GCMM models are also saved.
#' @export

	APDatPeak<-function(model1,model2,sample=1000, HDIprob=0.95) {

		m1PD<-GCMMsimsparams(model1, sample=sample) 
		m2PD<-GCMMsimsparams(model2, sample=sample) 
			m1apd<-c(); m2apd<-c(); m1p<-c(); m2p<-c()

			for (s in 1:sample) { 
			m1sims<-GCMMsims(PD=m1PD, s=s)
			m2sims<-GCMMsims(PD=m2PD, s=s)
				m1HDE<-activityHPDmean(m1sims,prob=0.5,scale="2pi",silent=TRUE)
					m1peak<-m1HDE$summary[1,1]
				m2HDE<-activityHPDmean(m2sims,prob=0.5,scale="2pi",silent=TRUE)
					m2peak<-m2HDE$summary[1,1]
				m1apd[s]<-calcAPD(m2peak,m1sims)
				m2apd[s]<-calcAPD(x=m1peak,curve=m2sims)
					m1p[s]<-m1peak; m2p[s]<-m2peak
				progress(s=s,sample=sample)
			}

		m1<-HDI(m1apd,prob=HDIprob)
		m2<-HDI(m2apd,prob=HDIprob)
		output<-rbind(m1,m2); colnames(output)<-c("Mean","Lower95","Upper95")
			rownames(output)<-c("Model 1", "Model 2")
				print(round(output,3))

		out<-list(HDI=output, PD=list(m1APD=m1apd, m2APD=m2apd, m1peak=m1p, m2peak=m2p))
			class(out)<-"GCMMestimate"

		return(out)

	}


###########################
#' Activity Probability Density Point Plot
#' 
#' @description Plot of GCMM activity with predicted activity probability density at particular time points
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param time Time point for which to plot predicted activity probability density
#' @param ymax Upper limit of y axis
#' @param scale Scale for which to plot the activity curve, either "2pi" for 0,2pi or "pi" for -pi, pi; default="2pi"
#' @param cex Size of plotted point; default=2.5
#' @param col Colour of plotted point
#' @param axisunits Units for xaxis 
#' @return No return value; prints plot of activity curve with activity probability density prediction at the specified time point and returns dataframe of time points and activity probability density
#' @export

	APDpointplot<-function(model,time,ymax="NULL",scale="2pi",cex=2.5,col="lightseagreen", axisunits="radians") {

	mixtureplot(model, scale=scale, ymax=ymax)
	APD<-c()
	for (i in 1:length(time) ) {
		a<-calcAPD(x=time[1],curve=model$GCMMmixture)
		APD<-c(APD,a)
			graphics::points(x=time[i],y=a,pch=16,cex=cex)
			graphics::points(x=time[i],y=a,pch=16,cex=cex*.55,col=col)
	}

	df<-data.frame(time, APD)
	return(df)
}


