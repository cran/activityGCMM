
#' @title GCMM Components Plot
#' @description Plot of activity curves for the separate components in the circular mixture model
#' @param model Model output from \code{\link{GCMM}} function, object of class \code{GCMM}
#' @param rug Logical argument for whether to plot a rug of the raw values. Plotting the rug for the separate components requires that saveclustIDs=TRUE when running GCMM. default=FALSE
#' @param scale Scale for the plot, either "pi" (-pi, pi) or "2pi" (0, 2pi); default is that recommended by the GCMM function
#' @param ruglwd Line width for rug plot
#' @param ymax Value for upper limit of y-axis
#' @param lwd Line width for activity curve lines
#' @param lty Line type for activity curve
#' @param col Character vector for colours for the activity curve lines and rug plot; must be of equal length to the number of components
#' @param axisunits Scale to use for the xaxis, either "radians", "time", "solar", "sun", or "none"; default="radians"
#' @param xlines Whether to include lines on the graph for the x axis labels; default=TRUE
#' @return Plot of the separate components of the circular mixture model
#' 
#' @examples
#' \donttest{ FoxActivityGCMM<-GCMM(data=redfoxsample$Radians, 
#'            RE1=redfoxsample$SamplingPeriod, family="vonmises", autorun=FALSE,
#'            adapt=0, sample=300, burnin=300, thin=1,n.chains=2)
#'         componentsplot(FoxActivityGCMM) }
#'  
#' @export

	componentsplot<-function(model,rug=FALSE,ruglwd=2,lwd=3,col=c("black","grey40","grey60","grey80"),
		scale="NULL", ymax="NULL", lty=1,
		axisunits=c("radians","sun","solar","time","none"),xlines=TRUE) {


	if (ymax=="NULL") { ymax<-yMax(models=list(model),type="components")*.4 }

	# Plot: 
	if(rug==TRUE){ yax="r" } else { yax="i" }
	if(scale=="NULL") { scale<-model$scale }
	if (scale=="pi") { lim1<-(-1*pi); lim2<-pi; xcenter<-c("midnight") }
		if (scale=="2pi") { lim1<-0; lim2<-2*pi; xcenter<-c("noon") }
		xx <- seq(lim1, lim2, length = 128)

		for(i in 1:length(model$GCMMcomponents)) {
				bw<-overlap::getBandWidth(model$GCMMcomponents[[i]], kmax=3)/5   
				dens<-overlap::densityFit(model$GCMMcomponents[[i]], xx, bw)
				graphics::plot(dens*as.numeric(model$P[i]),x=xx,type="l",col=col[i],yaxs=yax,xaxt="n",xaxs="i",lwd=lwd,lty=lty,xlim=c(lim1,lim2),
					ylim=c(0,ymax),main="",xlab="Time",ylab="Density")
				graphics::box(lwd=2); graphics::par(new=TRUE)
			}

	# Rug: 
	if(rug==TRUE) {
		if(scale=="pi") { data<- convertRad(model$data,to="pi") }
		if(scale=="2pi") { data<- convertRad(model$data,to="2pi") }
			df<-data.frame(Data=data,Clust=as.numeric(model$clust))
			df<-df[order(df$Clust),]
			data<-data
			cdata<-list()
			s<-1
			for(i in 1:length(model$GCMMcomponents)) { 
				N<-as.numeric(table(df$Clust)[i] )
				N2<-s+N-1 
				c<-data[s:N2] 
				s<-s+N
					rug(c,col=col[i],lwd=ruglwd,ticksize=0.035); graphics::box(lwd=2)
				}
			}
	xaxis(axisunits[1],lines=xlines); graphics::box(lwd=2)
	}


#' @title GCMM Mixture Plot
#' @description Plot of estimated activity curve from the circular mixture model
#' @param model Model output from \code{\link{GCMM}} function, object of class \code{GCMM}
#' @param rug Logical argument for whether to plot a rug of the raw values. Plotting the rug for the separate components requires that saveclustIDs=TRUE when running \code{\link{GCMM}} or \code{\link{updateGCMM}}. default=FALSE
#' @param scale Scale for the plot, either "pi" (-pi, pi) or "2pi" (0, 2pi); default is that recommended by the \code{\link{GCMM}} function
#' @param col Line colour for plot
#' @param ruglwd Line width for rug plot
#' @param lwd Line width for activity curve
#' @param lty Line type for activity curve
#' @param ymax Value to use as y-axis maximum
#' @param axisunits Scale to use for the xaxis, either "radians", "time", "solar", "sun", or "none"; default="radians"
#' @param xlines Whether to include lines on the graph for the x axis labels; default=TRUE
#' @return Prints mixture plot of the estimated activity curve from the circular mixture model
#' 
#' @examples
#' \donttest{ FoxActivityGCMM<-GCMM(data=redfoxsample$Radians, 
#'               RE1=redfoxsample$SamplingPeriod, family="vonmises", autorun=FALSE,
#'               adapt=0, sample=300, burnin=300, thin=1, n.chains=2)
#'            mixtureplot(FoxActivityGCMM) }
#'  
#' @export
	mixtureplot<-function(model,rug=FALSE,ruglwd=2,lwd=3,scale="NULL",ymax="NULL",col="black", lty=1,
		axisunits=c("radians","time","sun","solar","none"),xlines=TRUE) {

	if(rug==TRUE){yax="r"} else {yax="i" }
	if(scale=="NULL") { scale<-model$scale }
	if (scale=="pi") { lim1<-(-1*pi); lim2<-pi; xcenter<-c("midnight") }
		if (scale=="2pi") { lim1<-0; lim2<-2*pi; xcenter<-c("noon") }
	if (ymax=="NULL") { ymax<-yMax(models=list(model))*.8 }
		overlap::densityPlot(as.numeric(model$GCMMmixture),xcenter=xcenter,xscale=NA,yaxs=yax,xaxs="i",extend=NA,lwd=lwd,adjust=5,xaxt="n",xlim=c(lim1,lim2),ylim=c(0,ymax),main="",col=col,lty=lty)
 
	if(rug==TRUE) {
		if(scale=="pi") { data<- convertRad(model$data,to="pi") }
		if(scale=="2pi") { data<- convertRad(model$data,to="2pi") }
		rug(data,col="grey15",lwd=lwd,ticksize=0.035);graphics::box(lwd=2)
	}
	xaxis(axisunits[1],lines=xlines); graphics::box(lwd=2)

}


#' @title GCMM Combined Plot
#' @description Combined plot of estimated activity curve from mixture model and separate mixture components
#' @param model Model output from \code{\link{GCMM}} function, object of class \code{GCMM}
#' @param rug Logical argument for whether to plot a rug of the raw values. Plotting the rug for the separate components requires that saveclustIDs=TRUE when running GCMM. default=FALSE
#' @param scale Scale for the plot, either "pi" (-pi, pi) or "2pi" (0, 2pi); default is that recommended by the GCMM function
#' @param ymax Value for upper limit of y-axis
#' @param ruglwd Line width for rug plot
#' @param lwdc Line width for activity curve lines for components
#' @param lwdm Line width for activity curve lines from mixture
#' @param colc Character vector for colours for the activity curve lines and rug plot for components; must be of equal length to the number of components
#' @param ltyc Line type for activity curves for components
#' @param ltym Line type for activity curves for mixture
#' @param colm Character vector for colour of activity curve line from mixture model
#' @param axisunits Scale to use for the xaxis, either "radians", "time", "solar", "sun", or "none"; default="radians"
#' @param xlines Whether to include lines on the graph for the x axis labels; default=TRUE
#' @return Prints combined plot of estimated activity curve from mixture model and separate mixture components
#' 
#' @examples
#' \donttest{ FoxActivityGCMM<-GCMM(data=redfoxsample$Radians, 
#'               RE1=redfoxsample$SamplingPeriod, family="vonmises", autorun=FALSE,
#'               adapt=0, sample=300, burnin=300, thin=1,n.chains=2)
#'            comboplot(FoxActivityGCMM) }
#'  
#' @export

	comboplot<-function(model,rug=FALSE,ruglwd=2,ltyc=2,ltym=1,lwdc=3,lwdm=3,colc=c("grey40","grey55","grey70","grey85"), colm="black", scale="NULL", ymax="NULL",
		axisunits=c("radians","sun","solar","time","none"),xlines=TRUE) {

	if (ymax=="NULL") { ymax<-yMax(models=list(model),type="mixture") }

	# Plot: 
	if(rug==TRUE){ yax="r" } else { yax="i" }
	if(scale=="NULL") { scale<-model$scale }
	if (scale=="pi") { lim1<-(-1*pi); lim2<-pi; xcenter<-c("midnight") }
		if (scale=="2pi") { lim1<-0; lim2<-2*pi; xcenter<-c("noon") }
		xx <- seq(lim1, lim2, length = 128)

		for(i in 1:length(model$GCMMcomponents)) {
				bw<-overlap::getBandWidth(model$GCMMcomponents[[i]], kmax=3)/5   
				dens<-overlap::densityFit(model$GCMMcomponents[[i]], xx, bw)
				graphics::plot(dens*as.numeric(model$P[i]),x=xx,type="l",col=colc[i],xaxt="n",yaxs=yax,xaxs="i",lwd=lwdc,xlim=c(lim1,lim2),
					ylim=c(0,ymax),main="",xlab="",ylab="",lty=ltyc)
				graphics::box(lwd=2); graphics::par(new=TRUE)
			}
		dens<-overlap::densityFit(model$GCMMmixture, xx, bw)
		graphics::plot(dens*1.05,x=xx,type="l",col=colm,yaxs=yax,xaxs="i",xaxt="n",lwd=lwdm,lty=ltym,xlim=c(lim1,lim2),
				ylim=c(0,ymax),main="",xlab="Time",ylab="Density")

	# Rug: 
	if(rug==TRUE) {
		if(scale=="pi") { data<- convertRad(model$data,to="pi") }
		if(scale=="2pi") { data<- convertRad(model$data,to="2pi") }
			df<-data.frame(Data=data,Clust=as.numeric(model$clust))
			df<-df[order(df$Clust),]
			data<-data
			cdata<-list()
			s<-1
			for(i in 1:length(model$GCMMcomponents)) { 
				N<-as.numeric(table(df$Clust)[i] )
				N2<-s+N-1 
				c<-data[s:N2] 
				s<-s+N
					rug(c,col=colc[i],lwd=ruglwd,ticksize=0.035); graphics::box(lwd=2)
				}
			}
		xaxis(axisunits[1],lines=xlines); graphics::box(lwd=2)
	}

############################################
#' @title Plot multiple GCMM activity curves
#' @description Plot of multiple GCMM activity curves
#' @param models List of one or more objects of class \code{GCMM} containing output from the \code{\link{GCMM}} function
#' @param type Type of activity plots, either "mixture" for mixture plots (default) or "components" for components plots
#' @param scale Scale for the plot, either "pi" (-pi, pi) or "2pi" (0, 2pi); default is that recommended by the GCMM function
#' @param ymax Value for upper limit of y-axis
#' @param col Vector of colours to use for the activity curve lines in the plot
#' @param lty Vector of line types to use for the activity curves
#' @param lwd Value for the width of the lines for the activity curves
#' @param axisunits Scale to use for the xaxis, either "radians", "time", "solar", "sun", or "none"; default="radians"
#' @param xlines Whether to include lines on the graph for the x axis labels; default=TRUE
#' @return Prints plot
#' @examples
#' \donttest{ FoxActivityGCMM<-GCMM(data=redfoxsample$Radians, 
#'               RE1=redfoxsample$SamplingPeriod, family="vonmises", autorun=FALSE,
#'               adapt=0, sample=300, burnin=300, thin=1, n.chains=2)
#'            HumanActivityGCMM<-GCMM(data=humanssample$Radians, RE1=humanssample$SamplingPeriod, 
#'              family="vonmises", autorun=FALSE, adapt=0, sample=300, burnin=300, thin=1, n.chains=2)
#'            multiplot(models=list(FoxActivityGCMM,HumanActivityGCMM)) }
#' 
#' @export

	multiplot<-function(models,ymax="NULL",scale="2pi",lwd=3,type=c("mixture","components"),
		lty=c(1,2,3,4,5),col=c("grey15","grey40","grey55","grey70"),axisunits=c("radians","sun","solar","time","none"),xlines=TRUE) {

	if (ymax=="NULL") { ymax<-yMax(models,type="mixture") }
	if (type[1]=="components") { cols<-list(); ltys<-list()
		for (m in 1:length(models)) {
				cols[[m]]<-c(rep(col[m],length(models[[m]]$GCMMcomponents)))
				ltys[[m]]<-c(rep(lty[m],length(models[[m]]$GCMMcomponents)))
			}	}

	for (i in 1:length(models)) { 			
		if (type[1]=="mixture") {
			mixtureplot(models[[i]], scale=scale, ymax=ymax, axisunits="none", xlines=FALSE, lty=lty[i],col=col[i]); graphics::par(new=TRUE)
			}
		if (type[1]=="components") {
			componentsplot(models[[i]], scale=scale, ymax=ymax, axisunits="none", xlines=FALSE, lty=lty[i],col=cols[[i]]); graphics::par(new=TRUE)
			}
		}
		xaxis(axisunits[1],lines=xlines); graphics::box(lwd=2)
	}

##########################################

#' @title Circular plot of GCMM means
#' @description Circular plot of GCMM means (circular intercepts) 
#' @param models List of one or more objects of class \code{GCMM} containing output from the \code{\link{GCMM}} function
#' @param col Vector of colours to use in the plot
#' @param axisunits Units to be used for the axis, either "radians", "sun", or "time"
#' @return Prints plot
#' @examples
#' \donttest{ FoxActivityGCMM<-GCMM(data=redfoxsample$Radians, 
#'         RE1=redfoxsample$SamplingPeriod, family="vonmises", autorun=FALSE,
#'         adapt=0, sample=300, burnin=300, thin=1, n.chains=2)
#'      HumanActivityGCMM<-GCMM(data=humanssample$Radians, RE1=humanssample$SamplingPeriod, 
#'         family="vonmises", autorun=FALSE, adapt=0, sample=300, burnin=300, thin=1, n.chains=2)
#'     circplotmeans(models=list(FoxActivityGCMM,HumanActivityGCMM)) }
#' 
#' @export

	circplotmeans<-function(models, col=c("cyan3","orchid","deeppink","dodgerblue"), axisunits=c("radians","sun","time") )	{

	graphics::par(mfrow=c(1,1))

	graphics::par(new=FALSE)
	for (m in 1:length(models) ) { 
		for (i in 1:models[[m]]$Nclust) {
			graphics::plot(circular::as.circular( seq(models[[m]]$output[i,2],models[[m]]$output[i,3],length=1000),
					type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),
					cex=1.4,shrink=1.7,start.sep=.15*m,axes=FALSE,tol=0.501,main="")
				graphics::points(circular::as.circular(seq(models[[m]]$output[i,2],models[[m]]$output[i,3],length=1000),
					type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),
					shrink=1.7,start.sep=.15*m,cex=.8,col=col[m])
				graphics::points(circular::as.circular(models[[m]]$output[i,1],
					type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),
					shrink=1.7,start.sep=.15*m,cex=1.2)
				graphics::par(new=TRUE)
			}
		}
			circaxis(axisunits)
	}


