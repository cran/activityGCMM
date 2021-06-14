###############################################

#' Plot GCMM activity curve with random intercepts
#' @description Plot GCMM activity curve with random intercepts
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param RE1 Logical vector for whether to plot GCMM activity curve with random intercepts from RE1; default=TRUE
#' @param RE2 Logical vector for whether to plot GCMM activity curve with random intercepts from RE2; default=FALSE
#' @param scale Scale for plotting the activity curve, either "2pi" for 0,2pi or "pi" for -pi,pi; default="2pi"
#' @param ymax Value for upper limit of y axis
#' @param axisunits Units for x axis, either "radians", "time", "solar", "sun", or "none"; default="radians"
#' @return No return value; prints plot of GCMM activity curve with random intercepts
#' 
#' @examples
#' \donttest{ FoxGCMMREs<-GCMM(data=redfoxsample$Radians, RE1=redfoxsample$SamplingPeriod, 
#'        saveREs=TRUE, scale=c("2pi"), family="vonmises", autorun=FALSE,
#'        adapt=0, sample=300, burnin=300, thin=1)
#'     plotREs(FoxGCMMREs) }
#' @export

	plotREs<-function(model, RE1=TRUE, RE2=FALSE, scale="NULL", ymax="NULL", axisunits="radians") {

	if (RE1==TRUE) { if (RE2==TRUE) { graphics::par(mfrow=c(1,2)) }} else { 
		graphics::par(mfrow=c(1,1)) }

	if (scale=="NULL") { scale<-model$scale }
	if (scale=="pi") { lim1<-(-1*pi); lim2<-pi; xcenter<-c("midnight") }
	if (scale=="2pi") { lim1<-0; lim2<-2*pi; xcenter<-c("noon") }
	if (ymax=="NULL") { ymax<-yMax(list(model)) }


	if (RE1==TRUE) {
		a<-model$alpha1[,4]
		alphas<-list(); mus<-list()
		for (i in 1:model$Nclust) {
			alphas[[i]] <- a[seq(1,length(a),model$Nclust)]
			mus[[i]] <- model$Mu[i]+alphas[[i]]
			}

		graphics::par(new=FALSE)
		overlap::densityPlot(as.numeric(model$GCMMmixture), xcenter=xcenter, xscale=NA, yaxs="i", xaxs="i", extend=NA, lwd=2, adjust=5, 
			xlim = c(lim1, lim2), ylim=c(0,ymax), main = "",col="white");graphics::par(new=TRUE)
			graphics::abline(v = c(-pi, 0, pi, pi/2, 3 * pi/2, -pi/2), lty = 3, col = "gray80")

			gcmmcomponents<-c(); gcmmmixture<-c()
			for (j in 1:length(mus[[i]]) ) {
				for (i in 1:model$Nclust) {
						if (model$family=="wrappedcauchy") { 
						gcmmcomponents[[i]]<-circular::rwrappedcauchy(as.numeric(model$P[i])*10000, circular::circular(mus[[i]][j]), model$c[i])
						}
					if (model$family=="vonmises") { 
						gcmmcomponents[[i]]<-as.numeric(circular::rvonmises(as.numeric(model$P[i])*10000, circular::circular(mus[[i]][j]), model$c[i])  )
						}
					gcmmmixture<-c(gcmmmixture,gcmmcomponents[[i]])
					mix<-convertRad(gcmmmixture,to="scale") 
					}
				graphics::par(new=TRUE)
				overlap::densityPlot(as.numeric(mix), xcenter=xcenter, xscale=NA, yaxs="i", xaxs="i", extend=NA, lwd=1, adjust=5, 
					xlim = c(lim1, lim2), ylim=c(0,ymax), main = "",col="grey85")
					graphics::par(new=TRUE)
			gcmmmixture<-c(); gcmmcomponents<-c()
			}
			overlap::densityPlot(as.numeric(model$GCMMmixture), xcenter=xcenter, xscale=NA, yaxs="i", xaxs="i", extend=NA, lwd=2, adjust=5, 
				xlim = c(lim1, lim2), ylim=c(0,ymax), main = "")
				graphics::box(lwd=2)
		}

	if (RE2==TRUE) {
		a<-model$alpha2[,4]
		alphas<-list(); mus<-list()
		for (i in 1:model$Nclust) {
			alphas[[i]] <- a[seq(1,length(a),model$Nclust)]
			mus[[i]] <- model$Mu[i]+alphas[[i]]
			}

		graphics::par(new=FALSE)
		overlap::densityPlot(as.numeric(model$GCMMmixture), xcenter=xcenter, xscale=NA, yaxs="i", xaxs="i", extend=NA, lwd=2, adjust=5, 
			xlim = c(lim1, lim2), ylim=c(0,ymax), main = "",col="white");graphics::par(new=TRUE)
			graphics::abline(v = c(-pi, 0, pi, pi/2, 3 * pi/2, -pi/2), lty = 3, col = "gray80")
			gcmmcomponents<-c(); gcmmmixture<-c()
				for (j in 1:length(mus[[i]]) ) {
					for (i in 1:model$Nclust) {
					if (model$family=="wrappedcauchy") { 
						gcmmcomponents[[i]]<-circular::rwrappedcauchy(as.numeric(model$P[i])*10000, circular::circular(mus[[i]][j]), model$c[i])
						}
					if (model$family=="vonmises") { 
						gcmmcomponents[[i]]<-as.numeric(circular::rvonmises(as.numeric(model$P[i])*10000, circular::circular(mus[[i]][j]), model$c[i])  )
						}
					gcmmmixture<-c(gcmmmixture,gcmmcomponents[[i]])
					mix<-convertRad(gcmmmixture,to="scale") 
					}
				graphics::par(new=TRUE)
				overlap::densityPlot(as.numeric(mix), xcenter=xcenter, xscale=NA, yaxs="i", xaxs="i", extend=NA, lwd=1, adjust=5, 
					xlim = c(lim1, lim2), ylim=c(0,ymax), main = "",col="grey80")
					graphics::par(new=TRUE)
			gcmmmixture<-c(); gcmmcomponents<-c()
			}
			overlap::densityPlot(as.numeric(model$GCMMmixture), xcenter=xcenter, xscale=NA, yaxs="i", xaxs="i", extend=NA, lwd=2, adjust=5, 
				xlim = c(lim1, lim2), ylim=c(0,ymax), main = "")
				graphics::box(lwd=2)
		}
	}




#########################################################

#' Random Effects Circular Plot
#' @description Circular plot of GCMM random intercepts and 95% HDI
#' @param model Object of class \code{GCMM} with output from \code{\link{GCMM}} function
#' @param RE1 Logical vector for whether to plot GCMM activity curve with random intercepts from RE1; default=TRUE
#' @param RE2 Logical vector for whether to plot GCMM activity curve with random intercepts from RE2; default=FALSE
#' @param axisunits Units for x axis, either "radians", "time", "solar", "sun", or "none"; default="radians"
#' @return No return value; prints circle plot of GCMM random intercepts and 95% HDI
#' @examples
#' \donttest{ FoxGCMMREs<-GCMM(data=redfoxsample$Radians, RE1=redfoxsample$SamplingPeriod, 
#'        saveREs=TRUE, scale=c("2pi"), family="vonmises", autorun=FALSE,
#'        adapt=0, sample=300, burnin=300, thin=1)
#'     circplotREs(FoxGCMMREs, axisunits="sun") }
#' @export

	circplotREs<-function(model, RE1=TRUE, RE2=FALSE, axisunits=c("radians","sun","time","solar","none"))  {

	if(axisunits=="radians") { labels<-c(0,"pi/2","pi","3pi/2",0.5) }
	if(axisunits=="time") { labels<-c("24:00","6:00","12:00","18:00",0.4) }
	if(axisunits=="sun") { labels<-c("midnight","sunrise","noon","sunset",0.3) }
	if(axisunits=="solar") { labels<-c("solar midnight","sunrise","solar noon","sunset",0.3) }
	if(axisunits=="none") { labels<-c("",0.3) }
		graphics::par(mfrow=c(1,1))

	if (RE1==TRUE) {
		a<-model$alpha1[,4]
		alphas<-list(); mus<-list()
		for (i in 1:model$Nclust) {
			alphas[[i]] <- a[seq(1,length(a),by=model$Nclust)]
			mus[[i]] <- model$Mu[i]+alphas[[i]]
		sep<-seq(.2,5.9,length=length(a)+1)
			# windows() 			# only works in windows
			grDevices::dev.new() 		# might not work with RStudio
			graphics::plot(circular::as.circular(model$Mu[i], type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),
				shrink=2.5,start.sep=.07,axes=FALSE,tol=0.501,main="RE1")
				graphics::points(circular::as.circular(seq(model$output[i,2],model$output[i,3],length=1000),
					type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),shrink=3,start.sep=.04,cex=.4)
				C<-paste("C",i,sep="")
				graphics::text(x=0, y=0, C, font=2, cex=0.7)
					graphics::text(x=0, y=.7, labels[1], cex=as.numeric(labels[5]))
					graphics::text(x=.7, y=0, labels[2], cex=as.numeric(labels[5]))
					graphics::text(x=0, y=-.7, labels[3], cex=as.numeric(labels[5]))
					graphics::text(x=-.7, y=.0, labels[4], cex=as.numeric(labels[5]))
			for (j in 1:length(mus[[i]]) ) {
				graphics::points(circular::as.circular(seq(model$Mu[i]+model$alpha1[j,1],model$Mu[i]+model$alpha1[j,3],length=10000), 
					type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),start.sep=sep[j+1],cex=.2,col="grey75",type = "l")#pch=".")
				graphics::points(circular::as.circular(model$Mu[i]+model$alpha1[j,4], type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"), 
					start.sep=sep[j+1],col="grey15",cex=.4)
				}
			graphics::points(circular::as.circular(seq(model$output[i,2],model$output[i,3],length=1000),
				type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"), shrink=3,start.sep=.04,cex=.4)
			graphics::points(circular::as.circular(model$Mu[i],type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),shrink=3,start.sep=.03,cex=.5)
			}
		}

	if (RE2==TRUE) {
		a<-model$alpha2[,4]
		alphas<-list(); mus<-list()
		for (i in 1:model$Nclust) {
			alphas[[i]] <- a[seq(1,length(a),model$Nclust)]
			mus[[i]] <- model$Mu[i]+alphas[[i]]
				sep<-seq(.2,5.9,length=length(a)+1)
			# windows() 
			grDevices::dev.new()
			graphics::plot(circular::as.circular(model$Mu[i],type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),shrink=2.5,start.sep=.07,axes=FALSE,tol=0.501,main="RE2")
				graphics::points(circular::as.circular(seq(model$output[i,2],model$output[i,3],length=1000),type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),shrink=3,start.sep=.04,cex=.4)
				C<-paste("C",i,sep=""); graphics::text(x=0, y=0, C, font=2, cex=0.7)
					graphics::text(x=0, y=.7, labels[1], cex=as.numeric(labels[5]))
					graphics::text(x=.7, y=0, labels[2], cex=as.numeric(labels[5]))
					graphics::text(x=0, y=-.7, labels[3], cex=as.numeric(labels[5]))
					graphics::text(x=-.7, y=.0, labels[4], cex=as.numeric(labels[5]))
			for (j in 1:length(mus[[i]]) ) {
					graphics::points(circular::as.circular(seq(model$Mu[i]+model$alpha1[j,1],model$Mu[i]+model$alpha1[j,3],length=10000),
						type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),start.sep=sep[j+1],cex=.2,col="grey75",type = "l")#pch=".")
					graphics::points(circular::as.circular(model$Mu[i]+model$alpha1[j,4],
						type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),
						start.sep=sep[j+1],col="grey15",cex=.4)
					}
				graphics::points(circular::as.circular(seq(model$output[i,2],model$output[i,3],length=1000),type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),
					shrink=3,start.sep=.04,cex=.4)
				graphics::points(circular::as.circular(model$Mu[i],type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="clock"),shrink=3,start.sep=.03,cex=.5)
				}
			}


	}


