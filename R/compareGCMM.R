#' Compare GCMM parameters or estimates
#' 
#' @description Compare two posterior distributions of GCMM parameters or estimates
#' @param model1 Object either of class \code{GCMM} with output from \code{\link{GCMM}} function or of class \code{GCMMestimate} with output from other activityGCMM functions
#' @param model2 Object of either class \code{GCMM}, \code{GCMMestimate}, or a vector or single value for which to compare with model 1; if left blank, it is assumed
#'     that both parameters from arguments \code{p1} and \code{p2} are from model1
#' @param p1 Name of first parameter to compare, from model1
#' @param p2 Second parameter to compare, either name of a parameter, estimate, or a vector of values
#' @param sample Number of posterior samples from which to build HDI; default=1000
#' @param plot Logical argument for whether to plot histograms of p1, p2, and the difference between them; default=TRUE
#' @return Returns object of class \code{GCMMestimate} containing a list, including \code{PD} containing the posterior distributions of p1, p2 and the difference between them,
#'    and summary information (HDIs and PDS)
#' @export

		compareGCMM<-function(model1,p1,model2="NULL",p2,sample=1000,plot=TRUE) { 

	if (class(model1)=="GCMM") {
		MCMC<-rbind(model1$runjags$mcmc[[1]],model1$runjags$mcmc[[2]],model1$runjags$mcmc[[3]])
			MCMC<-samplerows(MCMC,n=sample)
	if(sum(as.numeric(p1==colnames(MCMC)))==1) {
		message(paste("Extracting posterior samples of parameter",p1,".....")) 
			for (i in 1:ncol(MCMC) ) {
				if (p1==colnames(MCMC)[i]) {
				PD1<-MCMC[,i] } } 
		} else {
			if (is.character(p1)==TRUE) {
				message("No parameter found matching p1") } else {
			if (is.numeric(p1)==TRUE) {
				message("No parameter found matching p1; assumed to be an external value...")
				PD1<-rep(p1,length(MCMC))
			} }
		}
	}

	if (class(model1)=="GCMMestimate") {
		if(sum(as.numeric(p1==names(model1$PD)))==1) {
			message(paste("Extracting posterior samples of estimate",p1,".....")) 
			for (i in 1:length(model1$PD) ) { if (p1==names(model1$PD)[i]) { PD1<-model1$PD[i] } } 
		} else {
			if (is.character(p1)==TRUE) {
				message("No parameter found matching p1") } else {
			if (is.numeric(p1)==TRUE) {
				message("No parameter found matching p1; assumed to be an external value...")
				PD1<-rep(p1,length(MCMC))
			} } }
	}

	if (length(model2)==1) {
		if (is.numeric(p2)==TRUE) { 
			message("No parameter found matching p2; assumed to be an external value...")	
			PD2<-rep(p2,sample)
			}
		if(sum(as.numeric(p2==colnames(MCMC)))==1) {		
			message(paste("Extracting from model1 posterior samples of parameter",p2,".....")) 
			for (i in 1:ncol(MCMC) ) { if (p2==colnames(MCMC)[i]) { PD2<-MCMC[,i] } }
			} else {
				if (is.character(p2)==TRUE) {
					message("No parameter found matching p2") } 
				}
		}


	if (class(model2)=="GCMM") {
		MCMC2<-rbind(model2$runjags$mcmc[[1]],model2$runjags$mcmc[[2]],model2$runjags$mcmc[[3]])
			MCMC2<-samplerows(MCMC2,n=sample)
		if(sum(as.numeric(p2==colnames(MCMC2)))==1) {
			message(paste("Extracting posterior samples of parameter",p2,".....")) 
				for (i in 1:ncol(MCMC2) ) {
					if (p2==colnames(MCMC2)[i]) {
					PD2<-MCMC2[,i] } } 
			} else {
				if (is.character(p2)==TRUE) {
					message("No parameter found matching p2") } 
		} }

	if (class(model2)=="GCMMestimate") {
		if(sum(as.numeric(p2==names(model2$PD)))==1) {
			message(paste("Extracting posterior samples of estimate",p2,".....")) 
			for (i in 1:length(model2$PD) ) { if (p2==names(model2$PD)[i]) { PD2<-model2$PD[i] } } 
		} else {
			if (is.character(p2)==TRUE) {
				message("No parameter found matching p2") } 
		}	}

	if (class(PD1)=="list") { PD1<-PD1[[1]] }
	if (class(PD2)=="list") { PD2<-PD2[[1]] }

	Difference<-c(PD1-PD2)
	dEst<-HDI(Difference)
	names(dEst)<-c("Difference Mean","Lower95","Upper95")
	pds<-PDS(dEst)

	p1Est<-HDI(PD1)
	names(p1Est)<-c("p1 Mean","Lower95","Upper95")

	p2Est<-HDI(PD2)
	names(p2Est)<-c("p2 Mean","Lower95","Upper95")

	PD<-list(Difference=Difference, PD1=PD1, PD2=PD2)
	out<-list(PD=PD, p1Est=p1Est, p2Est=p2Est, dEst=dEst, PDS=pds)

	message("")
	print(dEst); message("")
	print(pds); message("");message("")
	message(paste("Proportion of posterior of difference > 0:",as.numeric(pds[1])))
	message(paste("Proportion of posterior of difference < 0:",as.numeric(pds[2])))
	message(paste("Posterior distribution support: ",as.numeric(pds[3])*100,"%",sep="")); message("");message("")

	if(plot==TRUE) { 
		graphics::par(mfrow=c(3,1),mar=c(4.2,4.1,2,2))
		graphics::hist(PD1,breaks=100,prob=TRUE,main="",xlab="P1 Estimate", col="aquamarine3",border=col)
			graphics::segments(x0=p1Est[2], y0=-.001, x1=p1Est[3], y1=-.001, lwd=5,col="grey15"); graphics::box(lwd=2)
		graphics::hist(PD2,breaks=100,prob=TRUE,main="", xlab="P2 Estimate", col="aquamarine3",border=col)
			graphics::segments(x0=p2Est[2], y0=-.001, x1=p2Est[3], y1=-.001, lwd=5,col="grey15"); graphics::box(lwd=2)
		graphics::hist(Difference,breaks=100,prob=TRUE,main="",xlab="Difference Estimate",col="lightseagreen",border="lightseagreen")
			graphics::segments(x0=dEst[2], y0=.001, x1=dEst[3], y1=.001, lwd=5,col="red2");graphics::box(lwd=2)
			graphics::abline(v=0,lwd=5,col="red3")			
		}

	class(out)<-"GCMMestimate"
	return(out)

	}