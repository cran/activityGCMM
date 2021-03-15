#' @title Generalized circular mixed effect mixture models (GCMM) 
#' 
#' @description Bayesian parametric generalized circular mixed effect mixture models (GCMM) for estimating animal activity curves 
#'     from camera traps and other nested data structures using JAGS. Data distributions currently supported include 
#'     von Mises and wrapped Cauchy, with one or two random effects fit as random circular intercepts. The GCMM function
#'     automatically selects the number of components for the mixture model (supporting up to 4 mixture components)
#'     and runs the model in 'JAGS' through R. The number of clusters can also be manually selected. The function 
#'     returns the model summary and the activity curve estimated from the circular mixture model, with additional 
#'     information from the analysis provided in the output as a list of class \code{GCMM}.
#' 
#' Package: activityGCMM
#' Version: 1.0.0
#' Date: 2020-11-08 
#' Author: Liz AD campbell
#' 
#' @details The number of clusters is automatically selected based on a Bayesian linear finite normal mixture model via
#'     the \code{mclust} package. The Bayesian parametric GCMM is fit using 'JAGS' through R using the \code{runjags}
#'     package.   
#'
#' @param data Vector of observations in radians (0 to 2pi)
#' @param RE1 Vector identifying random effect for observations (e.g. camera trap ID)
#' @param RE2 Vector identifying second random effect for observations (e.g. study site, year, season, sampling period)
#' @param scale Scale of observations, either 0 to 2pi ("2pi") or -pi to pi ("pi")
#' @param family Probability distribution, either "vonmises" or "wrappedcauchy"
#' @param kmax Maximum number to test for vonmises kappa parameter; default=15
#' @param autojags Whether to use autorun.jags from runjags to automatically run model until convergence of MCMC chains and minimum effective sample size achieved; default=TRUE
#' @param maxtime Maximum time to let autorun.jags run for, can be input as minutes (e.g. "45m") or hours (e.g. "2h"); default="30m", increase for more complex models or larger datasets
#' @param thin Thinning rate for MCMC chains, i.e. how many samples are saved. For longer models, thin can be increased to reduce computer memory requirements
#' @param sample If autojags=FALSE, the number of MCMC samples per chain (which is multiplied by thin); default=10000
#' @param burnin If autojags=FALSE, the burnin for the MCMC chains which are not saved; default=5000
#' @param adapt adaptation to use for MCMC chains; default=1000
#' @param n.chains number of MCMC chains; default=3
#' @param Nclust Number of components for mixture models; if not provided, the function will estimate the number of clusters; if provided, values must be provided for clustmeans
#' @param clustmeans A vector equal in length to Nclust of the potential means for each component in the mixture models
#' @param saveREs Whether random intercepts are saved in output; recommended to save only one of saveREs, saveResids or saveYExp at one time due to memory limitations
#' @param saveclustIDs Whether to save component cluster identification for the data points; default=FALSE
#' @param saveResids Whether model residuals are saved in output; recommended to save only one of saveREs, saveResids or saveYExp at one time due to memory limitations
#' @param saveYExp Whether expected Y values based on model are saved in output; recommended to save only one of saveREs, saveResids or saveYExp at one time due to memory limitations
#' @param saveJAGS Logical argument of whether to save runjags output; default=FALSE
#' @param plot.type Which MCMC diagnostic plots from runjags; options: "trace" for MCMC traceplots, "hist" for histograms of posterior distributions; default="trace"
#' @param printvars Which of the saved variables to print in the output (the others will be saved in the summary); default does not include random intercepts, residuals, or expected Y
#' @param plotvars Which of the saved variables to plot in the diagnostic plots; recommended to avoid including random intercepts, residuals, or expected Y as each produces a separate plot; default: intercepts, concentration paramters ("K" for vonmises and "rho" for wrapped cauchy), cluster weights
#' 
#' @returns Returns object of class \code{GCMM} which is a list containing analysis results and details. A plot of the estimated activity curve from the mixed effect mixture model is printed.
#' @returns \code{output} GCMM model output summer
#' @returns \code{GCMMmixture} Vectors of simulated values from mixture model
#' @returns \code{GCMMcomponents} Vectors of simulated values from each component in the mixture model
#' @returns \code{runjags} GCMM model output from JAGS of class \code{runjags} from \code{runjags} package; see \code{\link[runjags]{run.jags}}
#' 
#' @keywords cameratrap; circular; vonmises; wrappedcauchy; mixture; activity; temporal; Bayesian
#' @author Liz AD Campbell
#' 
#' @importFrom mclust Mclust mclustBIC
#' 
#' @examples
#' data(redfoxsample)
#' \donttest{ FoxActivityGCMM<-GCMM(data=redfoxsample$Radians, RE1=redfoxsample$SamplingPeriod, 
#'     scale=c("2pi"), family="vonmises", autojags=FALSE,
#'     adapt=0, sample=300, burnin=300, thin=1,n.chains=2  ) }
#'  
#' @export
	GCMM<-function(data, RE1, RE2=NULL, scale="2pi", kmax=15, family=c("vonmises","wrappedcauchy"),
		autojags=TRUE, maxtime="30m", thin=2, burnin=5000, sample=5000, adapt=1000, n.chains=3,
		saveREs=FALSE, saveResids=FALSE, saveclustIDs=FALSE, saveYExp=FALSE, saveJAGS=TRUE, Nclust="NULL", clustmeans=NULL,
		plot.type=c("trace"), printvars=monitor, plotvars=c("CircularIntercept","ClustProb","sigma")  )  {

message(""); message("activityGCMM: Bayesian generalized circular mixed effect mixture models for estimating animal temporal activity curves")
message("Bugs or comments can be reported to lizadcampbell@atlasgoldenwolf.org"); message("")

GCMMout<-NULL
requireNamespace("mclust",quietly=TRUE)

if (Nclust=="NULL") { 
	message("Identifying clusters....")

	# Model-based k-clustering with mclust:
	# fit data (assuming 0-2pi)
		mclustfit <- mclust::Mclust(data)
			S2pi<-summary(mclustfit,parameters=TRUE)
			clustmeans2pi<-S2pi$mean
			Nclusts2pi<-length(clustmeans2pi)

	# check for circular data; converting to [-pi,pi] to check number of clusters:
		datpirad<-data; for (i in 1:length(data)) { if (data[i]>pi) { datpirad[i] <- data[i]-2*pi } }
		invisible(mclustfitpi <- mclust::Mclust(datpirad))
			Spi<-summary(mclustfitpi,parameters=TRUE) 
			clustmeanspi<-Spi$mean
			Nclustspi<-length(clustmeanspi)

	# selecting smaller: 
	Nclust<-min(Nclusts2pi, Nclustspi)
	if (Nclusts2pi <= Nclustspi) { 
		clustmeans<-clustmeans2pi; scale<-"2pi" } else {
		clustmeans<-clustmeanspi; scale<-"pi"; dat<-datpirad }

	message(paste("Number of suggested clusters:",Nclust))
	message(paste("Suggested scale for data:",scale)); message("")

	if (Nclusts2pi <= Nclustspi) {  mclustout<-mclustfit } else { mclustout<-mclustfitpi }

	if(Nclust>4) { 
	message("The suggested number of clusters is >4. This may result from small sample sizes. activityGCMM currently") 
	message("supports mixture models with up to 4 components. The next-best supported number of clusters are shown below.") 

	# Selecting next best number of clusters
		if (Nclusts2pi <= Nclustspi) {  
				BIC<-mclustfit$BIC } else {	 
				BIC<-mclustfitpi$BIC }		
		BICn<-BIC[-Nclust,]
			nextbestE<-as.numeric( names( sort(BICn[,1],decreasing=TRUE) )[1] )
			nextbestEBIC<- sort(BICn[,1],decreasing=TRUE) [1]
			nextbestV<-as.numeric( names( sort(BICn[,2],decreasing=TRUE) )[1] ) 
			nextbestVBIC<- sort(BICn[,2],decreasing=TRUE) [1]
		if(nextbestEBIC<nextbestVBIC) { nextbest<- nextbestE } else { nextbest<-nextbestV }

	print(BIC)
	message(""); message(paste("Next-best supported number of clusters: ",nextbest))

	fit <- stats::kmeans(data, nextbest) 
		cmeans<-stats::aggregate(data,by=list(fit$cluster),FUN=mean)
			names(cmeans)<-c("Cluster","Mean")
	Nclust<-nextbest
	clustmeans<-cmeans[,2]

		if(nextbest>4) {
			BICn<-BIC[1:4,]
				nextbestE<-as.numeric( names( sort(BICn[,1],decreasing=TRUE) )[1] )
				nextbestEBIC<- sort(BICn[,1],decreasing=TRUE) [1]
				nextbestV<-as.numeric( names( sort(BICn[,2],decreasing=TRUE) )[1] ) 
				nextbestVBIC<- sort(BICn[,2],decreasing=TRUE) [1]
				if(nextbestEBIC<nextbestVBIC) { nextbest<- nextbestE } else { nextbest<-nextbestV }
		message(""); message("The next-best supported number of clusters is also >4")
		message(paste("Best supported number of clusters <5: ",nextbest))
		message(paste("Recommended scale: ",as.character(scale))); message("")
			fit <- stats::kmeans(data, nextbest) 
			cmeans<-stats::aggregate(data,by=list(fit$cluster),FUN=mean)
				names(cmeans)<-c("Cluster","Mean")
		Nclust<-nextbest
		clustmeans<-cmeans[,2]
		}

	}


} else { clustmeans<-clustmeans; mclustout<-NULL }


	message(""); message("Running circular model...."); message("")

	# creating vector for candidate clusters:
	clust=rep(NA,length(data)) 
	for (c in 1:Nclust) {
		diff<-abs(data-as.numeric(clustmeans[c]))
			h1<-utils::head(sort(diff),2)
			for (i in 1:length(data)) { if (h1[1]==diff[i]) { clust[i] <- c } }
			for (i in 1:length(data)) { if (h1[2]==diff[i]) { clust[i] <- c } }
		}


	I0 <- function(kappa) besselI(kappa,0) # Bessel function of the first order
		kappas <- 1:kmax; I0_s <- sapply(kappas, I0)

	RE1<-as.numeric(as.factor(as.numeric(RE1)))
		if(length(RE2)>0) { RE2<-as.numeric(as.factor(as.numeric(RE2))) }

	if (scale=="pi") { lim1<-(-1*pi); lim2<-pi }
		if (scale=="2pi") { lim1<-0; lim2<-2*pi }

	if (family=="vonmises") {
		if (length(RE2)>0) {
			JAGSdata=list(y=data, N=length(data), C=1000000, z=rep(0,length(data)), pi=3.14159, 
				I0s=I0_s, kappas=kappas, p=rep(length=kmax, 1/kmax), onesRepNclust=rep(1,Nclust), lim1=lim1, lim2=lim2,	 
				RE1=as.numeric(RE1), nRE1=length(unique(RE1)),
				RE2=as.numeric(RE2), nRE2=length(unique(RE2)),
				Nclust=Nclust, clust=clust)
			} else {
			JAGSdata=list(y=data, N=length(data), C=1000000, z=rep(0,length(data)), pi=3.14159, 
				I0s=I0_s, kappas=kappas, p=rep(length=kmax, 1/kmax), onesRepNclust=rep(1,Nclust),	 
				RE1=as.numeric(RE1), nRE1=length(unique(RE1)), lim1=lim1, lim2=lim2,
				Nclust=Nclust, clust=clust)
			}
		}
	if (family=="wrappedcauchy") {
		if (length(RE2)>0) {
			JAGSdata=list(y=data, N=length(data), C=1000000, z=rep(0,length(data)), pi=3.14159, 
				onesRepNclust=rep(1,Nclust), lim1=lim1, lim2=lim2,	 
				RE1=as.numeric(RE1), nRE1=length(unique(RE1)),
				RE2=as.numeric(RE2), nRE2=length(unique(RE2)),
				Nclust=Nclust, clust=clust)
			} else {
			JAGSdata=list(y=data, N=length(data), C=1000000, z=rep(0,length(data)), pi=3.14159,  
				onesRepNclust=rep(1,Nclust),  lim1=lim1, lim2=lim2,	 
				RE1=as.numeric(RE1), nRE1=length(unique(RE1)),
				Nclust=Nclust, clust=clust)
			}
		}

VM2REs1C = " model {					   
    for(i in 1:N) {
       phi[i] <- -K*cos(y[i]-mu[i]) + log(2*pi*IO_hat) + C 
        	z[i] ~ dpois(phi[i])						   
		mu[i]<-CircularIntercept + alpha1[RE1[i]] + alpha2[RE2[i]]
 	}
      CircularIntercept ~ dunif(lim1,lim2)
		k ~ dcat(p[])		
   	      K <- kappas[k]
	      IO_hat <- I0s[k]
   	 for (r in 1:nRE1) { alpha1[r] ~  dnorm(0, tau1) } 
	 for (r in 1:nRE2) { alpha2[r] ~  dnorm(0, tau2) }
		sigma1 ~ dunif(0, 1.5) 
		    	tau1 <- 1/(sigma1*sigma1)
		sigma2 ~ dunif(0, 1.5)
		    	tau2 <- 1/(sigma2*sigma2)
	for(i in 1:N) {
		YExp[i]<-CircularIntercept + alpha1[RE1[i]] + alpha2[RE2[i]] 
		Resid[i]<-y[i]-YExp[i]
 		} 
}"

VM1RE1C = " model {					   
    for(i in 1:N) {
       phi[i] <- -K*cos(y[i]-mu[i]) + log(2*pi*IO_hat) + C 
        	z[i] ~ dpois(phi[i])						   
		mu[i]<-CircularIntercept + alpha1[RE1[i]]
 	}
      CircularIntercept ~ dunif(lim1,lim2)
		k ~ dcat(p[])		
   	      K <- kappas[k]
	      IO_hat <- I0s[k]
 	 for (r in 1:nRE1) { alpha1[r] ~  dnorm(0, tau1) } 
		sigma1 ~ dunif(0, 1.5) 
		    	tau1 <- 1/(sigma1*sigma1)
	for(i in 1:N) {
		YExp[i]<-CircularIntercept + alpha1[RE1[i]] 
		Resid[i]<-y[i]-YExp[i]
 		} 
   }"

WC2REs1C = " model {					   
    for(i in 1:N) { 
		phi[i] <- -log((1-pow(rhoC,2)) / (2*pi*(1+pow(rhoC,2)-2*rhoC*cos(y[i]-mu[i]))) ) + C
        	z[i] ~ dpois(phi[i])						   
		mu[i]<-CircularIntercept + alpha1[RE1[i]] + alpha2[RE2[i]]
 	}
      CircularIntercept ~ dunif(lim1,lim2)
		rhoC ~ dunif(0,1)
   	 for (r in 1:nRE1) { alpha1[r] ~  dnorm(0, tau1) } 
	 for (r in 1:nRE2) { alpha2[r] ~  dnorm(0, tau2) }
		sigma1 ~ dunif(0, 1.5) 
		    	tau1 <- 1/(sigma1*sigma1)
		sigma2 ~ dunif(0, 1.5)
		    	tau2 <- 1/(sigma2*sigma2) 
	for(i in 1:N) {
		YExp[i]<-CircularIntercept + alpha1[RE1[i]] + alpha2[RE2[i]] 
		Resid[i]<-y[i]-YExp[i]
 		} 
}"

WC1RE1C = " model {					   
    for(i in 1:N) {
		phi[i] <- -log((1-pow(rhoC,2)) / (2*pi*(1+pow(rhoC,2)-2*rhoC*cos(y[i]-mu[i]))) ) + C
        	z[i] ~ dpois(phi[i])						   
		mu[i]<-CircularIntercept + alpha1[RE1[i]]
 	}
      CircularIntercept ~ dunif(lim1,lim2)
		rhoC ~ dunif(0,1)
 	for (r in 1:nRE1) { alpha1[r] ~  dnorm(0, tau1) } 
		sigma1 ~ dunif(0, 1.5) 
		    	tau1 <- 1/(sigma1*sigma1)
	for(i in 1:N) {
		YExp[i]<-CircularIntercept + alpha1[RE1[i]]  
		Resid[i]<-y[i]-YExp[i]
 		} 
   }"


VM2REs2Cs = " model { 
	for(i in 1:N) {
      	phi[i] <- -kappa_hat[i]*cos(y[i]-mu[i]) + log(2*pi*IO_hat[i]) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		IO_hat[i]<-IOh[clust[i]]
		kappa_hat[i]<-K[clust[i]]
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
 	}
	for (clustIdx in 1:Nclust) {
		k[clustIdx] ~ dcat(p[])	
		K[clustIdx] <- kappas[k[clustIdx]]
		IOh[clustIdx] <- I0s[k[clustIdx]]
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1)
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx])
		sigma2[clustIdx] ~ dunif(0,1)
			tau2[clustIdx]  <- 1/(sigma2[clustIdx]*sigma2[clustIdx]) 
	}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2]) 	 }
 	for (r in 1:nRE2) { 
		alpha2[1, r] ~  dnorm(0, tau2[1]) 
		alpha2[2, r] ~  dnorm(0, tau2[2])  }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)
 
	for(i in 1:N) {
		YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
		Resid[i]<-y[i]-YExp[i]
	 	}
}"


VM1RE2Cs = " model {	
	for(i in 1:N) {
      	phi[i] <- -kappa_hat[i]*cos(y[i]-mu[i]) + log(2*pi*IO_hat[i]) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		IO_hat[i]<-IOh[clust[i]]
		kappa_hat[i]<-K[clust[i]]
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] 
 	}
	for (clustIdx in 1:Nclust) {
		k[clustIdx] ~ dcat(p[])	
		K[clustIdx] <- kappas[k[clustIdx]]
		IOh[clustIdx] <- I0s[k[clustIdx]]
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1) 
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx]) 
		}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2]) 	 }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)

	for(i in 1:N) {
	YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] 
	Resid[i]<-y[i]-YExp[i]
 	}
}"



VM2REs3Cs = " model {
	for(i in 1:N) {
      	phi[i] <- -kappa_hat[i]*cos(y[i]-mu[i]) + log(2*pi*IO_hat[i]) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		IO_hat[i]<-IOh[clust[i]]
		kappa_hat[i]<-K[clust[i]]
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
 	} 
	for (clustIdx in 1:Nclust) {
		k[clustIdx] ~ dcat(p[])	
		K[clustIdx] <- kappas[k[clustIdx]]
		IOh[clustIdx] <- I0s[k[clustIdx]]
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1) 
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx])
		sigma2[clustIdx] ~ dunif(0,1)
			tau2[clustIdx]  <- 1/(sigma2[clustIdx]*sigma2[clustIdx]) 
	}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2]) 
		alpha1[3, r] ~  dnorm(0, tau1[3])	 }
 	for (r in 1:nRE2) { 
		alpha2[1, r] ~  dnorm(0, tau2[1]) 
		alpha2[2, r] ~  dnorm(0, tau2[2])
		alpha2[3, r] ~  dnorm(0, tau2[3])  }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)

	for(i in 1:N) {
	YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
	Resid[i]<-y[i]-YExp[i]
 	}
}"

VM2REs4Cs = " model {
	for(i in 1:N) {
      	phi[i] <- -kappa_hat[i]*cos(y[i]-mu[i]) + log(2*pi*IO_hat[i]) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		IO_hat[i]<-IOh[clust[i]]
		kappa_hat[i]<-K[clust[i]]
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
 	}
	for (clustIdx in 1:Nclust) {
		k[clustIdx] ~ dcat(p[])	
		K[clustIdx] <- kappas[k[clustIdx]]
		IOh[clustIdx] <- I0s[k[clustIdx]]
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1)
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx])
		sigma2[clustIdx] ~ dunif(0,1)
			tau2[clustIdx]  <- 1/(sigma2[clustIdx]*sigma2[clustIdx]) 
	}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2]) 
		alpha1[3, r] ~  dnorm(0, tau1[3])
		alpha1[4, r] ~  dnorm(0, tau1[4])	 }
 	for (r in 1:nRE2) { 
		alpha2[1, r] ~  dnorm(0, tau2[1]) 
		alpha2[2, r] ~  dnorm(0, tau2[2])
		alpha2[3, r] ~  dnorm(0, tau2[3])
		alpha2[4, r] ~  dnorm(0, tau2[4])  }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)
 
	for(i in 1:N) {
	YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
	Resid[i]<-y[i]-YExp[i]
 	}
}"

VM1RE3Cs = " model {
	for(i in 1:N) {
      	phi[i] <- -kappa_hat[i]*cos(y[i]-mu[i]) + log(2*pi*IO_hat[i]) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		IO_hat[i]<-IOh[clust[i]]
		kappa_hat[i]<-K[clust[i]]
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] 
 	}
	for (clustIdx in 1:Nclust) {
		k[clustIdx] ~ dcat(p[])	
		K[clustIdx] <- kappas[k[clustIdx]]
		IOh[clustIdx] <- I0s[k[clustIdx]]
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1) 
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx]) 
		}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2])
		alpha1[3, r] ~  dnorm(0, tau1[3]) 	 }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)

	for(i in 1:N) {
	YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] 
	Resid[i]<-y[i]-YExp[i]
 	}
}"

VM1RE4Cs = " model {	
	for(i in 1:N) {
      	phi[i] <- -kappa_hat[i]*cos(y[i]-mu[i]) + log(2*pi*IO_hat[i]) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		IO_hat[i]<-IOh[clust[i]]
		kappa_hat[i]<-K[clust[i]]
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] 
 	}
	for (clustIdx in 1:Nclust) {
		k[clustIdx] ~ dcat(p[])	
		K[clustIdx] <- kappas[k[clustIdx]]
		IOh[clustIdx] <- I0s[k[clustIdx]]
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1) 
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx]) 
		}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2])
		alpha1[3, r] ~  dnorm(0, tau1[3])
		alpha1[4, r] ~  dnorm(0, tau1[4]) 	 }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)

	for(i in 1:N) {
	YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] 
	Resid[i]<-y[i]-YExp[i]
 	}
}"


WC2REs3Cs = " model {	
	for(i in 1:N) {
		phi[i] <- -log((1-pow(rho[i],2)) / (2*pi*(1+pow(rho[i],2)-2*rho[i]*cos(y[i]-mu[i]))) ) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		rho[i]<-rhoC[clust[i]]      
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
 	}
	for (clustIdx in 1:Nclust) {
		rhoC[clustIdx] ~ dunif(0,1)		
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1) 
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx])
		sigma2[clustIdx] ~ dunif(0,1) 
			tau2[clustIdx]  <- 1/(sigma2[clustIdx]*sigma2[clustIdx]) 
	}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2]) 
		alpha1[3, r] ~  dnorm(0, tau1[3])	 }
 	for (r in 1:nRE2) { 
		alpha2[1, r] ~  dnorm(0, tau2[1]) 
		alpha2[2, r] ~  dnorm(0, tau2[2])
		alpha2[3, r] ~  dnorm(0, tau2[3])  }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)

	for(i in 1:N) {
		YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
		Resid[i]<-y[i]-YExp[i]
 	}
}"

WC2REs4Cs = " model {	
	for(i in 1:N) {
		phi[i] <- -log((1-pow(rho[i],2)) / (2*pi*(1+pow(rho[i],2)-2*rho[i]*cos(y[i]-mu[i]))) ) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		rho[i]<-rhoC[clust[i]]      
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
 	}
	for (clustIdx in 1:Nclust) {
		rhoC[clustIdx] ~ dunif(0,1)		
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1) 
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx])
		sigma2[clustIdx] ~ dunif(0,1) 
			tau2[clustIdx]  <- 1/(sigma2[clustIdx]*sigma2[clustIdx]) 
	}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2]) 
		alpha1[3, r] ~  dnorm(0, tau1[3])
		alpha1[4, r] ~  dnorm(0, tau1[4])	 }
 	for (r in 1:nRE2) { 
		alpha2[1, r] ~  dnorm(0, tau2[1]) 
		alpha2[2, r] ~  dnorm(0, tau2[2])
		alpha2[3, r] ~  dnorm(0, tau2[3]) 
		alpha2[4, r] ~  dnorm(0, tau2[4])  }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)

	for(i in 1:N) {
		YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
		Resid[i]<-y[i]-YExp[i]
 	}
}"


WC2REs2Cs = " model {	
	for(i in 1:N) {
		phi[i] <- -log((1-pow(rho[i],2)) / (2*pi*(1+pow(rho[i],2)-2*rho[i]*cos(y[i]-mu[i]))) ) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		rho[i]<-rhoC[clust[i]]     
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
 	}
	for (clustIdx in 1:Nclust) {
		rhoC[clustIdx] ~ dunif(0,1)			
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1) 
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx])
		sigma2[clustIdx] ~ dunif(0,1) 
			tau2[clustIdx]  <- 1/(sigma2[clustIdx]*sigma2[clustIdx]) 
	}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2]) 	 }
 	for (r in 1:nRE2) { 
		alpha2[1, r] ~  dnorm(0, tau2[1]) 
		alpha2[2, r] ~  dnorm(0, tau2[2])  }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)

	for(i in 1:N) {
		YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] + alpha2[clust[i],RE2[i]] 
		Resid[i]<-y[i]-YExp[i]
 	}
}"



WC1RE3Cs = " model {	
	for(i in 1:N) {
		phi[i] <- -log((1-pow(rho[i],2)) / (2*pi*(1+pow(rho[i],2)-2*rho[i]*cos(y[i]-mu[i]))) ) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		rho[i]<-rhoC[clust[i]]
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]]  
 	}
	for (clustIdx in 1:Nclust) {
		rhoC[clustIdx] ~ dunif(0,1)			
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1) 
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx])
	}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2]) 
		alpha1[3, r] ~  dnorm(0, tau1[3])	 }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)

	for(i in 1:N) {
		YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]]
		Resid[i]<-y[i]-YExp[i]
 	}
}"


WC1RE4Cs = " model {	
	for(i in 1:N) {
		phi[i] <- -log((1-pow(rho[i],2)) / (2*pi*(1+pow(rho[i],2)-2*rho[i]*cos(y[i]-mu[i]))) ) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		rho[i]<-rhoC[clust[i]]      
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]]
 	}
	for (clustIdx in 1:Nclust) {
		rhoC[clustIdx] ~ dunif(0,1)			
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1) 
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx])
	}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2]) 
		alpha1[3, r] ~  dnorm(0, tau1[3])
		alpha1[4, r] ~  dnorm(0, tau1[4])	 }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)

	for(i in 1:N) {
		YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]]
		Resid[i]<-y[i]-YExp[i]
 	}
}"


WC1RE2Cs = " model {	
	for(i in 1:N) {
		phi[i] <- -log((1-pow(rho[i],2)) / (2*pi*(1+pow(rho[i],2)-2*rho[i]*cos(y[i]-mu[i]))) ) + C 
      	z[i] ~ dpois(phi[i])						   
		clust[i] ~ dcat(ClustProb[1:Nclust])
		rho[i]<-rhoC[clust[i]]      
		mu[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]] 
 	}
	for (clustIdx in 1:Nclust) {
		rhoC[clustIdx] ~ dunif(0,1)		
		CircularIntercept[clustIdx] ~ dunif(lim1,lim2)
		sigma1[clustIdx] ~ dunif(0,1) 
			tau1[clustIdx]  <- 1/(sigma1[clustIdx]*sigma1[clustIdx])
	}
	for (r in 1:nRE1) { 
		alpha1[1, r] ~  dnorm(0, tau1[1]) 
		alpha1[2, r] ~  dnorm(0, tau1[2]) 	 }
	ClustProb[1:Nclust] ~ ddirch(onesRepNclust)

	for(i in 1:N) {
		YExp[i]<-CircularIntercept[clust[i]] + alpha1[clust[i],RE1[i]]
		Resid[i]<-y[i]-YExp[i]
 	}
}"


	### MCMC chain initialization: 
	initsVM1RE <- function () { list( CircularIntercept=stats::rnorm(Nclust,clustmeans,.001), 
		ClustProb=stats::runif(Nclust,0,1), sigma1=stats::runif(Nclust,0,.5)	 ) }
	initsVM2REs <- function () { list( CircularIntercept=stats::rnorm(Nclust,clustmeans,.001), 
		ClustProb=stats::runif(Nclust,0,1), sigma1=stats::runif(Nclust,0,.5), sigma2=stats::runif(Nclust,0,.5)	 ) }
	initsWC1RE <- function () { list( CircularIntercept=stats::rnorm(Nclust,clustmeans,.001) ) } 
	initsWC2REs <- function () { list( CircularIntercept=stats::rnorm(Nclust,clustmeans,.001) ) } 				


	### Specifying model & inits: 
   	if(Nclust==1) { 												
   		if (family=="vonmises") {									
   	   		if (length(RE2)>0) { 	model<-VM2REs1C; inits<-initsVM2REs } else { 	
   							model<-VM1RE1C; inits<-initsVM1RE }	}		
   		if (family=="wrappedcauchy") {								
   			if (length(RE2)>0) { 	model<-WC2REs1C; inits<-initsWC2REs } else { 
   							model<-WC1RE1C; inits<-initsWC1RE }	}		
   		}													
	if(Nclust==2) { 												
		if (family=="vonmises") {
			if (length(RE2)>0) {	model<-VM2REs2Cs; inits<-initsVM2REs } else { 
							model<-VM1RE2Cs; inits<-initsVM1RE }	}
		if (family=="wrappedcauchy") {
			if (length(RE2)>0) { 	model<-WC2REs2Cs; inits<-initsWC2REs } else { 
							model<-WC1RE2Cs; inits<-initsWC1RE } }
		}													
	if(Nclust==3) { 												
		if (family=="vonmises") {									
			if (length(RE2)>0) {	model<-VM2REs3Cs; inits<-initsVM2REs } else { 	
							model<-VM1RE3Cs; inits<-initsVM1RE }	}	
		if (family=="wrappedcauchy") {								
			if (length(RE2)>0) { 	model<-WC2REs3Cs; inits<-initsWC2REs } else { 
							model<-WC1RE3Cs; inits<-initsWC1RE } }		
		}													
	if(Nclust==4) { 												
		if (family=="vonmises") {									
			if (length(RE2)>0) { 	model<-VM2REs4Cs; inits<-initsVM2REs } else { 	
							model<-VM1RE4Cs; inits<-initsVM1RE }	}	
		if (family=="wrappedcauchy") {								
			if (length(RE2)>0) { 	model<-WC2REs4Cs; inits<-initsWC2REs } else { 
							model<-WC1RE4Cs; inits<-initsWC1RE } }		
		}													


	### Monitored variables:
	if (family=="vonmises") { monitor<-c("CircularIntercept","K","sigma1") }
	if (family=="wrappedcauchy") { monitor<-c("CircularIntercept","rhoC","sigma1") }
		if(Nclust>1) { monitor<-c(monitor,"ClustProb") }
		if (length(RE2)>0) { monitor<-c(monitor,"sigma2") }
		if (saveREs==TRUE) { 
			if (length(RE2)>0) { monitor<-c(monitor,"alpha1","alpha2") } else {
			monitor<-c(monitor,"alpha1") }	}
		if (saveResids==TRUE) { monitor<-c(monitor,"Resid") }
		if (saveYExp==TRUE) { monitor<-c(monitor,"YExp") }
		if (saveclustIDs==TRUE) { monitor<-c(monitor,"clust") }


	if(autojags==TRUE) {outGCMM<-runjags::autorun.jags(model=model, method="parallel", monitor=monitor, 
					n.chains=3, thin=thin, adapt=adapt, data=JAGSdata, inits=inits, max.time=maxtime)
		} else {	outGCMM<-runjags::run.jags(model=model, method="parallel", monitor=monitor, 
					n.chains=n.chains, thin=thin, adapt=adapt, data=JAGSdata, inits=inits, 
					burnin=burnin, sample=sample) 
		}

	out<-summary(outGCMM)
	output<-cbind(out[,4],out[,1],out[,3],out[,5],out[,11],out[,9])
	colnames(output)<-c("Mean","Lower95","Upper95","SD","rhat","ESS")
	output<-round(output,3)

	mu<-summary(outGCMM,vars="CircularIntercept"); mu<-round(mu[,4],2)
	if(Nclust>1) { p<-summary(outGCMM,vars="ClustProb"); p<-round(p[,4],2)  } else { p<-1 }
	if (family=="vonmises") { k<-summary(outGCMM,vars="K"); k<-round(k[,4],2) }	
	if (family=="wrappedcauchy") { rho<-summary(outGCMM,vars="rhoC"); rho<-round(rho[,4],2) }

	if(saveclustIDs==TRUE) { clustIDs<-summary(outGCMM,vars="clust"); clustIDs<-clustIDs[,6] } else { clustIDs<-NULL } 
	if(saveResids==TRUE) { E<-summary(outGCMM,vars="Resid"); E<-E[,4] } else { E<-NULL } 
	if(saveYExp==TRUE) { YExp<-summary(outGCMM,vars="YExp") } else { YExp<-NULL } 
	if (saveREs==TRUE) { alpha1<-summary(outGCMM,vars="alpha1")
			if (length(RE2)>0) { alpha2<-summary(outGCMM,vars="alpha2") } else { alpha2<-NULL }
	   	 } else { alpha1<-NULL; alpha2<-NULL } 
	if (saveJAGS==TRUE) { runjags<-outGCMM } else { runjags<-NULL }

	gcmmcomponents<-c()
	gcmmmixture<-c()
		for (i in 1:Nclust) { 
			if (family=="vonmises") { 
				gcmmcomponents[[i]]<-circular::rvonmises(as.numeric(p[i])*10000, circular::circular(mu[i]), k[i])  
				gcmmmixture<-c(gcmmmixture,gcmmcomponents[[i]])
				 }
			if (family=="wrappedcauchy") { 
				gcmmcomponents[[i]]<-circular::rwrappedcauchy(p[i]*10000, circular::circular(mu[i]), rho[i] )  
				gcmmmixture<-c(gcmmmixture,gcmmcomponents[[i]])
				}
		}


	### Mixture Plot: 
	overlap::densityPlot(as.numeric(gcmmmixture),xscale=NA,yaxs="i",extend=NA,lwd=3,adjust=5,xlim=c(0,2*pi),main="")
		graphics::abline(v=c(pi,pi/2,3*pi/2),lty=2,col="gray60"); graphics::box(lwd=2)

	if(family=="vonmises"){ c<-k }
	if(family=="wrappedcauchy"){ c<-rho }

	GCMMout<<-new.env()
	GCMMout$JAGSdata<<-JAGSdata			
	GCMMoutput<-list(family=family,data=data,RE1=RE1,RE2=RE2,Nclust=Nclust,scale=scale,Mu=mu,P=p,c=c,model=outGCMM,
		GCMMmixture=gcmmmixture,GCMMcomponents=gcmmcomponents,output=output,#summary(outGCMM),
		saveclustIDs=saveclustIDs,saveResids=saveResids,saveYExp=saveYExp,saveREs=saveREs,saveJAGS=saveJAGS,
		clustIDs=clustIDs,E=E,YExp=YExp,alpha1=alpha1,alpha2=alpha2,runjags=runjags)
		class(GCMMoutput) <- "GCMM"

	message("");message("")
	print(output)

	minESS<-min(output[,6],na.rm=TRUE)
	maxrhat<-max(output[,5],na.rm=TRUE)
	message(""); message(paste("All MCMC chains reached convergence (rhat<1.05):",maxrhat<1.05))
	message(paste("    Maximum rhat value:",round(maxrhat,4))); message("")
	message(paste("Minimum effective sample size:",round(minESS,0)))
	if(minESS<10000) { message("    NOTE: A minimum ESS of 10000 is recommended, consider increasing the number of samples") }

message("");message("----------------------")
message("GCMM analysis complete");message("");message("")

		return(GCMMoutput)

}




#' @title Executable example of GCMM function
#' @description Example of applying generalized circular mixed effect mixture model with activityGCMM 
#'     using data included in the package
#' @return Prints message with example of GCMM function using data included in the package
#' @examples { exampleGCMM() }
#' @export
	exampleGCMM<- function() { message("Try example with sample data included in the package: data(redfoxsample)")
			message("Run a generalized circular mixed effect mixture vonmises model using:")
			message("    GCMM(data=redfoxsample$Radians, RE1=redfoxsample$CameraTrapID, family=\"vonmises\",") 
     			message("           autojags=FALSE, adapt=0, burnin=500, sample=500)")	}



#' @title GCMM components plot
#' @description Plot of activity curves for the separate components in the circular mixture model
#' 
#' @param model Model output from GCMM function, object of class \code{GCMM}
#' @param rug Logical argument for whether to plot a rug of the raw values. Plotting the rug for the separate components requires that saveclustIDs=TRUE when running GCMM. default=FALSE
#' @param ruglwd Line width for rug plot
#' @param lwd Line width for activity curve lines
#' @param col Character vector for colours for the activity curve lines and rug plot; must be of equal length to the number of components
#' 
#' @return Plot of the separate components of the circular mixture model
#' 
#' @examples
#' \donttest{ componentsplot(FoxActivityGCMM) }
#'  
#' @export
	componentsplot<-function(model,rug=FALSE,ruglwd=2,lwd=3,col=c("black","grey40","grey60","grey80")) {

	# Set max for ylim: 
		max<-c()
		for(i in 1:length(model$GCMMcomponents)) {
			xx <- seq(0, 2 * pi, length = 128)				
			bw<-overlap::getBandWidth(model$GCMMcomponents[[i]], kmax=3)/5   
			   dens<-overlap::densityFit(model$GCMMcomponents[[i]], xx, bw)
			 	max<-c(max,max(dens))#  toPlot<-cbind(x=xx, y=dens)
			}
		ymax<-max(max)

	# Plot: 
		for(i in 1:length(model$GCMMcomponents)) {
			overlap::densityPlot(as.numeric(model$GCMMcomponents[[i]]),
				col=col[i],xscale=NA,xaxs="i",extend=NA,lwd=lwd,adjust=5,xlim=c(0,2*pi),
				ylim=c(0,ymax),main="",xlab="Time")
				graphics::abline(h=0,col="white")
				graphics::abline(v=c(pi,pi/2,3*pi/2),lty=2,col="gray60")
				graphics::box(lwd=2); graphics::par(new=TRUE)
		}

	# Rug: 
	if(rug==TRUE) {
			df<-data.frame(Data=model$data,Clust=as.numeric(model$clust))
			df<-df[order(df$Clust),]
			data<-model$data
			cdata<-list()
			s<-1
			for(i in 1:length(model$GCMMcomponents)) { 
				N<-as.numeric(table(df$Clust)[i] )
				N2<-s+N-1 
				c<-data[s:N2] 
				s<-s+N
					rug(c,col=col[i],lwd=ruglwd); graphics::box(lwd=2)
				}
			}
	}




#' @title Extend GCMM analysis
#' @description Extend GCMM analysis using \code{\link[runjags]{extend.jags}} from package \code{runjags}
#' @seealso \code{\link{GCMM}} \code{\link[runjags]{extend.jags}} 
#'
#' @param model Object of class \code{GCMM} that is produced by the \code{\link{GCMM}} function
#' @param sample Number of iterations per MCMC chain
#' @param burnin Number of iterations per MCMC chain to be discarded as a burn-in
#' 
#' @return Returns an object of class \code{GCMM} with a list of analysis details and output; see \code{\link{GCMM}}. A mixture plot of the estimated activity curve is also printed.
#' 
#' @examples
#' \donttest{ updateFoxGCMM<-updateGCMM(FoxActivityGCMM) }
#' 
#' @export
	updateGCMM<-function(model, burnin=0, sample=10000) {

	outGCMM<-runjags::extend.jags(model$runjags,burnin=burnin,sample=sample)

	if(model$saveclustIDs==TRUE) { clustIDs<-summary(outGCMM,vars="clust"); clustIDs<-clustIDs[,6] } else { clustIDs<-NULL } 
	if(model$saveResids==TRUE) { E<-summary(outGCMM,vars="Resid"); E<-E[,4] } else { E<-NULL } 
	if(model$saveYExp==TRUE) { YExp<-summary(outGCMM,vars="YExp") } else { YExp<-NULL } 
	if (model$saveREs==TRUE) { 
		alpha1<-summary(outGCMM,vars="alpha1")
		if (length(model$RE2)>0) { alpha2<-summary(outGCMM,vars="alpha2") } else { alpha2<-NULL }
	    } else { alpha1<-NULL; alpha2<-NULL } 
	if (model$saveJAGS==TRUE) { runjags<-outGCMM } else { runjags<-NULL }

	out<-summary(outGCMM)
	output<-cbind(out[,4],out[,1],out[,3],out[,5],out[,11],out[,9])
	colnames(output)<-c("Mean","Lower95","Upper95","SD","rhat","ESS")
	output<-round(output,3)

	message("")
	print(output)

	mu<-summary(outGCMM,vars="CircularIntercept"); mu<-round(mu[,4],2)
	if(model$Nclust>1) { p<-summary(outGCMM,vars="ClustProb"); p<-round(p[,4],2)  } else { p<-1 }
	if (model$family=="vonmises") { k<-summary(outGCMM,vars="K"); k<-round(k[,4],2) }	
	if (model$family=="wrappedcauchy") { rho<-summary(outGCMM,vars="rhoC"); rho<-round(rho[,4],2) }
	if(model$family=="vonmises"){ c<-k }
	if(model$family=="wrappedcauchy"){ c<-rho }

	gcmmcomponents<-c()
	gcmmmixture<-c()
		for (i in 1:model$Nclust) { 
			if (model$family=="vonmises") { 
				gcmmcomponents[[i]]<-circular::rvonmises(as.numeric(p[i])*10000, circular::circular(mu[i]), kappa=c[i])  
				gcmmmixture<-c(gcmmmixture,gcmmcomponents[[i]])
				 }
			if (model$family=="wrappedcauchy") { 
				gcmmcomponents[[i]]<-circular::rwrappedcauchy(p[i]*10000, circular::circular(mu[i]), c[i] )  
				gcmmmixture<-c(gcmmmixture,gcmmcomponents[[i]])
				}
		}

	overlap::densityPlot(as.numeric(gcmmmixture),xscale=NA,yaxs="i",extend=NA,lwd=3,adjust=5,xlim=c(0,2*pi),main="")
		graphics::abline(v=c(pi,pi/2,3*pi/2),lty=2,col="gray60"); graphics::box(lwd=2)

	GCMMoutput<-list(family=model$family,data=model$data,RE1=model$RE1,RE2=model$RE2,
		Nclust=model$Nclust,scale=model$scale,Mu=mu,P=p,c=c,#model=outGCMM,
		GCMMmixture=gcmmmixture,GCMMcomponents=gcmmcomponents,output=output,#summary(outGCMM),
		saveclustIDs=model$saveclustIDs,saveResids=model$saveResids,saveYExp=model$saveYExp,saveREs=model$saveREs,saveJAGS=model$saveJAGS,
		clustIDs=clustIDs,E=E,YExp=YExp,alpha1=alpha1,alpha2=alpha2,runjags=runjags)
		class(output) <- "GCMM"
		class(GCMMoutput) <- "GCMM"

	minESS<-min(output[,6],na.rm=TRUE)
	maxrhat<-max(output[,5],na.rm=TRUE)
	message(""); message(paste("All MCMC chains reached convergence (rhat<1.05):",maxrhat<1.05))
	message(paste("    Maximum rhat value:",round(maxrhat,4))); message("")
	message(paste("Minimum effective sample size:",round(minESS,0)))
	if(minESS<10000) { message("    NOTE: A minimum ESS of 10000 is recommended, consider increasing the number of samples") }

	message("");message("----------------------")
	message("GCMM updated analysis complete");message("");message("")

		return(GCMMoutput)
	}


