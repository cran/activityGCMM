#' @title Generalized circular mixed effect mixture (GCMM) model 
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
#' Version: 1.0.1
#' Date: 2021-06-06 
#' Author: Liz AD Campbell
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
#' @param thin Thinning rate for MCMC chains, i.e. how many samples are saved. For longer models, thin can be increased to reduce computer memory requirements
#' @param sample If autojags=FALSE, the number of MCMC samples per chain (which is multiplied by thin); default=10000
#' @param burnin If autojags=FALSE, the burnin for the MCMC chains which are not saved; default=5000
#' @param adapt adaptation to use for MCMC chains; default=1000
#' @param n.chains number of MCMC chains; default=3
#' @param autorun Logical argument for whether to autmatically extend the analyses to achieve MCMC chain convergence and a specified minimum effective sample size (ESS) for all parameters; default=TRUE
#' @param minESS Minimum effective sample size (ESS) from the posterior distribution desired for all paramerers; default=5000, though a minimum ESS of 10000 is recommended
#' @param maxrep Maximum number of extensions of the analysis if \code{autorun=TRUE}; default=5
#' @param Nclust Number of components for mixture models; if not provided, the function will estimate the number of clusters; if provided, values must be provided for clustmeans
#' @param clustmeans A vector equal in length to Nclust of the potential means for each component in the mixture models
#' @param saveREs Whether random intercepts are saved in output; recommended to save only one of saveREs, saveResids or saveYExp at one time due to memory limitations
#' @param saveclustIDs Whether to save component cluster identification for the data points; default=FALSE
#' @param saveResids Whether model residuals are saved in output; recommended to save only one of saveREs, saveResids or saveYExp at one time due to memory limitations
#' @param saveYExp Whether expected Y values based on model are saved in output; recommended to save only one of saveREs, saveResids or saveYExp at one time due to memory limitations
#' @param saveJAGS Logical argument of whether to save runjags output; default=FALSE
#' 
#' @returns Returns object of class \code{GCMM} which is a list containing analysis results and details. A plot of the estimated activity curve from the mixed effect mixture model is printed.
#' @returns \code{output} GCMM model output summary
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
#'     scale=c("2pi"), family="vonmises", autorun=FALSE,
#'     adapt=0, sample=300, burnin=300, thin=1,n.chains=2  ) }
#'  
#' @export
	GCMM<-function(data, RE1, RE2=NULL, scale="2pi", kmax=15, family=c("vonmises","wrappedcauchy"),
		autorun=TRUE, minESS=5000, maxrep=5, thin=2, burnin=5000, sample=5000, adapt=1000, n.chains=3,
		saveREs=FALSE, saveResids=FALSE, saveclustIDs=FALSE, saveYExp=FALSE, saveJAGS=TRUE, Nclust="NULL", clustmeans=NULL )  {

message(""); message("activityGCMM: Bayesian generalized circular mixed effect mixture models for estimating animal temporal activity curves")
message("Bugs or comments can be sent to lizadcampbell@atlasgoldenwolf.org or posted on www.atlasgoldenwolf.org/ladcampbellstatsforum"); message("")

requireNamespace("mclust",quietly=TRUE)

	testjags<-runjags::testjags(silent=TRUE)
		if (testjags$JAGS.available==FALSE) { 
					stop("JAGS software is not available on your computer. To download JAGS, visit https://mcmc-jags.sourceforge.io/ ") 
					}

if (Nclust=="NULL") { 
	message("Identifying clusters....")

	# Model-based k-clustering with mclust (0-2pi):
		data<-convertRad(data,to="2pi")
		mclustfit <- mclust::Mclust(data)
			S2pi<-summary(mclustfit,parameters=TRUE)
			clustmeans2pi<-S2pi$mean
			Nclusts2pi<-length(clustmeans2pi)

	# check for circular data; converting to [-pi,pi] to check number of clusters:
		datpirad<-convertRad(data,to="pi")
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


	} else { clustmeans<-clustmeans; mclustout<-NULL; mclustfit<-NULL; mclustfitpi<-NULL }

	if(scale=="pi") { data<-convertRad(data,to="pi") }
	if(scale=="2pi") { data<-convertRad(data,to="2pi") }


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
	initsVM1RE <- function () { list( CircularIntercept=stats::rnorm(Nclust,clustmeans,.0001), 
		ClustProb=stats::runif(Nclust,0,1), sigma1=stats::runif(Nclust,0,.5)	 ) }
	initsVM2REs <- function () { list( CircularIntercept=stats::rnorm(Nclust,clustmeans,.0001), 
		ClustProb=stats::runif(Nclust,0,1), sigma1=stats::runif(Nclust,0,.5), sigma2=stats::runif(Nclust,0,.5)	 ) }
	initsWC1RE <- function () { list( CircularIntercept=stats::rnorm(Nclust,clustmeans,.0001), rhoC=stats::runif(Nclust,0.1,0.9))}#,
	##	alpha1=stats::rnorm(length(unique(RE1)),0,0.001) ) } 
	initsWC2REs <- function () { list( CircularIntercept=stats::rnorm(Nclust,clustmeans,.0001), rhoC=stats::runif(Nclust,0.1,0.9))}#,
	##	alpha1=stats::rnorm(length(unique(RE1)),0,0.001),alpha2=stats::rnorm(length(unique(RE2)),0,0.001) ) } 				


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
	
	outGCMM<-runjags::run.jags(model=model, method="parallel", monitor=monitor, 
					n.chains=n.chains, thin=thin, adapt=adapt, data=JAGSdata, inits=inits, 
					burnin=burnin, sample=sample) 

	out<-summary(outGCMM)
	output<-cbind(out[,4],out[,1],out[,3],out[,5],out[,11],out[,9])
	colnames(output)<-c("Mean","Lower95","Upper95","SD","rhat","ESS")
	output<-round(output,3)


	if(autorun==TRUE) {
		maxrhat<-max(output[,5],na.rm=TRUE)
		i<-1
		if(maxrhat>1.05) { repeat {
				i<-i+1 
				if(maxrhat<1.05) { break } #else {
  	      		if ( i==maxrep ) { message("Maximum attempts reached. Analysis stopped before all rhat < 1.05; max rhat =",round(maxrhat,3)," Try increasing the number of samples."); break }   
				message(""); message(""); message(paste("MCMC chains have not reached convergence (rhat<1.05); max rhat =",round(maxrhat,3)))
				message("    Increasing MCMC iterations....");message(""); message("")
				outGCMM<-runjags::extend.jags(outGCMM, burnin=burnin, sample=sample, combine=FALSE)
					out<-summary(outGCMM)
					output<-cbind(out[,4],out[,1],out[,3],out[,5],out[,11],out[,9])
					colnames(output)<-c("Mean","Lower95","Upper95","SD","rhat","ESS")
					output<-round(output,3)
					maxrhat<-max(output[,5],na.rm=TRUE)
			 }	 }
			
		ESSmin<-min(output[,6],na.rm=TRUE)
		i2<-1
			if(maxrhat<1.05) { if(ESSmin<minESS) { repeat {
				i2<-i2+1
				if(ESSmin>minESS) { break } 
  	    			if (i2==maxrep) { message(paste("Maximum attempts reached. Analysis stopped before minimum ESS achieved. minESS =",round(ESSmin,0))) 
					 message("   Try increasing the number of samples."); break }   
				message(""); message(""); message(paste("Minimum ESS not yet achieved; minESS =",round(ESSmin,0) ))
				message("    Increasing MCMC iterations...."); message(""); message("")
				outGCMM<-runjags::extend.jags(outGCMM, burnin=0, sample=sample)
					out<-summary(outGCMM)
					output<-cbind(out[,4],out[,1],out[,3],out[,5],out[,11],out[,9])
					colnames(output)<-c("Mean","Lower95","Upper95","SD","rhat","ESS")
					output<-round(output,3)
						ESSmin<-min(output[,6],na.rm=TRUE)
				 } } }
		}




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
				gcmmcomponents[[i]]<-circular::rvonmises(as.numeric(p[i])*10000, circular::circular(convertRad(mu[i],to="2pi")), k[i])  
				gcmmmixture<-c(gcmmmixture,gcmmcomponents[[i]])
				 }
			if (family=="wrappedcauchy") { 
				gcmmcomponents[[i]]<-circular::rwrappedcauchy(p[i]*10000, circular::circular(convertRad(mu[i],to="2pi")), rho[i] )  
				gcmmmixture<-c(gcmmmixture,gcmmcomponents[[i]])
				}
		}


	### Mixture Plot: 
	if (scale=="pi") { lim1<-(-1*pi); lim2<-pi; xcenter<-c("midnight") }
	if (scale=="2pi") { lim1<-0; lim2<-2*pi; xcenter<-c("noon") }
		overlap::densityPlot(as.numeric(gcmmmixture),xcenter=xcenter,xscale=NA,yaxs="i",xaxs="i",extend=NA,lwd=3,adjust=5,xlim=c(lim1,lim2),main="")
			graphics::abline(v=c(-pi,0,pi,pi/2,3*pi/2,-pi/2),lty=2,col="gray60"); graphics::box(lwd=2)
 

	if(family=="vonmises"){ c<-k }
	if(family=="wrappedcauchy"){ c<-rho }
		
	GCMMoutput<-list(family=family,data=data,RE1=RE1,RE2=RE2,Nclust=Nclust,clustmeans=clustmeans,scale=scale,Mu=mu,P=p,c=c,model=outGCMM,
		GCMMmixture=gcmmmixture,GCMMcomponents=gcmmcomponents,output=output,#summary(outGCMM),
		saveclustIDs=saveclustIDs,saveResids=saveResids,saveYExp=saveYExp,saveREs=saveREs,saveJAGS=saveJAGS,
			mclustfit=mclustfit,mclustfitpi=mclustfitpi, 
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


#######################################################################

#' @title Executable example of GCMM function
#' @description Example of applying generalized circular mixed effect mixture model with activityGCMM 
#'     using data included in the package
#' @return Prints message with example of GCMM function using data included in the package
#' @examples { exampleGCMM() }
#' @export
	exampleGCMM<- function() { message("Try example with sample data included in the package: data(redfoxsample)")
			message("Run a generalized circular mixed effect mixture vonmises model using:")
			message("    FoxGCMM<-GCMM(data=redfoxsample$Radians, RE1=redfoxsample$CameraTrapID, family=\"vonmises\",") 
     			message("           autorun=TRUE, burnin=2000, sample=5000, thin=2, minESS=2000)")	}


#####################################################

#' @title Extend GCMM analysis
#' @description Extend GCMM analysis using \code{\link[runjags]{extend.jags}} from package \code{runjags}
#' @seealso \code{\link{GCMM}} \code{\link[runjags]{extend.jags}} 
#'
#' @param model Object of class \code{GCMM} that is produced by the \code{\link{GCMM}} function
#' @param sample Number of iterations per MCMC chain
#' @param burnin Number of iterations per MCMC chain to be discarded as a burn-in
#' @param autorun Whether to automatically extend the analysis until MCMC chain convergence and minimum effective sample size (ESS) is achieved; default is TRUE
#' @param minESS Desired minimum effective sample size (MCMC) when automatically extending the analysis using \code{autorun=TRUE}; default is 5000
#' @param maxrep Maximum number of times to automatically extend the analysis if MCMC chains have not converged or the minimum effective sample size is not reached; default=5
#' @param saveREs Whether random intercepts are saved in output; recommended to save only one of saveREs, saveResids or saveYExp at one time due to memory limitations; default=FALSE
#' @param saveclustIDs Whether to save component cluster identification for the data points; default=FALSE
#' @param saveResids Whether model residuals are saved in output; recommended to save only one of saveREs, saveResids or saveYExp at one time due to memory limitations
#' @param drop.chain A number indicating which MCMC chain to drop from the updated analysis. This may be useful if one chain happens to converge on opposite clusters than the others. 
#' 
#' @return Returns an object of class \code{GCMM} with a list of analysis details and output; see \code{\link{GCMM}}. A mixture plot of the estimated activity curve is also printed.
#' 
#' @examples
#' \donttest{ FoxActivityGCMM<-GCMM(data=redfoxsample$Radians, 
#'               RE1=redfoxsample$SamplingPeriod, family="vonmises", autorun=FALSE,
#'               adapt=0, sample=300, burnin=300, thin=1, n.chains=2)
#'            updateFoxGCMM<-updateGCMM(FoxActivityGCMM, sample=300, autorun=FALSE) }
#' 
#' @export
	updateGCMM<-function(model, burnin=0, sample=10000, saveclustIDs=FALSE, saveREs=FALSE, saveResids=FALSE, 
		autorun=TRUE, minESS=5000, maxrep=5, drop.chain=0) {

	add<-FALSE; addmonitor<-c()
	if (model$saveclustIDs==FALSE) { if(saveclustIDs==TRUE) { add<-TRUE; addmonitor<-c(addmonitor,"clust") } }
	if (model$saveREs==FALSE) { if(saveREs==TRUE) { add<-TRUE; addmonitor<-c(addmonitor,"alpha1")
		if (length(model$RE2>0)) { add.monitor<-c(addmonitor, "alpha2") } } }
	if (model$saveResids==FALSE) { if(saveResids==TRUE) { add<-TRUE; addmonitor<-c(addmonitor,"Resid") } }
	if(add==TRUE) { 	
		outGCMM<-runjags::extend.jags(model$runjags,burnin=burnin,sample=sample,add.monitor=addmonitor) } else {
	if(drop.chain!=0) { 
		outGCMM<-runjags::extend.jags(model$runjags,burnin=burnin,sample=sample,drop.chain=drop.chain) } else {
		outGCMM<-runjags::extend.jags(model$runjags,burnin=burnin,sample=sample, combine=TRUE)
	  }	}

	if(saveclustIDs==TRUE) { clustIDs<-summary(outGCMM,vars="clust"); clustIDs<-clustIDs[,6] } else { clustIDs<-NULL } 
	if(model$saveResids==TRUE) { E<-summary(outGCMM,vars="Resid"); E<-E[,4] } else { E<-NULL } 
	if(model$saveYExp==TRUE) { YExp<-summary(outGCMM,vars="YExp") } else { YExp<-NULL } 
	if (saveREs==TRUE) { 
		alpha1<-summary(outGCMM,vars="alpha1")
			if (length(model$RE2)>0) { alpha2<-summary(outGCMM,vars="alpha2") } else { alpha2<-NULL }
			} else { 
				if (model$saveREs==TRUE) {
					alpha1<-summary(outGCMM,vars="alpha1")
						if (length(model$RE2)>0) { alpha2<-summary(outGCMM,vars="alpha2") } else { alpha2<-NULL }
	   						 } else { 
								alpha1<-NULL; alpha2<-NULL } }
	if (model$saveJAGS==TRUE) { runjags<-outGCMM } else { runjags<-NULL }

	out<-summary(outGCMM)
	output<-cbind(out[,4],out[,1],out[,3],out[,5],out[,11],out[,9])
	colnames(output)<-c("Mean","Lower95","Upper95","SD","rhat","ESS")
	output<-round(output,3)

	if(autorun==TRUE) {
		maxrhat<-max(output[,5],na.rm=TRUE)
		i<-1
		if(maxrhat>1.05) { repeat {
				i<-i+1 
				if(maxrhat<1.05) { break } 
  	      		if ( i==maxrep ) { message("Maximum attempts reached. Analysis stopped before all rhat < 1.05; max rhat =",round(maxrhat,3)," Try increasing the number of samples."); break }   
				message(""); message(""); message(paste("MCMC chains have not reached convergence (rhat<1.05); max rhat =",round(maxrhat,3)))
				message("    Increasing MCMC iterations....");message(""); message("")
				outGCMM<-runjags::extend.jags(outGCMM, burnin=burnin, sample=sample, combine=FALSE)
					out<-summary(outGCMM)
					output<-cbind(out[,4],out[,1],out[,3],out[,5],out[,11],out[,9])
					colnames(output)<-c("Mean","Lower95","Upper95","SD","rhat","ESS")
					output<-round(output,3)
					maxrhat<-max(output[,5],na.rm=TRUE)
				 }	 }	
		ESSmin<-min(output[,6],na.rm=TRUE)
		i2<-1
			if(maxrhat<1.05) { if(ESSmin<minESS) { repeat {
				i2<-i2+1
				if(ESSmin>minESS) { break } 
  	    			if (i2==maxrep) { message(paste("Maximum attempts reached. Analysis stopped before minimum ESS achieved. minESS =",round(ESSmin,0))) 
					 message("   Try increasing the number of samples."); break }   
				message(""); message(""); message(paste("Minimum ESS not yet achieved; minESS =",round(ESSmin,0) ))
				message("    Increasing MCMC iterations...."); message(""); message("")
				outGCMM<-runjags::extend.jags(outGCMM, burnin=0, sample=sample, combine=TRUE)
					out<-summary(outGCMM)
					output<-cbind(out[,4],out[,1],out[,3],out[,5],out[,11],out[,9])
					colnames(output)<-c("Mean","Lower95","Upper95","SD","rhat","ESS")
					output<-round(output,3)
						ESSmin<-min(output[,6],na.rm=TRUE)
				 } } }
		}

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
				gcmmcomponents[[i]]<-circular::rvonmises(as.numeric(p[i])*10000, circular::circular(convertRad(mu[i],to="2pi")), kappa=c[i])  
				gcmmmixture<-c(gcmmmixture,gcmmcomponents[[i]])
				 }
			if (model$family=="wrappedcauchy") { 
				gcmmcomponents[[i]]<-circular::rwrappedcauchy(p[i]*10000, circular::circular(convertRad(mu[i],to="2pi")), c[i] )  
				gcmmmixture<-c(gcmmmixture,gcmmcomponents[[i]])
				}
		}

	overlap::densityPlot(as.numeric(gcmmmixture),xscale=NA,yaxs="i",extend=NA,lwd=3,adjust=5,xlim=c(0,2*pi),main="")
		graphics::abline(v=c(pi,pi/2,3*pi/2),lty=2,col="gray60"); graphics::box(lwd=2)

	GCMMoutput<-list(family=model$family,data=model$data,RE1=model$RE1,RE2=model$RE2,
		Nclust=model$Nclust,scale=model$scale,Mu=mu,P=p,c=c,
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

	message(""); message("----------------------")
	message("GCMM updated analysis complete"); message(""); message("")

		return(GCMMoutput)
	}


