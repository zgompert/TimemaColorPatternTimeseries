## rate scaling models for T. cristinae stripe
## version trying to account for rate scaling
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


## read data
d<-read.table("Tcris_master_33_year.csv",header=TRUE,sep=",")

## t. cris 
k<-which(d$species=="T. cristinae")
dk<-d[k,]
sites<-unique(dk$location)
Ns<-length(sites)
Nys<-length(sites)
for(i in 1:Ns){
	dki<-dk[dk$location==sites[i],]
	ny<-length(unique(dki$year))
	Nys[i]<-ny
}

## focus on sites with 5 or more years
goodSites<-sites[Nys > 5]
Ngs<-length(goodSites)

## list of objects
yS<-vector("list",Ngs)
nS<-vector("list",Ngs)
yM<-vector("list",Ngs)
nM<-vector("list",Ngs)
yrs<-vector("list",Ngs)
pid<-vector("list",Ngs)
pS<-vector("list",Ngs)
pM<-vector("list",Ngs)
zS<-vector("list",Ngs)
zM<-vector("list",Ngs)
## loop over good sites
for(x in 1:Ngs){
	dkx<-dk[dk$location==goodSites[x],]
	yS[[x]]<-tapply(X=dkx$striped,INDEX=dkx$year,sum)
	nS[[x]]<-tapply(X=dkx$striped+dkx$unstriped,INDEX=dkx$year,sum)
	yrs[[x]]<-as.numeric(names(yS[[x]]))
	yM[[x]]<-tapply(X=dkx$melanistic,INDEX=dkx$year,sum)
	nM[[x]]<-tapply(X=dkx$striped+dkx$unstriped+dkx$melanistic,INDEX=dkx$year,sum)
	pid[[x]]<-rep(x,length(yrs[[x]]))
	pS[[x]]<-yS[[x]]/nS[[x]]
	pM[[x]]<-yM[[x]]/nM[[x]]
	zS[[x]]<-2*asin(sqrt(pS[[x]]))
	zM[[x]]<-2*asin(sqrt(pM[[x]]))
}

## now compute change and rates
dzS<-vector("list",Ngs)
dzM<-vector("list",Ngs)
dpS<-vector("list",Ngs)
dpM<-vector("list",Ngs)
rzS<-vector("list",Ngs)
rzM<-vector("list",Ngs)
rpS<-vector("list",Ngs)
rpM<-vector("list",Ngs)
dtS<-vector("list",Ngs)
dtM<-vector("list",Ngs)
mins<-40
minD<-1e-3
fixmin<-function(zz,mm){
	zz[zz<mm]<-mm
	zz
}

bndChng<-function(xx){
	yy<-xx
	for(i in 1:length(xx)){
		yy[i]<-max(xx[i],0)
	}
	return(yy)
}

for(x in 1:Ngs){
	nn<-length(yS[[x]])
	dzMat<-matrix(NA,nrow=nn,ncol=nn)
	dpMat<-matrix(NA,nrow=nn,ncol=nn)
	dtMat<-matrix(NA,nrow=nn,ncol=nn)
	for(a in 1:(nn-1)){for(b in (a+1):nn){
		if(nS[[x]][a] > mins & nS[[x]][b] > mins){
			pbar<-(pS[[x]][b]+pS[[x]][a])/2
			p1<-rbinom(n=1000,size=nS[[x]][a],prob=pbar)/nS[[x]][a]
			p2<-rbinom(n=1000,size=nS[[x]][b],prob=pbar)/nS[[x]][b]
			dpnull<-abs(p2-p1)
			dpMat[a,b]<-mean(bndChng(abs(pS[[x]][b]-pS[[x]][a])-dpnull))
			z1<-2*asin(sqrt(p1))
			z2<-2*asin(sqrt(p2))
			dznull<-abs(z2-z1)
			dzMat[a,b]<-mean(bndChng(abs(zS[[x]][b]-zS[[x]][a])-dpnull))
			dtMat[a,b]<-abs(yrs[[x]][b]-yrs[[x]][a])	
		}
	}}		
	dzS[[x]]<-fixmin(as.vector(dzMat)[is.na(as.vector(dzMat))==FALSE],minD)
	dpS[[x]]<-fixmin(as.vector(dpMat)[is.na(as.vector(dpMat))==FALSE],minD)
	dtS[[x]]<-as.vector(dtMat)[is.na(as.vector(dtMat))==FALSE]
	rpS[[x]]<-dpS[[x]]/dtS[[x]]
	rzS[[x]]<-dzS[[x]]/dtS[[x]]
	for(a in 1:(nn-1)){for(b in (a+1):nn){
		if(nM[[x]][a] > mins & nM[[x]][b] > mins){
			pbar<-(pM[[x]][b]+pM[[x]][a])/2
			p1<-rbinom(n=1000,size=nM[[x]][a],prob=pbar)/nM[[x]][a]
			p2<-rbinom(n=1000,size=nM[[x]][b],prob=pbar)/nM[[x]][b]
			dpnull<-abs(p2-p1)
			dpMat[a,b]<-mean(bndChng(abs(pM[[x]][b]-pM[[x]][a])-dpnull))
			z1<-2*asin(sqrt(p1))
			z2<-2*asin(sqrt(p2))
			dznull<-abs(z2-z1)
			dzMat[a,b]<-mean(bndChng(abs(zM[[x]][b]-zM[[x]][a])-dpnull))
			dtMat[a,b]<-abs(yrs[[x]][b]-yrs[[x]][a])	
		}
	}}		
	dzM[[x]]<-fixmin(as.vector(dzMat)[is.na(as.vector(dzMat))==FALSE],minD)
	dpM[[x]]<-fixmin(as.vector(dpMat)[is.na(as.vector(dpMat))==FALSE],minD)
	dtM[[x]]<-as.vector(dtMat)[is.na(as.vector(dtMat))==FALSE]
	rpM[[x]]<-dpM[[x]]/dtM[[x]]
	rzM[[x]]<-dzM[[x]]/dtM[[x]]
}


## Bayesian model, using log10 rate and time
## using transform, but doesn't really matter
### Stripe
dat<-list(X=log10(unlist(dtS)),Y=log10(unlist(rzS)),
	  pop=rep(1:28,unlist(lapply(rzS,length))),
	N=length(unlist(rzS)),J=Ngs)
fitzS<-stan(data=dat,"hlm2.stan",iter=8000)

betas<-extract(fitzS,"beta")
alphas<-extract(fitzS,"alpha")
mus<-extract(fitzS,"mu")

b_est<-apply(betas[[1]],2,quantile,probs=c(.5,.05,.95))
a_est<-apply(alphas[[1]],2,quantile,probs=c(.5,.05,.95))

xx<-seq(0,1.6,.01)
m_post<-matrix(nrow=16000,ncol=length(xx))
for(k in 1:16000){
        m_post[k,]<-mus[[1]][k,1] + xx * mus[[1]][k,2]
}

m_est<-apply(m_post,2,quantile,probs=c(.5,.05,.95))

a<-1:length(dat$X)
aord<-sample(a,length(a),replace=FALSE)

cs<-c(brewer.pal(n=9,"Set1"),brewer.pal(n=8,"Set2"),brewer.pal(n=11,"Set3"))
pid<-rep(1:28,unlist(lapply(rzS,length)))
pdf("ratesStripe-err.pdf",width=5,height=5)
par(mar=c(4.5,4.5,.5,.5))

ca<-1;cl<-1.3
plot(dat$X,dat$Y,type='n',xlab="Log10 interval (generations)",ylab="Log10 rate",cex.lab=cl,cex.axis=ca)
polygon(c(xx,rev(xx)),c(m_est[2,],rev(m_est[3,])),col=alpha("gray50",.5),border=NA)
points(dat$X[aord],dat$Y[aord],pch=19,col=alpha(cs,.5)[pid[aord]])
for(i in 1:Ngs){
abline(a=a_est[1,i],b=b_est[1,i],col=cs[i])}
lines(xx,m_est[1,],lwd=2)
apply(mus[[1]],2,quantile,probs=c(.5,.05,.95)) 
mtext(expression(paste(beta," = -0.77 (-0.88, -0.66)",sep="")),side=1,line=-2,cex=1.1)

dev.off()

### Melanic
dat<-list(X=log10(unlist(dtM)),Y=log10(unlist(rzM)),
	  pop=rep(1:28,unlist(lapply(rzM,length))),
	N=length(unlist(rzM)),J=Ngs)
fitzM<-stan(data=dat,"hlm2.stan",iter=8000)

betas<-extract(fitzM,"beta")
alphas<-extract(fitzM,"alpha")
mus<-extract(fitzM,"mu")

b_est<-apply(betas[[1]],2,quantile,probs=c(.5,.05,.95))
a_est<-apply(alphas[[1]],2,quantile,probs=c(.5,.05,.95))

xx<-seq(0,1.6,.01)
m_post<-matrix(nrow=16000,ncol=length(xx))
for(k in 1:16000){
        m_post[k,]<-mus[[1]][k,1] + xx * mus[[1]][k,2]
}

m_est<-apply(m_post,2,quantile,probs=c(.5,.05,.95))

a<-1:length(dat$X)
aord<-sample(a,length(a),replace=FALSE)

cs<-c(brewer.pal(n=9,"Set1"),brewer.pal(n=8,"Set2"),brewer.pal(n=11,"Set3"))
pid<-rep(1:28,unlist(lapply(rzS,length)))
pdf("ratesMelanic-err.pdf",width=5,height=5)
par(mar=c(4.5,4.5,.5,.5))

ca<-1;cl<-1.3
plot(dat$X,dat$Y,type='n',xlab="Log10 interval (generations)",ylab="Log10 rate",cex.lab=cl,cex.axis=ca)
polygon(c(xx,rev(xx)),c(m_est[2,],rev(m_est[3,])),col=alpha("gray50",.5),border=NA)
points(dat$X[aord],dat$Y[aord],pch=19,col=alpha(cs,.5)[pid[aord]])
for(i in 1:Ngs){
abline(a=a_est[1,i],b=b_est[1,i],col=cs[i])}
lines(xx,m_est[1,],lwd=2)
apply(mus[[1]],2,quantile,probs=c(.5,.05,.95)) ## number not righ
mtext(expression(paste(beta," = -0.83 (-0.94, -0.72)",sep="")),side=1,line=-2,cex=1.1)

dev.off()

## compare melanistic and stripe
betas<-extract(fitzS,"beta")
alphas<-extract(fitzS,"alpha")
mus<-extract(fitzS,"mu")

bS_est<-apply(betas[[1]],2,quantile,probs=c(.5,.05,.95))
aS_est<-apply(alphas[[1]],2,quantile,probs=c(.5,.05,.95))

betas<-extract(fitzM,"beta")
alphas<-extract(fitzM,"alpha")
mus<-extract(fitzM,"mu")

bM_est<-apply(betas[[1]],2,quantile,probs=c(.5,.05,.95))
aM_est<-apply(alphas[[1]],2,quantile,probs=c(.5,.05,.95))

oo<-data.frame(goodSites,bS_est[1,],bM_est[1,],10^aS_est[1,]/2,10^aM_est[1,]/2)

write.table(oo,file="RatesTcris-err.csv",row.names=FALSE,col.names=TRUE,quote=FALSE)

cor.test(bM_est[1,],bS_est[1,])
#data:  bM_est[1, ] and bS_est[1, ]
#t = 2.061, df = 26, p-value = 0.04944
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.001927136 0.656087374
#sample estimates:
#      cor 
#0.3747346 

cor.test(aM_est[1,],aS_est[1,])
#data:  aM_est[1, ] and aS_est[1, ]
#t = 1.1481, df = 26, p-value = 0.2614
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.1671008  0.5478492
#sample estimates:
#      cor 
#0.2196707 

cor.test(aM_est[1,],bM_est[1,])
#data:  aM_est[1, ] and bM_est[1, ]
#t = -0.97199, df = 26, p-value = 0.34
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.5237397  0.1997828
#sample estimates:
#       cor 
#-0.1872507

cor.test(aS_est[1,],bS_est[1,])
#t = 1.7496, df = 26, p-value = 0.09198
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.05520264  0.62228579
#sample estimates:
#      cor 
#0.3245585 


## direct comparison of change using same cutoff for stripe and M
## using signed change, unlike above!!
for(x in 1:Ngs){
	nn<-length(yS[[x]])
	dzMat<-matrix(NA,nrow=nn,ncol=nn)
	dpMat<-matrix(NA,nrow=nn,ncol=nn)
	dtMat<-matrix(NA,nrow=nn,ncol=nn)
	for(a in 1:(nn-1)){for(b in (a+1):nn){
		if(nS[[x]][a] > mins & nS[[x]][b] > mins){
			pbar<-(pS[[x]][b]+pS[[x]][a])/2
			p1<-rbinom(n=1000,size=nS[[x]][a],prob=pbar)/nS[[x]][a]
			p2<-rbinom(n=1000,size=nS[[x]][b],prob=pbar)/nS[[x]][b]
			dpnull<-(p2-p1)
			dpMat[a,b]<-mean(bndChng((pS[[x]][b]-pS[[x]][a])-dpnull))
			z1<-2*asin(sqrt(p1))
			z2<-2*asin(sqrt(p2))
			dznull<-(z2-z1)
			dzMat[a,b]<-mean(bndChng((zS[[x]][b]-zS[[x]][a])-dpnull))
			dtMat[a,b]<-abs(yrs[[x]][b]-yrs[[x]][a])	
		}
	}}		
	dzS[[x]]<-fixmin(as.vector(dzMat)[is.na(as.vector(dzMat))==FALSE],minD)
	dpS[[x]]<-fixmin(as.vector(dpMat)[is.na(as.vector(dpMat))==FALSE],minD)
	dtS[[x]]<-as.vector(dtMat)[is.na(as.vector(dtMat))==FALSE]
	rpS[[x]]<-dpS[[x]]/dtS[[x]]
	rzS[[x]]<-dzS[[x]]/dtS[[x]]
	for(a in 1:(nn-1)){for(b in (a+1):nn){
		if(nS[[x]][a] > mins & nS[[x]][b] > mins){
			pbar<-(pM[[x]][b]+pM[[x]][a])/2
			p1<-rbinom(n=1000,size=nM[[x]][a],prob=pbar)/nM[[x]][a]
			p2<-rbinom(n=1000,size=nM[[x]][b],prob=pbar)/nM[[x]][b]
			dpnull<-(p2-p1)
			dpMat[a,b]<-mean(bndChng((pM[[x]][b]-pM[[x]][a])-dpnull))
			z1<-2*asin(sqrt(p1))
			z2<-2*asin(sqrt(p2))
			dznull<-(z2-z1)
			dzMat[a,b]<-mean(bndChng((zM[[x]][b]-zM[[x]][a])-dpnull))
			dtMat[a,b]<-abs(yrs[[x]][b]-yrs[[x]][a])	
		}
	}}		
	dzM[[x]]<-fixmin(as.vector(dzMat)[is.na(as.vector(dzMat))==FALSE],minD)
	dpM[[x]]<-fixmin(as.vector(dpMat)[is.na(as.vector(dpMat))==FALSE],minD)
	dtM[[x]]<-as.vector(dtMat)[is.na(as.vector(dtMat))==FALSE]
	rpM[[x]]<-dpM[[x]]/dtM[[x]]
	rzM[[x]]<-dzM[[x]]/dtM[[x]]
}

cor.test(unlist(dzS),unlist(dzM))

#t = -1.3253, df = 978, p-value = 0.1854
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.10468607  0.02033577
#sample estimates:
#       cor 
#-0.0423409 

