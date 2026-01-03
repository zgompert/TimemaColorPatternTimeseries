library(data.table)

## relative to:
#240016H1 = VP green from C

## colinearity plots for all homologous chromsomes


## colinearity plots for all homologous chromsomes
## chrom number, gs, gsr1, 24_0016 H1
chtab<-matrix(c(1,8483,12,9,
	2,14640,6,6,
	3,42935,2,3,
	4,42912,1,2,
	5,18722,7,5,
	6,9928,8,8,
	7,10660,10,13,
	8,7748,11,1,
	9,16151,5,10,
	10,14160,4,12,
	11,12033,9,11,
	12,12380,13,7,
	13,14101,3,4),nrow=13,ncol=4,byrow=TRUE)

####################### 24_0028 ########################################
## read synteny dat
#########################################################################
## 0016h1 x 0028h1
##############
dat<-fread("out_synteny_TcrE240016H1_TcrE240028H1.psl",header=FALSE)
dfdat<-as.data.frame(dat)

## target = 0016h1
## query = 0028h1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
g1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
g2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% g1Ch) & (dfdat[,10] %in% g2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_g1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_g2<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

pdf("SynTcrisStripe_240016h1_240028h1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (24_0028 H1)",ylab="T. cristinae (24_0016 H1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_g1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_g2,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs, gsr1, 24_0016 H1, 24_0028 H1
chtab<-matrix(c(1,8483,12,9,8,
	2,14640,6,6,6,
	3,42935,2,3,2,
	4,42912,1,2,1,
	5,18722,7,5,13,
	6,9928,8,8,6,
	7,10660,10,13,12,
	8,7748,11,1,5,
	9,16151,5,10,10,
	10,14160,4,12,4,
	11,12033,9,11,11,
	12,12380,13,7,9,
	13,14101,3,4,3),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_240016H1_240028H1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        if(chtab[i,4] > 9){
                val<-paste(chtab[i,4],sep="")
        }else{
                val<-paste("0",chtab[i,4],sep="")
        }
        tcr_g1<-grep(x=subDfdat[,14],paste("hap1_Chr",val,sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g2<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
        cc<-tcr_g1[tcr_g1 %in% tcr_g2]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="24_0028 H1",ylab="24_0016 H1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

##############
## 0028h1 x 0028h2
##############
dat<-fread("out_synteny_TcrE240028H1_TcrE240028H2.psl",header=FALSE)
dfdat<-as.data.frame(dat)

## target = 0028h1
## query = 0028h2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
g1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
g2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% g1Ch) & (dfdat[,10] %in% g2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_g1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_g2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

pdf("SynTcrisStripe_240028h1_240028h2.pdf",width=6,height=6)## 1:1 chroms
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (24_0028 H2)",ylab="T. cristinae (24_0028 H1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_g1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_g2,las=2)
box()
dev.off()


pdf("AlnPlotTcris_240028H1_240028H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g1<-grep(x=subDfdat[,14],paste("hap1_Chr",val,sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
        cc<-tcr_g1[tcr_g1 %in% tcr_g2]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="24_0028 H2",ylab="24_0028 H1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

####################### 24_0029 ########################################
## read synteny dat
#########################################################################
## 0016h1 x 0029h1
##############
dat<-fread("out_synteny_TcrE240016H1_TcrE240029H1.psl",header=FALSE)
dfdat<-as.data.frame(dat)

## target = 0016h1
## query = 0029h1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
g1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
g2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% g1Ch) & (dfdat[,10] %in% g2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_g1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_g2<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

pdf("SynTcrisStripe_240016h1_240029h1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (24_0029 H1)",ylab="T. cristinae (24_0016 H1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_g1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_g2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs, gsr1, 24_0016 H1, 24_0029 H1
chtab<-matrix(c(1,8483,12,9,8,
	2,14640,6,6,7,
	3,42935,2,3,2,
	4,42912,1,2,1,
	5,18722,7,5,5,
	6,9928,8,8,6,
	7,10660,10,13,13,
	8,7748,11,1,4,
	9,16151,5,10,12,
	10,14160,4,12,10,
	11,12033,9,11,11,
	12,12380,13,7,9,
	13,14101,3,4,3),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_240016H1_240029H1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        if(chtab[i,4] > 9){
                val<-paste(chtab[i,4],sep="")
        }else{
                val<-paste("0",chtab[i,4],sep="")
        }
        tcr_g1<-grep(x=subDfdat[,14],paste("hap1_Chr",val,sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g2<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
        cc<-tcr_g1[tcr_g1 %in% tcr_g2]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="24_0029 H1",ylab="24_0016 H1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

##############
## 0029h1 x 0029h2
##############
dat<-fread("out_synteny_TcrE240029H1_TcrE240029H2.psl",header=FALSE)
dfdat<-as.data.frame(dat)

## target = 0029h1
## query = 0029h2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
g1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
g2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% g1Ch) & (dfdat[,10] %in% g2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_g1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_g2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

pdf("SynTcrisStripe_240029h1_240029h2.pdf",width=6,height=6)## 1:1 chroms
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (24_0029 H2)",ylab="T. cristinae (24_0029 H1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_g1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_g2,las=2)
box()
dev.off()

pdf("AlnPlotTcris_240029H1_240029H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g1<-grep(x=subDfdat[,14],paste("hap1_Chr",val,sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
        cc<-tcr_g1[tcr_g1 %in% tcr_g2]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="24_0029 H2",ylab="24_0029 H1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

####################### 24_0030 ########################################
## read synteny dat
#########################################################################
## 0016h1 x 0030h1
##############
dat<-fread("out_synteny_TcrE240016H1_TcrE240030H1.psl",header=FALSE)
dfdat<-as.data.frame(dat)

## target = 0016h1
## query = 0030h1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
g1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
g2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% g1Ch) & (dfdat[,10] %in% g2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_g1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_g2<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

pdf("SynTcrisStripe_240016h1_240030h1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (24_0030 H1)",ylab="T. cristinae (24_0016 H1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_g1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_g2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs, gsr1, 24_0016 H1, 24_0030 H1
chtab<-matrix(c(1,8483,12,9,7,
	2,14640,6,6,6,
	3,42935,2,3,1,
	4,42912,1,2,2,
	5,18722,7,5,13,
	6,9928,8,8,4,
	7,10660,10,13,12,
	8,7748,11,1,5,
	9,16151,5,10,10,
	10,14160,4,12,8,
	11,12033,9,11,11,
	12,12380,13,7,9,
	13,14101,3,4,3),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_240016H1_240030H1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        if(chtab[i,4] > 9){
                val<-paste(chtab[i,4],sep="")
        }else{
                val<-paste("0",chtab[i,4],sep="")
        }
        tcr_g1<-grep(x=subDfdat[,14],paste("hap1_Chr",val,sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g2<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
        cc<-tcr_g1[tcr_g1 %in% tcr_g2]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="24_0030 H1",ylab="24_0016 H1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

##############
## 0030h1 x 0030h2
##############
dat<-fread("out_synteny_TcrE240030H1_TcrE240030H2.psl",header=FALSE)
dfdat<-as.data.frame(dat)

## target = 0030h1
## query = 0030h2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
g1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
g2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% g1Ch) & (dfdat[,10] %in% g2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_g1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_g2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

pdf("SynTcrisStripe_240030h1_240030h2.pdf",width=6,height=6)## 1:1 chroms
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (24_0030 H2)",ylab="T. cristinae (24_0030 H1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_g1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_g2,las=2)
box()
dev.off()

pdf("AlnPlotTcris_240030H1_240030H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g1<-grep(x=subDfdat[,14],paste("hap1_Chr",val,sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
        cc<-tcr_g1[tcr_g1 %in% tcr_g2]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="24_0030 H2",ylab="24_0030 H1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

####################### 24_0089 ########################################
## read synteny dat
#########################################################################
## 0016h1 x 0089h1
##############
dat<-fread("out_synteny_TcrE240016H1_TcrE240089H1.psl",header=FALSE)
dfdat<-as.data.frame(dat)

## target = 0016h1
## query = 0089h1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
g1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
g2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% g1Ch) & (dfdat[,10] %in% g2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_g1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_g2<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

pdf("SynTcrisStripe_240016h1_240089h1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (24_0089 H1)",ylab="T. cristinae (24_0016 H1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_g1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_g2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs, gsr1, 24_0016 H1, 24_0089 H1
chtab<-matrix(c(1,8483,12,9,5,
	2,14640,6,6,13,
	3,42935,2,3,12,
	4,42912,1,2,8,
	5,18722,7,5,10,
	6,9928,8,8,7,
	7,10660,10,13,11,
	8,7748,11,1,3,
	9,16151,5,10,6,
	10,14160,4,12,2,
	11,12033,9,11,4,
	12,12380,13,7,9,
	13,14101,3,4,1),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_240016H1_240089H1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        if(chtab[i,4] > 9){
                val<-paste(chtab[i,4],sep="")
        }else{
                val<-paste("0",chtab[i,4],sep="")
        }
        tcr_g1<-grep(x=subDfdat[,14],paste("hap1_Chr",val,sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g2<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
        cc<-tcr_g1[tcr_g1 %in% tcr_g2]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="24_0089 H1",ylab="24_0016 H1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

##############
## 0089h1 x 0089h2
##############
dat<-fread("out_synteny_TcrE240089H1_TcrE240089H2.psl",header=FALSE)
dfdat<-as.data.frame(dat)

## target = 0089h1
## query = 0089h2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
g1Ch<-names(xx)[xx>200]
xx<-table(dfdat[,10])
g2Ch<-names(xx)[xx>200]
keep<-(dfdat[,14] %in% g1Ch) & (dfdat[,10] %in% g2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_g1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_g2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

pdf("SynTcrisStripe_240089h1_240089h2.pdf",width=6,height=6)## 1:1 chroms
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (24_0089 H2)",ylab="T. cristinae (24_0089 H1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_g1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_g2,las=2)
box()
dev.off()

pdf("AlnPlotTcris_240089H1_240089H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g1<-grep(x=subDfdat[,14],paste("hap1_Chr",val,sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
        cc<-tcr_g1[tcr_g1 %in% tcr_g2]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="24_0089 H2",ylab="24_0089 H1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

####################### 24_0179 ########################################
## read synteny dat
#########################################################################
## 0016h1 x 0179h1
##############
dat<-fread("out_synteny_TcrE240016H1_TcrE240179H1.psl",header=FALSE)
dfdat<-as.data.frame(dat)

## target = 0016h1
## query = 0179h1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
g1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
g2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% g1Ch) & (dfdat[,10] %in% g2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_g1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_g2<-as.numeric(gsub(pattern="hap1_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

pdf("SynTcrisStripe_240016h1_240179h1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (24_0179 H1)",ylab="T. cristinae (24_0016 H1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_g1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_g2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs, gsr1, 24_0016 H1, 24_0179 H1
chtab<-matrix(c(1,8483,12,9,10,
	2,14640,6,6,7,
	3,42935,2,3,2,
	4,42912,1,2,1,
	5,18722,7,5,6,
	6,9928,8,8,8,
	7,10660,10,13,13,
	8,7748,11,1,5,
	9,16151,5,10,11,
	10,14160,4,12,12,
	11,12033,9,11,3,
	12,12380,13,7,9,
	13,14101,3,4,4),nrow=13,ncol=5,byrow=TRUE)

pdf("AlnPlotTcris_240016H1_240179H1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        if(chtab[i,4] > 9){
                val<-paste(chtab[i,4],sep="")
        }else{
                val<-paste("0",chtab[i,4],sep="")
        }
        tcr_g1<-grep(x=subDfdat[,14],paste("hap1_Chr",val,sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g2<-grep(x=subDfdat[,10],pattern=paste("hap1_Chr",val,sep=""))
        cc<-tcr_g1[tcr_g1 %in% tcr_g2]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="24_0179 H1",ylab="24_0016 H1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()

##############
## 0179h1 x 0179h2
##############
dat<-fread("out_synteny_TcrE240179H1_TcrE240179H2.psl",header=FALSE)
dfdat<-as.data.frame(dat)

## target = 0179h1
## query = 0179h2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
g1Ch<-names(xx)[xx>200]
xx<-table(dfdat[,10])
g2Ch<-names(xx)[xx>200]
keep<-(dfdat[,14] %in% g1Ch) & (dfdat[,10] %in% g2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_g1<-as.numeric(gsub(pattern="hap1_Chr",x=colnames(tab),replacement=""))
tc_g2<-as.numeric(gsub(pattern="hap2_Chr",x=rownames(tab),replacement=""))

## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

pdf("SynTcrisStripe_240179h1_240179h2.pdf",width=6,height=6)## 1:1 chroms
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (24_0179 H2)",ylab="T. cristinae (24_0179 H1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_g1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_g2,las=2)
box()
dev.off()

pdf("AlnPlotTcris_240179H1_240179H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g1<-grep(x=subDfdat[,14],paste("hap1_Chr",val,sep=""))
        if(chtab[i,5] > 9){
                val<-paste(chtab[i,5],sep="")
        }else{
                val<-paste("0",chtab[i,5],sep="")
        }
        tcr_g2<-grep(x=subDfdat[,10],pattern=paste("hap2_Chr",val,sep=""))
        cc<-tcr_g1[tcr_g1 %in% tcr_g2]
        subd<-subDfdat[cc,]
        xub<-max(subd[,13]);yub<-max(subd[,17])

        plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="24_0179 H2",ylab="24_0179 H1")
        title(main=paste("Chromosome",i),cex.main=1.4)
        N<-dim(subd)[1]
        for(j in 2:N){
                if(subd[j,9]=="++"){
                        lines(subd[j,12:13],subd[j,16:17])
                }
                else{
                        lines(subd[j,12:13],yub-subd[j,16:17])
                }
        }
}
dev.off()
