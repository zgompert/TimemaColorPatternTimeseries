library(data.table)

aln<-fread("syn_blocks.txt",header=TRUE)

## 1 + number of genomes * 3 = 103
startID<-seq(3,103,3)
endID<-seq(4,103,3)
strandID<-seq(2,103,3)

## make distinct matrixes
aln<-as.data.frame(aln)
alnStart<-as.matrix(aln[,startID])
alnEnd<-as.matrix(aln[,endID])
alnStrand<-aln[,strandID]

sp<-gsub(pattern="_ch8.start",colnames(aln)[startID],replacement="")
sp8<-gsub(pattern=".start",colnames(aln)[startID],replacement="")

## dotplots
Nsp<-length(sp)
cs<-c("darkorange3","dimgray")
pdf("tcris_ch8_dotplots.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(4.5,4.5,1,1))
for(i in 1:(Nsp-1)){
	for(j in (i+1):Nsp){
		cat(i,",",j,"\n")
		midA<-(alnStart[,i] + alnEnd[,i])/2
		midB<-(alnStart[,j] + alnEnd[,j])/2
		cid<-(alnStrand[,i]==alnStrand[,j])+1
		plot(midA,midB,pch=19,xlab=sp[i],ylab=sp[j],col=cs[cid],cex.lab=1.2)
	}
}
dev.off()

## seems to work, try ribbon plot
chinfo<-read.table("chsizes.txt",header=TRUE)
## alternative order for plot
chinfo<-read.table("alt_chsizes.txt",header=TRUE)

## get sort, synteny is hash order
ord<-rep(NA,34)
for(i in 1:34){
	ord[i]<-which(sp8==chinfo[i,3])
}

# flip alignments as needed
RalnStart<-alnStart
RalnEnd<-alnEnd
flip<-c(3:12,21,22,27:29,31:32)
#flip<-c(1,6,10,14,17,19:22)
for(i in flip){
	a<-ord[i]
	RalnStart[,a]<-chinfo[i,2]-alnEnd[,a]
	RalnEnd[,a]<-chinfo[i,2]-alnStart[,a]
}

## green, melanic, stripe
cs<-c("chartreuse4","chocolate4","cadetblue")
cid<-as.numeric(as.factor(chinfo$color))
lid<-1+(chinfo$location=="r")

## plot
library(scales)
pdf("tcris_ch8_ribbon.pdf",width=8,height=11)
par(mar=c(4.5,7,.5,.5))
plot(c(0,1),xlim=c(1,max(chinfo[,2])),ylim=c(.5,34.5),type='n',
	     xlab="Position (bp)",ylab="",axes=FALSE)
for(i in 1:(Nsp-1)){
	j<-i+1
	a<-ord[i];b<-ord[j]
	for(k in 1:dim(alnStart)[1]){
		polygon(c(RalnStart[k,a],RalnEnd[k,a],
			  RalnEnd[k,b],RalnStart[k,b]),
			c(i,i,j,j),border=NA,col=alpha("gray",.5))
	}
}
for(i in 1:34){
	lines(c(1,chinfo[i,2]),c(i,i),lwd=2,col=cs[cid[i]],lty=lid[i])
}
axis(1)
par(xpd=TRUE)
text(rep(-1.2e7,34),1:34,chinfo$ID)
par(xpd=FALSE)
dev.off()

## dotplots sorted
ord2<-rep(NA,34)
for(i in 1:34){
        ord2[i]<-which(chinfo[,3]==sp8[i])
}
Nsp<-length(sp)
cs<-c("darkorange3","dimgray")
pdf("tcris_ch8_dotplots2.pdf",height=12,width=8.57)
par(mfrow=c(7,5))
par(mar=c(3,3,.5,.5))
for(i in ord){
        for(j in ord){
                cat(i,",",j,"\n")
                if(i != j){
                        midA<-(RalnStart[,i] + RalnEnd[,i])/2
                        midB<-(RalnStart[,j] + RalnEnd[,j])/2
                        cid<-(alnStrand[,i]==alnStrand[,j])+1
                        plot(midA,midB,pch=19,xlab="",ylab="",col=cs[cid],cex.lab=1.1,axes=FALSE)
                        mtext(chinfo$ID[ord2][i],side=1,line=.7)
                        mtext(chinfo$ID[ord2][j],side=2,line=.7)
                        box()
                } else{
                        plot(0,0,type='n',axes=FALSE,xlab="",ylab="")
                }
        }
        plot(0,0,type='n',axes=FALSE,xlab="",ylab="")
}
dev.off()
