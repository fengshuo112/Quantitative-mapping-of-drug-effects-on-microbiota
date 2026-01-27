
options(stringsAsFactors = F)
#setwd("/Volumes/Flower0501/16s/Qiime1")
arg1 <- commandArgs(T)
if(length(arg1) < 0){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg1)
library(stringr)
#arg1=c("/Volumes/Flower0501/2022_16s/CY/QC_Raw/*LengthCount.txt","/Volumes/Flower0501/2022_16s/CY/QC_Raw/Result_Table_Info_Raw.txt", "/Volumes/Flower0501/2022_16s/CY/QC_Merge/Result_Table_Info_Merge.txt", "/Volumes/Flower0501/2022_16s/CY/QC_End/Result_Length_End.txt", "/Volumes/Flower0501/2022_16s/CY/QC_End/Result_Table_Info_End.txt","/Volumes/Flower0501/2022_16s/CY/QC_End" )

dir=arg1[length(arg1)]
#dir = "QC_End"
############## write.table
d1=read.table(arg1[2],sep='\t',header=T)
#d1=read.table("Result_Table_Info_Raw.txt",sep = "\t",header = T)
r1=as.character(d1[,1])
#=gsub('\\_R1','',r1)
#r3=gsub('\\_R2','',r2)
r3=unlist(lapply(r1,function(x){str_split(x,"_")[[1]][1]}))
d1[,1]=r3
#d11=unique(d1[,1:2])
d11 <- aggregate(d1[,2],by=list(SampleID=d1$SampleID),max)
colnames(d11) <- colnames(d1)[1:2]
rownames(d11)=as.character(d11[,1])

#d2=read.table("Result_Table_Info_Merge.txt",sep='\t',header=T,row.names=1)
d2=read.table(arg1[3],sep='\t',header=T,row.names=1)
#r1=gsub('\\.join','',rownames(d2))
r1=unlist(lapply(rownames(d2),function(x){str_split(x,"_")[[1]][1]}))
rownames(d2)=r1
d2=d2[rownames(d11),]
#d3=read.table("QC_End/Result_Table_Info_End.txt",sep='\t',header=T,row.names=1)
d3=read.table(arg1[5],sep='\t',header=T,row.names=1)
#row_order = gsub("_","",rownames(d11))
r4 <- unlist(lapply(rownames(d3),function(x){str_split(x,"_")[[1]][1]}))
rownames(d3) <-r4
row_order <- rownames(d11)
d3=d3[row_order,]
dd_re=cbind(d11,d2,d3)
colnames(dd_re)=c('SampleID','PE Reads','Merge Tags','Length_Merge(bp)','AvgLen_Merge(bp)','GC_Merge(%)','Q20_Merge(%)','Q30_Merge(%)','SubSample Tags','AvgLen_Sub(bp)')

write.table(dd_re,paste0(dir,'/Table.txt'),sep='\t',quote=F,row.names=F,col.names=T)



##################### arg[1] Raw
arg <- Sys.glob(arg1[1])
#arg=Sys.glob("QC_Raw/*LengthCount.txt")
dd1=c()
dd2=c()
for(i in 1:length(arg)){
	data=read.table(arg[i],sep='\t',header=T)
	d1=as.character(data[,1])
	d2=as.numeric(data[,2])
	dd1=c(dd1,d1)
	dd2=c(dd2,d2)
}

q1=unique(dd1)
re1=c()
for(i in 1:length(q1)){
	g1=grep(paste0('^',q1[i],'$'),dd1)
	r1=sum(dd2[g1])
	re1=c(re1,r1)
}
names(re1)=q1

xy0=matrix(NA,nrow=length(q1),ncol=3)
colnames(xy0)=c('x','y','num')
for(i in 1:length(q1)){
	g1=unlist(strsplit(q1[i],'-'))
	if(length(g1)==1){
		xy0[i,1]=as.numeric(g1[1])
		xy0[i,2]=as.numeric(g1[1])
		xy0[i,3]=as.numeric(re1[i])
	}
	
	if(length(g1)==2){
		xy0[i,1]=as.numeric(g1[1])
		xy0[i,2]=as.numeric(g1[2])
		xy0[i,3]=as.numeric(re1[i])
	}
}

ord=order(xy0[,1],decreasing=F)
xy0=xy0[ord,]
xy=xy0
re=matrix(NA,nrow=0,ncol=3)
if(class(xy)=='numeric'){
	re=rbind(re,xy)
} else {
	for(i in 1:(nrow(xy)-1)){
		x1=xy[i,2]
		x2=xy[i+1,1]
		x3=x2-x1

		if(x3==1){
			re=rbind(re,xy[i,])
		} 
		if(x3!=1){
			x1=min(xy[i,1],xy[i+1,1])
			x2=max(xy[i,2],xy[i+1,2])
			xs=xy[i,3]+xy[i+1,3]
			c4=c(x1,x2,xs)
			xy[i+1,1]=c4[1]
			xy[i+1,2]=c4[2]
			xy[i+1,3]=c4[3]
			if(i == (nrow(xy)-1) ){
				re=rbind(re,xy[i,])
			}
		}	
	
	}
}
cc=c()
for(i in 1:nrow(re)){
	cc=c(cc,paste0(re[i,1],'-',re[i,2]))
}

r_text=cbind(re,cc)
#	write.table(r_text,paste0(dir,'/Pic_length.txt'),sep='\t',quote=F,col.names=T,row.names=F)

pdf(paste0(dir,'/01_Raw_Data_Length.pdf'))
par(mar=c(6,6,6,4))
d1=log10(as.numeric(re[,3]))
i1=which(d1<0)
d1[i1]=0
xx=barplot(d1,col='lightblue',ylab='',xlab='',main='Length Distribution',xaxt='n')
text(x=xx,y=par('usr')[3]-0.1,adj=1,labels=cc,xpd=T,srt=45)
mtext(side=1,adj=0.5,line=3.5,cex=1.5,text='Length(bp)')
mtext(side=2,adj=0.5,line=3,cex=1.5,text='log( Read Counts )')

mtext(side=3,adj=0.5,line=4,cex=1.5,text='Analysis of Raw Data',xpd=T)

dev.off()

############################ arg[2]

data=read.table(arg1[2],sep='\t',header=T,row.names=1)
#data=read.table("Result_Table_Info_Raw.txt",,sep='\t',header=T,row.names=1)

######## sequence

dd=as.numeric(data[,1])
ord=order(dd)
dd1=dd[ord]
#p<-ggplot(data, aes(x =Sequence))
#p + geom_density(color = "black", fill = "cornflowerblue",alpha=0.4)

pdf(paste0(dir,'/02_Raw_Data_Sequence.pdf'))
par(mar=c(6,6,6,4))
mm=mean(dd1)
plot(dd1,col='cornflowerblue',xlab='Sample',ylab='Number of Sequence',main="Distribution of Samples' Sequence",cex.lab=1.5,type='l',lwd=2)
abline(h = mm, col = "red",lwd=2,lty=2)
#text(x=mean(par('usr')[c(1,2)]),y=mm*1.2,labels=paste0('mean = ',mm),cex=1.3)
m1=round(mm,2)
text(x=mean(par('usr')[c(1,2)]),y=mm,labels=paste0('mean = ',m1),cex=1.3)
dev.off()
#plot(dd1,col='cornflowerblue',xlab='Sample',ylab='Sequence',main="Distribution of Samples' Sequence",cex.lab=1.5,type='b',cex=0.1)
######## Avglen 

d2=as.numeric(data[,3])
mm=round(mean(d2),2)
pdf(paste0(dir,'/03_Raw_Data_AveLen.pdf'))
par(mar=c(6,6,6,4))
stripchart(d2, vertical = TRUE,method = "jitter", pch = 8,col = 'lightgreen',cex=1.5,ylab='AveLen(bp)',xlab=paste0('N=',nrow(data)),main="Distribution of AvgLen",cex.lab=1.5)
abline(h = mm, col = "red",lwd=2,lty=2)
m1=round(mm,2)
text(x=mean(par('usr')[c(1,2)]),y=mm,labels=paste0('mean = ',m1),cex=1.3)
dev.off()
######## QC

q20=as.numeric(data[,5])*100
q30=as.numeric(data[,6])*100
dd=list(q20,q30)
names(dd)=c('Q20','Q30')
pdf(paste0(dir,'/04_Raw_Data_QC.pdf'))
par(mar=c(6,6,6,4))
stripchart(dd, vertical = TRUE,method = "jitter", pch = c(3,4),col = c('lightgreen','lightblue'),ylab='Percentage %',xlab="Quality Score",main="Distribution of Quality",cex.lab=1.5)#yaxt='n'
abline(h = c(90,95), col = c('springgreen4','lightgreen'),lwd=c(2,1),lty=2)
abline(h = c(80,85), col = c('steelblue','lightblue'),lwd=c(2,1),lty=2)
axis(side=4,at=c(80,90),labels=c(80,90),las=2)
dev.off()

png(paste0(dir,'/04_Raw_Data_QC.png'), bg="transparent")
par(mar=c(6,6,6,4))
stripchart(dd, vertical = TRUE,method = "jitter", pch = c(3,4),col = c('lightgreen','lightblue'),ylab='Percentage %',xlab="Quality Score",main="Distribution of Quality",cex.lab=1.5)#yaxt='n'
abline(h = c(90,95), col = c('springgreen4','lightgreen'),lwd=c(2,1),lty=2)
abline(h = c(80,85), col = c('steelblue','lightblue'),lwd=c(2,1),lty=2)
axis(side=4,at=c(80,90),labels=c(80,90),las=2)
dev.off()


###########################################
###### Subsample
##################### arg[1]
#dd=read.table("QC_End/Result_Length_End.txt",sep='\t',header=T)
dd=read.table(arg1[4],sep='\t',header=T)
ord=order(as.numeric(dd[,1]))
dd=dd[ord,]
l1=10
mm=ceiling((max(dd[,1])-min(dd[,1]))/l1)
while(mm > 30){
	l1=l1+5
	mm=ceiling((max(dd[,1])-min(dd[,1]))/l1)
}

re_sub=matrix(NA,nrow=0,ncol=3)
ll=c()
for(i in 1:mm){
	if(i==1){
		x1=dd[1,1]
		x2=x1+l1
	} else {
		x1=x2+1
		x2=x1+l1
	}
	if(x2>dd[nrow(dd),1]) x2=dd[nrow(dd),1]
	i1=which(dd[,1] >= x1)
	i2=which(dd[,1] <= x2)
	i3=intersect(i1,i2)
	if(length(i3)>0){
		c3=sum(dd[i3,2])
	} else {
		c3=0
	}
	ll=c(ll,paste0(x1,'-',x2))
	r1=c(x1,x2,c3)
	re_sub=rbind(re_sub,r1)
	if(x2==dd[nrow(dd),1]) break
}
colnames(re_sub)=c('x','y','num')



d1=log10(as.numeric(re_sub[,3]))
i1=which(d1<0)
d1[i1]=0
pdf(paste0(dir,'/05_SubSample_Length.pdf'))
par(mar=c(6,6,6,4))
xx=barplot(d1,col='lightblue',ylab='',xlab='',main='Length Distribution')
text(x=xx,y=par('usr')[3]-0.1,adj=1,labels=ll,xpd=T,srt=45)
mtext(side=1,adj=0.5,line=3.5,cex=1.5,text='Length(bp)')
mtext(side=2,adj=0.5,line=3,cex=1.5,text='log( Read Counts )')

mtext(side=3,adj=0.5,line=4,cex=1.5,text='Analysis of SubSample Data',xpd=T)

dev.off()

png(paste0(dir,'/05_SubSample_Length.png'), bg="transparent")
par(mar=c(6,6,6,4))
xx=barplot(d1,col='lightblue',ylab='',xlab='',main='Length Distribution')
text(x=xx,y=par('usr')[3]-0.1,adj=1,labels=ll,xpd=T,srt=45)
mtext(side=1,adj=0.5,line=3.5,cex=1.5,text='Length(bp)')
mtext(side=2,adj=0.5,line=3,cex=1.5,text='log( Read Counts )')

mtext(side=3,adj=0.5,line=4,cex=1.5,text='Analysis of SubSample Data',xpd=T)

dev.off()
############################ arg[5]

data=dd_re[,c(2,3,9)]


######## Avglen 

d2=as.numeric(dd_re[,10])
mm=mean(d2)
pdf(paste0(dir,'/06_SubSample_AveLen.pdf'))
par(mar=c(6,6,6,4))
stripchart(d2, vertical = TRUE,method = "jitter", pch = 8,col = 'lightgreen',cex=1.5,ylab='AveLen(bp)',xlab=paste0('N=',nrow(data)),main="Distribution of AvgLen",cex.lab=1.5)
abline(h = mm, col = "red",lwd=2,lty=2)
#text(x=mean(par('usr')[c(1,2)]),y=mm,labels=paste0('mean = ',mm),cex=1.3)
m1=round(mm,2)
text(x=mean(par('usr')[c(1,2)]),y=mm,labels=paste0('mean = ',m1),cex=1.3)
dev.off()
######## barplot

d1=as.numeric(data[,1])
d2=as.numeric(data[,2])
d3=as.numeric(data[,3])

tt1=max(d3)
if(is.na(tt1)){
	tt2=d3[!is.na(d3)]
	y1=max(tt2)*0.8
} else {
	y1=max(d3)*0.8
}

tt3=max(d1)
if(is.na(tt3)){
	tt4=d1[!is.na(d1)]
	y2=max(tt2)
} else {
	y2=max(d1)
}
pdf(paste0(dir,'/07_SubSample_Sequence.pdf'))
par(mar=c(6,6,6,2))
plot(d1,col='lightgreen',type='l',ylab='Sequence',xlab='Sample',main='Sequence Distribution',cex.lab=1.5,ylim=c(y1,y2),lwd=2)
lines(d2,col='lightcoral',type='l',lwd=2)
lines(d3,col='cornflowerblue',type='l',lwd=2)
legend('topright',legend=c("PE Reads","Merge Tags",'SubSample Tags'),lty=1,col=c('lightgreen','lightcoral','cornflowerblue'),border=NA,bty='n',xpd=T,lwd=2)
dev.off()
png(paste0(dir,'/07_SubSample_Sequence.png'), bg="transparent")
par(mar=c(6,6,6,2))
plot(d1,col='lightgreen',type='l',ylab='Sequence',xlab='Sample',main='Sequence Distribution',cex.lab=1.5,ylim=c(y1,y2),lwd=2)
lines(d2,col='lightcoral',type='l',lwd=2)
lines(d3,col='cornflowerblue',type='l',lwd=2)
legend('topright',legend=c("PE Reads","Merge Tags",'SubSample Tags'),lty=1,col=c('lightgreen','lightcoral','cornflowerblue'),border=NA,bty='n',xpd=T,lwd=2)
dev.off()
#############










