arg3 <- commandArgs(T)
if(length(arg3) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}
print(arg3)
#arg3 <- c("/home/wmj/02_Metagenome/16s/ZW_20190918/Analysis/15_RDA_and_Heatmap/Total", "A" ,"B" , "C" ,"D"  , "/home/wmj/02_Metagenome/16s/ZW_20190918/bin/Total-color.txt" )
setwd(arg3[1])
######################	参数
arg1 <- arg3[2:(length(arg3)-1)]


###########	自定义函数区

###Normaliztion,对行
F_normalize <- function(x) {
	center <- sweep(x, 1, apply(x, 1, min),'-') #在行的方向上减去最小值
	R <- apply(x, 1, max) - apply(x,1,min)   #算出极差，即列上的最大值-最小值
	x_star<- sweep(center, 1, R, "/")        #把减去均值后的矩阵在行的方向上除以极差向量
	x_star <- x_star
	return(x_star)
}

##########	去除重复的种属信息
F_unique <- function(data){
	t1 <- as.character(unique(data))
	if(sum(is.na(t1)) >0){
		tt <- t1[-which(is.na(t1))]
	} else {
		tt <- t1
	}
	
	if (length(tt)==0) return(NA)
	
	yy <- matrix(NA,ncol=3)
	for (i in 1:length(tt)){
		at <- which(data == tt[i])

		tar <- c()
		temp <- at
		while(TRUE){
			max <- length(temp) +temp[1]-1
			d1 <- c(temp[1]:max)
			d2 <- which(d1 != temp)
			if(length(d2) == 0){
				cc <- c(temp[1],temp[length(temp)],tt[i])
				yy <- rbind(yy,cc)
				break
			} else{
				cc <-c(temp[1],temp[d2[1]-1],tt[i]) 
				yy <- rbind(yy,cc)
				temp <- temp[d2]
			}

		}
		
		
	}
	return(yy)
	
	
}

#####################################################
grpInfo = read.table(file=arg3[length(arg3)],header=T,sep="\t",check.names = FALSE, comment.char = "")
rownames(grpInfo)=as.character(grpInfo[,1])


gut <- read.csv(file="04_gut.csv",header=T,row.names=1)
ii=intersect(rownames(grpInfo),colnames(gut))
grpInfo=grpInfo[ii,]
gut=gut[,ii]
gut <- F_normalize(gut)
abundance <- read.csv(file="04_abundance.csv",header=T,row.names=1)
n1 <- which(is.na(abundance[,"Phylum"]))
if(length(n1)>0){
	abundance <- abundance[-n1,]
}
ii <- rownames(abundance)
gut <- gut[ii,]

###########	HEATMAP 准备
data <- t(gut)
data_row <- rownames(data)


######	HEATMAP 颜色
colorsChoice<- colorRampPalette(c("DarkSlateGray","WHITE","red"))
pal=colorsChoice(100)
breaks<-seq(0,1,length.out=101)


#####	HEATMAP 布局
len1 <- length(arg1)+1
t1 <- 2*len1-1
c1 <- c(0:t1)
c2 <- rep(0,len1)
cc <- c(c1,c2)
lay_heat <- matrix(cc,nrow=3,ncol=len1,byrow=T)

if (len1 %% 2 == 0){ ##求余
	lay_heat[3,(len1/2):(len1/2+1)] <- t1+1
} else {
	tt <- (len1 %/% 2) +1 ##求整
	lay_heat[3,tt] <- t1+1
}




############	计算布局
######	GRID 布局
lay_grid <- matrix(c(1:3),ncol=1,nrow=3)
######	种属分布布局
lay_text <- matrix(c(1:3),ncol=1,nrow=3)

####	总布局
lay1 <- max(lay_heat) +lay_grid
lay2 <- cbind(lay_heat,lay1)
lay3 <- max(lay2) +lay_text
lay_final <- lay2



####### 输出文件长宽比
Hei <- 1.1*length(gut[,1])
wid <- rep(1,length(lay_final[1,-1]))
ww <- c(wid,sum(wid)-1)
Wid <- sum(ww) * length(gut[1,])/(ncol(lay_final)-2)
if(Wid < 100) Wid=100
ff <- paste("Final.",nrow(gut),".pdf",sep="")

pdf(ff,height=Hei,width=Wid)
#pdf(ff,height=150,width=100)
layout(mat=lay_final,,heights=c(1,10,1),width=ww)

################	HEATMAP
###################	第一行


for (i in 1:length(arg1)){
	par(mar= c(0,1,0,1),mex=0.5,lwd=10,cex=8,font=2)
	plot(0,0,xlab="",ylab="",axes=F,xlim=c(1,12),ylim=c(0,0.5),col="white")
	segments(x0=1,y0=0,x1=11,y1=0) ### 线段
	m1 <- arg1[i]
	text(x=7, y=0.1,srt = 0, adj = 0.5, labels =m1,xpd = TRUE) ###标签
}

################	第二行
###### 括号
nco <- length(abundance[1,])-4
len_down <- length(which(abundance[,nco] == 0)) 	###下
len_up <- length(which(abundance[,nco] == 1))	###上
len_all <- length(abundance[,1]) 

par(mar = c(0,2,2,1),mex=0.5,cex=6,lwd=10,font=2,xaxs = "i", yaxs = "i",adj=1)
yy1 <- c(1:len_all)
plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,4),ylim=c(0,len_all),col="white")

if(len_down==0){
	y.top=0
} else {
	y.top=yy1[len_down]

}

rect(xleft = 0, ybottom = 0, xright = 4, ytop = y.top,col="#E1FFFF",border="#E1FFFF") ### 0-down
rect(xleft = 0, ybottom = y.top, xright = 4, ytop = yy1[len_all],col="#FFF0F5",border="#FFF0F5") ### 1-up
text(x=4, y=yy1-0.5,srt = 0, labels = colnames(data),xpd = TRUE) ###X





###  Heatmap
for (i in 1:length(arg1)){
	par(mar = c(0,0,2,1),mex=0.5,cex=6,lwd=1,font=2)
	m1 <- paste("^",arg1[i],"$*",sep="")
	c0=grep(m1,as.character(grpInfo[,2]))
	c1=as.character(grpInfo[c0,1])

	d2 <- data[c1,]
	x2 <- c(1:nrow(d2))
	y2 <- c(1:ncol(d2))
	lab <- length(d2[,1])
	image(x=x2,y=y2,z=d2,xlab="",ylab="",breaks=breaks,col=pal,axes=FALSE)
	text(x=x2+0.2, y=par("usr")[4] +0.5,srt =0, adj = 1, labels = c(1:lab),xpd = TRUE) ###X
	abline(h=c(0:ncol(d2))+0.5,v=c(0:nrow(d2))+0.5,col="black",xpd=F) ###网格

}

######## 第三行
### 图例
breaks2<-breaks[-length(breaks)]
par(mar = c(3,0,1,1),mex=0.5,lwd=1,font=2,cex.axis=2,cex=4)
#image(x=0:length(breaks2), y=1,z=matrix(breaks2)*1.001,col=pal[1:length(breaks)-1],axes=FALSE,breaks=breaks,xlab="",ylab="",xaxt="n")
image(x=0:length(breaks2), y=1,z=matrix(breaks2)*1.001,col=pal[1:length(breaks)-1],axes=FALSE,breaks=breaks,xlab="",ylab="",xaxt="n")
at_br <-seq(0,length(breaks2),length.out=3)
#la_br <- breaks[at_br+1]
la_br <- breaks[at_br+1]
abline(v=0:(length(breaks2)-1),col=pal[1:length(breaks)-1],xpd=F) ###网格
#axis(1,at=at_br,labels=la_br,col="white",las=1)
axis(1,at=at_br,labels=la_br,col="white",las=1)




########## TEXT
###########第一行
par(mar= c(0,1,0,1),mex=0.5,lwd=10,cex=8,font=2,xaxs = "i", yaxs = "i",adj=0.5)
	
#plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,4),ylim=c(0,1),col="white")
#segments(x0=0.1,y0=0,x1=0.9,y1=0)
#text(x=0.5, y=0.2,srt =0,, labels ="Species",xpd = TRUE) ###标签
#segments(x0=1.1,y0=0,x1=1.9,y1=0)
#text(x=1.5, y=0.2,srt =0, labels ="Genus",xpd = TRUE) 
#segments(x0=2.1,y0=0,x1=2.9,y1=0)
#text(x=2.5, y=0.2,srt =0, labels ="Family",xpd = TRUE) 
#segments(x0=3.1,y0=0,x1=3.9,y1=0)
#text(x=3.5, y=0.2,srt =0, labels ="Phylum",xpd = TRUE) 


plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,3),ylim=c(0,0.5),col="white")
segments(x0=0.1,y0=0,x1=0.9,y1=0)
text(x=0.5, y=0.1,srt =0,, labels ="Genus",xpd = TRUE) ###标签
segments(x0=1.1,y0=0,x1=1.9,y1=0)
text(x=1.5, y=0.1,srt =0, labels ="Family",xpd = TRUE) 
segments(x0=2.1,y0=0,x1=2.9,y1=0)
text(x=2.5, y=0.1,srt =0, labels ="Phylum",xpd = TRUE) 

#########第二行

c1 <- length(abundance[1,])
c2 <- c1-3
tt <-abundance[,c2:c1]
len_tt <- length(tt[,1])

######计算坐标与标签
lab <- tt
yy_t <-list("Species","Genus","Family","Phylum")
for(i in 2:length(tt[1,])){
	
	yy_t[[i]] <-  F_unique(tt[,i])
}




par(mar = c(0,1,2,1),mex=0.5,cex=6,lwd=10,font=4,xaxs = "i", yaxs = "i",adj=0)
yy3 <- c(1:len_tt)
#plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,4),ylim=c(0,len_tt),col="white")
#rect(xleft = 0, ybottom = 0, xright = 4, ytop = yy3[len_down],col="#E1FFFF",border="#E1FFFF") ### 0-down
#rect(xleft = 0, ybottom = yy3[len_down], xright = 4, ytop = yy3[len_tt],col="#FFF0F5",border="#FFF0F5") ### 1-up

plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,3),ylim=c(0,len_tt),col="white")
rect(xleft = 0, ybottom = 0, xright = 3, ytop = yy3[len_down],col="#E1FFFF",border="#E1FFFF") ### 0-down
rect(xleft = 0, ybottom = yy3[len_down], xright = 3, ytop = yy3[len_tt],col="#FFF0F5",border="#FFF0F5") ### 1-up

xx=0.4
for(i in 2:length(tt[1,])){	
	if(length(yy_t[[i]])==1){
		xx =xx+1
		next
	}
	
	yy0 <- as.numeric(yy_t[[i]][,1])-0.5
	yy1 <- as.numeric(yy_t[[i]][,2])-0.5

	yy2 <- (yy0+yy1)/2
	

	xx0 <- rep(xx,length(yy2))

	text(x=xx0,y=yy2,srt =0, labels =yy_t[[i]][,3],xpd = TRUE)
	segments(x0=xx-0.05,y0=yy0,x1=xx-0.05,y1=yy1)
	segments(x0=xx-0.1,y0=yy0,x1=xx-0.05,y1=yy0)
	segments(x0=xx-0.1,y0=yy1,x1=xx-0.05,y1=yy1)
	
	xx =xx+1
	
}

###############第三行
plot(0,0,xlab="",ylab="",axes=F,xlim=c(1,5.3),ylim=c(0,0.2),col="white")

dev.off()

###############################################
###############################################
###############################################
ff0 <- paste("Final.",nrow(gut),".png",sep="")

png(ff0,height=Hei*480/7,width=Wid*480/7)
#pdf(ff,height=150,width=100)
layout(mat=lay_final,,heights=c(1,10,1),width=ww)

################	HEATMAP
###################	第一行


for (i in 1:length(arg1)){
	par(mar= c(0,1,0,1),mex=0.5,lwd=10,cex=8,font=2)
	plot(0,0,xlab="",ylab="",axes=F,xlim=c(1,12),ylim=c(0,0.5),col="white")
	segments(x0=1,y0=0,x1=11,y1=0) ### 线段
	m1 <- arg1[i]
	text(x=7, y=0.1,srt = 0, adj = 0.5, labels =m1,xpd = TRUE) ###标签
}

################	第二行
###### 括号
nco <- length(abundance[1,])-4
len_down <- length(which(abundance[,nco] == 0)) 	###下
len_up <- length(which(abundance[,nco] == 1))	###上
len_all <- length(abundance[,1]) 

par(mar = c(0,2,2,1),mex=0.5,cex=6,lwd=10,font=2,xaxs = "i", yaxs = "i",adj=1)
yy1 <- c(1:len_all)
plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,4),ylim=c(0,len_all),col="white")

if(len_down==0){
	y.top=0
} else {
	y.top=yy1[len_down]

}

rect(xleft = 0, ybottom = 0, xright = 4, ytop = y.top,col="#E1FFFF",border="#E1FFFF") ### 0-down
rect(xleft = 0, ybottom = y.top, xright = 4, ytop = yy1[len_all],col="#FFF0F5",border="#FFF0F5") ### 1-up
text(x=4, y=yy1-0.5,srt = 0, labels = colnames(data),xpd = TRUE) ###X





###  Heatmap
for (i in 1:length(arg1)){
	par(mar = c(0,0,2,1),mex=0.5,cex=6,lwd=1,font=2)
	m1 <- paste("^",arg1[i],"$*",sep="")
	c0=grep(m1,as.character(grpInfo[,2]))
	c1=as.character(grpInfo[c0,1])

	d2 <- data[c1,]
	x2 <- c(1:nrow(d2))
	y2 <- c(1:ncol(d2))
	lab <- length(d2[,1])
	image(x=x2,y=y2,z=d2,xlab="",ylab="",breaks=breaks,col=pal,axes=FALSE)
	text(x=x2+0.2, y=par("usr")[4] +0.5,srt =0, adj = 1, labels = c(1:lab),xpd = TRUE) ###X
	abline(h=c(0:ncol(d2))+0.5,v=c(0:nrow(d2))+0.5,col="black",xpd=F) ###网格

}

######## 第三行
### 图例
breaks2<-breaks[-length(breaks)]
par(mar = c(3,0,1,1),mex=0.5,lwd=1,font=2,cex.axis=2,cex=4)
#image(x=0:length(breaks2), y=1,z=matrix(breaks2)*1.001,col=pal[1:length(breaks)-1],axes=FALSE,breaks=breaks,xlab="",ylab="",xaxt="n")
image(x=0:length(breaks2), y=1,z=matrix(breaks2)*1.001,col=pal[1:length(breaks)-1],axes=FALSE,breaks=breaks,xlab="",ylab="",xaxt="n")
at_br <-seq(0,length(breaks2),length.out=3)
#la_br <- breaks[at_br+1]
la_br <- breaks[at_br+1]
abline(v=0:(length(breaks2)-1),col=pal[1:length(breaks)-1],xpd=F) ###网格
#axis(1,at=at_br,labels=la_br,col="white",las=1)
axis(1,at=at_br,labels=la_br,col="white",las=1)




########## TEXT
###########第一行
par(mar= c(0,1,0,1),mex=0.5,lwd=10,cex=8,font=2,xaxs = "i", yaxs = "i",adj=0.5)
	
#plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,4),ylim=c(0,1),col="white")
#segments(x0=0.1,y0=0,x1=0.9,y1=0)
#text(x=0.5, y=0.2,srt =0,, labels ="Species",xpd = TRUE) ###标签
#segments(x0=1.1,y0=0,x1=1.9,y1=0)
#text(x=1.5, y=0.2,srt =0, labels ="Genus",xpd = TRUE) 
#segments(x0=2.1,y0=0,x1=2.9,y1=0)
#text(x=2.5, y=0.2,srt =0, labels ="Family",xpd = TRUE) 
#segments(x0=3.1,y0=0,x1=3.9,y1=0)
#text(x=3.5, y=0.2,srt =0, labels ="Phylum",xpd = TRUE) 


plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,3),ylim=c(0,0.5),col="white")
segments(x0=0.1,y0=0,x1=0.9,y1=0)
text(x=0.5, y=0.1,srt =0,, labels ="Genus",xpd = TRUE) ###标签
segments(x0=1.1,y0=0,x1=1.9,y1=0)
text(x=1.5, y=0.1,srt =0, labels ="Family",xpd = TRUE) 
segments(x0=2.1,y0=0,x1=2.9,y1=0)
text(x=2.5, y=0.1,srt =0, labels ="Phylum",xpd = TRUE) 

#########第二行

c1 <- length(abundance[1,])
c2 <- c1-3
tt <-abundance[,c2:c1]
len_tt <- length(tt[,1])

######计算坐标与标签
lab <- tt
yy_t <-list("Species","Genus","Family","Phylum")
for(i in 2:length(tt[1,])){
	
	yy_t[[i]] <-  F_unique(tt[,i])
}




par(mar = c(0,1,2,1),mex=0.5,cex=6,lwd=10,font=4,xaxs = "i", yaxs = "i",adj=0)
yy3 <- c(1:len_tt)
#plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,4),ylim=c(0,len_tt),col="white")
#rect(xleft = 0, ybottom = 0, xright = 4, ytop = yy3[len_down],col="#E1FFFF",border="#E1FFFF") ### 0-down
#rect(xleft = 0, ybottom = yy3[len_down], xright = 4, ytop = yy3[len_tt],col="#FFF0F5",border="#FFF0F5") ### 1-up

plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,3),ylim=c(0,len_tt),col="white")
rect(xleft = 0, ybottom = 0, xright = 3, ytop = yy3[len_down],col="#E1FFFF",border="#E1FFFF") ### 0-down
rect(xleft = 0, ybottom = yy3[len_down], xright = 3, ytop = yy3[len_tt],col="#FFF0F5",border="#FFF0F5") ### 1-up

xx=0.4
for(i in 2:length(tt[1,])){	
	if(length(yy_t[[i]])==1){
		xx =xx+1
		next
	}
	
	yy0 <- as.numeric(yy_t[[i]][,1])-0.5
	yy1 <- as.numeric(yy_t[[i]][,2])-0.5

	yy2 <- (yy0+yy1)/2
	

	xx0 <- rep(xx,length(yy2))

	text(x=xx0,y=yy2,srt =0, labels =yy_t[[i]][,3],xpd = TRUE)
	segments(x0=xx-0.05,y0=yy0,x1=xx-0.05,y1=yy1)
	segments(x0=xx-0.1,y0=yy0,x1=xx-0.05,y1=yy0)
	segments(x0=xx-0.1,y0=yy1,x1=xx-0.05,y1=yy1)
	
	xx =xx+1
	
}

###############第三行
plot(0,0,xlab="",ylab="",axes=F,xlim=c(1,5.3),ylim=c(0,0.2),col="white")

dev.off()

###############################################




