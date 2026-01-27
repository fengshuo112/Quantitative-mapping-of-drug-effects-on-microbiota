arg3 <- commandArgs(T)
if(length(arg3) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}
print(arg3)
#arg3=c("/home/wmj/视频/16s/qiime1/TT/Analysis/15_RDA_and_Heatmap/Total","T1Control"  ,"T1Classical" ,"T1T1H"  , "T1T1L", "/home/wmj/视频/16s/qiime1/TT/bin/Total-color.txt"  , "T1Control"   , "T1Classical" , "T1T1H" ,'0.05')
setwd(arg3[1])
P=as.numeric(arg3[length(arg3)])
######################	参数
gg=grep('*color.txt',arg3)

arg1 <- arg3[2:(gg-1)]	###Heatmap 
arg2 <-arg3[(gg+1):(length(arg3)-1)]	## Grid;第一位必须是基准;第二位必须是正常组，后面是给药组（或者相反）
#########	P值
###########	自定义函数区
### 判断形状
F_type <- function(data){
	circle <- which(data[1] < data) -1 # 1
	diamond <- which(data[1] > data) -1 #18
	len <- length(data)
	star <- which( (data[1] > data[2] & data[1] < data[3:len]) | (data[1] < data[2] & data[1] > data[3:len])) +len-1 #2
	all <- len+len-2-1
	Type <- rep(-1,all)
	Type[circle] <- 1
	Type[diamond] <- 18
	Type[star] <- 5
	return(Type)
}

###Normaliztion,对行
F_normalize <- function(x) {
	center <- sweep(x, 1, apply(x, 1, min),'-') #在行的方向上减去最小值
	R <- apply(x, 1, max) - apply(x,1,min)   #算出极差，即列上的最大值-最小值
	x_star<- sweep(center, 1, R, "/")        #把减去均值后的矩阵在行的方向上除以极差向量
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

grpInfo = read.table(file=arg3[gg],header=T,sep="\t",check.names = FALSE, comment.char = "")
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
data <- t(F_normalize(gut))
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


########	Grid 准备
len2 <- length(arg2)
gut_row <- rownames(gut)
gut_col <- colnames(gut)
##### 存储均值
abundance_grid <- matrix(NA,ncol=len2,nrow=length(gut_row))
colnames(abundance_grid) <- arg2
rownames(abundance_grid) <- gut_row

##### 存储形状
len_star <- len2-2
type <- matrix(NA,ncol=len2+len_star-1,nrow=length(gut_row))
star_name <- paste(arg2[-(1:2)],"Star",sep="-")

rownames(type) <- gut_row

######### 计算均值
for (i in 1:len2){
	c1 <- grep(paste0('^',arg2[i],'$'),as.character(grpInfo[,2]),value=F,fixed = F,ignore.case = F)
	c2=as.character(grpInfo[c1,1])
	tt <- gut[,c2]
	abundance_grid[,i] <- apply(tt,1,mean)
}


######### 计算形状
# 1:空心圆 +
# 18:实心菱形 -
# * : * -为方便，记为5
 type <- t(apply(abundance_grid,1,F_type))
colnames(type) <- c(arg2[-1],star_name)

## 挑选所选中的OTU的P值,显著
all_p <- read.csv(file="04_all_p.csv",header=T,row.names=1)
all_p_col <- colnames(all_p)
all_cc <- c()
for (i in 2:len2){
	m1 <- paste(arg2[1],".*",arg2[i],sep="")
	c1 <- grep(m1,all_p_col,value=F,fixed = F,ignore.case = F)
	all_cc <- c(all_cc,c1)
}
all_pp <- all_p[,all_cc]
for (i in 1:length(all_pp[1,])){
	type[which(all_pp[,i] > P ) ,i] <- -1
}
write.csv(type,file="04_type.csv",row.names=T)


############	计算布局
######	GRID 布局
lay_grid <- matrix(c(1:3),ncol=1,nrow=3)
######	种属分布布局
lay_text <- matrix(c(1:3),ncol=1,nrow=3)

####	总布局
lay1 <- max(lay_heat) +lay_grid
lay2 <- cbind(lay_heat,lay1)
lay3 <- max(lay2) +lay_text
lay_final <- cbind(lay2,lay3)



####### 输出文件长宽比
Hei <- 1.1*length(gut[,1])
wid <- rep(1,length(lay_final[1,-1]))
ww <- c(wid,sum(wid)-2)
Wid <- sum(ww) * length(gut[1,])/(length(lay_final[1,-1])-2)
if(Wid < 100) Wid=100
pdf("Final.pdf",height=Hei,width=Wid)
#pdf("Final.pdf",height=150,width=200)
layout(mat=lay_final,heights=c(1,20,1),width=ww)

################	HEATMAP
###################	第一行


for (i in 1:length(arg1)){
	#par(mai= c(0,0.1,0,0.1),mex=0.5,lwd=1,cex=3)
	par(mar= c(0,1,0,1),mex=0.5,lwd=10,cex=8,font=2)
	plot(0,0,xlab="",ylab="",axes=F,xlim=c(1,12),ylim=c(0,0.5),col="white")
	segments(x0=1,y0=0,x1=11,y1=0) ### 线段
	m1 <- arg1[i]
	text(x=7, y=0.1,srt = 0, adj = 1, labels =m1,xpd = TRUE) ###标签
}

################	第二行
###### 括号
#par(mai= c(0,0.1,0,0),mex=0.5,cex=0.3,lwd=1,xaxs = "i", yaxs = "i")
#par(mai= c(0,0.1,0.1,0),mex=0.5,cex=1,lwd=1,xaxs = "r", yaxs = "r")
nco <- length(abundance[1,])-4
len_down <- length(which(abundance[,nco] == 0)) 	###下
len_up <- length(which(abundance[,nco] == 1))	###上
len_all <- length(abundance[,1]) 

par(mar = c(0,2,2,1),mex=0.5,cex=6,lwd=10,font=2,xaxs = "i", yaxs = "i",adj=1)
yy1 <- c(1:len_all)
plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,4),ylim=c(0,len_all),col="white")

rect(xleft = 0, ybottom = 0, xright = 4, ytop = yy1[len_down],col="#E1FFFF",border="#E1FFFF") ### 0-down
rect(xleft = 0, ybottom = yy1[len_down], xright = 4, ytop = yy1[len_all],col="#FFF0F5",border="#FFF0F5") ### 1-up
text(x=4, y=yy1-0.5,srt = 0, labels = colnames(data),xpd = TRUE) ###X
#segments(x0=0.1,y0=yy1[1]-0.5,x1=0.1,y1=yy1[len_down]-0.5) ###down |
#segments(x0=0.1,y0=yy1[1]-0.5,x1=0.5,y1=yy1[1]-0.5)	## 下 -
#segments(x0=0.1,y0=yy1[len_down]-0.5,x1=0.5,y1=yy1[len_down]-0.5)	## 上 -
#segments(x0=0.1,y0=yy1[len_down+1]-0.5,x1=0.1,y1=yy1[len_all]-0.5)###up 
#segments(x0=0.1,y0=yy1[len_down+1]-0.5,x1=0.5,y1=yy1[len_down+1]-0.5)
#segments(x0=0.1,y0=yy1[len_all]-0.5,x1=0.5,y1=yy1[len_all]-0.5)







###  Heatmap
for (i in 1:length(arg1)){
	par(mar = c(0,0,2,1),mex=0.5,cex=6,lwd=1,font=2)
	m1 <- paste("^",arg1[i],"[0-9]*",sep="")
	c1 <- grep(m1,data_row,value=F,fixed = F,ignore.case = F)
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
par(mar = c(5,0,5,1),mex=0.5,lwd=1,font=2,cex.axis=2,cex=4)
#image(x=0:length(breaks2), y=1,z=matrix(breaks2)*1.001,col=pal[1:length(breaks)-1],axes=FALSE,breaks=breaks,xlab="",ylab="",xaxt="n")
image(x=0:length(breaks2), y=1,z=matrix(breaks2)*1.001,col=pal[1:length(breaks)-1],axes=FALSE,breaks=breaks,xlab="",ylab="",xaxt="n")
at_br <-seq(0,length(breaks2),length.out=3)
la_br <- breaks[at_br+1]
abline(v=0:(length(breaks2)-1),col=pal[1:length(breaks)-1],xpd=F) ###网格
axis(1,at=at_br,labels=la_br,col="white",las=1)


##############	GRID
#######第一行
par(mar= c(0,1,0,1),mex=0.5,lwd=10,cex=6,font=2,adj=0,xaxs = "i", yaxs = "i",adj=0.5)
plot(0,0,xlab="",ylab="",axes=F,xlim=c(1,len2),ylim=c(0,1),col="white")
text(x=c(2:len2)-0.5, y=0.1,srt =90, labels =arg2[-1],xpd = TRUE) ###标签

#######第二行
par(mar = c(0,1,2,1),mex=0.5,cex=6,lwd=10,font=2,xaxs = "i", yaxs = "i",mgp=c(1,0,0))
#mgp=c(3, 1, 0) 三个坐标成分的位置.第一个参数是轴标签相对轴位置的距离,以文本行作为参照单位的.第二个参数表示刻度标记的距离,最后一个参数是轴位置到轴线的距离(常常是0).正值表示在图形外,负值表示在图形内.
len_row <- length(gut_row)
yy2 <- c(1:len_row)
plot(0,0,xlab="",ylab="",axes=F,xlim=c(1,len2),ylim=c(0,len_row),col="white")
rect(xleft = 1, ybottom = 0, xright = len2, ytop = yy2[len_down],col="#E1FFFF",border="#E1FFFF") ### 0-down
rect(xleft = 1, ybottom = yy2[len_down], xright = len2, ytop = yy2[len_all],col="#FFF0F5",border="#FFF0F5") ### 1-up

abline(h=c(0:len_row),v=c(1:len2),col="black",xpd=F) ###网格
for(i in 1:(len2-1)){
	xx <- rep(i+0.5,len_row)
	points(x=xx,y=c(0:(len_row-1))+0.5,pch=type[,i])
}

cc <- c(len2:length(type[1,]))
tt <- apply(type[,cc],1,mean)
aa <- which(tt > 0)
axis(side=2,at=aa+0.2,labels=rep("*",length(aa)),col="black",las=1,tick=F)
#######第三行
plot(0,0,xlab="",ylab="",axes=F,xlim=c(1,5.3),ylim=c(0,0.2),col="white")

########## TEXT
###########第一行
par(mar= c(0,1,0,1),mex=0.5,lwd=5,cex=8,font=2,xaxs = "i", yaxs = "i",adj=0.5)
plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,4),ylim=c(0,1),col="white")
segments(x0=0.1,y0=0,x1=0.9,y1=0)
text(x=0.5, y=0.2,srt =0,, labels ="Species",xpd = TRUE) ###标签
segments(x0=1.1,y0=0,x1=1.9,y1=0)
text(x=1.5, y=0.2,srt =0, labels ="Genus",xpd = TRUE) 
segments(x0=2.1,y0=0,x1=2.9,y1=0)
text(x=2.5, y=0.2,srt =0, labels ="Family",xpd = TRUE) 
segments(x0=3.1,y0=0,x1=3.9,y1=0)
text(x=3.5, y=0.2,srt =0, labels ="Phylum",xpd = TRUE) 


#########第二行

c1 <- length(abundance[1,])
c2 <- c1-3
tt <-abundance[,c2:c1]
len_tt <- length(tt[,1])

######计算坐标与标签
lab <- tt
yy_t <-list("Species","Genus","Family","Phylum")
for(i in 1:length(tt[1,])){
	
	yy_t[[i]] <-  F_unique(tt[,i])
}




par(mar = c(0,1,2,1),mex=0.5,cex=6,lwd=10,font=4,xaxs = "i", yaxs = "i",adj=0)
yy3 <- c(1:len_tt)
plot(0,0,xlab="",ylab="",axes=F,xlim=c(0,4),ylim=c(0,len_tt),col="white")

rect(xleft = 0, ybottom = 0, xright = 4, ytop = yy3[len_down],col="#E1FFFF",border="#E1FFFF") ### 0-down
rect(xleft = 0, ybottom = yy3[len_down], xright = 4, ytop = yy3[len_tt],col="#FFF0F5",border="#FFF0F5") ### 1-up

xx=0.4
for(i in 1:length(tt[1,])){	
	
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


