arg1 <- commandArgs(T)
if(length(arg1) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}

print(arg1)
###########	测试
#arg1 <- c( "/home/wmj/02_Metagenome/16s/20191017_LY/林杨16s原始数据/Analysis/15_RDA_and_Heatmap/Total","NC-DM", "HNK-DM", "2" ,"/home/wmj/02_Metagenome/16s/20191017_LY/林杨16s原始数据/bin/Total-color.txt"          )
#arg1 <- c("/Volumes/Flower0501/2022_16s/CY/15_RDA_and_Heatmap","Z1-Z2",2,"/Volumes/Flower0501/2022_16s/CY/Total-color.txt")
library(vegan) 
setwd(arg1[1])
arg.tt <- arg1[2:(length(arg1)-2)]
arg <- arg.tt

###########	自定义函数区
#### 判断 Factor
F_Factor<- function(data){
	if(length(data)==1) return(as.numeric(data))
	t1 <- as.character(unique(data))
	if(length(t1) == 1){
		return(as.numeric(t1))
	} else {
		return(0.5)
	}
}

###	挑选种属信息
F_ifno <- function(data){
	tt <-unlist(strsplit( as.character(data),split=";"))
	cc <- c("s__","g__","f__","p__")
	mm <- c()
	for(i in 1:length(cc)){
		s <- paste("",cc[i],".*",sep="")
		s1 <- grep(s,tt,value=F,fixed = F,ignore.case = F)
		if(length(s1) == 0) {
			mm <- c(mm,NA)			
		} else {
			hh <- unlist(strsplit( tt[s1],split="__"))
			#nn <- paste("-",hh[2],sep="")
			mm <- c(mm,hh[2])
		}
		
	}
	return(mm)
} 
#####################################################
pvalue <- read.csv(file="03_p_value.csv",header=T,row.names=1)
all_p<- read.csv(file="03_all_p_value.csv",header=T,row.names=1)

grpInfo = read.table(file=arg1[length(arg1)],header=T,sep="\t",check.names = FALSE, comment.char = "")
rownames(grpInfo)=as.character(grpInfo[,1])
gut <- read.csv(file="03_all_gut.csv",header=T,row.names=1)
ii=intersect(rownames(grpInfo),colnames(gut))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}
grpInfo=grpInfo[ii,]
gut=gut[,ii]



gut_col <- colnames(gut)
grp <- length(arg)

###### UP DOWN
mm <- c()
for (i in 1:grp){
	group <- unlist(strsplit(arg[i],split="-"))
	mm <- c(mm,group[1],group[2])
}


ab_cc= c(mm,arg,"Factor")

abundance <- matrix(NA,ncol=length(ab_cc),nrow=length(rownames(gut)))
colnames(abundance) <- ab_cc ### 1:UP ,0:down
rownames(abundance) <- rownames(gut)

########## 计算均值
for (i in 1:length(mm)){
	c1 <- grep(paste0('^',mm[i],'$'),as.character(grpInfo[,2]),value=F,fixed = F,ignore.case = F)
	c2=as.character(grpInfo[c1,1])
	abundance[,i] <- apply(gut[,c2],1,mean)
}


########分组判断   ##### b高为1,为up
for(i in 1:length(arg))
{
	b=2*i
	a=b-1
	c=length(mm)+i
	d1=as.numeric(abundance[,a])
	d2=as.numeric(abundance[,b])

	len_up <- which(d2 >= d1 )


	abundance[,c] <- 0
	if(length(len_up)>0){
		abundance[len_up,c] <- 1 
	}	
}

######## 判断 Factor
d=b+length(arg)
len_F <- apply(as.matrix(abundance[,(b+1):d]),1,F_Factor)
abundance[,length(abundance[1,])] <- len_F


###########	处理 种属信息
###########	处理 种属信息
info <- read.csv(file="03_info.csv",header=T,row.names=1)

tt <- apply(info,1,F_ifno)
temp <- t(tt)
temp <- cbind(temp,rownames(info))
colnames(temp) <- c("Species","Genus","Family","Phylum","ID")

d1 <- cbind(abundance,ID=rownames(abundance))


d2 <- merge(d1,temp,by.x="ID",by.y="ID",sort=F)

d_2 <- d2[,-1]
rownames(d_2) <- d2[,1]
ord2 <- order(d_2$Factor,d_2$Phylum,d_2$Family,d_2$Genus,d_2$Species,decreasing = F)
write.csv(d_2[ord2,],file="03_all_abundance.csv",row.names=T)
############排序
#########删除冲突部分
l_ct <- which(d2$Factor == 0.5)
if(length(l_ct)!=0){
	d3 <- d2[-l_ct,]
	p3 <- pvalue[-l_ct,]
	g3 <- gut[-l_ct,]
	all_p3 <- all_p[-l_ct,]
} else {
	d3 <- d2
	p3 <- pvalue
	g3 <- gut
	all_p3 <- all_p
}
ord <- order(d3$Factor,d3$Phylum,d3$Family,d3$Genus,d3$Species,decreasing = F)


d4 <- d3[,-1]
rownames(d4) <- d3[,1]


############
write.csv(all_p3[ord,],file="04_all_p.csv",row.names=T)
write.csv(p3[ord,],file="04_p_value.csv",row.names=T)
write.csv(d4[ord,],file="04_abundance.csv",row.names=T)
write.csv(g3[ord,],file="04_gut.csv",row.names=T)
write.csv(temp[,-length(temp[1,])],file="04_all_info.csv",row.name=T)

########## rda
dd=read.csv("01_all_gut.csv",header=T,row.names=1)
dd1=dd[,ii]
dd2=t(dd1)
gg=grpInfo
gr=unique(as.character(gg[,2]))
info=matrix(0,ncol=length(gr),nrow=nrow(gg))
colnames(info)=gr
rownames(info)=rownames(gg)
for(i in 1:length(gr)){
	i1=as.character(gg[,2]) %in% gr[i]
	info[i1,i]=1
}


cc=decorana(dd2)
ss=sum(cc$evals)
#根据看分析结果中Axis Lengths的第一轴的大小,即ss
#如果大于4.0,就应选CCA（基于单峰模型，典范对应分析）
#如果在3.0-4.0之间，选RDA和CCA均可
#如果小于3.0, RDA的结果会更合理（基于线性模型，冗余分析）
if(ss <= 4){
	print(paste0("SS:",ss,",Choose RDA"))
	pca<- rda(as.data.frame(dd2)~.,as.data.frame(info),scale=T) 
} else {
	print(paste0("SS:",ss,",Choose CCA"))
	pca<- cca(as.data.frame(dd2)~.,as.data.frame(info),scale=T) 
}

pdf('RDA.Raw.pdf')
plot(pca)
dev.off()
########### rda or cca

ax <- envfit(pca,info,permu=999)    #######每个环境因子显著性检验
xy_gr <- ax$vectors$arrows
pp <- ax$vectors$pvals

#perm=permutest(pca,permu=999) ###########检验环境因子相关显著性（Monte Carlo permutation test）
pp0=anova(pca)
pp=pp0$`Pr(>F)`[1]

##### PLOT
pdf(file="RDA.Genus.pdf")
d_pp=pca$CCA$v[,c(1,ncol(pca$CCA$v))]
pre <- data.frame(d_pp)  #提取样本得分
rr=rep('grey',nrow(d_pp))
pre=cbind(pre,rr,rr,rr)
colnames(pre)=c('RDA1','RDA2','Color','Size','Label')
pre[,'Label']=NA
pre[rownames(d4),'Label']=as.character(d4[,'Genus'])
pre[,'Size']=apply(dd1[rownames(pre),],1,sum)
pre[,'Size']=( (pre[,'Size']-min(pre[,'Size']))/(max(pre[,'Size'])-min(pre[,'Size'])) +1 ) *1.2
pre[,'Color']='lightgray'
p_ii=rownames(d4)[!is.na(d4[,'Genus'])]
pre[p_ii,'Color']='firebrick1'
rda1 =round(pca$CCA$eig[1]/sum(pca$CCA$eig)*100,2) #第一轴标签
rda2 =round(pca$CCA$eig[2]/sum(pca$CCA$eig)*100,2) #第二轴标签
#biplot(pca,scaling=3,xlab="RDA1",ylab="RDA2",display="sp",col="blue",xlim=c(-1,1),type="text",main='Genus') #sp:物种，si:样方，bp:环境因子
tt=which(is.na(pre[,'Label']))
par(mar=c(6,6,6,6))
plot(x=as.numeric(pre[tt,1]),y=as.numeric(pre[tt,2]),
	xlab=paste0('RDA1(',rda1,"%)"),ylab=paste0('RDA2(',rda2,"%)"),main='RDA Plot(Genus Level)',
	col=as.character(pre[tt,'Color']),cex=as.numeric(pre[tt,'Size']),pch=16)
txt <- paste0("Permutation Test:P=",pp)
mtext(text=txt,side=3,adj=0)
abline(h=0,col='gray',lwd=2,lty=2)
abline(v=0,col='gray',lwd=2,lty=2)
points(x=as.numeric(pre[-tt,1]),y=as.numeric(pre[-tt,2]),pch=16,col=as.character(pre[-tt,'Color']),cex=as.numeric(pre[-tt,'Size']),xpd=T)
for(i in 1:nrow(pre)){
	ll=pre[i,'Label']
	if(is.na(ll)){
		next
	} else {
		x1=as.numeric(pre[i,1])
		y1=as.numeric(pre[i,2])
		text(x=x1,y=y1,labels=ll,adj=0.5,col='black',cex=0.5,xpd=T)	
	}
}

min_x=min(as.numeric(pre[,1]))
max_x=max(as.numeric(pre[,1]))
min_y=min(as.numeric(pre[,2]))
max_y=max(as.numeric(pre[,2]))
for(i in 1:nrow(xy_gr)){
	x1=as.numeric(xy_gr[i,1])
	y1=as.numeric(xy_gr[i,2])
	if(x1<(min_x)*1.2) x1=par('usr')[1]-(max_x-min_x)*0.1
	if(x1>(max_x)*1.2) x1=par('usr')[2]+(max_x-min_x)*0.1
	if(y1<(min_y)*1.2) y1=par('usr')[3]-(max_y-min_y)*0.1
	if(y1>(max_y)*1.2) y1=par('usr')[4]+(max_y-min_y)*0.1
	arrows(x0=0,y0=0,x1=x1,y1=y1,code=2,col="steelblue",length=0.1,xpd=T) ##YES
	text(x=x1,y=y1,labels=paste0(gr[i],'(',round(as.numeric(xy_gr[i,1]),2),',',round(as.numeric(xy_gr[i,2]),2),')'),col="steelblue",font=2,xpd=T,cex=1.2)
}

dev.off()


png("RDA.Genus.png",height=560,width=560)
par(mar=c(6,6,6,6))
plot(x=as.numeric(pre[tt,1]),y=as.numeric(pre[tt,2]),
	xlab=paste0('RDA1(',rda1,"%)"),ylab=paste0('RDA2(',rda2,"%)"),main='RDA Plot(Genus Level)',
	col=as.character(pre[tt,'Color']),cex=as.numeric(pre[tt,'Size']),pch=16)
txt <- paste0("Permutation Test:P=",pp)
mtext(text=txt,side=3,adj=0)
abline(h=0,col='gray',lwd=2,lty=2)
abline(v=0,col='gray',lwd=2,lty=2)
points(x=as.numeric(pre[-tt,1]),y=as.numeric(pre[-tt,2]),pch=16,col=as.character(pre[-tt,'Color']),cex=as.numeric(pre[-tt,'Size']),xpd=T)
for(i in 1:nrow(pre)){
	ll=pre[i,'Label']
	if(is.na(ll)){
		next
	} else {
		x1=as.numeric(pre[i,1])
		y1=as.numeric(pre[i,2])
		text(x=x1,y=y1,labels=ll,adj=0.5,col='black',cex=0.5,xpd=T)	
	}
}

min_x=min(as.numeric(pre[,1]))
max_x=max(as.numeric(pre[,1]))
min_y=min(as.numeric(pre[,2]))
max_y=max(as.numeric(pre[,2]))
for(i in 1:nrow(xy_gr)){
	x1=as.numeric(xy_gr[i,1])
	y1=as.numeric(xy_gr[i,2])
	if(x1<(min_x)*1.2) x1=par('usr')[1]-(max_x-min_x)*0.1
	if(x1>(max_x)*1.2) x1=par('usr')[2]+(max_x-min_x)*0.1
	if(y1<(min_y)*1.2) y1=par('usr')[3]-(max_y-min_y)*0.1
	if(y1>(max_y)*1.2) y1=par('usr')[4]+(max_y-min_y)*0.1
	arrows(x0=0,y0=0,x1=x1,y1=y1,code=2,col="steelblue",length=0.1,xpd=T) ##YES
	text(x=x1,y=y1,labels=paste0(gr[i],'(',round(as.numeric(xy_gr[i,1]),2),',',round(as.numeric(xy_gr[i,2]),2),')'),col="steelblue",font=2,xpd=T,cex=1.2)
}

dev.off()



