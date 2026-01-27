arg1 <- commandArgs(T)
if(length(arg1) < 5){
	cat("Argument: Data_File TaxFile Out_Dir Group  Color\n")
	quit('no')
}
print(arg1)
###########	测试
#arg1 <- c("Total.gut.sub1W.filter_4_noChimera_rmAlignFail.otu_table.txt"  , "/home/wmj/视频/16s/qiime1/TT/Analysis/15_RDA_and_Heatmap/Total","T1Control-T1Classical" ,"T1Control-T1T1H" ,"T1Control-T1T1L", "/home/wmj/视频/16s/qiime1/TT/bin/Total-color.txt" )
#arg1 <- c("Total.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.txt","/Volumes/Flower0501/16s/Qiime1/Analysis/4_common_OTUs/Total/","Con-Treat","/Volumes/Flower0501/16s/Qiime1/Total-color.txt")
#arg1 <- c("/Volumes/Flower0501/2022_16s/CY/0_input_data/Reotu_table.tsv","/Volumes/Flower0501/2022_16s/CY/0_input_data/taxonomy.tsv","/Volumes/Flower0501/2022_16s/CY/15_RDA_and_Heatmap/","Z1-Z2","/Volumes/Flower0501/2022_16s/CY/Total-color.txt")
# 计算p值

setwd(arg1[3])
arg <- arg1[4:(length(arg1)-1)]
grp <- length(arg)

#########	P值
library(vegan) 
####################自定义函数区######################
### Tukey's HSD
F_Tukey <- function(data) {
	ss <- names(data)
	model <- aov(data ~ ss)
	posthoc <- TukeyHSD(model)
	return(posthoc$ss[4])
} 


#####################################################


#####数据读取

gut <- read.table(file=arg1[1],header=T,row.names=1,sep="\t",check.names = FALSE)
tax <- read.table(file=arg1[2],header=T,sep="\t",check.names = FALSE)
gut$tax <-lapply(rownames(gut),function(x){tax[tax$`Feature ID` == x,"Taxon"]})
#gut <- readr::read_table2(arg1[1],col_names = TRUE)
grpInfo = read.table(file=arg1[length(arg1)],header=T,sep="\t",check.names = FALSE, comment.char = "")
rownames(grpInfo)=as.character(grpInfo[,1])
ii=intersect(rownames(grpInfo),colnames(gut))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}
grpInfo=grpInfo[ii,]
gut1=gut[,ii]
gut=cbind(gut1,unlist(gut$tax))
colnames(gut)[ncol(gut)] <- "tax"
## 数据预处理,保留含10%数据的OTU
cc.na <-c()
for(i in 1:nrow(gut)){
	dd <- gut[i,-ncol(gut)]
	cc.na <- c(cc.na,length(which(dd!=0)))

}
n <- floor((ncol(gut)-1)*0.1)
gut <- gut[which(cc.na >n),]


gut_row <- rownames(gut)
otu <- sub(pattern="ASV", replacement="OTU",gut_row)  ###替换
rownames(gut) <- otu
gut_col <- colnames(gut)
#p_value <- vector('list',grp)
#names(p_value) <- arg
p_value <- matrix(NA,ncol=grp,nrow=length(gut_row))
print(paste0("pvalue Dim:",dim(p_value)))
rownames(p_value) <- rownames(gut)
colnames(p_value) <- arg

write.csv(gut,file="01_all_gut.csv",row.names=T)
write.csv(gut[,-length(gut[1,])],file="02_all_gut.csv",row.names=T)

info <- gut[,length(gut[1,])]
names(info) <- otu
write.csv(info,file="03_info.csv",row.names=T)

#############数据重组并Tukey检验
for (i in 1:grp){
	group <- unlist(strsplit(arg[i],split="-"))
	m1 <- paste("^",group[1],"$",sep="")
	c1 <- grep(m1,as.character(grpInfo[,2]),value=F,fixed = F,ignore.case = F)
	c11=as.character(grpInfo[c1,1])
	m2 <- paste("^",group[2],"$",sep="")
	c2 <- grep(m2,as.character(grpInfo[,2]),value=F,fixed = F,ignore.case = F)
	c22=as.character(grpInfo[c2,1])
	data1 <- gut[,c11]
	data2 <- gut[,c22]
	D <- rep(group[1],length(c11))
	H <- rep(group[2],length(c22))
	temp <- as.matrix(cbind(data1,data2))
	colnames(temp) <- c(D,H)	
	p_value[,i] <- apply(temp,1,F_Tukey)
}


write.csv(p_value,file="02_p_value.csv",row.names=T)




## abundance
dd <- gut[,-length(gut[1,])]
cc <- c()
for(i in 1:length(arg)){
	g1 <- unlist(strsplit(arg[i],"-"))
	cc<-c(cc,g1)
}

c1 <- unique(cc)
re <- matrix(NA,ncol=length(c1),nrow=nrow(dd))
colnames(re)<-c1
rownames(re)<-rownames(dd)
for(i in 1:length(c1)){
	g1 <- grep(paste0("^",c1[i],"$"),as.character(grpInfo[,2]))
	g2=as.character(grpInfo[g1,1])
	d1 <- dd[,g2]
	d2 <-apply(d1,1,mean)
	re[,i]<-d2
}

write.csv(re,"01_all_abundance.csv",row.names=T)


