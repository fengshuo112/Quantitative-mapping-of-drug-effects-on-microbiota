arg <- commandArgs(T)
if(length(arg) != 3){
        cat("Argument: Data_File Group_File Out_Dir\n")
        quit('no')
}
print(arg)
#arg <- c("/home/Project/16S/LSZ_20190908/Analysis/2_Diversity/Total.Genus.abundance.txt", "/home/Project/16S/LSZ_20190908/bin/Total-color.txt", "/home/Project/16S/LSZ_20190908/Analysis/2_Diversity" )
dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
setwd(outDir)
grp_col <- 2;

SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1] ####### NOTE



library(vegan);
library(MASS);


data <- read.table(dataFile,header=T,row.names=1,check.names = FALSE);
grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo) <- grpInfo[,1]
grp <- basename(outDir)
#nn <- grep(paste("^",grp,".*",sep=""),colnames(grpInfo))
#nn1 <- colnames(grpInfo)[nn]
#nn2 <- strsplit(nn1,".",fixed=T)[[1]][2]
#nn3 <- grep(paste("^",nn2,".*",sep=""),colnames(grpInfo))
#ii <- intersect(colnames(data),grpInfo[,1])
ii <- intersect(rownames(grpInfo),colnames(data))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),colnames(data))
	print(ii1)
}

grpInfo <- grpInfo[ii,]
data <- data[,ii]
Abundance <- t(data);
Alpha <- as.matrix(diversity(Abundance,index="shannon",MARGIN=1,base=exp(1)));

groupname <- c();
for(i in 1:length(rownames(Alpha)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste0('^',names(Alpha[i,]),'$'),grpInfo[,1]),grp_col]));
}

rownames(Alpha) <- groupname;
Groups <- as.character(levels(as.factor(grpInfo[,grp_col])))
Group <- c();
Group.mean  <- c();
Group.sd  <- c();
for(i in 1:length(Groups))
{
	Group[[i]] <- as.matrix(Alpha[grep(paste("^",Groups[i],"$",sep=""),rownames(Alpha)),]);
	Group.mean[i] <- mean (Group[[i]]);
	Group.sd[i] <- sd (Group[[i]]);
}

GG.Utest <- c();
GG.Asterisk <- c();
for(m in 1:(length(Groups)-1))
{
	for(n in (m+1):length(Groups))
	{
		GG=paste(Groups[m],"-",Groups[n],sep="");
		GG.Utest[[GG]]=wilcox.test(Group[[m]],Group[[n]]);
		if(GG.Utest[[GG]]$p.value < 0.001){
			GG.Asterisk[[GG]] <- '***' ;
		}
		else if(GG.Utest[[GG]]$p.value < 0.01){
                        GG.Asterisk[[GG]] <- ' **' ;
		}
		else if(GG.Utest[[GG]]$p.value < 0.05){
                        GG.Asterisk[[GG]] <- '  *' ;
		}
		else{
			GG.Asterisk[[GG]] <- '   ' ;
		}
	}
}


col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
pdf(paste(SamID,".barplot.pdf",sep=""),width=4,height=4)

y.lim <- max(Group.mean+Group.sd)*1.3;
err.wid <- 0.1
sig.tick <- 0.03
gap <- 0.008
step <- 0.03

bar.pos <- barplot(Group.mean,border=NA,col=rainbow(length(Groups)),names.arg=Groups,ylim=c(0,y.lim),ylab='Volume',cex.names=1.0,cex.lab=1.0,cex=1.2)
rect(par('usr')[1],0,par('usr')[2],y.lim,border=NA)
segments(bar.pos,Group.mean-Group.sd,bar.pos,Group.mean+Group.sd)
segments(c(bar.pos-err.wid,bar.pos-err.wid),c(Group.mean-Group.sd,Group.mean+Group.sd),c(bar.pos+err.wid,bar.pos+err.wid),c(Group.mean-Group.sd,Group.mean+Group.sd),lwd=0.6)

for(m in 1:(length(Groups)-1))
{
        for(n in (m+1):length(Groups))
       {
		GG=paste(Groups[m],"-",Groups[n],sep="");
		if(GG.Utest[[GG]]$p.value < 0.05)
		{
			lines(rep(c(bar.pos[m],bar.pos[n]-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
			text(median(c(bar.pos[m],bar.pos[n]-gap)),y.lim*(1-step*(n+m)+0.02),labels=GG.Asterisk[[GG]]);
		}
	}
}
		

dev.off()
