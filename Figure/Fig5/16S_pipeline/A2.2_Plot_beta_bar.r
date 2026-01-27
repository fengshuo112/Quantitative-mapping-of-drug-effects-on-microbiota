arg <- commandArgs(T)
if(length(arg) != 3){
        cat("Argument: Data_File Group_File Out_Dir\n")
        quit('no')
}
print(arg)
#arg <-c("/home/wmj/Metagenome-2/16s/Work/WXD/WXD-3/Analysis/2_Diversity/Total/Total.All.Genus.abundance.txt", "/home/wmj/Metagenome-2/16s/Work/WXD/WXD-3/bin/gut-all.info" ,"/home/wmj/Metagenome-2/16s/Work/WXD/WXD-3/Analysis/2_Diversity/Total"    )
dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
setwd(outDir)
grp_col <- 3;

SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1] ####### NOTE


library(vegan);
library(MASS);


data <- read.table(dataFile,header=T,row.names=1,check.names = FALSE);
grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo) <- grpInfo[,1]
#grp <- basename(outDir)
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

groupname <- c();
for(i in 1:length(rownames(Abundance)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste0('^',rownames(Abundance)[i],'$'),grpInfo[,1]),grp_col]));
}

rownames(Abundance) <- groupname;
Groups <- as.character(levels(as.factor(grpInfo[,grp_col])))
Abundance.Group <- c();
Beta.Group <- c();
Beta.Group.mean  <- c();
Beta.Group.sd  <- c();

for(i in 1:length(Groups))
{
        Abundance.Group[[i]] <- as.matrix(Abundance[grep(paste("^",Groups[i],"$",sep=""),rownames(Abundance)),]);
	Beta.Group[[i]] <- vegdist(Abundance.Group[[i]],method="bray");
	Beta.Group.mean[i] <- mean(Beta.Group[[i]]);
	Beta.Group.sd[i] <- sd(Beta.Group[[i]]);
}

Beta.GG.Utest <- c();
Beta.GG.Asterisk <- c();
for(m in 1:(length(Groups)-1))
{
	for(n in (m+1):length(Groups))
	{
		GG=paste(Groups[m],"-",Groups[n],sep="");
		Beta.GG.Utest[[GG]]=wilcox.test(Beta.Group[[m]],Beta.Group[[n]]);
		if(Beta.GG.Utest[[GG]]$p.value < 0.001){
			Beta.GG.Asterisk[[GG]] <- '***' ;
		}
		else if(Beta.GG.Utest[[GG]]$p.value < 0.01){
                        Beta.GG.Asterisk[[GG]] <- '**' ;
		}
		else if(Beta.GG.Utest[[GG]]$p.value < 0.05){
			Beta.GG.Asterisk[[GG]] <- '*' ;
		}
		else{
			Beta.GG.Asterisk[[GG]] <- '   ' ;
		}
	}
}



pdf(paste(SamID,".beta.barplot.pdf",sep=""),width=4,height=4)

y.lim <- max(Beta.Group.mean+Beta.Group.sd)*1.3;
err.wid <- 0.1
sig.tick <- 0.03
gap <- 0.008
step <- 0.03

bar.pos <- barplot(Beta.Group.mean,border=NA,col=rainbow(length(Groups)),names.arg=Groups,ylim=c(0,y.lim),ylab='Sorensen index',cex.names=1.0,cex.lab=1.0,cex=1.2)
rect(par('usr')[1],0,par('usr')[2],y.lim,border=NA)
segments(bar.pos,Beta.Group.mean-Beta.Group.sd,bar.pos,Beta.Group.mean+Beta.Group.sd)
segments(c(bar.pos-err.wid,bar.pos-err.wid),c(Beta.Group.mean-Beta.Group.sd,Beta.Group.mean+Beta.Group.sd),c(bar.pos+err.wid,bar.pos+err.wid),c(Beta.Group.mean-Beta.Group.sd,Beta.Group.mean+Beta.Group.sd),lwd=0.6)

for(m in 1:(length(Groups)-1))
{
        for(n in (m+1):length(Groups))
       {
		GG=paste(Groups[m],"-",Groups[n],sep="");
		if(Beta.GG.Utest[[GG]]$p.value < 0.05)
		{
			lines(rep(c(bar.pos[m],bar.pos[n]-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
			text(median(c(bar.pos[m],bar.pos[n]-gap)),y.lim*(1-step*(n+m)+0.02),labels=Beta.GG.Asterisk[[GG]]);
		}
	}
}


dev.off()
