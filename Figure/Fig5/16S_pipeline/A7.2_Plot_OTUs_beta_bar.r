arg <- commandArgs(T)
if(length(arg) < 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column color_column\n")
        quit('no')
}
print(arg)
#arg <- c("/Volumes/Flower0501/2022_16s/Mala/Group2/Analysis/8_OTUs_beta/bray-curtis_distance_matrix.tsv","/Volumes/Flower0501/2022_16s/Mala/Group2/Analysis/0_input_data/sample-metadata.txt","/Volumes/Flower0501/2022_16s/Mala/Group2/Analysis/8_OTUs_beta","4",5)
dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);
color_col <- as.numeric(arg[5])

grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo) <- as.character(grpInfo[,1])
FileName <- strsplit(basename(dataFile),'_',fixed=T)[[1]];
SamID <- "Total";
MetID <- FileName[1];
GID <- colnames(grpInfo)[grp_col];
Len <- ncol(grpInfo)-1


Beta <- as.matrix(data.frame(read.table(dataFile,header=T,row.names = 1,check.names = FALSE)));

ii <- intersect(rownames(grpInfo),rownames(Beta))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}





Beta <- as.matrix(Beta[ii,ii])
grpInfo <- grpInfo[ii,]



groupname <- c();
for(i in 1:length(rownames(Beta)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste("^",rownames(Beta)[i],"$",sep=""),grpInfo[,1]),grp_col]));
}

rownames(Beta) <- groupname;
colnames(Beta) <- groupname;

Groups <- as.character(unique(grpInfo[,grp_col]))
G_cc=as.character(unique(grpInfo[,color_col]))
names(G_cc)=Groups
#gg <- read.table(arg[5],header=T)
#rownames(gg) <- as.character(gg[,1])
#ii <- intersect(rownames(grpInfo),rownames(gg))
#gg <- gg[ii,]
gg <- grpInfo
Groups.1 <-as.character(unique(gg[,2])) 


Beta.Group <- c();
Beta.Group.mean  <- c();
Beta.Group.sd  <- c();

for(i in 1:length(Groups))
{
	Beta.Group[[i]] <- Beta[grep(paste("^",Groups[i],"$",sep=""),rownames(Beta)),grep(paste("^",Groups[i],"$",sep=""),colnames(Beta))];
	Beta.Group.mean[[i]] <- mean(Beta.Group[[i]]);
	Beta.Group.sd[[i]] <- sd(Beta.Group[[i]]);
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
                        Beta.GG.Asterisk[[GG]] <- ' **' ;
		}
		else if(Beta.GG.Utest[[GG]]$p.value < 0.05){
			Beta.GG.Asterisk[[GG]] <- '  *' ;
		}
		else{
			Beta.GG.Asterisk[[GG]] <- '   ' ;
		}
	}
}

pdf(paste(outDir,"/",SamID,".",GID,".",MetID,".resample.otu.barplot.pdf",sep=""),width=2*Len,height=2*Len)
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=G_cc
y.lim <- max(mapply("+",Beta.Group.mean,Beta.Group.sd))*1.3;
err.wid <- 0.1
sig.tick <- 0.03
gap <- 0.008
step <- 0.03
par(mgp=c(4,1,0),tck=-0.02,mar=c(5,5,1,1),cex.lab=1.2,cex.axis=1.2)
bar.pos <- barplot(as.vector(unlist(Beta.Group.mean)),border=NA,col=col.line[Groups],ylim=c(0,y.lim),ylab="",las=1,space=0.5,xaxt="n")
mtext(at=bar.pos,side=1,text=Groups,line=0.5,cex=1,adj=0.5)
mtext(side=2,text='Sorensen index',line=3.5,cex=1.8,adj=0.5)
#mtext(side=1,text=Groups.1,at=bar.pos,line=0.5,cex=1.5)
rect(par('usr')[1],0,par('usr')[2],y.lim,border=NA)
segments(bar.pos,mapply("-",Beta.Group.mean,Beta.Group.sd),bar.pos,mapply("+",Beta.Group.mean,Beta.Group.sd))
segments(c(bar.pos-err.wid,bar.pos-err.wid),c(mapply("-",Beta.Group.mean,Beta.Group.sd),mapply("+",Beta.Group.mean,Beta.Group.sd)),c(bar.pos+err.wid,bar.pos+err.wid),c(mapply("-",Beta.Group.mean,Beta.Group.sd),mapply("+",Beta.Group.mean,Beta.Group.sd)),lwd=0.6)

for(m in 1:(length(Groups)-1))
{
        for(n in (m+1):length(Groups))
       {
		GG=paste(Groups[m],"-",Groups[n],sep="");
		if(Beta.GG.Utest[[GG]]$p.value < 0.05)
		{
			par(lty=1,lwd=2,cex=2)
			lines(rep(c(bar.pos[m],bar.pos[n]-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
			text(median(c(bar.pos[m],bar.pos[n]-gap)),y.lim*(1-step*(n+m)+0.02),labels=Beta.GG.Asterisk[[GG]]);
		}
	}
}


dev.off()

png(paste(outDir,"/",SamID,".",GID,".",MetID,".resample.otu.barplot.png",sep=""),width=2*Len*480/6,height=2*Len*480/6, bg="transparent")
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=G_cc
y.lim <- max(mapply("+",Beta.Group.mean,Beta.Group.sd))*1.3;
err.wid <- 0.1
sig.tick <- 0.03
gap <- 0.008
step <- 0.03
par(mgp=c(4,1,0),tck=-0.02,mar=c(5,5,1,1),cex.lab=1.2,cex.axis=1.2)
bar.pos <- barplot(as.vector(unlist(Beta.Group.mean)),border=NA,col=col.line[Groups],ylim=c(0,y.lim),ylab="",las=1,space=0.5,xaxt="n")
mtext(at=bar.pos,side=1,text=Groups,line=0.5,cex=1,adj=0.5)
mtext(side=2,text='Sorensen index',line=3.5,cex=1.8,adj=0.5)
#mtext(side=1,text=Groups.1,at=bar.pos,line=0.5,cex=1.5)
rect(par('usr')[1],0,par('usr')[2],y.lim,border=NA)
segments(bar.pos,mapply("-",Beta.Group.mean,Beta.Group.sd),bar.pos,mapply("+",Beta.Group.mean,Beta.Group.sd))
segments(c(bar.pos-err.wid,bar.pos-err.wid),c(mapply("-",Beta.Group.mean,Beta.Group.sd),mapply("+",Beta.Group.mean,Beta.Group.sd)),c(bar.pos+err.wid,bar.pos+err.wid),c(mapply("-",Beta.Group.mean,Beta.Group.sd),mapply("+",Beta.Group.mean,Beta.Group.sd)),lwd=0.6)

for(m in 1:(length(Groups)-1))
{
  for(n in (m+1):length(Groups))
  {
    GG=paste(Groups[m],"-",Groups[n],sep="");
    if(Beta.GG.Utest[[GG]]$p.value < 0.05)
    {
      par(lty=1,lwd=2,cex=2)
      lines(rep(c(bar.pos[m],bar.pos[n]-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
      text(median(c(bar.pos[m],bar.pos[n]-gap)),y.lim*(1-step*(n+m)+0.02),labels=Beta.GG.Asterisk[[GG]]);
    }
  }
}


dev.off()


#############boxplot

pdf(paste(outDir,"/",SamID,".",GID,".",MetID,".resample.otu.boxplot.pdf",sep=""),width=2*Len,height=2*Len)

par(mgp=c(4,1,0),tck=-0.02,mar=c(2,6,2,2),cex.lab=1.2,cex.axis=1.2,bty='l')

xx=boxplot(Beta~groupname,border=col.line[Groups],ylab="",las=1,outline=F,xaxt="n")
mtext(side=2,text='Sorensen index',line=3.5,cex=1.8,adj=0.5)
xx=Groups
#mtext(side=1,text=Groups.1,at=xx,line=0.5,cex=1.5)
text(x=1:length(Groups),y=par("usr")[3],labels=Groups,xpd=T,srt=45,pos=2,font=3,cex=1,offset=0)

y.lim=par('usr')[4]
for(m in 1:(length(Groups)-1))
{
        for(n in (m+1):length(Groups))
       {
		GG=paste(Groups[m],"-",Groups[n],sep="");
		if(Beta.GG.Utest[[GG]]$p.value < 0.05)
		{
			par(lty=1,lwd=2,cex=2)
			lines(rep(c(m,n-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
			text(median(c(m,n-gap)),y.lim*(1-step*(n+m)+0.02),labels=Beta.GG.Asterisk[[GG]]);
		}
	}
}




dev.off()

