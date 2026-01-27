arg <- commandArgs(T)
if(length(arg) < 5){
  cat("Argument: Input_path Group_File Out_path Group_col Color_col\n")
  quit('no')
}
pacman::p_load(tidyverse,MicrobiotaProcess)
#arg <- c("/Volumes/Flower0501/2022_16s/Mala/Group2/Analysis/0_input_data/","/Volumes/Flower0501/2022_16s/Mala/Group2/Analysis/0_input_data/sample-metadata.txt","/Volumes/Flower0501/2022_16s/Mala/Group2/Analysis/7_OTUs_alpha/",4,5)
input_path <- arg[1]
otu <- paste0(input_path,"otu_table-rename.qza")
rep <- paste0(input_path,"rep-seqs-rename.qza")
tree <- paste0(input_path,"rooted-tree.qza")
tax <- paste0(input_path,"taxonomy.qza")
sample <- paste0(input_path,"sample-metadata.txt")
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeqza=tree)

ps <- ps_dada2
Alpha <- get_alphaindex(ps) %>% as.data.frame()

groupFile <- arg[2]
SamID <- "Total"
MetIDs <- c("Observe", "Chao1", "Shannon", "Simpson")
outDir <- arg[3]

grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo)<-as.character(grpInfo[,1])
# FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
# SamID <-FileName[1];
# MetID <- FileName[4];


#Alpha <- as.matrix(data.frame(read.table(dataFile,header=T,row.names = 1,check.names = FALSE)));
#MetID <- "Observe"
for(MetID in MetIDs){
  print(MetID)
  ii <- intersect(rownames(grpInfo),rownames(Alpha))
  if(length(ii)!=nrow(grpInfo)){
    print('Data Lost')
    ii1=setdiff(rownames(grpInfo),ii)
    print(ii1)
  }
  
  
  Alphaj <- as.matrix(Alpha[ii,MetID])
  grpInfo <- grpInfo[ii,]
  grp_col <- as.numeric(arg[4]);
  groupname <- c();
  for(i in 1:length(rownames(Alpha)))
  {
    groupname <- append(groupname,as.character(grpInfo[grep(paste("^",rownames(Alpha)[i],"$",sep=""),as.character(grpInfo[,1])),grp_col]));
  }
  
  rownames(Alphaj) <- groupname;
  #Groups <- as.character(levels(as.factor(grpInfo[,grp_col])))
  Groups <- as.character(unique(grpInfo[,grp_col]))
  color_col <- as.numeric(arg[5])
  G_cc=as.character(unique(grpInfo[,color_col]))
  names(G_cc)=Groups
  Alpha.Group <- c();
  Alpha.Group.mean  <- c();
  Alpha.Group.sd  <- c();
  for(i in 1:length(Groups))
  {
    Alpha.Group[[i]] <- as.matrix(Alphaj[grep(paste("^",Groups[i],"$",sep=""),rownames(Alphaj)),]);
    Alpha.Group.mean[[i]] <- mean (Alpha.Group[[i]]);
    Alpha.Group.sd[[i]] <- sd (Alpha.Group[[i]]);
  }
  
  
  Alpha.GG.Utest <- c();
  Alpha.GG.Asterisk <- c();
  for(m in 1:(length(Groups)-1))
  {
    for(n in (m+1):length(Groups))
    {
      GG=paste(Groups[m],"-",Groups[n],sep="");
      Alpha.GG.Utest[[GG]]=wilcox.test(Alpha.Group[[m]],Alpha.Group[[n]],exact = FALSE);
      if(Alpha.GG.Utest[[GG]]$p.value < 0.001){
        Alpha.GG.Asterisk[[GG]] <- '***' ;
      }
      else if(Alpha.GG.Utest[[GG]]$p.value < 0.01){
        Alpha.GG.Asterisk[[GG]] <- ' **' ;
      }
      else if(Alpha.GG.Utest[[GG]]$p.value < 0.05){
        Alpha.GG.Asterisk[[GG]] <- '  *' ;
      }
      else{
        Alpha.GG.Asterisk[[GG]] <- '   ' ;
      }
    }
  }
  
  png(paste(outDir,"/",SamID,".",MetID,".final.otu.barplot.png",sep=""),width=4*480/7,height=4*480/7, bg="transparent")
  #col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
  col.line=G_cc
  y.lim <- max(mapply(`+`,Alpha.Group.mean,Alpha.Group.sd))*1.3;
  err.wid <- 0.1
  sig.tick <- 0.03
  gap <- 0.008
  
  
  step <- 0.03
  
  bar.pos <- barplot( as.vector(unlist(Alpha.Group.mean)),border=NA,col=col.line[Groups],names.arg="",ylim=c(0,y.lim),ylab=paste(MetID,"Index",sep=" "),cex.names=1.0,cex.lab=1.0,cex=1.2)
  #text(x=seq((num_groups-1),(num_groups*top-1),num_groups),y=par("usr")[3]-0.5,labels=Top.name,xpd=T,srt=45,pos=2,font=3,cex=1.1,offset=0)
  text(x=bar.pos,y=par("usr")[3]-0.5,labels=Groups,xpd=T,srt=45,pos=2,font=3,cex=1,offset=0)
  #mtext(at=bar.pos,text=Groups,srt=-45)
  rect(par('usr')[1],0,par('usr')[2],y.lim,border=NA)
  segments(bar.pos,mapply("-",Alpha.Group.mean,Alpha.Group.sd),bar.pos,mapply("+",Alpha.Group.mean,Alpha.Group.sd))
  segments(c(bar.pos-err.wid,bar.pos-err.wid),c(mapply("-",Alpha.Group.mean,Alpha.Group.sd),mapply("+",Alpha.Group.mean,Alpha.Group.sd)),c(bar.pos+err.wid,bar.pos+err.wid),c(mapply("-",Alpha.Group.mean,Alpha.Group.sd),mapply("+",Alpha.Group.mean,Alpha.Group.sd)),lwd=0.6)
  
  for(m in 1:(length(Groups)-1))
  {
    for(n in (m+1):length(Groups))
    {
      GG=paste(Groups[m],"-",Groups[n],sep="");
      if(Alpha.GG.Utest[[GG]]$p.value < 0.05)
      {
        lines(rep(c(bar.pos[m],bar.pos[n]-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
        text(median(c(bar.pos[m],bar.pos[n]-gap)),y.lim*(1-step*(n+m)+0.02),labels=Alpha.GG.Asterisk[[GG]]);
      }
    }
  }
  
  dev.off()
  
  
  
  pdf(paste(outDir,"/",SamID,".",MetID,".final.otu.barplot.pdf",sep=""),width=4,height=4)
  #col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
  col.line=G_cc
  y.lim <- max(mapply(`+`,Alpha.Group.mean,Alpha.Group.sd))*1.3;
  err.wid <- 0.1
  sig.tick <- 0.03
  gap <- 0.008
  
  
  step <- 0.03
  
  bar.pos <- barplot( as.vector(unlist(Alpha.Group.mean)),border=NA,col=col.line[Groups],names.arg="",ylim=c(0,y.lim),ylab=paste(MetID,"Index",sep=" "),cex.names=1.0,cex.lab=1.0,cex=1.2)
  #text(x=seq((num_groups-1),(num_groups*top-1),num_groups),y=par("usr")[3]-0.5,labels=Top.name,xpd=T,srt=45,pos=2,font=3,cex=1.1,offset=0)
  text(x=bar.pos,y=par("usr")[3]-0.5,labels=Groups,xpd=T,srt=45,pos=2,font=3,cex=1,offset=0)
  #mtext(at=bar.pos,text=Groups,srt=-45)
  rect(par('usr')[1],0,par('usr')[2],y.lim,border=NA)
  segments(bar.pos,mapply("-",Alpha.Group.mean,Alpha.Group.sd),bar.pos,mapply("+",Alpha.Group.mean,Alpha.Group.sd))
  segments(c(bar.pos-err.wid,bar.pos-err.wid),c(mapply("-",Alpha.Group.mean,Alpha.Group.sd),mapply("+",Alpha.Group.mean,Alpha.Group.sd)),c(bar.pos+err.wid,bar.pos+err.wid),c(mapply("-",Alpha.Group.mean,Alpha.Group.sd),mapply("+",Alpha.Group.mean,Alpha.Group.sd)),lwd=0.6)
  
  for(m in 1:(length(Groups)-1))
  {
    for(n in (m+1):length(Groups))
    {
      GG=paste(Groups[m],"-",Groups[n],sep="");
      if(Alpha.GG.Utest[[GG]]$p.value < 0.05)
      {
        lines(rep(c(bar.pos[m],bar.pos[n]-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
        text(median(c(bar.pos[m],bar.pos[n]-gap)),y.lim*(1-step*(n+m)+0.02),labels=Alpha.GG.Asterisk[[GG]]);
      }
    }
  }
  
  
  dev.off()
  
  
}
