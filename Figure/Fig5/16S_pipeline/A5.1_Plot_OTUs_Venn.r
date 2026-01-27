library('Vennerable')

arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <-c("/home/Metagenome-2/16S/ZX-2/Analysis/5_Sharing_OTUs/Total/Total.gut.otu.counts.txt","/home/Metagenome-2/16S/ZX-2/Analysis/5_Sharing_OTUs/Total/Total.group.txt" ,"/home/Metagenome-2/16S/ZX-2/Analysis/5_Sharing_OTUs/Total","2"    )

dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1];


data <- read.table(dataFile,header=T,row.names=1,check.names = FALSE);
grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo)=as.character(grpInfo[,1])
ii <- intersect(rownames(grpInfo),colnames(data))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}
data=data[,ii]
grpInfo=grpInfo[ii,]

GID <- colnames(grpInfo)[grp_col];

Counts <- data.frame(data[,1:ncol(data)],check.names = FALSE);
rownames(Counts) <- rownames(data);

groupname <- c();
for(i in 1:length(colnames(Counts)))
{
        groupname <- append(groupname,as.character(grpInfo[grep(paste0('^',colnames(Counts)[i],'$'),grpInfo[,1]),grp_col]));
}

colnames(Counts) <- groupname;
Groups <- unique(as.character(grpInfo[,grp_col]))
Counts.Group <- c();
Type.Counts.Group <- c();
OTUs.Exists.Group <- c();
for(i in 1:length(Groups))
{
        Counts.Group[[i]] <- as.matrix(Counts[,grep(paste("^",Groups[i],"$",sep=""),colnames(Counts))]);
        rownames(Counts.Group[[i]]) <- rownames(Counts);
	Type.Counts.Group[[i]] <- as.matrix(data.frame(Counts.Group[[i]][,1:ncol(Counts.Group[[i]])],sum=rowSums(Counts.Group[[i]][,1:ncol(Counts.Group[[i]])])));

	Type.Counts.Group[[i]] <- Type.Counts.Group[[i]][,-(1:(ncol(Type.Counts.Group[[i]])-1))];
	Type.Counts.Group[[i]][Type.Counts.Group[[i]]==0] <- NA;
	OTUs.Exists.Group[[i]]  <- names(na.omit(Type.Counts.Group[[i]]));
}

names(OTUs.Exists.Group) <- Groups;

pdf(paste(outDir,"/",SamID,".",GID,".venn.pdf",sep=""),w=5,h=5);
data <- Venn(OTUs.Exists.Group);
plot(data,doWeight=T);
