arg <- commandArgs(T)
if(length(arg) < 2){
  cat("Argument: Input_path Out_path\n")
  quit('no')
}
pacman::p_load(tidyverse,MicrobiotaProcess)

outDir <- arg[2]
input_path <- arg[1]
otu <- paste0(input_path,"otu_table-rename.qza")
rep <- paste0(input_path,"rep-seqs-rename.qza")
tree <- paste0(input_path,"rooted-tree.qza")
tax <- paste0(input_path,"taxonomy.qza")
sample <- paste0(input_path,"sample-metadata.txt")
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeqza=tree)

ps <- ps_dada2

pcoares <- get_pcoa(obj = ps_dada2,distmethod="euclidean", method="hellinger")
pcoaplot <- ggordpoint(obj=pcoares, biplot=TRUE, speciesannot=TRUE,
                       pc=c(1,2),factorNames=c("Group"), ellipse=TRUE) + 
  scale_color_manual(values=c("#2874C5", "#EABF00","#2874C5", "#EABF00","#2874C5", "#EABF00"))
pcoaplot
png(paste0(outDir,"bray_curtis_Total.pcoa.circle.png"),width=4*480/3,height=4*480/4, bg="transparent")
pcoaplot
dev.off()

pdf(paste0(outDir,"bray_curtis_Total.pcoa.circle.pdf"),width = 8,height = 6)
pcoaplot
dev.off()

pcares <- get_pca(obj=ps_dada2, method="hellinger")
pcaplot <- ggordpoint(obj=pcares, biplot=TRUE, speciesannot=TRUE,
                      pc=c(1,2),factorNames=c("Group"), ellipse=TRUE) + 
  scale_color_manual(values=c("#2874C5", "#EABF00","#2874C5", "#EABF00","#2874C5", "#EABF00"))

png(paste0(outDir,"bray_curtis_Total.pca.circle.png"),width=4*480/3,height=4*480/4, bg="transparent")
pcaplot
dev.off()

pdf(paste0(outDir,"bray_curtis_Total.pca.circle.pdf"),width = 8,height = 6)
pcaplot
dev.off()
