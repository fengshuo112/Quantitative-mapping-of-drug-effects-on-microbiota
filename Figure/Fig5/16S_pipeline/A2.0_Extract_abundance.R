arg <- commandArgs(T)
if(length(arg) < 3){
  cat("Argument: Input_path taxlevel Out_path\n")
  quit('no')
}
print(arg)
pacman::p_load(tidyverse,MicrobiotaProcess,magrittr)

outDir=arg[3]
input_path <- arg[1]
otu <- paste0(input_path,"otu_table-rename.qza")
rep <- paste0(input_path,"rep-seqs-rename.qza")
tree <- paste0(input_path,"rooted-tree.qza")
tax <- paste0(input_path,"taxonomy.qza")
sample <- paste0(input_path,"sample-metadata.txt")
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeqza=tree)


# 数据的分布
ps <- ps_dada2
if(arg[2]==1){
  tax="Phylum"
  phytax <- get_taxadf(obj=ps_dada2, taxlevel=2)
  abundance <- phytax@otu_table %>% as.data.frame()
  
}else{
  tax="Genus"
  gentax <- get_taxadf(obj=ps_dada2, taxlevel=6)
  abundance <- gentax@otu_table %>% as.data.frame()
}


norm = as.data.frame(t(t(abundance)/colSums(abundance,na=T))) %>% round(digits = 5)
norm$ID <- rownames(norm)
norm <- select(norm,ID,everything())


write.table(norm,paste0(outDir,"Total.",tax,".abundance.txt"),row.names = FALSE,quote = FALSE)