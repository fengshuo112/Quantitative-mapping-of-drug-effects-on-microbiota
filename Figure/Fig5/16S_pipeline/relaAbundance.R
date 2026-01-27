setwd('E:/Lab_Work/FS/16S/FS1108/')

library(tidyverse)
library(MicrobiotaProcess)
library(showtext)
library(patchwork)

# 16s
rename <- function(df,label){
  df <- df[!grepl("__un_",rownames(df)),]
  df$name <- unlist(lapply(rownames(df),function(x) str_match_all(x,paste0(label,"__([A-Za-z0-9\\[\\]-]+)"))[[1]][2]))
  dat_final <- aggregate(df[,c(1:(ncol(df)-1))],by=list(name=df$name),sum)
  dat_final <- column_to_rownames(dat_final,"name")
  return(dat_final)
}

# silva
input_path <-'Analysis/0_input_data/'
otu <- paste0(input_path,"otu_table-rename.qza")
tax <- paste0(input_path,"taxonomy_silva.qza")
sample <- paste0(input_path,"sample-metadata.txt")
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax, mapfilename=sample)

# greengene
tax <- paste0(input_path,"taxonomy_gg.qza")
ps_dada2 <- import_qiime2(otuqza = otu,taxaqza = tax,mapfilename = sample)

# genus,silva
leveltax <- get_taxadf(obj=ps_dada2, taxlevel=6)
df_genus <- leveltax@otu_table %>% as.data.frame()
df_genus <- rename(df_genus,"g__g")
df_genus_re <- apply(df_genus, 2, function(x) x/sum(x)) %>% as.data.frame()

write.csv(df_genus, 'Analysis/Abundance_silva/Genus_reads.csv')
write.csv(df_genus_re, 'Analysis/Abundance_silva/Genus_relaAbundance.csv')

# family,silva
leveltax <- get_taxadf(obj=ps_dada2, taxlevel=5)
df_family <- leveltax@otu_table %>% as.data.frame()
df_family <- rename(df_family,"f__f")
df_family_re <- apply(df_family, 2, function(x) x/sum(x)) %>% as.data.frame()
write.csv(df_family, 'Analysis/Abundance_silva/Family_reads.csv')
write.csv(df_family_re, 'Analysis/Abundance_silva/Family_relaAbundance.csv')

# order,silva
leveltax <- get_taxadf(obj=ps_dada2, taxlevel=4)
df_order <- leveltax@otu_table %>% as.data.frame()
df_order <- rename(df_order,"o__o")
df_order_re <- apply(df_order, 2, function(x) x/sum(x)) %>% as.data.frame()
write.csv(df_order, 'Analysis/Abundance_silva/Order_reads.csv')
write.csv(df_order_re, 'Analysis/Abundance_silva/Order_relaAbundance.csv')

# class,silva
leveltax <- get_taxadf(obj=ps_dada2, taxlevel=3)
df_class <- leveltax@otu_table %>% as.data.frame()
df_class <- rename(df_class,"c__c")
df_class_re <- apply(df_class, 2, function(x) x/sum(x)) %>% as.data.frame()
write.csv(df_class, 'Analysis/Abundance_silva/Class_reads.csv')
write.csv(df_class_re, 'Analysis/Abundance_silva/Class_relaAbundance.csv')

# phylum,silva
leveltax <- get_taxadf(obj=ps_dada2, taxlevel=2)
df_phylum <- leveltax@otu_table %>% as.data.frame()
df_phylum <- rename(df_phylum,"p__p")
df_phylum_re <- apply(df_phylum, 2, function(x) x/sum(x)) %>% as.data.frame()
write.csv(df_phylum, 'Analysis/Abundance_silva/Phylum_reads.csv')
write.csv(df_phylum_re, 'Analysis/Abundance_silva/Phylum_relaAbundance.csv')

# genus,gg
leveltax <- get_taxadf(obj=ps_dada2, taxlevel=6)
df_genus <- leveltax@otu_table %>% as.data.frame()
rownames(df_genus)[1]
df_genus <- rename(df_genus,"g")
df_genus_re <- apply(df_genus, 2, function(x) x/sum(x)) %>% as.data.frame()

write.csv(df_genus, 'Analysis/Abundance_gg/Genus_reads.csv')
write.csv(df_genus_re, 'Analysis/Abundance_gg/Genus_relaAbundance.csv')

# family,gg
leveltax <- get_taxadf(obj=ps_dada2, taxlevel=5)
df_family <- leveltax@otu_table %>% as.data.frame()
rownames(df_family)[1]
df_family <- rename(df_family,"f")
df_family_re <- apply(df_family, 2, function(x) x/sum(x)) %>% as.data.frame()
write.csv(df_family, 'Analysis/Abundance_gg/Family_reads.csv')
write.csv(df_family_re, 'Analysis/Abundance_gg/Family_relaAbundance.csv')

# order,gg
leveltax <- get_taxadf(obj = ps_dada2,taxlevel=4)
df_order <- leveltax@otu_table %>% as.data.frame()
rownames(df_order)[1]
df_order <- rename(df_order,"o")
df_order_re <- apply(df_order, 2, function(x) x/sum(x)) %>% as.data.frame()
write.csv(df_order, 'Analysis/Abundance_gg/Order_reads.csv')
write.csv(df_order_re, 'Analysis/Abundance_gg/Order_relaAbundance.csv')

# class,gg
leveltax <- get_taxadf(obj = ps_dada2,taxlevel=3)
df_class <- leveltax@otu_table %>% as.data.frame()
rownames(df_class)[1]
df_class <- rename(df_class,"c")
df_class_re <- apply(df_class, 2, function(x) x/sum(x)) %>% as.data.frame()
write.csv(df_class, 'Analysis/Abundance_gg/Class_reads.csv')
write.csv(df_class_re, 'Analysis/Abundance_gg/Class_relaAbundance.csv')

# phylum,gg
leveltax <- get_taxadf(obj=ps_dada2, taxlevel=2)
df_phylum <- leveltax@otu_table %>% as.data.frame()
df_phylum <- rename(df_phylum,"p")
df_phylum_re <- apply(df_phylum, 2, function(x) x/sum(x)) %>% as.data.frame()
write.csv(df_phylum, 'Analysis/Abundance_gg/Phylum_reads.csv')
write.csv(df_phylum_re, 'Analysis/Abundance_gg/Phylum_relaAbundance.csv')

# df_genus_new <- df_genus %>% 
#   select(colnames(df_genus)[grep('-',colnames(df_genus))])
# df_genus_old <- read.csv('E:/Lab_Work/FS/16S/FS0801/Analysis/Abundance_gg/Genus_reads.csv')
# df_genus_none <- df_genus_new[!rownames(df_genus_new) %in% df_genus_old$X,]
# 
# df_genus_yu <- df_genus %>% 
#   select(colnames(df_genus)[!grepl('-',colnames(df_genus))])
# df_genus_yu_none <- df_genus_yu[rownames(df_genus_none),]
