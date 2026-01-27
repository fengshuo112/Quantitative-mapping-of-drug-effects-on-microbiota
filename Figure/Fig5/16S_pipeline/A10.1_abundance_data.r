arg <- commandArgs(T)
if(length(arg) < 1){
  cat("Argument: Input_path Out_path\n")
  quit('no')
}
library(tidyverse)

rename <- function(df,label){
  df <- df[!grepl("__un_",rownames(df)),]
  df$name <- unlist(lapply(rownames(df),function(x) str_match_all(x,paste0(label,"__([A-Za-z1-9\\[\\]-]+)"))[[1]][2]))
  dat_final <- aggregate(df[,c(1:(ncol(df)-1))],by=list(name=df$name),sum)
  dat_final <- column_to_rownames(dat_final,"name")
  return(dat_final)
}
print(arg)

#arg <- c("/Volumes/Flower0501/2022_16s/CY/CY20220706/Analysis/0_input_data/","/Volumes/Flower0501/2022_16s/CY/CY20220706/Analysis/10_LEfSe/Total/","/Volumes/Flower0501/2022_16s/CY/CY20220706/Total-color.txt",2)
input_path <- arg[1]
output_path <- arg[2]

df_otu <- read.table(paste0(input_path,"Reotu_table.tsv"),header = T,row.names = 1,check.names = FALSE)
df_abund <-  matrix(ncol = ncol(df_otu),nrow = nrow(df_otu))%>% as.data.frame()
for(i in 1:ncol(df_otu)){
  Sum <- sum(df_otu[,i])
  for(n in 1:nrow(df_otu)){
    df_abund[n,i] <- round(df_otu[n,i]/Sum,5)
  } 
}

colnames(df_abund) <- colnames(df_otu)
rownames(df_abund) <- rownames(df_otu)

df_abund <- df_abund[rowSums(df_abund)>0,colSums(df_otu)>0] 





tax <- c("p","c","o","f","g","s")
All_nam <- c("Phylum","Class","Order","Family","Genus","Species")
df_tax_all <- read.table(paste0(input_path,"taxonomy.tsv"),header=T,sep = "\t",row.names = 1)

df_abund <- cbind(df_abund,df_tax_all[rownames(df_abund),"Taxon"])
colnames(df_abund)[ncol(df_abund)] <- "Lineage"
df_abund <- df_abund %>% rownames_to_column("#OTU ID")
df_abund[,1] <- sub("ASV","",df_abund[,1])
write.table(df_abund,paste0(output_path,"Total.normolized.out.tab"),quote = FALSE,row.names = FALSE,sep = "\t")

##### Genus abundance#####
df_abund[,"genus"] <- sapply(df_abund$Lineage,function(x) str_split(x,";",simplify = T)[6])
df_abund[,"genus"] <- sapply(unlist(df_abund["genus"]),function(x) str_match(x,"g__([A-Za-z0-9\\[\\]\\-]+)")[2])
df_genus <- df_abund[-which(is.na(df_abund["genus"])),] %>% dplyr::select(-c("Lineage","#OTU ID"))
df_genus <- aggregate(df_genus[1:(ncol(df_genus)-1)],by = list(ID=df_genus$genus),FUN=sum,drop=FALSE)

write.table(df_genus,paste0(output_path,"Total.Genus.abundance.txt"),quote = FALSE,row.names = FALSE,sep = "\t")

##### KEGG Path abundance######
KeggPath <- paste0(input_path,"Total.kegg.path.csv")
df_kegg <- read.csv(KeggPath,header = T,check.names = FALSE)
df_kegg <- df_kegg[,-c(1,(ncol(df_kegg)-2):ncol(df_kegg))]
df_kegg_abund <- matrix(nrow = nrow(df_kegg),ncol = ncol(df_kegg)) %>% as.data.frame()
for(i in 1:nrow(df_kegg)){
  Sum <- sum(df_kegg[i,2:ncol(df_kegg)])
  for(n in 2:ncol(df_kegg)){
    df_kegg_abund[i,n] <- round(df_kegg[i,n]/Sum,5)
  }
}
df_kegg_abund[,1] <- df_kegg[,1]
colnames(df_kegg_abund) <- colnames(df_kegg)
colnames(df_kegg_abund)[1] <- "ID"
write.table(df_kegg_abund,paste0(output_path,"Total.kegg.L3.abundance.txt"),quote = FALSE,row.names = FALSE,sep = "\t")
