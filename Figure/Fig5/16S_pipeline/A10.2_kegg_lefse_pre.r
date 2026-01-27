arg <- commandArgs(T)
if(length(arg) < 2){
  cat("Argument: Input_path Out_path Group_File GrouCol\n")
  quit('no')
}


pacman::p_load(tidyverse)

make_lefse <- function(dat,group,arg){
  if(is.na(arg[5])){
    grp_col <- as.numeric(arg[4])
    group_name <- c("Class")
    for(i in colnames(dat)[2:ncol(dat)]){
      group_name <- c(group_name,group[grep(i,group[,1]),grp_col])
    }
    dat <- rbind(colnames(dat),dat)
    colnames(dat) <- group_name
    dat <- rbind(colnames(dat),dat)
  }else{
    grp_col <- as.numeric(arg[4])
    sub_col <- as.numeric(arg[5])
    group_name <- c("Class")
    sub_group_name <-  c("subClass")
    for(i in colnames(dat)[2:ncol(dat)]){
      group_name <- c(group_name,group[grep(i,group[,1]),grp_col])
    }
    for(i in colnames(dat)[2:ncol(dat)]){
      sub_group_name <- c(sub_group_name,group[grep(i,group$sample.id),grp_col])
    }
    anno <- rbind(data.frame(t(group_name)),data.frame(t(sub_group_name)))
    colnames(anno) <- colnames(dat)
    dat <- rbind(anno,colnames(dat),dat)
  }
  return(dat)
}

#arg <- c("/Volumes/Flower0501/2022_16s/CY/脚本/Total.kegg.path.txt","/Volumes/Flower0501/2022_16s/CY/脚本/","/Volumes/Flower0501/2022_16s/CY/脚本/Total-color.txt",2)
InFile <- arg[1]
output_path <- arg[2]
Group <- arg[3]

df <- read.table(InFile,header = T)
group <- read.table(Group,comment.char = "",header = T)
df_lefse <- make_lefse(df,group,arg)
write_tsv(df_lefse,paste0(output_path,"Total.kegg.4lefse.txt"),quote = "none",col_names = FALSE)
