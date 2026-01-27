library(pacman)
arg <- commandArgs(T)
if(length(arg) != 1){
cat("Argument:You should inout a path!\n")
quit('no')}
input_file = paste(arg[1],"otu_table.tsv",sep = "/")
output_file = paste(arg[1],"Reotu_table.tsv",sep = "/")
pacman::p_load(tidyverse,magrittr,stringr)
otu <- input_file %>% read.table(check.names = FALSE,header = T,sep="\t")

rown <- paste0("ASV",seq_len(nrow(otu)))
otu[,1] <- rown
colnames(otu)[1] <- paste0("ASV",colnames(data)[1])
write.table (otu,file =output_file,sep ="\t",row.names = F,quote = F)
print("Success!")
