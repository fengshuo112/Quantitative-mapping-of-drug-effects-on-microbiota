library(pacman)
arg <- commandArgs(T)
if(length(arg) != 1){
  cat("Argument: You should input a path!\n")
  quit('no')
}
Input_file <- paste(arg[1],"rep.fa",sep = "/")
Output_file <- paste(arg[1],"rep.xls",sep = "/")
pacman::p_load(tidyverse,magrittr,stringr)
rep <-  Input_file %>% read.delim(check.names = FALSE, row.names = 1) %>%
  set_rownames(paste0(">ASV", seq_len(nrow(.))))
write.table (rep,file =Output_file, sep ="\t", row.names = T,quote = F,col.names = F)
print("You are success!")
