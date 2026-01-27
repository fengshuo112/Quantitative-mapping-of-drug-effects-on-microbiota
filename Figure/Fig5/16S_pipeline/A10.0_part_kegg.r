df <- read.delim('../Analysis/0_input_data/Reotu_table.tsv')
group <- read.delim('Total-color.txt')
df <- df[,c('ASV',group$sample.id)]
write.table(df, '../Analysis/0_input_data/Reotu_table.tsv',sep = '\t', quote = F,row.names = F)

