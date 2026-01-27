arg <- commandArgs(T)
if(length(arg) < 1){
  cat("Argument: Out_Dir \n")
  quit('no')
}
library(tidyverse)

#arg <- c( "/Volumes/Flower0501/2022_16s/Mala/Group0/Analysis/10_LEfSe/Total/Total.Genus.res","/Volumes/Flower0501/2022_16s/Mala/Group0/Analysis/17_kegg_L3/Total" , "/Volumes/Flower0501/2022_16s/Mala/Group0/Total-color.txt" )
print(arg)
setwd(arg[2])
data <- read.table(arg[1],row.names=1,sep="\t")
data <- na.omit(data)
filename <- basename(arg[1])

grp <- read.table(arg[3],header=T,sep="\t", comment.char = "")
rownames(grp)<-as.character(grp[,1])

gg <- unique(as.character(grp[,2]))
col <-unique(as.character(grp[,3]))
names(col)=gg
dd <- data[,2:3]
dd <-dd[order(as.character(dd[,1]),as.numeric(dd[,2])),]
dd <- dd %>% rownames_to_column("ID")
dd$ID <- factor(dd$ID,levels =dd$ID )
colnames(dd) <- c("ID","group","value")
for(i in 1:length(gg)){
  dd[dd$group == gg[i],"color"]  <- col[i]
}


ggplot(dd,aes(ID,as.numeric(value),fill=group))+
  geom_histogram(stat = "identity",color="black",binwidth = 5)+
  scale_fill_manual(values = col,name="")+
  scale_y_continuous(expand = c(0,0))+
  labs(y="LDA Score (log 10 )",x="")+
  coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_text(colour = dd$color,face = "bold"),
        legend.position = "top",
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

ggsave(paste0(filename,".pdf"),width = max(dd$value)*2.5,height = 8)
ggsave(paste0(filename,".png"),width = max(dd$value)*5,height = 8,units ="cm",dpi=300)
