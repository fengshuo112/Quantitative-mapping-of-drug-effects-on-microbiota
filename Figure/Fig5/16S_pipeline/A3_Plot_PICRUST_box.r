arg <- commandArgs(T)
if(length(arg) < 3){
  cat("Argument: Input_File Out_path Color\n")
  quit('no')
}
pacman::p_load(ggplot2,reshape2)
#arg <- c("/Volumes/Flower0501/2022_16s/CY/Total.catagrized.kegg.csv","/Volumes/Flower0501/2022_16s/CY","/Volumes/Flower0501/2022_16s/CY/Total-color.txt")
outDir <- arg[2]
data <- read.csv(arg[1],row.names = 1,header = T)
Group <- read.table(arg[3],header = T,comment.char = "")
G_cc <- unique(Group[,2])

data_plot <- data[,(colnames(data)[2:(ncol(data)-3)])]
data_norm <- data_plot
for(i in 1:ncol(data_plot)){
  sample_sum <- sum(data_plot[,i])
  for(j in 1:nrow(data_plot)){
    data_norm[j,i] = round(data_plot[j,i]/sample_sum,3)
  }
}
data_pathway <- data[order(data$relati.abun,decreasing = T),][1:50,]
data_plot2 <- cbind(data$name,data$description,data$relati.abun,data_norm)
#colnames(data_plot2) <- c("name","description","sample","value")
colnames(data_plot2) <- c("name","description","relati.abun",colnames(data_plot))
dd <- data_plot2[data_plot2$description %in% data_pathway$description,]
#abundance <- data_pathway$relati.abun


for(i in 1:nrow(dd)){
  y <- dd[i,colnames(dd) %in% Group[Group[,2]==G_cc[1],][[1]]]
  x <- dd[i,colnames(dd) %in% Group[Group[,2]==G_cc[2],][[1]]]
  temp <- try(t.test(x,y)$p.value, silent = T)
  dd[i,"P_value"] <- ifelse(class(temp) == "try-error", NA, temp)
}

dd_p <- na.omit(dd)
dd_p$p.adj <- p.adjust(dd_p$P_value,method = "BH",n = nrow(dd_p))

dd_box <- melt(cbind(data[ncol(data)],data[2],data_norm))

for(i in 1:nrow(dd_p)){
  dd_p[i,"Shapen"]
}


png(paste0(outDir,"Total.kegg.boxplot.png"),height = 4*480/3,width = 4*480/2,bg = "transparent")
y <- ggplot(dd_p, aes(x=relati.abun, y=description)) + 
  geom_point(aes(color=-1*log10(p.adj),size = relati.abun,shape = factor(name,levels = c(unique(name)))))+
  scale_colour_gradient(low="blue",high="red") + 
  labs(color=expression(-log10),shape = "KEGG Path",size="Relative Abundance", x=" ",y="Pathway name",title="Pathway enrichment")+
  theme_classic()
# 画图确定x、y轴并根据通路富集基因数改变点的大小，根据p值大小改变点的颜色 + 自定义p值的渐变颜色 + 设置标题及标注名称
y                 #显示气泡图
dev.off()

pdf(paste0(outDir,"Total.kegg.boxplot.pdf"),width = 10,height = 6)
y <- ggplot(dd_p, aes(x=relati.abun, y=description)) + 
  geom_point(aes(color=-1*log10(p.adj),size = relati.abun,shape = factor(name,levels = c(unique(name)))))+
  scale_colour_gradient(low="blue",high="red") + 
  labs(color=expression(-log10),shape = "KEGG Path",size="Relative Abundance", x=" ",y="Pathway name",title="Pathway enrichment")+
  theme_classic()
# 画图确定x、y轴并根据通路富集基因数改变点的大小，根据p值大小改变点的颜色 + 自定义p值的渐变颜色 + 设置标题及标注名称
y 
dev.off()
