library(tidyverse)
library(fmsb)
setwd("D:/Desktop/supplementary/Figs/Fig2")
iden <- read_tsv("radar/d0.tab")
abun <- read_tsv("radar/d1.tab")
pre <- read_tsv("radar/d2.tab")
num<- read_tsv("radar/d5.tab")
id <- read.table("radar/id.txt",sep='\t',header=T)
data <- full_join(iden[,c(1,3)],abun[,c(1,3)])
data <- full_join(data,pre[,c(1,3)]) 
num <- inner_join(num,id,by=c('drugname'='Drug'))
data <- inner_join(data, num, by = c("protein" = "Target"))
df<- data[,c(1:6)]
df <- unique(df)

df$name<- paste(df$drugname,df$protein,  sep = "-")
df <- df[,c(2,3,4,6,7)]
colnames(df) <- c('averiden','averabun','averpre','number','name')

min <- data.frame(averiden=min(iden$averiden),averabun=min(abun$averabun),averpre=min(pre$averpre),number=min(num$N),name='min')
max <- data.frame(averiden=max(iden$averiden),averabun=max(abun$averabun),averpre=max(pre$averpre),number=max(num$N),name='max')

df <- rbind(min,df)
df <- rbind(max,df)
rownames(df) <- df$name
df <- df[,c(1:4)]


# plot_drug --------------------------------------------------------------------
library(RColorBrewer)
colors <- brewer.pal(10, "Set3")

titles <- rownames(df)[-c(1,2)]
colnames(df) <- c("1","2","3","4")
create_beautiful_radarchart <- function(data, color = "#00AFBB", axistype = 0,cex = 1.5,
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
   radarchart(
      data, axistype = axistype,
      # Customize the polygon
      pcol = color, 
      pfcol = scales::alpha(color, 0.5),
      plwd = 2, plty = 1,
      # Customize the grid
      cglcol = "grey", cglty = 1, cglwd = 0.8,
      # Customize the axis
      axislabcol = "grey", 
      # Variable labels
      vlcex = vlcex, vlabels = vlabels,
      caxislabels = caxislabels, title = title, ...
   )
}

pdf(file = "radar/radar.pdf",height = 4,width = 8)
op <- par(mar = c(2, 2, 2, 2),xpd=T)
par(mfrow = c(2,5))


for(i in 1:10){
   create_beautiful_radarchart(
      data = df[c(1, 2, i+2), ], caxislabels = NULL,
      color =colors, 
      #title = titles[i],
      title = paste0(strsplit(titles[i],split = "-")[[1]][1],
                     "\n",
                     strsplit(titles[i],split = "-")[[1]][2]),
      xpd=T,
      cex = 1,
      vlcex = 1.7,
      axistype =0
     # vlabels
   )
}
dev.off()
