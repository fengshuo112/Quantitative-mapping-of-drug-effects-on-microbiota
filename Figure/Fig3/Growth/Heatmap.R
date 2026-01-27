library(ggnewscale)
library(tidyverse)
library(ggsci)
library(tidyverse)
library(ggplot2)
library(scatterpie)
library(patchwork)

data <- read_tsv("result_auc.txt")



tax <- read_tsv("Strains.txt")
order <- read.table(header = F,text = "Bacteroides
Prevotella
Alistipes
Eubacterium
Bifidobacterium
Ruminococcus
Parabacteroides
Faecalibacterium
Clostridium
Akkermansia
Escherichia
Fusobacterium
")


ord <- c("Bacteroides" , "Prevotella" ,"Alistipes",  
         "Eubacterium", "Bifidobacterium",  "Ruminococcus",
         "Parabacteroides", "Faecalibacterium", "Clostridium",     
         "Akkermansia", "Escherichia","Fusobacterium")

pal_colors <- c(pal_npg("nrc")(7),pal_jama("default")(7))
pal_colors <- pal_colors[c(1:10,12:13)]
names(pal_colors) <- order$V1


data <- data %>% mutate(flag = case_when(
  auc_ic25 <= Conc ~ "yes",
  TRUE ~ "no"
))
tax$Genus <- factor(tax$Genus,levels = ord)
tax <- tax %>% arrange(Genus)
data$strain <- factor(data$strain,levels = tax$strain,labels = tax$label3)
orders <- read_tsv('drug_data.txt')

data <- data %>% mutate(Drug=factor(Drug,levels = rev(unique(orders$Drug))))
data <- data %>% mutate(y = as.integer(strain),
                        x = as.integer(Drug))

lab <- data %>% filter(flag == "yes")

tax$y <- 1:nrow(tax)
tax$Genus <- factor(tax$Genus,levels = unique(order$V1))
xx <- data %>% select(x,Drug) %>% distinct()

rect <- data %>% group_by(Drug) %>% 
  summarise(xmin = min(x)-0.5,xmax = max(x)+0.5)





pheat <- ggplot()+
  annotate("rect", xmin = 0+0.5, xmax = 10+0.5, 
           ymin = 0+0.5, ymax = 25+0.5, 
           size = 1,
           fill = "white", color = "gray40") +
  geom_tile(data = data,aes(x = x, y = y ,fill = auc_ic25),
            color="white",size=0.1)+
  #geom_text(data = lab,aes(x = x, y = y ,label = "▲"),size = 3,family = "calibri")+
  geom_text(data = lab,aes(x = x, y = y ,label = "+"),size = 2)+
  coord_equal(clip = "off") + 
  geom_text(data = tax,aes(x = -1, y = y, label = label3),
            hjust = 1, fontface = "italic",size = 2.5, angle = 45)+
  geom_text(data = xx,aes(x = x, y = -0.1, label = Drug),hjust = 1, angle = 0,size = 2.5)+
  coord_flip(clip = "off")+
  scale_fill_gradientn(trans = "log10",
                       colours = c('#550501','#AA0A02','#AA0A02','#FF5F57','#FF948F','#FFB7B4','#FFCFCD','#FFDFDD'),
                       #limits=c(-100,300),
                       na.value = "white",
                       breaks=c(0,0.5,1,4,13,40,121,364))+
  ggnewscale::new_scale_fill()+
  geom_rect(data = tax, aes(xmin = -0.5,xmax = -0.1,ymin = y-0.5,ymax = y+0.5,fill = Genus))+  # 添加矩形
  scale_fill_manual(values = pal_colors)+
  theme_void()+
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        plot.margin = margin(0,1,5,5,unit = "cm"))




data <- read_tsv("drug_data.txt")

data$Drug <- factor(data$Drug,levels = unique(data$Drug))


p1 <- ggplot()+
  geom_scatterpie(data=data,aes(x=x,y=y,r=0.4,group=Drug),
                  color = "white",
                  size = 0.2,
                  cols = c('Fraction excreted in feces',
                           'Other'))+
  geom_text(data=data,aes(x=1,y=y,label=paste0(100*data$`Fraction excreted in feces`,"%")),
            size = 1.5,hjust =0)+
  scale_fill_manual(values = c("#B997D3","gray80"))+
  ggtitle('Fecal elimination')+
  ylab('')+
  coord_equal(clip = "off")+
  theme_void()+
  theme(legend.position = 'none',
        plot.margin = margin(0,100,0,0),
        plot.title = element_text(size = 4,hjust = 0.5)
  )


p1
p2 <- ggplot()+
  geom_tile(data=data,aes(x=x,y=y,fill=Concentration),width=0.8,height=0.8)+
  coord_fixed(1)+
  geom_text(data=data,aes(x=1,y=y,label=data$Concentration),
            size=1.5,hjust =0.5)+
  scale_fill_gradientn(trans = "log10",
                       colours = c('#550501','#AA0A02','#AA0A02','#FF5F57','#FF948F','#FFB7B4','#FFCFCD','#FFDFDD'),
                       #limits=c(0,370),
                       na.value = "white",
                       breaks=c(0,0.5,1,4,13,40,121,364))+
  ggtitle("Estimated maximum\nconcentraiton in gut")+
  xlab('')+
  ylab('')+
  coord_equal(clip = "off")+
  theme_void()+
  theme(legend.position = 'none',
        plot.margin = margin(0,0,0,3),
        plot.title = element_text(size = 4,hjust = 0.5),
        
        axis.text = element_blank())


p2

p3=ggplot()+
  geom_text(data=data,aes(x=x,y=y,label=sub("[.] ","\n",ATC)),
            hjust = 0,size=1.5)+
  xlab('')+
  ylab('')+
  ggtitle("ATC classification")+
  theme_void()+
  theme(plot.margin = margin(2,2,0,2),
        legend.position = 'top',
        plot.title = element_text(size = 4,hjust = 0.5),
        axis.text = element_blank())

p3


pheat+p1+p2+p3+
  plot_layout(ncol = 4, byrow = TRUE)

ggsave("heatmap.pdf",width = 12,height = 4)








