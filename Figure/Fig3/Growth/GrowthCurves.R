library(openxlsx)
library(tidyverse)
library(growthcurver)
library(gridExtra)
library(ggsci)



strain <- read_tsv("Strains.txt")
plot_d_adj <- read_tsv("growth_data.txt")
p <- list()
plot_d_adj$strain <- gsub('Phocaeicola','Bacteroides',plot_d_adj$strain)
strain$strain <- gsub('Phocaeicola','Bacteroides',strain$strain)
plot_d_adj$con <- factor(plot_d_adj$con,levels = unique(plot_d_adj$con))
plot_d_adj <- plot_d_adj %>% arrange(strain)

for (strain_n in unique(plot_d_adj$strain)) {
  df <- plot_d_adj %>% filter(strain==strain_n)
    p[[strain_n]] <- 
    ggplot(df,aes(group=con,color=con)) +
    #geom_point(aes(t, n)) +
    scale_x_continuous(expand = c(0,0))+
    geom_smooth(aes(nortime, .fitted), se = F,size = 0.6) +
    ylab("OD")+
    #scale_color_hue()+
    #scale_colour_brewer(palette = "Set3")+
    #scale_color_npg()+
    scale_color_manual(values = c("#9A4BA2","#4DBBD5FF","#00A087FF","#3C5488FF",
                                    "#F39B7FFF","#8491B4FF","#91D1C2FF","#E64B35FF"))+
    ggtitle(strain$strain_reocrd[match(unique(df$strain),strain$strain)])+
    facet_grid(.~drug)+
    theme_bw()+
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = 'black'),
      strip.text.x = element_blank(),
      #strip.text.y = element_text(size = 12, angle = 0,color = "red", face = "bold.italic"),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank())
}


ptotal <- gridExtra::grid.arrange(grobs = p,nrow = 25)

ggsave(plot = ptotal,filename = "growthcurve.pdf",width = 8,height = 25)



ggplot(df,aes(group=con,color=con)) +
  #geom_point(aes(t, n)) +
  scale_x_continuous(expand = c(0,0))+
  geom_smooth(aes(nortime, .fitted), se = F,size = 0.6) +
  ylab("OD")+
  #scale_color_hue()+
  #scale_colour_brewer(palette = "Set3")+
  #scale_color_npg()+
  labs(color = "Concentration(μg/mL)")+
  scale_color_manual(values = c("#9A4BA2","#4DBBD5FF","#00A087FF","#3C5488FF",
                                "#F39B7FFF","#8491B4FF","#91D1C2FF","#E64B35FF"))+
  ggtitle(strain$strain_reocrd[match(unique(df$strain),strain$strain)])+
  facet_grid(.~drug)+
  theme_bw()+
  theme(
    #legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    #strip.text.y = element_text(size = 12, angle = 0,color = "red", face = "bold.italic"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank())

ggsave(filename = "growth_legend.pdf",width = 8,height = 2)




data <- read_tsv("result_auc.txt")
data <- data %>% mutate(flag = case_when(
  auc_ic25 <= Conc ~ "yes",
  TRUE ~ "no"
))
data$strain <- gsub('Phocaeicola','Bacteroides',data$strain)

orders <- read_tsv('drug_data.txt')

data <- data %>% mutate(Drug=factor(Drug,levels = rev(unique(orders$Drug))))

data <- data %>% filter(drug == "MTX")



df <- plot_d_adj %>% filter(drug =="MTX")

ll <- df %>% group_by(strain) %>%
  summarise(max = max(.fitted),
            xmax = max(t),
            con = )
ll <- ll %>% left_join(data %>% select(strain,auc_ic25,flag))

df$strain <- factor(df$strain,levels = unique(df$strain))
df$strain <- fct_relevel(df$strain,"Bifidobacterium animalis",after = Inf)

p <- 
  ggplot(df,aes(group=con,color=con)) +
  #ggplot(df) +
  scale_color_manual(values = c("0"="#E64B35FF","0.5"="#91D1C2FF","1"="#8491B4FF","4"="#F39B7FFF",
                                "12"="#3C5488FF","40"="#00A087FF","120"="#4DBBD5FF","360"="#9A4BA2"
                                ))+
  scale_x_continuous(expand = c(0,0),#breaks = pretty(df$t, n = 3),
                     breaks = pretty_breaks(3),
                     )+
  scale_y_continuous(#breaks = pretty(df$t, n = 3),
                     breaks = pretty_breaks(3),
  )+
  #scale_y_continuous(expand = c(0,0),limits = c(0,0.5))+
  geom_smooth(data=df,aes(t, .fitted), se = F) +
  ylab("OD")+
  #ggtitle(unique(df$strain))+
  #facet_grid(strain~drug)+
  facet_wrap(~ strain,  ncol = 8,scales = "free",
             labeller = as_labeller(function(x) strain$label4[match(x,strain$strain)])) +  
  theme_bw()+
  guides(color = "none")+
  theme(
    legend.position = "none",
    strip.text = element_text(face = "italic",size = 7),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 8,color = "black"),
    plot.title = element_text(size = 18),
    axis.title.x = element_blank())

library(ggnewscale)

ll$con <- as.character(ll$con)
ll$strain <- factor(ll$strain,levels = levels(df$strain))
p+
  new_scale_color()+
  new_scale_fill()+
geom_rect(data = ll,
          aes(xmin = 0, xmax = xmax/6, ymin = max*0.9,
              ymax = max*1.1,fill = auc_ic25),size = 0.2,color = "gray60")+
  geom_text(data = ll %>% filter(flag == "yes"),
            aes(x = xmax/12, y = max  ,
                label = "+"),size = 3,hjust = 0.5,vjust = 0.5)+
  
scale_color_gradientn(trans = "log10",
                       colours = c('#550501','#AA0A02','#AA0A02','#FF5F57','#FF948F','#FFB7B4','#FFCFCD','#FFDFDD'),
                       #limits=c(-100,300),
                       na.value = "white",
                       breaks=c(0,0.5,1,4,12,40,120,360))+

scale_fill_gradientn(trans = "log10",
                       colours = c('#550501','#AA0A02','#AA0A02','#FF5F57','#FF948F','#FFB7B4','#FFCFCD','#FFDFDD'),
                       #limits=c(-100,300),
                       na.value = "white",
                       breaks=c(0,0.5,1,4,12,40,120,360),
                     labels = c(0,"0.5","1","4","12","40","120","360"))+
  theme(legend.direction = 'horizontal', 
        legend.position = c(0.5,0.1),
        legend.text = element_text(size = 10))
 

ggsave(filename = "MTX.pdf",width = 7,height = 4)
