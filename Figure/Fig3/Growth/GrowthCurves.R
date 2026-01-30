library(openxlsx)
library(tidyverse)
library(growthcurver)
library(gridExtra)
library(ggsci)

setwd("./Growth curve")

strain <- read_tsv("Strains.txt")
plot_d_adj <- read_tsv("growth_data.txt")
p <- list()
plot_d_adj$con <- factor(plot_d_adj$con,levels = unique(plot_d_adj$con))
plot_d_adj <- plot_d_adj %>% arrange(strain)

for (strain_n in unique(plot_d_adj$strain)) {
  df <- plot_d_adj %>% filter(strain==strain_n)
    p[[strain_n]] <- 
    ggplot(df,aes(group=con,color=con)) +
    scale_x_continuous(expand = c(0,0))+
    geom_smooth(aes(nortime, .fitted), se = F,size = 0.6) +
    ylab("OD")+
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
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank())
}


ptotal <- gridExtra::grid.arrange(grobs = p,nrow = 25)

ggsave(plot = ptotal,filename = "growthcurve.pdf",width = 8,height = 25)



ggplot(df,aes(group=con,color=con)) +
  scale_x_continuous(expand = c(0,0))+
  geom_smooth(aes(nortime, .fitted), se = F,size = 0.6) +
  ylab("OD")+
  labs(color = "Concentration(μg/mL)")+
  scale_color_manual(values = c("#9A4BA2","#4DBBD5FF","#00A087FF","#3C5488FF",
                                "#F39B7FFF","#8491B4FF","#91D1C2FF","#E64B35FF"))+
  ggtitle(strain$strain_reocrd[match(unique(df$strain),strain$strain)])+
  facet_grid(.~drug)+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank())

ggsave(filename = "growth_legend.pdf",width = 8,height = 2)

