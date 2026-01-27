library(ggpubr)
library(tidyverse)
# fun ---------------------------------------------------------------------


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
data <- read_tsv("qpcr_data.txt")


df <- data %>% group_by(strain,gene,drug,time) %>% 
  mutate(mean=mean(value))
df <- df %>% mutate(deltadelta =value-df$mean[which(df$gene == "16S"& 
            df$strain==strain &df$drug == drug & df$time == time )])
df <- df %>% mutate(log =2^-deltadelta)


p <- list()
df$time <- factor(df$time,levels = c("0h","2h","24h"),
                  labels = c("0 h","2 h","24 h"))  #



i=0
for (strains in c("E. fergusonii" ,"B. ovatus"))  {  
   plot_d<- df %>% filter(strain == strains,gene != "16S")
  library(dplyr)
  
  
  sigs <- compare_means(log ~ drug, data = plot_d,group.by = c("time","gene"),
                        method = "t.test")
  
  sta=summarySE(plot_d, measurevar="log", groupvars=c("drug","gene","time"))
  sta <- sta %>% left_join(plot_d %>% select(drug,group) %>% distinct())
 
  library(dplyr)
  detach(package:plyr)
  y <- plot_d %>% group_by(time,gene) %>% summarise(ylog=max(log))
  sigs <- sigs %>% left_join(y)
  plot_d$time <- gsub("h"," h",plot_d$time)
  
  for (Gene in unique(plot_d$gene)) {
    plot_dd <- plot_d %>% filter(gene == Gene)
    plot_dd$drug<- factor(plot_dd$drug,levels = c('DMSO',"Canagliflozin"))
    sta$drug<- factor(sta$drug,levels = c('DMSO',"Canagliflozin"))
    mean_1 <- sta %>%filter(drug=='DMSO')
    mean_1 <- mean_1 %>% select(gene,time,log,group,strain) %>% 
    rename(log_2 = log)
    sta <- left_join(sta,mean_1)
    sigs <- left_join(sigs ,mean_1 %>% select(gene,time,log_2))
    i = i+1
    p[[i]] <- 
      ggplot(data = plot_dd ,aes(x = drug, y =1))+
      geom_bar(data= sta %>% filter(gene == Gene),aes(x = time, y =  log/log_2, fill = drug),
               width = 0.6,stat = "identity",position = "dodge")+
      scale_y_continuous(expand = c(0,0))+#,breaks = seq(0,20,by= 5)
      scale_fill_manual(values = c("gray30","#F8B6A7"),name = NULL)+
      geom_errorbar(data = sta %>% filter(gene == Gene),
                    aes(x = time,ymin = (log - sd)/log_2, ymax = (log + sd)/log_2,group = drug),
                    position = position_dodge(0.6),
                    size = 0.2,
                    width = 0.2)+
      
      geom_hline(yintercept = 1, linetype="dashed",color="gray80")+
      geom_text(data = sigs%>% filter(gene == Gene), aes(x = time, y = (1.05 *ylog)/log_2,
                                 label = paste0("p = ",ifelse(p>0.05,round(p,2),formatC(p, format = "e", digits = 1)))), 
                position = position_dodge(0.6), size = 1.5, vjust = 0) +  # 添加显著性标注
      geom_segment(data = sigs %>% filter(gene == Gene), 
                   aes(x = as.numeric(time) - 0.1, xend = as.numeric(time) + 0.2, y = (1.02 * ylog) / log_2, yend = (1.02 * ylog) / log_2), 
                   position = position_dodge(0.6), 
                   size = 0.1, color = "black") + 
       coord_cartesian(ylim =c(0,5))+
      ggtitle(paste0(unique(plot_dd$gene)),subtitle = unique(plot_dd$strain))+
      ylab('Relative mRNA level') +#Log[2]~fold~change))+
      xlab("")+
      theme_linedraw()+
      guides(fill = guide_legend(title = NULL, 
                                 label.theme = element_text(size = 4),  
                                 keywidth = 0.4,  
                                 keyheight = 0.4)  
      ) +  
      theme(text = element_text(colour = "black"))+
      theme(panel.grid = element_blank(),
            plot.subtitle = element_text(size = 8,face = "italic",hjust = 0.5),
            strip.text = element_text(color = 'black',size = 10),
            strip.background = element_blank(),
            legend.position = c(1, 1),
            legend.justification = c(1, 1) ,
            legend.background = element_blank() ,
            axis.text.x = element_text(size = 8),
            axis.title = element_text(size = 8),
            #legend.position = 'none',
            # axis.text.x = element_text(face = "italic"),
            plot.title = element_text(hjust = 0.5,size = 8)
      )
     }
  
  
}


for (strains in c("B. catenulatum" ,"C. sporogenes"))  {  # {  #"Canagliflozin","Methotrexate" )) {
   plot_d<- df %>% filter(strain == strains,gene != "16S")
   library(dplyr)
   
   sigs <- compare_means(log ~ drug, data = plot_d,group.by = c("time","gene"),
                         method = "t.test")
   
   sta=summarySE(plot_d, measurevar="log", groupvars=c("drug","gene","time"))
   sta <- sta %>% left_join(plot_d %>% select(drug,group) %>% distinct())
   
   library(dplyr)
   detach(package:plyr)
   y <- plot_d %>% group_by(time,gene) %>% summarise(ylog=max(log))
   sigs <- sigs %>% left_join(y)
   plot_d$time <- gsub("h"," h",plot_d$time)
   
   for (Gene in unique(plot_d$gene)) {
      plot_dd <- plot_d %>% filter(gene == Gene)
      plot_dd$drug<- factor(plot_dd$drug,levels = c('DMSO',"Methotrexate"))
      sta$drug<- factor(sta$drug,levels = c('DMSO',"Methotrexate"))
      mean_1 <- sta %>%filter(drug=='DMSO')
      mean_1 <- mean_1 %>% select(gene,time,log,group,strain) %>% 
        rename(log_2 = log)
      
      sta <- left_join(sta,mean_1)
      sigs <- left_join(sigs ,mean_1 %>% select(gene,time,log_2))
      i = i+1
      p[[i]] <- 
         ggplot(data = plot_dd ,aes(x = drug, y =1))+
         geom_bar(data=sta%>% filter(gene == Gene),aes(x = time, y =  log/log_2, fill = drug),
                  width = 0.6,stat = "identity",position = "dodge")+
         scale_y_continuous(expand = c(0,0),breaks = c(0,1,2,3))+#,breaks = seq(0,20,by= 5)
         scale_fill_manual(values = c("gray30","#909CCF"),name = NULL)+
         geom_errorbar(data = sta %>% filter(gene == Gene),
                       aes(x = time,ymin = (log - sd)/log_2, ymax = (log + sd)/log_2,group = drug),
                       position = position_dodge(0.6),
                       size = 0.2,
                       width = 0.2)+
        geom_segment(data = sigs %>% filter(gene == Gene), 
                     aes(x = as.numeric(time) - 0.1, xend = as.numeric(time) + 0.2, y = (1.02 * ylog) / log_2, yend = (1.02 * ylog) / log_2), 
                     position = position_dodge(0.6), 
                     size = 0.1, color = "black") + 
         geom_hline(yintercept = 1, linetype="dashed",color="gray80")+
         geom_text(data = sigs%>% filter(gene == Gene), aes(x = time, y = (1.05 *ylog)/log_2,
                                                            label = paste0("p = ",ifelse(p>0.05,round(p,2),formatC(p, format = "e", digits = 1)))), 
                   position = position_dodge(0.6), size = 1.5, vjust = 0) +  # 添加显著性标注
         coord_cartesian(ylim =c(0,3))+
         ggtitle(paste0(unique(plot_dd$gene)),subtitle = unique(plot_dd$strain))+
         ylab('Relative mRNA level') +#Log[2]~fold~change))+
         xlab("")+
         theme_linedraw()+
         guides(fill = guide_legend(title = NULL, 
                                    label.theme = element_text(size = 4),  
                                    keywidth = 0.4,  
                                    keyheight = 0.4)  
         ) +  
         theme(text = element_text(colour = "black"))+
         theme(panel.grid = element_blank(),
               plot.subtitle = element_text(size = 8,face = "italic",hjust = 0.5),
               strip.text = element_text(color = 'black',size = 10),
               strip.background = element_blank(),
               legend.position = c(1, 1),
               legend.justification = c(1, 1) ,
               legend.background = element_blank() ,
               axis.text.x = element_text(size = 8),
               axis.title = element_text(size = 8),
               #legend.position = 'none',
               # axis.text.x = element_text(face = "italic"),
               plot.title = element_text(hjust = 0.5,size = 8)
         )
      }
   
   
}
library(patchwork)
(p[[1]]+p[[2]])/(p[[3]]+p[[4]])
ggsave(filename = "expression_sep.pdf",width = 3,height = 3)
