setwd("D:/Desktop/supplementary/Figs/Quantitative-mapping-of-drug-effects-on-microbiota/Figure/Fig4/Feces_growth/")


library(tidyverse)
library(growthcurver)
library(ggpubr)
library(scales)
library(openxlsx)
library(ggsci)


sel <- read.xlsx("tableS3.xlsx",sheet = 1,startRow = 2)
sig_data1 <- read_tsv("sig_data.txt")
mean_sum <- read_tsv("mean_sum.txt")

# ATC ---------------------------------------------------------------------

atc <- read_tsv("extract_level_extend.txt")

sig_data1 <- sig_data1 %>% filter(group2 %in% sel$Drug)
sig_data1 <- sig_data1 %>% left_join(atc,by=c('group2'='drugname'))


unique(sig_data1$group2[which(sig_data1$p<=0.05)])

unique(sig_data1$l1[which(sig_data1$p<=0.05)])




# plot --------------------------------------------------------------------

bar <- mean_sum %>% 
  filter(Name %in% c('Control',sel$Drug)) %>% 
  left_join(sig_data1 %>% select(group2,p.signif,DrugID,l1,l2,par,p) %>% distinct(),
                              by = c('par'='par','Name'='group2'),multiple = "first")
###single levels

bar <- bar %>% filter(!is.na(Name),Name!="Control")
bar$l1[is.na(bar$l1)] <- "Unclassified"
bar$l2[is.na(bar$l2)] <- "Unclassified"
bar <- bar %>% arrange(l1)
bar$p.signif[which(bar$p.signif == "ns")] <- ''



heat <- bar %>% filter(par == unique(bar$par)[1])
heat$l2[is.na(heat$l2)] <- "Unclassified"

order <- heat %>% ungroup()%>% select(Name,l1,l2) %>% distinct()
order <- order %>% group_by(l1,l2) %>% mutate(N=n())


heat <- heat %>% arrange(desc(l1),desc(l2),desc(ratio))
heat$y <- 1:nrow(heat)

heat_data <- bar %>% left_join(heat %>% ungroup() %>% select(Name,y,l1,l2) %>% distinct())
heat_data <- heat_data %>% filter(par != "auc_e")
heat_data <- heat_data %>% mutate(par = fct_inorder(par),
                                x = as.integer(par))
xy <- heat_data %>% mutate(l1 = factor(l1,levels = unique(heat_data$l1)),
                           newy = y + (length(unique(heat_data$l1))-as.integer(l1))) 

X <- xy %>% ungroup() %>% select(x, par) %>% unique()
Y <- xy %>% ungroup() %>% select(newy, Name) %>% distinct()

YY <- xy %>% group_by(l1) %>% 
  summarise(y1=min(newy),y2=max(newy),.groups = "drop",x=6)
YY2 <- xy %>% group_by(l2) %>% 
  summarise(y1=min(newy),y2=max(newy),.groups = "drop",x=6.5)


mlog10_trans = function() trans_new("mlog10", function(x) -log10(x), function(y) 10**(-y))

YY2$l2 <- str_to_sentence(YY2$l2)

xy <- xy %>% mutate(size = case_when(
  p < 0.05 & p > 5e-5 ~ "small",
  p <= 5e-5 & p > 5e-8 ~ "mid",
  p <= 5e-8 ~ "large"
))
xy$size <- factor(xy$size)



ggplot(xy, aes(x, newy, fill = ratio)) + 
  geom_rect(aes(xmin = x - 0.5, xmax = x + 0.5, ymin = newy - 0.5, ymax = newy + 0.5), color = "gray80", size = 0.1) +
  geom_text(data = xy,aes(label = p.signif),size = 1)+
  geom_text(data = YY2, aes(x = x + 0.3, y = (y1 + y2) / 2, label = l2), inherit.aes = FALSE, hjust = "left", size = 1, angle = 0) +
  geom_segment(data = YY2 %>% dplyr::filter(y1 != y2),
               aes(x = x, xend = x, y = y1 - 0.3, yend = y2 + 0.3), inherit.aes = FALSE, color = "grey", size = 0.15) +
  scale_fill_gradient2(transform = "log2", low = "blue", mid = "white", high = "red",
                       breaks = c(0.5, 1, 2),
                       limits = c(0.5, 2),
                       labels = c(-1, 0, 1),
                       midpoint = 1, na.value = "#dddddd",
                       name = "log2fold change") +
  scale_x_continuous(position = "top", name = "", breaks = X$x, 
                     labels = c('Auc', 'Carrier ability', 'Growth rate', 'Time'), expand = c(0, 0)) +
  scale_y_reverse(name = "", breaks = Y$newy, labels = Y$Name, expand = c(0, 0)) +
  theme_minimal() + coord_fixed(clip = "off") +
  theme(strip.text.y = element_text(angle = 0, hjust = 0),
        panel.grid = element_blank(),
        axis.text.x.top = element_text(angle = 90, hjust = 0, color = "black", size = 3), 
        legend.position = "bottom", legend.box = "vertical",
        axis.text.y = element_text(size = 3, color = "black"),
        axis.title.y = element_text(size = 7, margin = margin(0, 0, 0, 0, unit = "cm")),
        strip.text.x = element_blank(),
        plot.margin = margin(0, 3, 0, 0, unit = "cm"),
        legend.key.height = unit(0.35, "line"),
        legend.margin = margin(0, 3, 0, 2, unit = "line"),
        legend.title = element_text(size = 3),
        legend.text = element_text(size = 3))



ggsave("./heatmap.pdf",width = 8, height =10, units = "cm")



# statistic ---------------------------------------------------------------

data_sig <- heat_data %>% filter(p.signif != '') %>%
  select(Name,par)

data_sig1 <- data_sig %>% group_by(par) %>% 
  summarise(N=n(),.groups = "drop")

data_sig_total <- data.frame(par="Total",N=length(unique(data_sig$Name)))

data_sig <- data_sig1 %>% bind_rows(data_sig_total)
data_sig <- data_sig %>% mutate(ratio = N/41)
data_sig <- data_sig %>% mutate(Other = 1-ratio)
data_sig <- data_sig %>%  pivot_longer(cols = c(ratio,Other), names_to = "condition", values_to = "value")

#write_tsv(data_sig,"./data_sig.txt")

ggplot(data = data_sig %>% filter(condition == "ratio"), aes(x = par, y = value,fill = par))+
  geom_bar(stat = "identity",color = "black",size = 0.1,width = 0.75)+
  scale_x_discrete(labels = c('Auc','Carrier ability','Growth rate','Time to \nmid exponential','Total'))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1)) +
  scale_fill_jama()+
  ylab("Fraction of drugs\n significantly impacting microbiota")+
  xlab("")+
  theme_classic() + 
  theme(
        axis.line = element_line(linewidth = .2),
        axis.ticks = element_line(linewidth = .2),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1,color = "black",
                                       size=15), 
        axis.text.y = element_text(size= 15,color = "black"),
        axis.title.y = element_text(size= 15,color = "black"),
        legend.position = "none",
        legend.key.height = unit(0.5, "line"),
        legend.key.width = unit(0.5, "line"),
        legend.margin = margin(0,3,0,2,unit="line"),
        legend.text = element_text(size=4))

  


ggsave("./par-bar.pdf",width = 4, height = 4)
