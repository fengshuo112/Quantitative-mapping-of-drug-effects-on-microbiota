library(tidyverse)
library(ggsci)
library(ggpubr)
library(patchwork)

main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=0.5, colour="black"),
                   axis.line.y=element_line(size=0.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=12),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=12),
                   text=element_text(family="sans", size=12),
                   plot.title=element_text(hjust = 0.5,vjust=0.5,size=12),
                   plot.subtitle=element_text(size=12))


# Statistic ---------------------------------------------------------------
data <- read_tsv("Dock_growth.txt")

data <- data %>% 
  complete(Drug,target,prediction,affect) %>% 
  group_by(Drug,affect,target,prediction) %>% 
  summarise(.,num=n_distinct(Strain,na.rm = T))
data <- data %>% filter(!(target == "without"& prediction == "decrease"))
data$group <- paste0(data$affect,"_",data$target,"_",data$prediction)

df <- data %>% ungroup() %>% group_by(Drug,affect) %>% 
  mutate(sum = sum(num),
         percent = num/sum)
df_median <- df %>% group_by(group) %>% 
  summarise(median = median(percent))

p1 <- ggplot(data = df %>% filter(affect == "inhibit"),
       mapping = aes(x = group,y = percent)) +
  geom_crossbar(data = df_median %>% filter(grepl("^inhibit",group)),
                size = .5,
                width = 0.2,
                mapping = aes(y = median, ymin = median, ymax = median,color = group), width = 0.5) +
  scale_fill_nejm()+
  scale_color_nejm()+
  geom_jitter(data = df %>% filter(affect == "inhibit"),
              mapping = aes(x = group,y = percent,fill = group),
              size=4,alpha = 0.6,
              shape = 21,  color = "transparent",
              height = 0,width = 0.1,show.legend = F) +
  geom_signif(data = df %>% filter(affect == "inhibit"),mapping = aes(group,percent),
              comparisons = list(c("inhibit_with_decrease","inhibit_with_no-decrease"),
                                 c("inhibit_with_no-decrease","inhibit_without_no-decrease"),
                                 c("inhibit_with_decrease","inhibit_without_no-decrease")),
              map_signif_level = function(p) sprintf("P = %.2g", p),
              tip_length = c(0,0,0,0,0,0),
              textsize = 6,
              y_position = c(1.15,1,1.3),
              test = "wilcox.test") +
    scale_y_continuous(#limits = c(0,1.2),
    labels = scales::percent_format(accuracy = 1),
                     breaks = c(0,0.25,0.5,0.75,1),
                     expand = c(0.05,0.05)) +
  scale_x_discrete(labels = c('With specific HA-MDTPH',
                              'With specific LA-MDTPH',
                              'Without specific MDTPH')) +
  theme_classic(base_line_size = 0.8) +
  coord_cartesian(ylim = c(0,1.4))+
  labs(title = "Sensitive",x="",y="") +
  main_theme+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "none",
    plot.title = element_text(size = 18))


p1


p2 <- 
  ggplot(data = df %>% filter(affect == "no-inhibit"),
         mapping = aes(x = group,y = percent)) +
  geom_crossbar(data = df_median %>% filter(grepl("no-inhibit",group)),
                size = .5,
                width = 0.2,
                mapping = aes(y = median, ymin = median, ymax = median,color = group), width = 0.5) +
  scale_fill_nejm()+
  scale_color_nejm()+
  geom_jitter(data = df %>% filter(affect == "no-inhibit"),
              mapping = aes(x = group,y = percent,fill = group),
              size=4,alpha = 0.6,
              shape = 21,  color = "transparent",
              height = 0,width = 0.1,show.legend = F) +
  geom_signif(data = df %>% filter(affect == "no-inhibit"),
              mapping = aes(group,percent),
              comparisons = list(c("no-inhibit_with_decrease","no-inhibit_with_no-decrease"),
                                 c("no-inhibit_with_no-decrease","no-inhibit_without_no-decrease"),
                                 c("no-inhibit_with_decrease","no-inhibit_without_no-decrease")),
              map_signif_level = function(p) sprintf("P = %.2g", p),
              tip_length = c(0,0,0,0,0,0),
              y_position = c(1.15,1,1.3),
              textsize = 6,
              test = "wilcox.test") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                    limits = c(0,1.5),
                     breaks = c(0,0.25,0.5,0.75,1),
                     expand = c(0.05,0.05)) +
  scale_x_discrete(labels = c('With specific HA-MDTPH',
                              'With specific LA-MDTPH',
                              'Without specific MDTPH')) +
  theme_classic(base_line_size = 0.8) +
  labs(title = "Resistant",x="",y="") +
  main_theme+
  theme(
    axis.text.x = element_text(angle = 15,size = 18,hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "none",
    plot.title = element_text(size = 18))
p2
p1/p2
ggsave(filename = "point.pdf",width = 5,height = 6.5)








# pie ---------------------------------------------------------------------

data_decrease <- df %>% 
  filter(prediction == "decrease") %>% 
  group_by(affect) %>%
  summarise(.,number = sum(num)) 


mylabs3 <- with(data_decrease,paste0(c("Affected","Unaffected"),"\n",round(number/sum(number),2)*100,"%","\n","N=",number))
data_decrease$mylabs <- mylabs3

p3 <- ggplot(data_decrease,aes(x = "",y = number,fill = as.character(affect))) +
  geom_bar(width = 1,stat = "identity",show.legend = F) +
  coord_polar("y",start = 0) +
  geom_text(aes(y =  number/2,
                label = mylabs),size = 3.5) +
  scale_fill_manual(values = c("#DD5555","#E6B4B5")) +
  labs(title = "High affinity score") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 15,
                              face = "bold",
                              hjust = 0.5)
  )

p3
data_nodec <- df %>% 
  filter(prediction == "no-decrease",target == "with") %>% 
  group_by(affect) %>%
  summarise(.,number = sum(num)) 

mylabs <- with(data_nodec,paste0(c("Affected","Unaffected"),"\n",round(number/sum(number),2)*100,"%","\n","N=",number))
p4 <- ggplot(data_nodec,aes(x = "",y = number,fill = as.character(affect))) +
  geom_bar(width = 1,stat = "identity",show.legend = F) +
  coord_polar("y",start = 0,direction = -1) +
  geom_text(aes(y = number/4 + c(0,cumsum(number)[-length(number)]),
                label = mylabs),
            size = 3.5) +
  scale_fill_manual(values = c("#5280A8","#C6DFF5")) +
  labs(title = "Low affinity score") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 15,
                              face = "bold",
                              hjust = 0.5)
  )


data_notar <- df %>% 
  filter(target == "without") %>% 
  group_by(affect) %>%
  summarise(.,number = sum(num))
mylabs2 <- with(data_notar,paste0(c("Affected","Unaffected"),"\n",round(number/sum(number),2)*100,"%","\n","N=",number))
p5 <- ggplot(data_notar,aes(x = "",y = number,fill = as.character(affect))) +
  geom_bar(width = 1,stat = "identity",show.legend = F) +
  coord_polar("y",start = 0) +
  geom_text(aes(y = number/3 + c(0,cumsum(number)[-length(number)]),
                label = mylabs2),size = 3.5) +
  scale_fill_manual(values = c("#E08627","#ECB67D")) +
  labs(title = "Without target") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 15,
                              face = "bold",
                              hjust = 0.5)
  )
p5

p3+p4+p5

ggsave(filename = "pieplot.pdf",width = 6,height = 6)




