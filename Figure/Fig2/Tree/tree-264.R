library(ggtree)
library(treeio)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(openxlsx)
library(dplyr)
library(tidyverse)
setwd('./Figure/Fig2/Tree')

tree <- read.tree("RAxML_bestTree-264.out")


dat3 <- read.xlsx("tableS1.xlsx",
                  startRow = 2,sheet = 1)
dat3 <- dat3 %>% filter(Target.protein == "DHFR")

group <- read_tsv("id.txt")

dat3 <- group %>% left_join(dat3, by = c("seq" = "Representative.sequence.ID")) 


tax <- read_tsv("supp_dock_DHFR.tab")
dat4 <- dat3 %>% left_join(tax,by=c("seq" = "rep"))

dat4 <- dat4 %>% mutate(Phylum = case_when(
  Phylum == "Bacteroidetes" ~ "Bacteroidetes",
  Phylum == "Firmicutes" ~ "Firmicutes",
  Phylum == "Proteobacteria" ~ "Proteobacteria",
  Phylum == "Actinobacteria" ~ "Actinobacteria",
  Phylum == "Verrucomicrobia" ~ "Verrucomicrobia",
  TRUE ~ "Others"
))
unique(dat4$seq[which(dat4$diff_zscore<=0)])
p <- ggtree(tree, layout="fan", open.angle=10, size=0.5)


 p1 <- rotate(p,node = 378)#379
 p1 <- rotate(p1,node = 516)#517
  p1+
   geom_tiplab(size = 0.5)


tree_df <- fortify(tree) 
dat1 <- p1[["data"]] %>% filter(isTip) %>% 
  select(label,y) %>% rename(node = y)
dat1 <- dat1 %>% mutate(nodes = ifelse(node %% 20 == 0,node,""))


p2 <- 
  p1 + 
  geom_fruit(
    data=dat1,
    geom=geom_text,
    mapping=aes(y=label, label = nodes),
    size = 2,
    offset=0.1,
    pwidth=0.04
  ) +
  geom_fruit(
    data=dat4,
    geom=geom_tile,
    mapping=aes(y=label, fill = newPhylum),
    size = 0.0005,
    offset=0.1,
    #pwidth=0.0001
    pwidth=0.3
    )+
  scale_fill_manual(values = c("#F1D69E", "#C75B4C", "#7A8DC6" ,"#B2E47A","#006666","#AA92C3"))
  
p2


dat2 <- dat4 %>% select(label,diff_zscore)
dat2$diff_zscore[which(dat2$diff_zscore> 3)] <- 3
dat2$diff_zscore[which(dat2$label=="Human")] <- 0
p3 <- p2 + 
  new_scale_fill() + 
  geom_fruit(
    data=dat2,
    geom=geom_tile,
    mapping=aes(y=label, fill=-diff_zscore),
    color = "gray80",
    size = 0.005,
    offset=0.1,   # The distance between external layers, default is 0.03 times of x range of tree.
    pwidth=0.3 # width of the external layer, default is 0.2 times of x range of tree.
  ) +
  scale_fill_gradient2(low="blue",mid = "white", high="red",
                      midpoint = 0,
                      limits = c(-1.3,1.3),
                      na.value = "blue",
                      name = "Affinity score"
                      #guide="none"
                      )


p3


dat4$abun <- 7 + log10(dat4$Average.relative.abundance)
p4 <- p3 + 
  new_scale_fill() +
  geom_fruit(
    data=dat4,
    geom=geom_col,
    mapping=aes(y=label, x=abun),  #The 'Abundance' of 'dat1' will be mapped to x
    fill = "#90BEA6",
    pwidth=0.2,
    offset=0.1,
    axis.params=list(
      axis="x", # add axis text of the layer.
      text.angle=-45, # the text size of axis.
      title = "Abundance",
      size = 1,
      title.height = 0,
      breaks = c(1,2,3),
      nbreak = 3, 
      labels = c(10^-8,10^-7,10^-6),
      hjust=0  # adjust the horizontal position of text of axis.
    ),
    grid.params=list() # add the grid line of the external bar plot.
  ) 

p4


p5 <- 
  p4 + 
  new_scale_fill() +
  geom_fruit(
    data = dat3,
    geom = geom_col,
    mapping = aes(y=label, x=Individual.count),
    pwidth=0.2,
    fill = "#90C2E3",
    offset=0.02,
    axis.params=list(
      axis="x", # add axis text of the layer.
      text.angle=-45, # the text size of axis.
      title = "Individuals",
      size = 10,
      title.height = 0,
      hjust=0  # adjust the horizontal position of text of axis.
    ),
    grid.params=list() # add the grid line of the external bar plot.
  ) 

p5


ggsave("DHFR-264.pdf", width = 10, height = 10)
 
