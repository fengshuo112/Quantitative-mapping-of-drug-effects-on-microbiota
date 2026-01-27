arg <- commandArgs(T)
if(length(arg) < 3){
  cat("Argument: Input_path Group_File Out_path\n")
  quit('no')
}
pacman::p_load(tidyverse,MicrobiotaProcess,
               magrittr,ggsci,stringr,patchwork,reshape2,showtext)
showtext_auto()
font_add("Times New Roman","Times New Roman.ttf")
color <- c("#fbe183", "#f4c40f", "#fe9b00", "#d8443c", "#9b3441", "#de597c", 
           "#e87b89", "#e6a2a6", "#aa7aa1", "#9f5691", "#633372", "#1f6e9c", 
           "#2b9b81", "#92c051", "#bd3106", "#d9700e", "#e9a00e", "#eebe04", 
           "#5b7314", "#c3d6ce", "#89a6bb", "#454b87")


rename <- function(df,label){
  df <- df[!grepl("__un_",rownames(df)),]
  df$name <- unlist(lapply(rownames(df),function(x) str_match_all(x,paste0(label,"__([A-Za-z1-9\\[\\]-]+)"))[[1]][2]))
  dat_final <- aggregate(df[,c(1:(ncol(df)-1))],by=list(name=df$name),sum)
  dat_final <- column_to_rownames(dat_final,"name")
  return(dat_final)
}


#arg <- c("/Volumes/Flower0501/2022_16s/Mala/Group2/Analysis/0_input_data/","Pre_treatment-Post_treatment1-Post_treatment3-Post_treatment6","/Volumes/Flower0501/2022_16s/Mala/Group2/Analysis/1_Abundance/")
input_path <- arg[1]
otu <- paste0(input_path,"otu_table-rename.qza")
rep <- paste0(input_path,"rep-seqs-rename.qza")
tree <- paste0(input_path,"rooted-tree.qza")
tax <- paste0(input_path,"taxonomy.qza")
sample <- paste0(input_path,"sample-metadata.txt")
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeqza=tree)

# 数据的分布
ps <- ps_dada2
group <- ps_dada2@sam_data %>% as.data.frame()
#group <- read.table(arg[2],header = T,sep="\t",check.names = FALSE,comment.char = "")
#group$Group <- sub("Treat","Ros",group$Group)%>%sub("Con","Nor",.)
Gr <- unlist(str_split(arg[2],"-"))

# 分组门水平物种组成图
phytax <- get_taxadf(obj=ps_dada2, taxlevel=2)
df_phylum <- phytax@otu_table %>% as.data.frame()
df_phylum <- rename(df_phylum,"p")
df_Ph_plot <- df_phylum[order(rowSums(df_phylum),decreasing = T),][1:5,]
norm_ph <-  as.data.frame(t(df_Ph_plot)/colSums(df_Ph_plot,na=T))
norm_ph$sample <- rownames(norm_ph)
dd_ph <- cbind(norm_ph,group[rownames(norm_ph),"Group"])
dd_Ph2 <- melt(dd_ph)
colnames(dd_Ph2) <- c("sample","group","Phylum","value")
dd_Ph2$sample <- factor(dd_Ph2$sample,levels = str_sort(unique(dd_Ph2$sample),numeric = TRUE))
dd_Ph2$Phylum <- factor(dd_Ph2$Phylum,levels = dput(rownames(df_Ph_plot)))
dd_Ph2$group <- factor(dd_Ph2$group,levels = Gr)
phybar <- ggplot(dd_Ph2, aes(x=sample, y=value*100, fill=Phylum)) +
  labs(x='', y='Relative Abundance (%)',title = "Top 5 Phylum")+
  geom_bar(stat="identity") +
  facet_wrap(.~group,scales="free_x")+
  guides(fill=guide_legend(reverse=F)) +
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = c("#a40000", "#16317d", "#007e2f", "#ffcd12", "#b86092", 
                               "#721b3e", "#00b7a7"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(vjust = .6,hjust = .5,family = "Times New Roman",face= "bold",size=16),
        legend.text = element_text(family = "Times New Roman",size = 10),
        plot.margin = unit(rep(1,4),'lines'),
        legend.title = element_blank())
phybar
# 分组属水平物种组成图
genustax <- get_taxadf(obj=ps_dada2, taxlevel=6)
df_genus <- genustax@otu_table %>% as.data.frame()
df_genus <- rename(df_genus,"g")


df_ge_plot <- df_genus[order(rowSums(df_genus),decreasing = T),][1:15,]
norm_ge <- as.data.frame(t(df_ge_plot)/colSums(df_ge_plot,na=T))
norm_ge$sample <- rownames(norm_ge)
dd_g <- cbind(norm_ge,group[rownames(norm_ge),"Group"])
dd_g2 <- melt(dd_g)
colnames(dd_g2) <- c("sample","group","genus","value")
dd_g2$sample <- factor(dd_g2$sample,levels = str_sort(unique(dd_g2$sample),numeric = TRUE))
dd_g2$genus <- factor(dd_g2$genus,levels = dput(rownames(df_ge_plot)))
dd_g2$group <- factor(dd_g2$group,levels = Gr)
genusbar <- ggplot(dd_g2, aes(x=sample, y=value*100, fill=genus)) +
  labs(x='', y='Relative Abundance (%)',title = "Top 15 Genus")+
  geom_bar(stat="identity") +
  facet_wrap(.~group,scales="free_x")+
  guides(fill=guide_legend(reverse=F)) +
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = rev(color))+
  guides(fill = guide_legend(ncol = 2))+ 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=8, color="black",
                                   face="plain"),
        plot.title = element_text(vjust = .6,hjust = .5,family = "Times New Roman",face= "bold",size=16),
        legend.text = element_text(family = "Times New Roman",size = 10),
        plot.margin = unit(rep(1,4),'lines'),
        legend.title = element_blank())
genusbar
png(paste(arg[3],"Total.abundance.barplot.png",sep="/"),height = 4*480/2,width = 4*480/2,bg = "transparent")
phybar+genusbar+plot_layout(ncol = 1)
dev.off()

pdf(paste(arg[3],"Total.abundance.barplot.pdf",sep="/"),height = 10,width = 12,family = 'GB1')
phybar+genusbar+plot_layout(ncol = 1)
dev.off()
