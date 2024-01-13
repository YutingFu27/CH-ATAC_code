setwd('~/CH/CH-cross/Fig3/liftover/')

library(data.table)
#mm
mm9_peak_align=fread('./motif/mm9_align.bed')
mm10_peak=fread('./motif/mm10_peak.bed')

trans_mm9_mm10=data.frame('mm9'=paste(mm9_peak_align$V1,mm9_peak_align$V2,mm9_peak_align$V3,sep = '-'),'mm10'=paste(mm10_peak$V1,mm10_peak$V2,mm10_peak$V3,sep = '-'))
rownames(trans_mm9_mm10)=trans_mm9_mm10$mm10

mm10_peak_align=fread('./motif/mm10_align.bed')
mm10_peak_align=paste(mm10_peak_align$V1,mm10_peak_align$V2,mm10_peak_align$V3,sep = '-')
mm10_peak_align[1:5]
trans_mm9_mm10[mm10_peak_align,]$mm9->mm9_align

hg38_peak=fread('./motif/mm10toHg38_peak.bed')
hg38_peak=paste(hg38_peak$V1,hg38_peak$V2,hg38_peak$V3,sep = '-')
trans_mm10_hg38=data.frame(mm10=mm10_peak_align,hg38=hg38_peak)
#trans_mm10_hg38=trans_mm10_hg38[!duplicated(trans_mm10_hg38$hg38),]
#rownames(trans_mm10_hg38)=trans_mm10_hg38$hg38

mm9_peak_total=fread('./mm9_peak.bed')
head(mm9_peak_total)
mm9_peak_total$V4=paste(mm9_peak_total$V1,mm9_peak_total$V2,mm9_peak_total$V3,sep = '-')

mm9_peak_total$tag='unalign'
mm9_peak_total[mm9_peak_total$V4%in%mm9_align,]$tag='align'

table(mm9_peak_total$tag)

mm10_conservehuman=fread('./motif/mm10_conservehuman.bed')
mm10_conservehuman=paste(mm10_conservehuman$V1,mm10_conservehuman$V2,mm10_conservehuman$V3,sep = '-')
mm10_conservehuman=unique(mm10_conservehuman)

library(dplyr)
intersect(rownames(trans_mm10_hg38),mm10_conservehuman)%>%length()
mm10_conservehuman_1=trans_mm10_hg38[trans_mm10_hg38$hg38%in%mm10_conservehuman,]$mm10%>%unique()

intersect(mm10_conservehuman_1,rownames(trans_mm9_mm10))%>%length()
mm9_conservehuman=trans_mm9_mm10[mm10_conservehuman_1,]$mm9


mm9_peak_total[mm9_peak_total$V4%in%mm9_conservehuman,]$tag='conserve'
table(mm9_peak_total$tag)
mm9_peak_unalign=mm9_peak_total[mm9_peak_total$tag=='unalign',]
mm9_peak_align=mm9_peak_total[mm9_peak_total$tag=='align',]
mm9_peak_conserve=mm9_peak_total[mm9_peak_total$tag=='conserve',]
fwrite(mm9_peak_unalign[,1:4],quote = F,row.names = F,col.names = F,sep = '\t',file = './motif/mouse_phy/mm9_unalign.bed')
fwrite(mm9_peak_align[,1:4],quote = F,row.names = F,col.names = F,sep = '\t',file = './motif/mouse_phy/mm9_align.bed')
fwrite(mm9_peak_conserve[,1:4],quote = F,row.names = F,col.names = F,sep = '\t',file = './motif/mouse_phy/mm9_conserve.bed')

#phy
library(readr)
setwd('~/CH/CH-cross/Fig3/liftover/motif/mouse_phy/')
phy.ls=list.files(pattern = 'phyres')
unalign_phy.ls='./mm10_unalign.phyres'
align_phy.ls='./mm10_align.phyres'
conserve_phy.ls='./mm10_conserve.phyres'

phy_unalign <- read_delim('./mm10_unalign.phyres',delim = " ",col_names = F)
colnames(phy_unalign)=c('peak','phy_score')
phy_unalign$class='N'

phy_align <- read_delim('./mm10_align.phyres',delim = " ",col_names = F)
colnames(phy_align)=c('peak','phy_score')
phy_align$class='S'

phy_conserve <- read_delim('./mm10_conserve.phyres',delim = " ",col_names = F)
colnames(phy_conserve)=c('peak','phy_score')
phy_conserve$class='SF'



library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
library(ggsignif)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(ggpubr)
all=rbind(phy_unalign,phy_align,phy_conserve)
# head(all)
df_p_val1 <- all %>%
  t_test(formula = phy_score ~ class) %>%
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>%
  add_xy_position(x='class')

p <- ggboxplot(all, x = "class", y = "phy_score",  color = "black", fill = 'class',
               width = 0.5, show.legend = F, bxp.errorbar = T,outlier.shape = NA) +theme_classic()+geom_signif(
                 comparisons = list(c('N','S'),c('S','SF'),c('N','SF')),y_position = c(1,1.5,1.8),map_signif_level = T)+ylim(0,2.5)+scale_fill_manual(values = c('#00AFBB', '#E7B800','#6495ED'))
p


p1 <- p +xlab('')+ylab('phyloP60 score')
p1
ggsave(p1,filename = './mouse_phy.pdf',height = 15,width = 15)

