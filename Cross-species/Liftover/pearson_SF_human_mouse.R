setwd('~/CH/CH-cross/Fig3/liftover/')

wo=fread('./motif/mm10_conservehuman_SF_wo.bed')
wo$hpeak=paste(wo$V4,wo$V5,wo$V6,sep = '-')
wo$mpeak=paste(wo$V1,wo$V2,wo$V3,sep = '-')
wo=wo[!duplicated(wo$hpeak),]
wo=wo[!duplicated(wo$mpeak),]
# SF_h=fread('./motif/human_phy/hg38_conserve.bed')
# SF_h=SF_h$V4
SF_h=wo$hpeak

library(data.table)
#mm
mm9_peak_align=fread('./motif/mm9_align.bed')
mm10_peak=fread('./motif/mm10_peak.bed')

trans_mm9_mm10=data.frame('mm9'=paste(mm9_peak_align$V1,mm9_peak_align$V2,mm9_peak_align$V3,sep = '-'),'mm10'=paste(mm10_peak$V1,mm10_peak$V2,mm10_peak$V3,sep = '-'))
rownames(trans_mm9_mm10)=trans_mm9_mm10$mm10

mm10_peak_align=fread('./motif/mm10_align.bed')
mm10_peak_align=paste(mm10_peak_align$V1,mm10_peak_align$V2,mm10_peak_align$V3,sep = '-')
mm10_peak_align[1:5]


hg38_peak=fread('./motif/mm10toHg38_peak.bed')
hg38_peak=paste(hg38_peak$V1,hg38_peak$V2,hg38_peak$V3,sep = '-')
trans_mm10_hg38=data.frame(mm10=mm10_peak_align,hg38=hg38_peak)
trans_mm10_hg38=trans_mm10_hg38[!duplicated(trans_mm10_hg38$hg38),]
rownames(trans_mm10_hg38)=trans_mm10_hg38$hg38


mm10_conservehuman_1=trans_mm10_hg38[trans_mm10_hg38$hg38%in%wo$mpeak,]$mm10%>%unique()

intersect(mm10_conservehuman_1,rownames(trans_mm9_mm10))%>%length()
mm9_conservehuman=trans_mm9_mm10[mm10_conservehuman_1,]$mm9
SF_m=mm9_conservehuman


pmat_m=readRDS('~/CH/CH-cross/Fig3/species/mm9/pmat.rds')
dim(pmat_m)


pmat_h=readRDS('~/CH/CH-cross/Fig3/species/intergrate_hg38/pmat_0306.rds')
dim(pmat_h)


pmat_h[1:5,1:5]

library(ArchR)
proj_h <- readRDS('~/CH/CH-cross/Fig3/species/intergrate_hg38/hg38_intergrate_bigFigure_0306.rds')
proj_mo <- readRDS('~/CH/CH-cross/Fig3/species/mm9/mm9_bigFigure.rds')


SF_h=gsub('-','_',SF_h)
SF_m=gsub('-','_',SF_m)

pmat_h_use=pmat_h[SF_h,]
dim(pmat_h_use)


# SF_m_old=fread('./motif/mouse_phy/mm9_conserve.bed')
# SF_m_old=SF_m_old$V4
# SF_m_old=gsub('-','_',SF_m_old)
# table(SF_m%in%SF_m_old)

pmat_m_use=pmat_m[SF_m,]
dim(pmat_m_use)

#seob
meta_h=proj_h@cellColData%>%as.data.frame()
meta_h=meta_h[colnames(pmat_h_use),]

meta_m=proj_mo@cellColData%>%as.data.frame()
meta_m=meta_m[colnames(pmat_m_use),]

library(Seurat)
seob_h=CreateSeuratObject(pmat_h_use,meta.data = meta_h)
avpr_h = as.data.frame(AverageExpression(seob_h,group.by = 'main_ct')[[1]])
dim(avpr_h)
avpr_h[1:5,1:5]


seob_m=CreateSeuratObject(pmat_m_use,meta.data = meta_m)
avpr_m = as.data.frame(AverageExpression(seob_m,group.by = 'main_ct')[[1]])
dim(avpr_m)
avpr_m[1:5,1:5]
ct=c('Cardiomyocyte','Endothelial','Enterocyte','Hepatocyte','Oligodendrocyte','Proximal tuble','Stromal cell','T cell')
intersect(colnames(avpr_h),c('Cardiomyocyte','Endothelial','Enterocyte','Hepatocyte','Oligodendrocyte','Proximal tuble','Stromal cell','T cell','Nephron cell'))
intersect(colnames(avpr_m),c('Cardiomyocyte','Endothelial','Enterocyte','Hepatocyte','Oligodendrocyte','Proximal tuble','Stromal cell','T cell','Nephron cell'))

colnames(avpr_h)[1:5]

avpr_h_use=avpr_h[,ct]
avpr_m_use=avpr_m[,ct]

#Cardiomyocyte
h_counts=avpr_h_use[,'Cardiomyocyte']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

m_counts=avpr_m_use[,'Cardiomyocyte']%>%as.data.frame()
m_norm = m_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "mouse" = m_norm
) %>%
  filter(human != 0, mouse != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)
library(ggplot2)
source('./plot_utils.R')
p1 = ggplot(tpm, aes(human, mouse)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Cardiomyocyte')
p1
ggsave(p1,filename = './Cardiomyocyte_pearson_SF.pdf',height = 15,width = 15)

#Endothelial
h_counts=avpr_h_use[,'Endothelial']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

m_counts=avpr_m_use[,'Endothelial']%>%as.data.frame()
m_norm = m_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "mouse" = m_norm
) %>%
  filter(human != 0, mouse != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)
library(ggplot2)
source('./plot_utils.R')
p1 = ggplot(tpm, aes(human, mouse)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Endothelial')
p1
ggsave(p1,filename = './Endothelial_pearson_SF.pdf',height = 15,width = 15)

#Enterocyte
h_counts=avpr_h_use[,'Enterocyte']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

m_counts=avpr_m_use[,'Enterocyte']%>%as.data.frame()
m_norm = m_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "mouse" = m_norm
) %>%
  filter(human != 0, mouse != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))


p1 = ggplot(tpm, aes(human, mouse)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Enterocyte')
p1
ggsave(p1,filename = './Enterocyte_pearson_SF.pdf',height = 15,width = 15)

#Hepatocyte
h_counts=avpr_h_use[,'Hepatocyte']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

m_counts=avpr_m_use[,'Hepatocyte']%>%as.data.frame()
m_norm = m_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "mouse" = m_norm
) %>%
  filter(human != 0, mouse != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))


p1 = ggplot(tpm, aes(human, mouse)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Hepatocyte')
p1
ggsave(p1,filename = './Hepatocyte_pearson_SF.pdf',height = 15,width = 15)

#Oligodendrocyte
h_counts=avpr_h_use[,'Oligodendrocyte']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

m_counts=avpr_m_use[,'Oligodendrocyte']%>%as.data.frame()
m_norm = m_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "mouse" = m_norm
) %>%
  filter(human != 0, mouse != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))


p1 = ggplot(tpm, aes(human, mouse)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Oligodendrocyte')
p1
ggsave(p1,filename = './Oligodendrocyte_pearson_SF.pdf',height = 15,width = 15)


#Proximal tuble
h_counts=avpr_h_use[,'Proximal tuble']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

m_counts=avpr_m_use[,'Proximal tuble']%>%as.data.frame()
m_norm = m_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "mouse" = m_norm
) %>%
  filter(human != 0, mouse != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))


p1 = ggplot(tpm, aes(human, mouse)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Proximal tuble')
p1
ggsave(p1,filename = './Proximal_tuble_pearson_SF.pdf',height = 15,width = 15)

#T cell
h_counts=avpr_h_use[,'T cell']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

m_counts=avpr_m_use[,'T cell']%>%as.data.frame()
m_norm = m_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "mouse" = m_norm
) %>%
  filter(human != 0, mouse != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))


p1 = ggplot(tpm, aes(human, mouse)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'T cell')
p1
ggsave(p1,filename = './T_cell_pearson_SF.pdf',height = 15,width = 15)

