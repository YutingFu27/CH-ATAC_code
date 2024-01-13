#human based on zf
setwd('~/CH/CH-cross/Fig3/liftover/')
hg38_peak_align=fread('./motif/hg38_align_zf.bed')
hg38_peak_align=paste(hg38_peak_align$V1,hg38_peak_align$V2,hg38_peak_align$V3,sep = '-')
hg38_peak_align[1:5]


dr11_peak=fread('./motif/hg38Todr11_peak.bed')
dr11_peak=paste(dr11_peak$V1,dr11_peak$V2,dr11_peak$V3,sep = '-')
trans_hg38_dr11=data.frame(hg38=hg38_peak_align,dr11=dr11_peak)


hg38_peak_total=fread('./hg38_peak_0918.bed')
head(hg38_peak_total)
hg38_peak_total$V4=paste(hg38_peak_total$V1,hg38_peak_total$V2,hg38_peak_total$V3,sep = '-')

hg38_peak_total$tag='unalign'
hg38_peak_total[hg38_peak_total$V4%in%hg38_peak_align,]$tag='align'

table(hg38_peak_total$tag)

hg38_conservezf=fread('./motif/hg38_conservezf.bed')
hg38_conservezf=paste(hg38_conservezf$V1,hg38_conservezf$V2,hg38_conservezf$V3,sep = '-')
hg38_conservezf=unique(hg38_conservezf)

library(dplyr)

hg38_conservezf_1=trans_hg38_dr11[trans_hg38_dr11$dr11%in%hg38_conservezf,]$hg38%>%unique()


hg38_peak_total[hg38_peak_total$V4%in%hg38_conservezf_1,]$tag='conserve'
table(hg38_peak_total$tag)

hg38_peak_unalign=hg38_peak_total[hg38_peak_total$tag=='unalign',]
hg38_peak_align=hg38_peak_total[hg38_peak_total$tag=='align',]
hg38_peak_conserve=hg38_peak_total[hg38_peak_total$tag=='conserve',]

fwrite(hg38_peak_unalign[,1:4],quote = F,row.names = F,col.names = F,sep = '\t',file = './motif/human_phy_basedzf/hg38_unalign.bed')
fwrite(hg38_peak_align[,1:4],quote = F,row.names = F,col.names = F,sep = '\t',file = './motif/human_phy_basedzf/hg38_align.bed')
fwrite(hg38_peak_conserve[,1:4],quote = F,row.names = F,col.names = F,sep = '\t',file = './motif/human_phy_basedzf/hg38_conserve.bed')

#cor
wo=fread('./motif/dr11_conservehuman_SF_wo.bed')
wo$hpeak=paste(wo$V4,wo$V5,wo$V6,sep = '-')
wo$zpeak=paste(wo$V1,wo$V2,wo$V3,sep = '-')
wo=wo[!duplicated(wo$hpeak),]
wo=wo[!duplicated(wo$zpeak),]
# SF_h=fread('./motif/human_phy/hg38_conserve.bed')
# SF_h=SF_h$V4
SF_h=wo$hpeak

dr11_peak_align=fread('./motif/dr11_align.bed')
dr11_peak_align=paste(dr11_peak_align$V1,dr11_peak_align$V2,dr11_peak_align$V3,sep = '-')
dr11_peak_align[1:5]


hg38_peak=fread('./motif/dr11toHg38_peak.bed')
hg38_peak=paste(hg38_peak$V1,hg38_peak$V2,hg38_peak$V3,sep = '-')
trans_dr11_hg38=data.frame(dr11=dr11_peak_align,hg38=hg38_peak)
trans_dr11_hg38=trans_dr11_hg38[!duplicated(trans_dr11_hg38$dr11),]
trans_dr11_hg38=trans_dr11_hg38[!duplicated(trans_dr11_hg38$hg38),]

dr11_conservehuman_1=trans_dr11_hg38[trans_dr11_hg38$hg38%in%wo$zpeak,]$dr11%>%unique()

SF_z=unique(dr11_conservehuman_1)



pmat_h=readRDS('~/CH/CH-cross/Fig3/species/intergrate_hg38/pmat_0306.rds')
dim(pmat_h)

pmat_z=readRDS('~/CH/CH-cross/Fig3/species/dr11/pmat_0307.rds')
dim(pmat_z)


library(ArchR)
proj_h <- readRDS('~/CH/CH-cross/Fig3/species/intergrate_hg38/hg38_intergrate_bigFigure_0306.rds')
proj_z <- readRDS('~/CH/CH-cross/Fig3/species/dr11/dr11_bigFigure_20230302.rds')


SF_h=gsub('-','_',SF_h)
SF_z=gsub('-','_',SF_z)

pmat_h_use=pmat_h[SF_h,]
dim(pmat_h_use)

pmat_z_use=pmat_z[SF_z,]
dim(pmat_z_use)

#seob
meta_h=proj_h@cellColData%>%as.data.frame()
meta_h=meta_h[colnames(pmat_h_use),]

meta_z=proj_z@cellColData%>%as.data.frame()
meta_z=meta_z[colnames(pmat_z_use),]

library(Seurat)
seob_h=CreateSeuratObject(pmat_h_use,meta.data = meta_h)
avpr_h = as.data.frame(AverageExpression(seob_h,group.by = 'main_ct')[[1]])
dim(avpr_h)
avpr_h[1:5,1:5]


seob_z=CreateSeuratObject(pmat_z_use,meta.data = meta_z)
avpr_z = as.data.frame(AverageExpression(seob_z,group.by = 'main_ct')[[1]])
dim(avpr_z)
avpr_z[1:5,1:5]
ct=c('Cardiomyocyte','Endothelial','Enterocyte','Hepatocyte','Oligodendrocyte','Proximal tuble','Stromal cell','T cell')




avpr_h_use=avpr_h
avpr_z_use=avpr_z

#Cardiomyocyte
h_counts=avpr_h_use[,'Cardiomyocyte']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

z_counts=avpr_z_use[,'Cardiomyocyte']%>%as.data.frame()
z_norm = z_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "zebrafish" = z_norm
) %>%
  filter(human != 0, zebrafish != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)
library(ggplot2)
#source('./plot_utils.R')
p1 = ggplot(tpm, aes(human, zebrafish)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Cardiomyocyte')
p1
ggsave(p1,filename = './Cardiomyocyte_pearson_SF.pdf',height = 15,width = 15)
#Endothelial
h_counts=avpr_h_use[,'Endothelial']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

z_counts=avpr_z_use[,'Endothelial']%>%as.data.frame()
z_norm = z_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "zebrafish" = z_norm
) %>%
  filter(human != 0, zebrafish != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

p1 = ggplot(tpm, aes(human, zebrafish)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Endothelial')
p1
ggsave(p1,filename = './Endothelial_pearson_SF.pdf',height = 15,width = 15)

#Enterocyte
h_counts=avpr_h_use[,'Enterocyte']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

z_counts=avpr_z_use[,'Enterocyte']%>%as.data.frame()
z_norm = z_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "zebrafish" = z_norm
) %>%
  filter(human != 0, zebrafish != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

p1 = ggplot(tpm, aes(human, zebrafish)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Enterocyte')
p1
ggsave(p1,filename = './Enterocyte_pearson_SF.pdf',height = 15,width = 15)

#Hepatocyte
h_counts=avpr_h_use[,'Hepatocyte']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

z_counts=avpr_z_use[,'Hepatocyte']%>%as.data.frame()
z_norm = z_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "zebrafish" = z_norm
) %>%
  filter(human != 0, zebrafish != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

p1 = ggplot(tpm, aes(human, zebrafish)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Hepatocyte')
p1
ggsave(p1,filename = './Hepatocyte_pearson_SF.pdf',height = 15,width = 15)

#Oligodendrocyte
h_counts=avpr_h_use[,'Oligodendrocyte']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

z_counts=avpr_z_use[,'Oligodendrocyte']%>%as.data.frame()
z_norm = z_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "zebrafish" = z_norm
) %>%
  filter(human != 0, zebrafish != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

p1 = ggplot(tpm, aes(human, zebrafish)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Oligodendrocyte')
p1
ggsave(p1,filename = './Oligodendrocyte_pearson_SF.pdf',height = 15,width = 15)

#Proximal tuble
h_counts=avpr_h_use[,'Proximal tuble']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

z_counts=avpr_z_use[,'Nephron cell']%>%as.data.frame()
z_norm = z_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "zebrafish" = z_norm
) %>%
  filter(human != 0, zebrafish != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

p1 = ggplot(tpm, aes(human, zebrafish)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Proximal tuble')
p1
ggsave(p1,filename = './Proximal tuble_pearson_SF.pdf',height = 15,width = 15)

#Stromal cell
h_counts=avpr_h_use[,'Stromal cell']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

z_counts=avpr_z_use[,'Stromal cell']%>%as.data.frame()
z_norm = z_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "zebrafish" = z_norm
) %>%
  filter(human != 0, zebrafish != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

p1 = ggplot(tpm, aes(human, zebrafish)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'Stromal cell')
p1
ggsave(p1,filename = './Stromal_cell_pearson_SF.pdf',height = 15,width = 15)

#T cell
h_counts=avpr_h_use[,'T cell']%>%as.data.frame()
h_norm = h_counts %>% rowSums() %>% {1e6 *. / sum(.)}

z_counts=avpr_z_use[,'T cell']%>%as.data.frame()
z_norm = z_counts %>% rowSums() %>% {1e6 *. / sum(.)}


library(tidyverse)
tpm = tibble(
  "human" = h_norm,
  "zebrafish" = z_norm
) %>%
  filter(human != 0, zebrafish != 0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

p1 = ggplot(tpm, aes(human, zebrafish)) + stat_scatter_density() + scale_color_viridis_c(option = 'D') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~"))) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(title = 'T cell')
p1
ggsave(p1,filename = './T_cell_pearson_SF.pdf',height = 15,width = 15)
