#compare uq
setwd('~/CH/ch-mouse/compare_kidney/compare_20230907')
library(ArchR)
load('./intergrated_kidney_proj_peak.rda')
proj2
meta=proj2@cellColData%>%as.data.frame()
meta$platform=''
meta[meta$Sample=='ch_kidney',]$platform='CHATAC'
meta[meta$Sample=='sci_Kidney',]$platform='Sci1'
table(meta$platform)

library(ggpubr)
library(ggsci)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggtext)
uq_df=meta[,c('nFrags','platform')]
colnames(uq_df)[1]='uq'

df_p_val1 <- uq_df %>% 
  t_test(formula = uq ~ platform) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='platform')
df_p_val1
p=ggboxplot(uq_df, x = "platform", y = "uq", color = "black", fill = 'platform',title = '',
            width = 0.5, show.legend = F, bxp.errorbar = T,outlier.shape = NA)+scale_fill_manual(values = c('#00AFBB', '#E7B800'))+ylab('Unique fragment')+ylim(0,40000)

p
ggsave(p,filename = './uq_sci_ch_kidney.pdf',height = 15,width = 15)

tss_df=meta[,c('TSSEnrichment','platform')]
colnames(tss_df)[1]='tss'

df_p_val1 <- tss_df %>% 
  t_test(formula = tss ~ platform) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='platform')
df_p_val1
p=ggboxplot(tss_df, x = "platform", y = "tss", color = "black", fill = 'platform',title = '',
            width = 0.5, show.legend = F, bxp.errorbar = T,outlier.shape = NA)+scale_fill_manual(values = c('#00AFBB', '#E7B800'))+ylab('TSS enrichment score')+stat_pvalue_manual(df_p_val1,label = '{p.signif}',tip.length = 0.01,y.position = 30)

p
ggsave(p,filename = './tss_sci_ch_kidney.pdf',height = 15,width = 15)


p1+p2

#IDR
#idr --samples ch_peaks_sort.narrowPeak sci_peak_sort.narrowPeak --input-file-type narrowPeak --rank p.value --output-file compare_idr --plot
