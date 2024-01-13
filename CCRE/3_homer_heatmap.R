rm(list = ls())

load('~/CH/CH-zf/datu/Fig1/0_datu/proj2_datu_anno_pmat.rda')

setwd('~/CH/CH-zf/datu/Fig2/homer/da/')
peak=proj2@peakSet
lineage.ls=sort(unique(peak@ranges@NAMES))

motif_list=list()

for (i in lineage.ls) {
  if(!grepl('unknown',i)){
    i=gsub(' ','_',i)
    i=gsub('/','-',i)
    motif=fread(paste0(i,'/knownResults.txt'))
    motif$lineage=i
    motif_list[[i]]=motif
  }
}



for (i in names(motif_list)) {
  tmp=motif_list[[i]]
  tmp$motif_name=reshape2::colsplit(tmp$`Motif Name`,'[(]',names = c('c1','c2'))$c1
  motif_list[[i]]=tmp
}

df=data.frame()
for (i in names(motif_list)) {
  tmp=motif_list[[i]]
  tmp=tmp[1:3,]
  tmp$score=-tmp$`Log P-value`
  tmp=tmp[,c('motif_name','score')]
  tmp$ct=i
  df=rbind(df,tmp)
}

df$motif_name=reshape2::colsplit(df$motif_name,'/',names = c('c1','c2'))$c1
df$score=as.numeric(df$score)
d <- df%>%  select("motif_name", "ct", "score")%>%dcast(motif_name~ct,value.var="score",fun.aggregate = sum)

d <- as.data.frame(d)
rownames(d)=d$motif_name
d=d[,-1]

colnames(d)=gsub('_',' ',colnames(d))
ct=c('B cell','T cell','Granulocyte','Macrophage','Erythroid','Rod photoreceptor cell','Retina bipolar cell','Radial glia','Cone photoreceptor cell','GABAergic neuron','Granule cell','Olfactory neuron','Oligodendrocyte','Neural cell','Enterocyte','Intestinal bulb cell','Hepatocyte','Nephron cell','Epithelial cell','Keratinocyte','Ionocyte','Pancreatic cell','Smooth muscle cell','Cardiomyocyte','Stromal cell','Spermatocyte','Primordial germ cell','Endothelial cell','Granulosa cell','Zona radiata cell')
unique(ct)

df=as.data.frame(df)
df$ct=gsub('_',' ',df$ct)
df$ct=factor(df$ct,levels = ct)

df=df[order(df$ct),]

use.motif=unique(df$motif_name)
d=d[,ct]
d=d[use.motif,]
d[is.na(d)]=0

library(pheatmap)
library(viridis)
p=pheatmap::pheatmap(d,
                     color = viridis(100,option = 'F'),
                     scale = 'none',
                     border_color = 'NA',
                     cluster_rows = F,
                     cluster_cols = F,
                     fontsize = 15,
                     width = 3,
                     height = 3)
p