rm(list = ls())
setwd('~/CH/CH-ew/bed_total/')
library(Seurat)
library(ArchR)
load('~/CH/CH-ew/res/ew_archR.rda')
setwd('~/CH/CH-ew/bed_total/')
addArchRThreads(threads =12)

proj2=readRDS('~/CH/CH-ew/Fig1/0_datu/proj2_datu_lineage_20230304.rds')
table(proj2$lineage)

gmat <- getMatrixFromProject(proj2,useMatrix = 'GeneScoreMatrix')
gene=getFeatures(proj2)
gene[1:5]

gmat=gmat@assays@data$GeneScoreMatrix
rownames(gmat)=gene

df=gmat

marker  <- c(
  'evm.TU.Chr02.2009',#Muscle cell
  'evm.TU.Chr10.1480',#epidermis_cell
  'evm.TU.Chr04.1255', #Cuticle_epidermis_cell
  'evm.TU.Chr01.933',#Mechanosensory cell
  'evm.TU.Chr02.1450',#Digestive_gland_cell
  "evm.TU.Chr04.2499", #Coelomocyte
  'evm.TU.Chr01.1402',#Erythrocruorin
  'evm.TU.Chr10.1273',#neural cell 
  'evm.TU.Chr04.1542',#Sperm_stem_cell
  'evm.TU.Chr01.3449',#Spermatids-Spermatogenesis
  'evm.TU.Chr03.1056'#Spermatogonia
)


marker=marker[marker%in%rownames(df)]
marker

meta=as.data.frame(proj2@cellColData)
meta$anno=gsub('_',' ',meta$anno)
table(meta$anno)

cell.types=c('Muscle cell','Epidermal cell','Cuticle epidermis cell','Mechanosensory cell','Digestive gland cell','Coelomocyte','Erythrocruorin','Neuron cell','Sperm stem cell','Spermatids-Spermatogenesis','Spermatogonia')
cell.use=rownames(meta[meta$anno%in%cell.types,])

meta1=meta[cell.use,]
df1=df[marker,]
df1=df1[,cell.use]

library(Seurat)
seob_use <- CreateSeuratObject(df1,meta.data = meta1)
Idents(seob_use) <- 'anno'
seob_use@active.ident <- factor(seob_use@active.ident,levels = rev(cell.types))
seob_use=NormalizeData(seob_use)
seob_use@assays$RNA@data[1:5,1:5]

dir.create('~/CH/CH-ew/Fig1/2_dotplot')
pdf(file = '~/CH/CH-ew/Fig1/2_dotplot/marker_dotplot.pdf',width = 9,height = 6)

DotPlot(seob_use,features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order = 3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
dev.off()