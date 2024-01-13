setwd('~/CH/CH-zf/datu/bed/')

proj2 <- readRDS('~/CH/CH-zf/datu/Fig1/0_datu/proj2_datu_anno_new.rds')

gmat <- getMatrixFromProject(proj2,useMatrix = 'GeneScoreMatrix')
gene=getFeatures(proj2)
gene[1:5]

gmat=gmat@assays@data$GeneScoreMatrix
rownames(gmat)=gene

df=gmat

marker<-c('mylz3',#Cardiomyocyte
          'acta2',#Smooth muscle cell
          'col1a2',#Stromal cell
          
          'cd79a',#B cell
          'cd8a',#T cell
          'ccr7',#naive T
          'fcer1g',#Macrophage
          
          'opn1mw1',#Cone photoreceptor cell
          'gad1b',#GABAergic neuron
          'zic2a',#Granule cell 
          'snap25a',#Neural cell
          'gnal',#Olfactory neuron
          'olig1',#Oligodendrocyte
          #   'rdh5',#Pigment cell
          'her4.1',#Radial glia
          'cabp5a',#Retina bipolar cell
          'grk1a',#Rod photoreceptor cell
          
          'hbba1',#Erythroid
          
          'ddx4',#Primordial germ cell
          'pimr52',#Spermatocyte
          
          'kdr',#Endothelial cell
          
          'fabp6',#Enterocyte
          'chia.1',#Intestinal bulb cell
          'apoa1b',#Hepatocyte
          'slc22a2',#Nephron cell
          'cldna',#Epithelial cell
          'atp1a1a.2',#Ionocyte
          'ptgdsb.1',#Keratinocyte
          #  'chia.6',#Tuft cell,
          'ctrl',#Pancreatic cell
          
          'yif1b',#Granulosa cell
          'mmp13a'#Zona radiata cell
)

marker=marker[marker%in%rownames(df)]
marker

cell.types=unique(proj2$anno)
cell.types[!grepl('unknown',cell.types)]
cell.types=c('Cardiomyocyte','Smooth muscle cell','Stromal cell','B cell','T cell','Naive T cell','Macrophage','Cone photoreceptor cell','GABAergic neuron','Granule cell','Neural cell','Olfactory neuron','Oligodendrocyte','Radial glia','Retina bipolar cell','Rod photoreceptor cell','Erythroid','Primordial germ cell','Spermatocyte','Endothelial cell','Enterocyte','Intestinal bulb cell','Hepatocyte','Nephron cell','Epithelial cell','Ionocyte','Keratinocyte','Pancreatic cell','Granulosa cell','Zona radiata cell')

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
pdf(file = '~/CH/CH-zf/datu/Fig1/4_dotplot/marker_dotplot.pdf',width = 9,height = 6)

DotPlot(seob_use,features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order = 3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
dev.off()

