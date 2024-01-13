rm(list = ls())
library(data.table)
library(Seurat)
library(ArchR)
addArchRChrPrefix(chrPrefix = FALSE)
#library(LSD)
library(reshape2)
addArchRThreads(threads = 1)
library(BSgenome.Dmelanogaster.UCSC.dm6)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Dmelanogaster.UCSC.dm6)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Dmelanogaster.UCSC.dm6.ensGene, OrgDb = org.Dm.eg.db)
proj2 <- readRDS('~/CH/CH-dm/Fig1/0_datu/proj2_datu_anno.rds')

gmat <- getMatrixFromProject(proj2,useMatrix = 'GeneScoreMatrix')
gene=getFeatures(proj2)
gene[1:5]

gmat=gmat@assays@data$GeneScoreMatrix
rownames(gmat)=gene

df=gmat

marker<-c('alphaTry',#Anterior enterocyte
          'lambdaTry',#Posterior Enterocyte
          'lab',#Middle Enterocyte
          'E(spl)m3-HLH',#Enteroblast
          'Mur18B',#Malpighian tubule
          'srp',#Fatbody/Hemocyte
          'Spn43Ad',#Fatbody
          'trh',#Trachel cell
          'grh',#Epithelial
          'LysP',#Salivery gland
          'Act87E',#Intersegmental muscle
          'Act88F',#Indirect flight muscle
          'Act79B',#Jump muscle
          
          'repo',#Glial cell
          'pros',#Central nerve cell
          'SoxN',#T4/5 neuron
          'brp',#Lamina cortex neuron
          
          'BicC',#Germline cell
          
          'Sfp96F',#Male accessory gland main cell 
          'tj',#Follicle stem cell
          'yellow-g',#Follicle cell
          # 'Ilp8',#Late-staged follicle cell
          
          'FASN3'#Oenocyte
          #   'CCKLR-17D3'#Fatbody
)

marker=marker[marker%in%rownames(df)]
marker

cell.types=c('Anterior enterocyte','Posterior Enterocyte','Middle Enterocyte','Enteroblast','Malpighian tubule','Fatbody/Hemocyte','Fatbody','Tracheal cell','Epithelial cell','Salivary gland cell','Intersegmental muscle','Indirect flight muscle','Jump muscle','Glial cell','Central nerve cell','T4/T5 neuron','Lamina cortex neuron','Germline cell','Male accessory gland main cell','Follicle stem cell','Follicle cell','Oenocyte')

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
pdf(file = '~/CH/CH-dm/Fig1/4_dotplot/marker_dotplot.pdf',width = 9,height = 6)

DotPlot(seob_use,features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order = 3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
dev.off()
