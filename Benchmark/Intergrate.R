#intergrating CH and sci kidney data using signac 
pmat=readRDS('./intergrated_kidney_pmat.rds')
pmat[1:5,1:5]
colnames(pmat)=reshape2::colsplit(colnames(pmat),'#',names = c('c1','c2'))$c2
library(Signac)
library(Seurat)
chrom_assay <- CreateChromatinAssay(
  counts = pmat,
  sep = c("-", "-"),fragments = './sort_inter_Kidney.tsv.gz'
)
load('intergrated_kidney_proj_peak.rda')#load proj
meta=proj2@cellColData%>%as.data.frame()
rownames(meta)=reshape2::colsplit(rownames(meta),'#',names = c('c1','c2'))$c2
obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = meta,
  names.field = 1)
table(obj$anno)
table(obj$Sample)

obj=RunTFIDF(obj)
obj=FindTopFeatures(obj,min.cutoff = 'q0')
obj=RunSVD(obj)
DepthCor(obj,n = 50)

obj=FindNeighbors(obj,reduction = 'lsi',dims = 2:15)
obj=FindClusters(obj,resolution = 0.8,verbose = F)
table(obj$seurat_clusters)

obj=RunUMAP(obj,reduction = 'lsi',dims = 2:15)

library(RColorBrewer)
library(ggsci)
DimPlot(obj,reduction = 'umap',label = T,cols = colorRampPalette(pal_npg(alpha = 0.6)(10))(35),pt.size = 1)+NoLegend()
DimPlot(obj,group.by = 'anno',label = T,reduction = 'umap',cols = colorRampPalette(pal_npg(alpha = 0.6)(10))(15),pt.size = 0.6)
DimPlot(obj,group.by = 'Sample',label = T,reduction = 'umap',cols = colorRampPalette(pal_npg(alpha = 0.6)(10))(15),pt.size = 0.6)


library(harmony)
table(obj$Sample)
meta_obj=obj@meta.data
meta_obj$platform=''
meta_obj[meta_obj$Sample=='ch_kidney',]$platform='CHATAC'
meta_obj[meta_obj$Sample=='sci_Kidney',]$platform='Sci1'
table(meta_obj$platform)
obj$platform=meta_obj$platform

obj=RunHarmony(obj,group.by.vars = 'Sample',reduction = 'lsi',assay.use = 'peaks',project.dim = F,lambda = 1.5)
obj=RunUMAP(obj,reduction = 'harmony',dims = 2:15)
obj=FindNeighbors(obj,reduction = 'harmony',dims = 2:15)
obj=FindClusters(obj,resolution = 0.8)
p=DimPlot(obj,group.by = 'platform',label = F,reduction = 'umap',cols = c('#20B2AA','#6495ED'),pt.size = 0.8)
p
ggsave(p,filename = 'intergrated_sci_ch_kidney.pdf',height = 10,width = 10)

saveRDS(obj,file = 'signac_intergrate_obj_new.rds')

p1=DimPlot(obj,reduction = 'umap',label = T,pt.size = 0.8)
p2=DimPlot(obj,reduction = 'umap',label = T,pt.size = 0.8,group.by = 'anno')
p1+p2
p2
