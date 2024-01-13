#gmat
setwd('~/CH-cross/cluster/Renbing/bed/')
proj <- readRDS('~/CH-cross/Fig3/species/hg38/hg38_proj_sample.rds')

gmat <- getMatrixFromProject(proj,useMatrix = 'GeneScoreMatrix')
gmat=gmat@assays@data$GeneScoreMatrix
gene=getFeatures(proj)
gene[1:5]
gmat[1:5,1:5]
rownames(gmat)=gene
gmat[1:5,1:5]
saveRDS(gmat,file = '~/CH-cross/Fig3/species/hg38/hg38_gmat.rds')



setwd('/media/ggj/Elements/sci-mouseATAC/sci1/')
proj <- readRDS('~/CH-cross/Fig3/species/mm9/mm9_proj_sample.rds')

gmat <- getMatrixFromProject(proj,useMatrix = 'GeneScoreMatrix')
gmat=gmat@assays@data$GeneScoreMatrix
gene=getFeatures(proj)
gene[1:5]
gmat[1:5,1:5]
rownames(gmat)=gene
gmat[1:5,1:5]
saveRDS(gmat,file = '~/CH-cross/Fig3/species/mm9/mm9_gmat.rds')


setwd('/media/ggj/Fhg38YT/CH/CH-cross/cluster/monkey/tissue/')
proj <- readRDS('~/CH-cross/Fig3/species/macaca/macaca_proj_sample.rds')

gmat <- getMatrixFromProject(proj,useMatrix = 'GeneScoreMatrix')
gmat=gmat@assays@data$GeneScoreMatrix
gene=getFeatures(proj)
gene=reshape2::colsplit(gene,':',names = c('c1','c2'))$c2
gene[1:5]
gmat[1:5,1:5]
rownames(gmat)=gene
gmat[1:5,1:5]
saveRDS(gmat,file = '~/CH-cross/Fig3/species/macaca/macaca_gmat.rds')


setwd('~/CH-zf/datu/bed/')
proj <- readRDS('~/CH-cross/Fig3/species/dr11/dr11_proj_sample.rds')

gmat <- getMatrixFromProject(proj,useMatrix = 'GeneScoreMatrix')
gmat=gmat@assays@data$GeneScoreMatrix
gene=getFeatures(proj)
gene=reshape2::colsplit(gene,':',names = c('c1','c2'))$c2
gene[1:5]
gmat[1:5,1:5]
rownames(gmat)=gene
gmat[1:5,1:5]
saveRDS(gmat,file = '~/CH-cross/Fig3/species/dr11/dr11_gmat.rds')


setwd('~/CH-dm/arrow/')
proj <- readRDS('~/CH-cross/Fig3/species/dm6/dm6_proj_sample.rds')

gmat <- getMatrixFromProject(proj,useMatrix = 'GeneScoreMatrix')
gmat=gmat@assays@data$GeneScoreMatrix
gene=getFeatures(proj)
gene=reshape2::colsplit(gene,':',names = c('c1','c2'))$c2
gene[1:5]
gmat[1:5,1:5]
rownames(gmat)=gene
gmat[1:5,1:5]
saveRDS(gmat,file = '~/CH-cross/Fig3/species/dm6/dm6_gmat.rds')



setwd('~/CH-ew/bed_total/')
proj <- readRDS('~/CH-cross/Fig3/species/ew/ew_proj_sample_20230304.rds')

gmat <- getMatrixFromProject(proj,useMatrix = 'GeneScoreMatrix')
gmat=gmat@assays@data$GeneScoreMatrix
gene=getFeatures(proj)

gene[1:5]
gmat[1:5,1:5]
rownames(gmat)=gene
gmat[1:5,1:5]
saveRDS(gmat,file = '~/CH-cross/Fig3/species/ew/ew_gmat_20230304.rds')


#h5ad
library(Seurat)
library(SeuratDisk)

#human
setwd('~/CH-cross/cluster/Renbing/bed/')
proj <- readRDS('~/CH-cross/Fig3/species/hg38/hg38_proj_sample.rds')
gmat=readRDS('~/CH-cross/Fig3/species/hg38/hg38_gmat.rds')
meta=proj@cellColData%>%as.data.frame()
meta=meta[,c('main_ct','lineage','cell_type')]
seob <-CreateSeuratObject(counts = gmat) 
SaveH5Seurat(seob,filename = '~/CH-cross/Fig3/species/hg38/hg38.h5Seurat')
Convert('~/CH-cross/Fig3/species/hg38/hg38.h5Seurat',dest = 'h5ad')
head(meta)
meta$cell=rownames(meta)
fwrite(meta,file = '~/CH-cross/Fig3/species/hg38/hg38_meta.csv',quote = F,row.names = F,sep = '\t')



#mouse
proj <- readRDS('~/CH-cross/Fig3/species/mm9/mm9_proj_sample.rds')
gmat=readRDS('~/CH-cross/Fig3/species/mm9/mm9_gmat.rds')
meta=proj@cellColData%>%as.data.frame()
meta=meta[,c('main_ct','lineage','cell_type')]
seob <-CreateSeuratObject(counts = gmat) 
SaveH5Seurat(seob,filename = '~/CH-cross/Fig3/species/mm9/mm9.h5Seurat')
Convert('~/CH-cross/Fig3/species/mm9/mm9.h5Seurat',dest = 'h5ad')
head(meta)
meta$cell=rownames(meta)
fwrite(meta,file = '~/CH-cross/Fig3/species/mm9/mm9_meta.csv',quote = F,row.names = F,sep = '\t')


#macaca
proj <- readRDS('~/CH-cross/Fig3/species/macaca/macaca_proj_sample.rds')

gmat=readRDS('~/CH-cross/Fig3/species/macaca/macaca_gmat.rds')
gmat[1:5,1:5]
rownames(gmat)=gsub('NA_','',rownames(gmat))

meta=proj@cellColData%>%as.data.frame()
meta=meta[,c('main_ct','lineage','cell_type')]
seob <-CreateSeuratObject(counts = gmat) 
SaveH5Seurat(seob,filename = '~/CH-cross/Fig3/species/macaca/macaca.h5Seurat')
Convert('~/CH-cross/Fig3/species/macaca/macaca.h5Seurat',dest = 'h5ad')
head(meta)
meta$cell=rownames(meta)
fwrite(meta,file = '~/CH-cross/Fig3/species/macaca/macaca_meta.csv',quote = F,row.names = F,sep = '\t')

#zebrafish
proj <- readRDS('~/CH-cross/Fig3/species/dr11/dr11_proj_sample.rds')

gmat=readRDS('~/CH-cross/Fig3/species/dr11/dr11_gmat.rds')
gmat[1:5,1:5]


meta=proj@cellColData%>%as.data.frame()
meta=meta[,c('main_ct','lineage','cell_type')]
seob <-CreateSeuratObject(counts = gmat) 
SaveH5Seurat(seob,filename = '~/CH-cross/Fig3/species/dr11/dr11.h5Seurat')
Convert('~/CH-cross/Fig3/species/dr11/dr11.h5Seurat',dest = 'h5ad')
head(meta)
meta$cell=rownames(meta)
fwrite(meta,file = '~/CH-cross/Fig3/species/dr11/dr11_meta.csv',quote = F,row.names = F,sep = '\t')


#dm
proj <- readRDS('~/CH-cross/Fig3/species/dm6/dm6_proj_sample.rds')

gmat=readRDS('~/CH-cross/Fig3/species/dm6/dm6_gmat.rds')
gmat[1:5,1:5]


meta=proj@cellColData%>%as.data.frame()
meta=meta[,c('main_ct','lineage','cell_type')]
seob <-CreateSeuratObject(counts = gmat) 
SaveH5Seurat(seob,filename = '~/CH-cross/Fig3/species/dm6/dm6.h5Seurat')
Convert('~/CH-cross/Fig3/species/dm6/dm6.h5Seurat',dest = 'h5ad')
head(meta)
meta$cell=rownames(meta)
fwrite(meta,file = '~/CH-cross/Fig3/species/dm6/dm6_meta.csv',quote = F,row.names = F,sep = '\t')


#ew
proj <- readRDS('~/CH-cross/Fig3/species/ew/ew_proj_sample_20230304.rds')

gmat=readRDS('~/CH-cross/Fig3/species/ew/ew_gmat.rds')
gmat[1:5,1:5]


meta=proj@cellColData%>%as.data.frame()
meta=meta[,c('main_ct','lineage','cell_type')]
seob <-CreateSeuratObject(counts = gmat) 
SaveH5Seurat(seob,filename = '~/CH-cross/Fig3/species/ew/ew.h5Seurat')
Convert('~/CH-cross/Fig3/species/ew/ew.h5Seurat',dest = 'h5ad')
head(meta)
meta$cell=rownames(meta)
fwrite(meta,file = '~/CH-cross/Fig3/species/ew/ew_meta.csv',quote = F,row.names = F,sep = '\t')
