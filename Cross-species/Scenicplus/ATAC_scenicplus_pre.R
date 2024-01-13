rm(list = ls())

library(data.table)
library(dplyr)
library(Matrix)
dir.create('~/CH/CH-cross/Fig3/Scenicplus/hg38')
setwd('~/CH/CH-cross/Fig3/Scenicplus/hg38/')
pmat <- readRDS('~/CH/CH-cross/Fig3/species/intergrate_hg38/pmat_0306.rds')
pmat[1:5,1:5]

savepath='~/CH/CH-cross/Fig3/Scenicplus/hg38/'

dim(pmat)
chr=reshape2::colsplit(rownames(pmat),'_',names=c('c1','c2','c3'))$c1
start=reshape2::colsplit(rownames(pmat),'_',names=c('c1','c2','c3'))$c2
end=reshape2::colsplit(rownames(pmat),'_',names=c('c1','c2','c3'))$c3

features = paste0(chr,':',start,'-',end)
rownames(pmat)=features
barcodes = colnames(pmat)
pmat[1:5,1:5]
dim(pmat)

#dir.create(savepath)
Matrix::writeMM(pmat, file = paste0(savepath, 'matrix.mtx'))
write.table(features, file = paste0(savepath, 'features.csv'), row.names = F, col.names = F)
write.table(barcodes, file = paste0(savepath, 'barcodes.csv'), row.names = F, col.names = F)

proj <- readRDS('~/CH/CH-cross/Fig3/species/intergrate_hg38/hg38_intergrate_bigFigure_0311.rds')
meta = as.data.frame(proj@cellColData)
meta.use=meta[colnames(pmat),]
write.csv(meta.use, file = paste0(savepath, 'metadata.csv'), row.names = F)