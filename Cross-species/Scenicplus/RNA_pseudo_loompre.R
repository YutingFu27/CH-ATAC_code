setwd('~/CH/CH-cross/Fig5/')
seob=readRDS('human_rna.rds')

DGE_tab_data=seob@assays$RNA@counts%>%as.data.frame()
phenotype_tab_data=data.frame('CellID'=colnames(seob), 'Celltype'=seob$Lineage)
#pseudo
## file IO
read_tabfile <- function(fname, sep='\t'){
  message(paste("loading", fname, sep='\t'))
  out <- read.table(fname, sep=sep, header = T)
  message(paste("loaded", fname, sep='\t'))
  return(out)
}


## function Pseudocell_analysis_pipeline
# required: DGE_tab_data, phenotype_tab_data          Can be path string or data.frame
# DGE_tab_data: colnames should be CellID
# phenotype_tab_data: should contain ordered CellID, Celltype
# 
# pseudocell_size:  size of pseudocell                DEFAULT: 20
#
# colname_str: set as the colname of pseudocell DGE,  DEFAULT: "Pseudocell_"
# save_out_file: save result or not,                  DEFAULT: TRUE
# out_file_name_str: name of output file,             DEFAULT: "OUT"
#
# return list[result_dge, new_phe] if not save_out_file
# result_dge: Pseudocell dge; new_phe: Pseudocell phenotype
#
Pseudocell_analysis_pipeline <- function(DGE_tab_data, phenotype_tab_data, sep=',',
                                         pseudocell_size = 20, colname_str = 'Pseudocell_',
                                         save_out_file = TRUE, out_file_name_str = 'OUT'){
  # io
  if(is.character(DGE_tab_data) && is.character(phenotype_tab_data)){
    # if input is fname
    # input file io
    message("loading DGE")
    data <- read_tabfile(DGE_tab_data, sep = sep)
    message("loading phenotype.tab")
    id <- read_tabfile(phenotype_tab_data, sep = sep)
    message("input file io done, processing begin")
    
  }else if(is.data.frame(DGE_tab_data) && is.data.frame(phenotype_tab_data)){
    # if input is data.frame
    data <- DGE_tab_data
    id <- phenotype_tab_data
  }else{
    stop('typeof DGE_tab_data and phenotype_tab_data should be string or data.frame')
  }
  
  # check and proc with table cols and rows name
  Inter.id<-id
  if (all(c("CellID", "Celltype") %in% colnames(Inter.id))){
    message("phenotype_tab_data has CellID, Celltype")
  }else if(length(colnames(Inter.id)) == 3){
    colnames(Inter.id)<-c("CellID", "Tissue", "Celltype")
  }else if(length(colnames(Inter.id)) == 2) {
    colnames(Inter.id)<-c("CellID", "Celltype")
  }else{
    message("colnames of phenotype_tab_data more than 3,
            set as CellID, Tissue, Celltype, others")
    colnames(Inter.id)[1:3] <- c("CellID", "Tissue", "Celltype")
  }
  rownames(Inter.id)<-Inter.id$CellID
  # intersection
  Inter<-data
  Inter1<-Inter[,as.character(Inter.id$CellID)]
  if (length(rownames(Inter1)) == 0){
    stop("No intersection in CellID of DGE_tab_data, phenotype_tab_data")
  }
  
  ## CORE
  #normalized to 10e6
  raw.10e6<-t(t(Inter1)/colSums(Inter1))*1000000
  Inter<-raw.10e6
  Inter<-as.data.frame(Inter)
  
  #pseudocell.size was set to 20 by default
  pseudocell.size = pseudocell_size 
  new_ids_list = list()
  for (i in 1:length(unique(Inter.id$Celltype))) {
    cluster_id = unique(Inter.id$Celltype)[i]
    cluster_cells <- rownames(Inter.id[Inter.id$Celltype == cluster_id,])
    cluster_size <- length(cluster_cells)		
    pseudo_ids <- floor((seq_along(cluster_cells)-1)/pseudocell.size)
    pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
    names(pseudo_ids) <- sample(cluster_cells)	
    new_ids_list[[i]] <- pseudo_ids		
  }
  
  #assign new cell ids
  new_ids <- unlist(new_ids_list)
  new_ids <- as.data.frame(new_ids)
  new_ids_length <- table(new_ids)
  
  new_colnames <- rownames(new_ids)  ###add
  all.data<-Inter[,as.character(new_colnames)] ###add
  all.data <- t(all.data)###add
  
  #aggregate
  new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                      list(name=new_ids[,1]),FUN=mean)
  rownames(new.data)<-new.data$name
  new.data<-new.data[,-1]
  
  new_ids_length<-as.matrix(new_ids_length)##
  short<-which(new_ids_length < floor(pseudocell_size * 0.8)) ## short if not 0
  if (length(short)>0){
    new_good_ids<-as.matrix(new_ids_length[-short,])##
  }else{
    new_good_ids<-as.matrix(new_ids_length)}
  
  new_good_ids <- as.matrix(new_ids_length[-short,])##
  new_good_ids <- as.matrix(new_ids_length)
  
  name_str <- colname_str
  result<-t(new.data)[,rownames(new_good_ids)]
  colnames(result)<-paste(name_str,colnames(result),sep="")
  rownames(result)<-rownames(Inter)
  
  cc<-gsub("[_]Cell.*$","",colnames(result))
  new.phe<-cbind(colnames(result), name_str, cc)
  colnames(new.phe)<-c("Sample_ID","Study_ID","Celltype")
  
  if(save_out_file){
    #output file name
    out_rds_str <- paste(out_file_name_str, "pseudocell", pseudocell_size, ".Rds", 
                         sep="_", collapse = NULL)
    out_pheno_str <- paste(out_file_name_str, "pseudocell", pseudocell_size, ".pheno.csv", 
                           sep="_", collapse = NULL)
    
    #output file io
    message("process done, saving data")
    saveRDS(result,file=out_rds_str)
    write.table(new.phe, file=out_pheno_str, sep=',', quote=F, row.names=F)
    message("pipeline exit successfully")
  }else{
    return(list(result, new.phe))
  }
}

Pseudocell_analysis_pipeline(
  DGE_tab_data, 
  phenotype_tab_data,
  pseudocell_size = 50,
  colname_str = 'Pseudocell_',
  save_out_file = TRUE, 
  out_file_name_str = 'hg38')

dge_pseudo=readRDS('./hg38_pseudocell_20_.Rds')
dge_pseudo[1:5,1:5]
meta_pseudo=fread('./hg38_pseudocell_20_.pheno.csv')
dge_pseudo=dge_pseudo[,meta_pseudo$Sample_ID]

seob_pseudo=CreateSeuratObject(dge_pseudo)
seob_pseudo$lineage=meta_pseudo$Celltype

seob=seob_pseudo
seob <- FindVariableFeatures(seob, selection.method = "vst", nfeatures = 3000)

seob <- ScaleData(seob)

seob <- RunPCA(seob)
ElbowPlot(seob, ndims = 50)

seob <- RunUMAP(seob, dims = 1:30)

seob <- FindNeighbors(seob,dims = 1:30)
#1.2聚类
#1.2.1选择cluster数目不变的最小resolution
seob <- FindClusters(seob,
                     resolution = 0.8, # 值越大，cluster 越多
                     random.seed = 1) 

DimPlot(seob,
        reduction = "umap",
        group.by  ="lineage",
        label = T,pt.size =0.2,raster =F) + NoLegend()

table(seob$lineage)
seob$lineage=gsub('Pseudocell_','',seob$lineage)



SaveLoom(object = seob,filename = '~/CH/CH-cross/Fig3/Scenicplus/hg38/hg38_rna_lineage_pseudo_20.loom',overwrite = T)
saveRDS(seob,file = '~/CH/CH-cross/Fig3/Scenicplus/hg38/hg38_rna_pseudo_20.rds')