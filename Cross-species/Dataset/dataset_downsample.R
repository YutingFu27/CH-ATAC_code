#sample 300cells per subcluster
#hg38
proj2 <- readRDS('~/CH/CH-cross/Fig3/species/hg38/hg38_proj_all.rds')

df=as.data.frame(proj2@cellColData)
table(df$subanno)
table(df$lineage)

df=df[!grepl('unknown',df$subanno),]
df=df[!grepl('Other',df$lineage),]
stat=table(df$subanno,df$lineage)%>%as.data.frame()
stat=stat[stat$Freq>0,]
summary(stat$Freq)

all.bc=c()
for (i in unique(stat$Var1)) {
  tmp=df[df$subanno==i,]
  num=dim(tmp)[1]
  if(num>300){
    bc=rownames(tmp[sample(num,300),])
  }
  else{
    bc=rownames(tmp)
  }
  all.bc=c(all.bc,bc)
}

proj=proj2[proj2$cellNames%in%all.bc,]
proj
table(proj$lineage)
table(proj$subanno)
table(proj$main_ct)
saveRDS(proj,file = '~/CH/CH-cross/Fig3/species/hg38/hg38_proj_sample.rds')
proj#31464



#mm9
proj2 <- readRDS('~/CH/CH-cross/Fig3/species/mm9/mm9_proj_all.rds')

df=as.data.frame(proj2@cellColData)
table(df$subanno)
table(df$lineage)

df=df[!grepl('unknown',df$subanno),]
df=df[!grepl('Other',df$lineage),]
stat=table(df$subanno,df$lineage)%>%as.data.frame()
stat=stat[stat$Freq>0,]
summary(stat$Freq)

all.bc=c()
for (i in unique(stat$Var1)) {
  tmp=df[df$subanno==i,]
  num=dim(tmp)[1]
  if(num>300){
    bc=rownames(tmp[sample(num,300),])
  }
  else{
    bc=rownames(tmp)
  }
  all.bc=c(all.bc,bc)
}

proj=proj2[proj2$cellNames%in%all.bc,]
proj
table(proj$lineage)
table(proj$subanno)
table(proj$main_ct)
saveRDS(proj,file = '~/CH/CH-cross/Fig3/species/mm9/mm9_proj_sample.rds')
proj #20138


#macaca
proj2 <- readRDS('~/CH/CH-cross/Fig3/species/macaca/macaca_proj_all.rds')

df=as.data.frame(proj2@cellColData)
table(df$subanno)
table(df$lineage)

df=df[!grepl('Unknown',df$subanno),]
stat=table(df$subanno,df$lineage)%>%as.data.frame()
stat=stat[stat$Freq>0,]
summary(stat$Freq)

all.bc=c()
for (i in unique(stat$Var1)) {
  tmp=df[df$subanno==i,]
  num=dim(tmp)[1]
  if(num>300){
    bc=rownames(tmp[sample(num,300),])
  }
  else{
    bc=rownames(tmp)
  }
  all.bc=c(all.bc,bc)
}

proj=proj2[proj2$cellNames%in%all.bc,]
proj
table(proj$lineage)
table(proj$subanno)
table(proj$main_ct)
saveRDS(proj,file = '~/CH/CH-cross/Fig3/species/macaca/macaca_proj_sample.rds')
proj 


#dr11
proj2 <- readRDS('~/CH/CH-cross/Fig3/species/dr11/dr11_proj_all.rds')

df=as.data.frame(proj2@cellColData)
table(df$subanno)
table(df$lineage)

df=df[!grepl('-unknown',df$subanno),]
stat=table(df$subanno,df$lineage)%>%as.data.frame()
stat=stat[stat$Freq>0,]
summary(stat$Freq)

all.bc=c()
for (i in unique(stat$Var1)) {
  tmp=df[df$subanno==i,]
  num=dim(tmp)[1]
  if(num>300){
    bc=rownames(tmp[sample(num,300),])
  }
  else{
    bc=rownames(tmp)
  }
  all.bc=c(all.bc,bc)
}

proj=proj2[proj2$cellNames%in%all.bc,]
proj
table(proj$lineage)
table(proj$subanno)
table(proj$main_ct)
saveRDS(proj,file = '~/CH/CH-cross/Fig3/species/dr11/dr11_proj_sample_20230302.rds')
proj 


#dm6
proj2 <- readRDS('~/CH/CH-cross/Fig3/species/dm6/dm6_proj_all.rds')

df=as.data.frame(proj2@cellColData)
table(df$subanno)
table(df$lineage)

df=df[!grepl('-unknown',df$subanno),]
df=df[!grepl('Other',df$lineage),]
stat=table(df$subanno,df$lineage)%>%as.data.frame()
stat=stat[stat$Freq>0,]
table(stat$Var1)
summary(stat$Freq)

all.bc=c()
for (i in unique(stat$Var1)) {
  tmp=df[df$subanno==i,]
  num=dim(tmp)[1]
  if(num>300){
    bc=rownames(tmp[sample(num,300),])
  }
  else{
    bc=rownames(tmp)
  }
  all.bc=c(all.bc,bc)
}

proj=proj2[proj2$cellNames%in%all.bc,]
proj
table(proj$lineage)
table(proj$subanno)
table(proj$main_ct)
saveRDS(proj,file = '~/CH/CH-cross/Fig3/species/dm6/dm6_proj_sample.rds')
proj 


#ew
proj2 <- readRDS('~/CH/CH-cross/Fig3/species/ew/ew_proj_all_20230304.rds')

df=as.data.frame(proj2@cellColData)
df$subanno=df$cell_type
table(df$subanno)
table(df$lineage)

df=df[!grepl('unknown',df$subanno),]

stat=table(df$subanno,df$lineage)%>%as.data.frame()
stat=stat[stat$Freq>0,]
table(stat$Var1)
summary(stat$Freq)

all.bc=c()
for (i in unique(stat$Var1)) {
  tmp=df[df$subanno==i,]
  num=dim(tmp)[1]
  if(num>300){
    bc=rownames(tmp[sample(num,300),])
  }
  else{
    bc=rownames(tmp)
  }
  all.bc=c(all.bc,bc)
}

proj=proj2[proj2$cellNames%in%all.bc,]
proj
proj$subanno=proj$cell_type
table(proj$lineage)
table(proj$subanno)
table(proj$main_ct)
saveRDS(proj,file = '~/CH/CH-cross/Fig3/species/ew/ew_proj_sample_20230304.rds')
proj
