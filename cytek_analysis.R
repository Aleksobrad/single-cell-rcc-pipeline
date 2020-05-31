###Sample IDs for RCC CyTEK. 
#TvsAdjNormal: 1==Tumor, 2==AdjNormal
#PatientID= Patient number 1-10 as in excel
#patientdate MonthDayYear

###cytek_analysis ALL PATIENTS LYMPHOID-- vs VIPER and Gene Expression
library(flowCore)
library(ggcyto)
set.seed(1234)

lymphoid=read.FCS("cytek_data/export_REDO_concat_1_T_Cells_CD45+_Live.fcs",transformation=F,truncate_max_range=F)
summary(lymphoid)
table(lymphoid@exprs[,29])
a=cbind(markernames(lymphoid),colnames(lymphoid@exprs)[5:27])
autoplot(lymphoid, "FSC-A", "FSC-H")
autoplot(lymphoid,"FJComp-BV480-A","FJComp-Pacific Blue-A",bins=100)
translist <- estimateLogicle(lymphoid, as.character(colnames(lymphoid@exprs)[1:27]))
after <- transform(lymphoid, translist)
summary(after)
autoplot(after, "FSC-A", "FSC-H")
autoplot(after,"FJComp-BV480-A","FJComp-Pacific Blue-A",bins=100)

i=sample(which(lymphoid@exprs[,29] %in% c(1,2,3,4,6,8,9)),250000)
dat=as.data.frame(t(after@exprs[i,c(5:18,20:26)]))
rownames(dat)=toupper(a[c(1:14,16:22),1])
cytek_seurat <- CreateSeuratObject(counts = dat)
cytek_seurat$patient=as.factor(lymphoid@exprs[i,29])
cytek_seurat$patient=mapvalues(cytek_seurat$patient, from = c("1","2","3","4","6","8","9"), to = c("Patient1","Patient2", "Patient3","Patient4","Patient5","Patient6","Patient7"))
cytek_seurat$tissue=as.factor(lymphoid@exprs[i,31])
cytek_seurat$tissue=mapvalues(cytek_seurat$tissue, from = c("1","2"), to = c("Tumor", "Normal"))
cytek_seurat@assays$RNA@scale.data=as.matrix(cytek_seurat@assays$RNA@data)
cytek_seurat <- RunPCA(cytek_seurat,features=rownames(cytek_seurat))
cytek_seurat <- RunUMAP(cytek_seurat, dims = 1:20, verbose = FALSE,umap.method="umap-learn",metric="correlation")
cytek_seurat <- FindNeighbors(cytek_seurat, dims = 1:20, verbose = FALSE)
cytek_seurat <- FindClusters(cytek_seurat, resolution=seq(0.01,0.1,by=0.01), verbose = FALSE,algorithm=1) 
clust=cytek_seurat@meta.data[,which(grepl("RNA_snn_res.",colnames(cytek_seurat@meta.data)))]
mat=cytek_seurat@assays$RNA@scale.data
out=sil_subsample_v2(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,0.1,by=0.01)
tiff("cytek_lymphoid_resolution.tiff", width = 8, height = 6, units = 'in', res = 100)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
dev.off()
best=tail(x[which(means==max(means))],n=1)
cytek_seurat$seurat_clusters=cytek_seurat@meta.data[,which(colnames(cytek_seurat@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(cytek_seurat) <- "seurat_clusters"
tiff("cytek_lymphoid_louvain.tiff", width = 8, height = 6, units = 'in', res = 100)
DimPlot(cytek_seurat, reduction = "umap",label = TRUE) + NoLegend()
dev.off()
cytek_seurat$gene_clustering=cytek_seurat$seurat_clusters
cytek_seurat.markers <- FindAllMarkers(cytek_seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.2,test.use = "wilcox")
top10 <- cytek_seurat.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tiff("cytek_lymphoid_heatmap.tiff", width = 8, height = 6, units = 'in', res = 100)
geneHeatmap_viper(cytek_seurat,cytek_seurat$seurat_clusters,unique(cytek_seurat.markers$gene))
dev.off()
tiff("cytek_lymphoid_louvain_split.tiff", width = 8, height = 6, units = 'in', res = 100)
DimPlot(cytek_seurat, reduction = "umap",label = TRUE,split.by = "tissue") + NoLegend()
dev.off()
tiff("cytek_lymphoid_louvain_split_patient.tiff", width = 8, height = 6, units = 'in', res = 100)
DimPlot(cytek_seurat, reduction = "umap",label = TRUE,split.by = "patient",ncol=3) + NoLegend()
dev.off()
saveRDS(cytek_seurat,file="cytek_seurat_lymphoid.rds")
rm(cytek_seurat,after,lymphoid)

myeloid <- read.FCS("cytek_data/export_RCC_Myeloid_LiveCD45.fcs", transformation=F,truncate_max_range = F)
summary(myeloid)
table(myeloid@exprs[,29])
b=cbind(markernames(myeloid),colnames(myeloid@exprs)[5:27])
autoplot(myeloid, "FSC-A", "FSC-H")
autoplot(myeloid,"FJComp-BV480-A","FJComp-Pacific Blue-A",bins=100)
translist <- estimateLogicle(myeloid, as.character(colnames(myeloid@exprs)[1:27]))
after_m <- transform(myeloid, translist)
summary(after_m)
autoplot(after_m, "FSC-A", "FSC-H")
autoplot(after_m,"FJComp-BV480-A","FJComp-Pacific Blue-A",bins=100)

i=sample(which(myeloid@exprs[,29] %in% c(1,2,3,4,6,8,9)),250000)
dat=as.data.frame(t(after_m@exprs[i,c(5:18,20:26)]))
rownames(dat)=toupper(b[c(1:14,16:22),1])
cytek_seurat_m <- CreateSeuratObject(counts = dat)
cytek_seurat_m$patient=myeloid@exprs[i,29]
cytek_seurat_m$tissue=as.factor(myeloid@exprs[i,31])
cytek_seurat_m$patient=mapvalues(cytek_seurat_m$patient, from = c("1","2","3","4","6","8","9"), to = c("Patient1","Patient2", "Patient3","Patient4","Patient5","Patient6","Patient7"))
cytek_seurat_m$tissue=mapvalues(cytek_seurat_m$tissue, from = c("1","2"), to = c("Tumor", "Normal"))
cytek_seurat_m@assays$RNA@scale.data=as.matrix(cytek_seurat_m@assays$RNA@data)
cytek_seurat_m <- RunPCA(cytek_seurat_m,features=rownames(cytek_seurat_m))
cytek_seurat_m <- RunUMAP(cytek_seurat_m, dims = 1:20, verbose = FALSE,umap.method="umap-learn",metric="correlation")
cytek_seurat_m <- FindNeighbors(cytek_seurat_m, dims = 1:20, verbose = FALSE)
cytek_seurat_m <- FindClusters(cytek_seurat_m, resolution=seq(0.01,0.1,by=0.01), verbose = FALSE,algorithm=1) 
clust=cytek_seurat_m@meta.data[,which(grepl("RNA_snn_res.",colnames(cytek_seurat_m@meta.data)))]
mat=cytek_seurat_m@assays$RNA@scale.data
out=sil_subsample_v2(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,0.1,by=0.01)
tiff("cytek_myeloid_resolution.tiff", width = 8, height = 6, units = 'in', res = 100)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
dev.off()
best=tail(x[which(means==max(means))],n=1)
cytek_seurat_m$seurat_clusters=cytek_seurat_m@meta.data[,which(colnames(cytek_seurat_m@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(cytek_seurat_m) <- "seurat_clusters"
tiff("cytek_myeloid_louvain.tiff", width = 8, height = 6, units = 'in', res = 100)
DimPlot(cytek_seurat_m, reduction = "umap",label = TRUE) + NoLegend()
dev.off()
cytek_seurat_m$gene_clustering=cytek_seurat_m$seurat_clusters
cytek_seurat_m.markers <- FindAllMarkers(cytek_seurat_m, only.pos = T, min.pct = 0.25, logfc.threshold = 0.2,test.use = "wilcox")
top10 <- cytek_seurat_m.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tiff("cytek_myeloid_heatmap.tiff", width = 8, height = 6, units = 'in', res = 100)
geneHeatmap_viper(cytek_seurat_m,cytek_seurat_m$seurat_clusters,unique(cytek_seurat_m.markers$gene))
dev.off()
tiff("cytek_myeloid_louvain_split.tiff", width = 8, height = 6, units = 'in', res = 100)
DimPlot(cytek_seurat_m, reduction = "umap",label = TRUE,split.by = "tissue") + NoLegend()
dev.off()
tiff("cytek_myeloid_louvain_split_patient.tiff", width = 8, height = 6, units = 'in', res = 100)
DimPlot(cytek_seurat_m, reduction = "umap",label = TRUE,split.by = "patient",ncol=3) + NoLegend()
dev.off()
saveRDS(cytek_seurat_m,file="cytek_seurat_myeloid.rds")
rm(cytek_seurat_m,after_m,myeloid)


########plots with same proteins in same order for CyTEK, gene exp, and VIPER:: LYMPHOID
gene_to_protein=read.csv("cytek_data/cytek_protein_to_gene.csv")
gene_to_protein=gene_to_protein[which(gene_to_protein[,2]!=""),]
gene_to_protein[,1]=toupper(gene_to_protein[,1])
cytek_seurat=readRDS("cytek_seurat_lymphoid.rds")
cytek_seurat=cytek_seurat[which(rownames(cytek_seurat) %in% gene_to_protein[,1]),]
Idents(cytek_seurat) <- "seurat_clusters"
cytek_seurat.markers <- FindAllMarkers(cytek_seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.2,test.use = "wilcox")
top10 <- cytek_seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
geneHeatmap_viper(cytek_seurat,cytek_seurat$seurat_clusters,c(unique(cytek_seurat.markers$gene),setdiff(rownames(cytek_seurat),unique(cytek_seurat.markers$gene))))
DimPlot(cytek_seurat, reduction = "umap",label = TRUE) + NoLegend()
DimPlot(cytek_seurat, reduction = "umap",label = TRUE,split.by = "patient",ncol=3) + NoLegend()
lymphoid_proteins=c(unique(cytek_seurat.markers$gene),setdiff(rownames(cytek_seurat),unique(cytek_seurat.markers$gene)))
rownames(gene_to_protein)=gene_to_protein[,1]
lymphoid_genes=toupper(unlist(strsplit(as.character(gene_to_protein[lymphoid_proteins,2]),", ")))

gene=readRDS("cd45pos_merged_seurat.rds")
Idents(gene)="seurat_clusters"
gene=subset(gene,subset=patient %in% c("Patient1","Patient2","Patient3","Patient4","Patient5","Patient6","Patient7"))
dat=gene@assays$RNA@data
dat=log10(apply(dat,2,function(x){x/sum(x)*1000000})+1)
dat=dat[lymphoid_genes,]
dat=t(apply(dat,1,function(x){(x-mean(x))/sd(x)}))
gene_a=CreateSeuratObject(counts = dat)
gene_a@assays$RNA@scale.data=as.matrix(gene_a@assays$RNA@data)
gene_a <- RunPCA(gene_a,features=rownames(gene_a))
gene_a <- RunUMAP(gene_a, dims = 1:(nrow(gene_a)-1), verbose = FALSE,umap.method="umap-learn",metric="correlation")
gene_a <- FindNeighbors(gene_a, dims = 1:(nrow(gene_a)-1), verbose = FALSE)
gene_a <- FindClusters(gene_a, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=gene_a@meta.data[,which(grepl("RNA_snn_res.",colnames(gene_a@meta.data)))]
mat=gene_a@assays$RNA@scale.data
out=sil_subsample_v2(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
gene_a$seurat_clusters=gene_a@meta.data[,which(colnames(gene_a@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(gene_a) <- "seurat_clusters"
DimPlot(gene_a, reduction = "umap",label = TRUE) + NoLegend()
gene_a$patient=gene$patient
gene_a$tissue=gene$tissue
gene_a$gene_clustering=gene$seurat_clusters
geneHeatmap_viper(gene_a,gene_a$seurat_clusters,lymphoid_genes)
DimPlot(gene_a, reduction = "umap",label = TRUE,split.by = "patient",ncol=3) + NoLegend()
saveRDS(gene_a,file="gene_cytek_comparison_lymphoid.rds")
rm(gene_a,gene)


list_cd45pos=list()
cd45pos_patients=c("Patient1","Patient2","Patient3","Patient4","Patient5","Patient6","Patient7")
for(pt in cd45pos_patients){
  print(paste(pt,"CD45pos VIPER",sep=" "))
  s1_vp_meta=readRDS(paste(pt,"_cd45pos_vp.rds",sep=""))
  s1=readRDS(paste(pt,"_CD45pos.rds",sep=""))
  s1_vp_meta <- s1_vp_meta[intersect(lymphoid_genes,rownames(s1_vp_meta)),]
  s1_vp_meta_seurat <- CreateSeuratObject(counts = s1_vp_meta)
  s1_vp_meta_seurat@assays$RNA@scale.data=as.matrix(s1_vp_meta_seurat@assays$RNA@data)
  s1_vp_meta_seurat$tissue=s1$tissue[colnames(s1_vp_meta_seurat)]
  s1_vp_meta_seurat$patient=s1$patient[colnames(s1_vp_meta_seurat)]
  s1_vp_meta_seurat$gene_clustering=s1$seurat_clusters[colnames(s1_vp_meta_seurat)]
  s1_vp_meta_seurat$l=s1$l[colnames(s1_vp_meta_seurat)]
  list_cd45pos=c(list_cd45pos,s1_vp_meta_seurat)
  rm(s1,s1_vp_meta_seurat,s1_vp_meta)
}
pos_genes=lapply(list_cd45pos,rownames)
pos_genes_intersect=pos_genes[[1]]
for(i in 2:length(pos_genes)){
  pos_genes_intersect=intersect(pos_genes_intersect,pos_genes[[i]])
}
list_cd45pos=lapply(list_cd45pos,function(x){x[pos_genes_intersect,]})

vp=readRDS("cd45pos_merged_seurat_viper.rds")
vp_a=merge(list_cd45pos[[1]],y=list_cd45pos[2:length(list_cd45pos)],project = "RCC_SC_VIPER")
rm(list_cd45pos)
x=as.matrix(vp_a@assays$RNA@data)
x["CD4",]=x["IL7R",]
vp_a@assays$RNA@scale.data=x
vp_a$patient=vp$patient[colnames(vp_a)]
vp_a$tissue=vp$tissue[colnames(vp_a)]
vp_a$gene_clustering=vp$gene_clustering[colnames(vp_a)]
vp_a=vp_a[,which(vp_a$gene_clustering!=9)]
vp_a <- RunPCA(vp_a,features=rownames(vp_a))
vp_a <- RunUMAP(vp_a, dims = 1:(nrow(vp_a)-1), verbose = FALSE,umap.method="umap-learn",metric="correlation")
vp_a <- FindNeighbors(vp_a, dims = 1:(nrow(vp_a)-1), verbose = FALSE)
vp_a <- FindClusters(vp_a, resolution=seq(0.05,0.2,by=0.01), verbose = FALSE,algorithm=1) 
clust=vp_a@meta.data[,which(grepl("RNA_snn_res.",colnames(vp_a@meta.data)))]
mat=vp_a@assays$RNA@scale.data
out=sil_subsample_v2(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.05,0.2,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
vp_a$seurat_clusters=vp_a@meta.data[,which(colnames(vp_a@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(vp_a) <- "seurat_clusters"
DimPlot(vp_a, reduction = "umap",label = TRUE) + NoLegend()
geneHeatmap_viper(vp_a,vp_a$seurat_clusters,lymphoid_genes[lymphoid_genes %in% pos_genes_intersect])
DimPlot(vp_a, reduction = "umap",label = TRUE,split.by = "patient",ncol=3) + NoLegend()
saveRDS(vp_a,file="viper_cytek_comparison_lymphoid.rds")
rm(vp_a,vp)


########plots with same proteins in same order for CyTEK, gene exp, and VIPER-- MYELOID
gene_to_protein=read.csv("cytek_data/cytek_protein_to_gene.csv")
gene_to_protein=gene_to_protein[which(gene_to_protein[,2]!=""),]
gene_to_protein[,1]=toupper(gene_to_protein[,1])
cytek_seurat=readRDS("cytek_seurat_myeloid.rds")
cytek_seurat=cytek_seurat[which(rownames(cytek_seurat) %in% gene_to_protein[,1]),]
Idents(cytek_seurat) <- "seurat_clusters"
cytek_seurat.markers <- FindAllMarkers(cytek_seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.2,test.use = "wilcox")
top10 <- cytek_seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
geneHeatmap_viper(cytek_seurat,cytek_seurat$seurat_clusters,c(unique(cytek_seurat.markers$gene),setdiff(rownames(cytek_seurat),unique(cytek_seurat.markers$gene))))
DimPlot(cytek_seurat, reduction = "umap",label = TRUE) + NoLegend()
DimPlot(cytek_seurat, reduction = "umap",label = TRUE,split.by = "patient",ncol=3) + NoLegend()
myeloid_proteins=c(unique(cytek_seurat.markers$gene),setdiff(rownames(cytek_seurat),unique(cytek_seurat.markers$gene)))
rownames(gene_to_protein)=gene_to_protein[,1]
myeloid_genes=toupper(unlist(strsplit(as.character(gene_to_protein[myeloid_proteins,2]),", ")))


gene=readRDS("cd45pos_merged_seurat.rds")
Idents(gene)="seurat_clusters"
gene=subset(gene,subset=patient %in% c("Patient1","Patient2","Patient3","Patient4","Patient5","Patient6","Patient7"))
dat=gene@assays$RNA@data
dat=log10(apply(dat,2,function(x){x/sum(x)*1000000})+1)
myeloid_genes=intersect(myeloid_genes,rownames(dat))
dat=dat[myeloid_genes,]
dat=t(apply(dat,1,function(x){(x-mean(x))/sd(x)}))
gene_a=CreateSeuratObject(counts = dat)
gene_a@assays$RNA@scale.data=as.matrix(gene_a@assays$RNA@data)
gene_a <- RunPCA(gene_a,features=rownames(gene_a))
gene_a <- RunUMAP(gene_a, dims = 1:(nrow(gene_a)-1), verbose = FALSE,umap.method="umap-learn",metric="correlation")
gene_a <- FindNeighbors(gene_a, dims = 1:(nrow(gene_a)-1), verbose = FALSE)
gene_a <- FindClusters(gene_a, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=gene_a@meta.data[,which(grepl("RNA_snn_res.",colnames(gene_a@meta.data)))]
mat=gene_a@assays$RNA@scale.data
out=sil_subsample_v2(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
gene_a$seurat_clusters=gene_a@meta.data[,which(colnames(gene_a@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(gene_a) <- "seurat_clusters"
DimPlot(gene_a, reduction = "umap",label = TRUE) + NoLegend()
gene_a$patient=gene$patient
gene_a$tissue=gene$tissue
gene_a$gene_clustering=gene$seurat_clusters
geneHeatmap_viper(gene_a,gene_a$seurat_clusters,myeloid_genes)
DimPlot(gene_a, reduction = "umap",label = TRUE,split.by = "patient",ncol=3) + NoLegend()
saveRDS(gene_a,file="gene_cytek_comparison_myeloid.rds")
rm(gene_a,gene)


list_cd45pos=list()
cd45pos_patients=c("Patient1","Patient2","Patient3","Patient4","Patient5","Patient6","Patient7")
for(pt in cd45pos_patients){
  print(paste(pt,"CD45pos VIPER",sep=" "))
  s1_vp_meta=readRDS(paste(pt,"_cd45pos_vp.rds",sep=""))
  s1=readRDS(paste(pt,"_CD45pos.rds",sep=""))
  s1_vp_meta <- s1_vp_meta[intersect(myeloid_genes,rownames(s1_vp_meta)),]
  s1_vp_meta_seurat <- CreateSeuratObject(counts = s1_vp_meta)
  s1_vp_meta_seurat@assays$RNA@scale.data=as.matrix(s1_vp_meta_seurat@assays$RNA@data)
  s1_vp_meta_seurat$tissue=s1$tissue[colnames(s1_vp_meta_seurat)]
  s1_vp_meta_seurat$patient=s1$patient[colnames(s1_vp_meta_seurat)]
  s1_vp_meta_seurat$gene_clustering=s1$seurat_clusters[colnames(s1_vp_meta_seurat)]
  s1_vp_meta_seurat$l=s1$l[colnames(s1_vp_meta_seurat)]
  list_cd45pos=c(list_cd45pos,s1_vp_meta_seurat)
  rm(s1,s1_vp_meta_seurat,s1_vp_meta)
}
pos_genes=lapply(list_cd45pos,rownames)
pos_genes_intersect=pos_genes[[1]]
for(i in 2:length(pos_genes)){
  pos_genes_intersect=intersect(pos_genes_intersect,pos_genes[[i]])
}
list_cd45pos=lapply(list_cd45pos,function(x){x[pos_genes_intersect,]})

vp=readRDS("cd45pos_merged_seurat_viper.rds")
vp_a=merge(list_cd45pos[[1]],y=list_cd45pos[2:length(list_cd45pos)],project = "RCC_SC_VIPER")
rm(list_cd45pos)
x=as.matrix(vp_a@assays$RNA@data)
vp_a@assays$RNA@scale.data=x
vp_a$patient=vp$patient[colnames(vp_a)]
vp_a$tissue=vp$tissue[colnames(vp_a)]
vp_a$gene_clustering=vp$gene_clustering[colnames(vp_a)]
vp_a=vp_a[,which(vp_a$gene_clustering!=9)]
vp_a <- RunPCA(vp_a,features=rownames(vp_a))
vp_a <- RunUMAP(vp_a, dims = 1:(nrow(vp_a)-1), verbose = FALSE,umap.method="umap-learn",metric="correlation")
vp_a <- FindNeighbors(vp_a, dims = 1:(nrow(vp_a)-1), verbose = FALSE)
vp_a <- FindClusters(vp_a, resolution=seq(0.05,0.2,by=0.01), verbose = FALSE,algorithm=1) 
clust=vp_a@meta.data[,which(grepl("RNA_snn_res.",colnames(vp_a@meta.data)))]
mat=vp_a@assays$RNA@scale.data
out=sil_subsample_v2(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.05,0.2,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
vp_a$seurat_clusters=vp_a@meta.data[,which(colnames(vp_a@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(vp_a) <- "seurat_clusters"
DimPlot(vp_a, reduction = "umap",label = TRUE) + NoLegend()
geneHeatmap_viper(vp_a,vp_a$seurat_clusters,myeloid_genes[myeloid_genes %in% pos_genes_intersect])
DimPlot(vp_a, reduction = "umap",label = TRUE,split.by = "patient",ncol=3) + NoLegend()
saveRDS(vp_a,file="viper_cytek_comparison_myeloid.rds")
rm(vp_a,vp)





