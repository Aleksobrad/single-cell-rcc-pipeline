#kaplan-meier and test for enrichment of tumor mac signature (GSEA) with recurrence/ time-to-recurrence
cd45pos=readRDS("cd45pos_merged_seurat.rds")
cd45neg=readRDS("cd45neg_merged_seurat.rds")
Idents(cd45pos) <- "seurat_clusters"
Idents(cd45neg) <- "seurat_clusters"
cd45pos_vp=readRDS("cd45pos_merged_seurat_viper.rds")
cd45neg_vp=readRDS("cd45neg_merged_seurat_viper.rds")
Idents(cd45pos_vp) <- "seurat_clusters"
Idents(cd45neg_vp) <- "seurat_clusters"
cd45pos=cd45pos[,colnames(cd45pos_vp)]
cd45neg=cd45neg[,colnames(cd45neg_vp)]
cd45pos$seurat_clusters=cd45pos_vp$seurat_clusters
cd45neg$seurat_clusters=cd45neg_vp$seurat_clusters
Idents(cd45pos) <- "seurat_clusters"
Idents(cd45neg) <- "seurat_clusters"
markers.cd45pos.viper <- FindAllMarkers(cd45pos_vp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,test.use = "t")
markers.cd45neg.viper <- FindAllMarkers(cd45neg_vp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,test.use = "t")
markers.cd45pos <- FindAllMarkers(cd45pos, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
markers.cd45neg <- FindAllMarkers(cd45neg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
write.csv(markers.cd45pos,file="CD45POS_marker_genes.csv")
write.csv(markers.cd45neg,file="CD45NEG_marker_genes.csv")
write.csv(markers.cd45pos.viper,file="CD45POS_marker_proteins.csv")
write.csv(markers.cd45neg.viper,file="CD45NEG_marker_proteins.csv")

Mac_signature_genexp=markers.cd45pos$gene[which(markers.cd45pos$cluster=="Macrophage")]
Mac_signature_vp=markers.cd45pos.viper$gene[which(markers.cd45pos.viper$cluster=="Macrophage")]

recurrence_annotation=data.frame(sample=c("VC006","VC007","VC008","VC009","VC010","VC011","VC001","VC002","VC003","VC004","VC005"),
                                 recurrence=c(T,T,T,T,T,T,F,F,F,F,F),
                                 MRN=c(5201663,5485822,3557750,4601078,5574669,5697201,1488960,5332167,5436080,3475752,5799053),
                                 Gender=c("F","F","F","M","M","F","F","M","M","F","F"),
                                 Tumor_size_cm=c(8,8.3,3.5,10,15,10,8.3,9,4,5.2,11.3),
                                 Grade=c("3","2","2","2/3","3","3","1","2","2/3","2","3"),
                                 Margins=rep("Negative",11),
                                 Stage=c("pT2","pT2","pT3b","pT3a","pT3a","pT3a","pT2","pT3b","pT3b","pT3b","pT3a"),
                                 Time_to_recurrence_months=c(12,5,15,82,8,12,110,86,113,35,64))
rownames(recurrence_annotation)=recurrence_annotation[,1]

#DESeq2 analysis
data=read.table("Vinson_est_counts_genes_kallisto.txt",sep="\t",header=T)
data=data[which(rowSums(data)>0),]
group=recurrence_annotation[colnames(data),2]
data_normalized=log10(apply(data,2,function(x){x/sum(x)*1000000})+1)
data_normalized_scaled=data_normalized[which(apply(data_normalized,1,sd)>0),]
data_normalized_scaled=t(apply(data_normalized_scaled,1,function(x){(x-mean(x))/sd(x)}))
pca=prcomp(t(data_normalized))
plot(pca$x[,1],pca$x[,2],pch=19,cex=1,col=as.numeric(group)+1,xlab="PC1",ylab="PC2")
text(pca$x[,1], pca$x[,2]-.5, colnames(data_normalized),cex=0.5)

data=data[,setdiff(colnames(data),c("VC005","VC007","VC008"))]
group=recurrence_annotation[colnames(data),2]
data_normalized=log10(apply(data,2,function(x){x/sum(x)*1000000})+1)
pca=prcomp(t(data_normalized))
plot(pca$x[,1],pca$x[,2],pch=19,cex=1,col=as.numeric(group)+1,xlab="PC1",ylab="PC2")
text(pca$x[,1], pca$x[,2]-.5, colnames(data_normalized),cex=0.5)

library(DESeq2)
colData <- matrix(NA,nrow=ncol(data),ncol=2)
rownames(colData) <- colnames(data)
colnames(colData) <- c("condition","type")
colData[1:4,"condition"] <- "no_recurrence"
colData[5:8,"condition"] <- "recurrence"
colData[,"type"] <- "paired-end"
colData <- as.data.frame(colData)

dds <- DESeqDataSetFromMatrix(countData = data,colData = colData,design= ~ condition)
dds <- DESeq(dds)
res <- results(dds, name="condition_recurrence_vs_no_recurrence")
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered),file="res.csv")
res <- read.csv("res.csv", header=TRUE, colClasses=c("character",rep("numeric",6)))

##TUMOR MACROPHAGE SIGNATURE
library(atools)
set.seed(1234)
res=res[which(!is.na(res$pvalue)),]
x=res$log2FoldChange
names(x)=res$X
x=x[which(!is.na(x))]
g=rep(100,length(Mac_signature_genexp))
names(g)=Mac_signature_genexp
gsea(signature = x, geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
ledge=gsea(signature = x, geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))$ledge

means=apply(data_normalized,1,mean)
sds=apply(data_normalized,1,sd)
data_normalized_scaled=(data_normalized-means)/sds
data_normalized_scaled=data_normalized_scaled[which(sds>0),]
gsea(signature = data_normalized_scaled[,1], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = data_normalized_scaled[,2], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = data_normalized_scaled[,3], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = data_normalized_scaled[,4], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = data_normalized_scaled[,5], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = data_normalized_scaled[,6], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = data_normalized_scaled[,7], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = data_normalized_scaled[,8], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))

clinical_dat=recurrence_annotation[colnames(data_normalized_scaled),]
clinical_dat$TumorMacrophageGSEA=c(-5.23,-1.42,-5.05,-0.77,4.09,-.38,4.11,2.03)
library(survival)
library(survminer)
clinical_dat$SurvObj=with(clinical_dat,Surv(Time_to_recurrence_months,recurrence))
res.cox1 <- coxph(SurvObj ~ ., data =  clinical_dat[,c(10,11)])
summary(res.cox1)
ggforest(res.cox1, data = clinical_dat,fontsize = 1,main = "Tumor Macrophage GSEA Hazard Ratio")
clinical_dat$TumorMacrophage=clinical_dat$TumorMacrophageGSEA
clinical_dat$TumorMacrophage[which(clinical_dat$TumorMacrophageGSEA<0)]="low"
clinical_dat$TumorMacrophage[which(clinical_dat$TumorMacrophageGSEA>0)]="high"
clinical_dat$TumorMacrophage=as.factor(clinical_dat$TumorMacrophage)
km.by.resp <- survfit(SurvObj ~ TumorMacrophage, data = clinical_dat, conf.type = "log-log")
km.by.resp
ggsurvplot(km.by.resp, data = clinical_dat, pval = T,ylab="Fraction Without Recurrence",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(90, 0.9))

library("pheatmap")
select <- names(ledge)
select=unique(c("APOE","C1QA","C1QB","C1QC",names(ledge)))
df <- as.data.frame(colData(dds)[,c("condition")])
colnames(df)="condition"
pheatmap(data_normalized_scaled[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df,fontsize_row = 8)



##
##VIPER TUMOR MACROPHAGE SIGNATURE-- SINGLE-CELL NETWORKS!
## 
library(viper)
set.seed(1234)
data_normalized_scaled=data_normalized[which(apply(data_normalized,1,sd)>0),]
data_normalized_scaled=t(apply(data_normalized_scaled,1,function(x){(x-mean(x))/sd(x)}))
filenames <- list.files("sc_nets/", pattern="Patient*", full.names=TRUE)
filenames=filenames[grepl("pos",filenames)]
nets=lapply(filenames,readRDS)
library(atools)
g=rep(100,50)
names(g)=Mac_signature_vp[1:50]
vp=viper(data_normalized_scaled,nets,method="none")
vp[,2]=-vp[,2]
x=apply(vp,1,function(x){mean(x[5:8])-mean(x[1:4])})
gsea(signature = x, geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
ledge=gsea(signature = x, geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))$ledge
vp=viper(data_normalized_scaled,nets,method="none")
gsea(signature = vp[,1], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = vp[,2], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = vp[,3], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = vp[,4], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = vp[,5], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = vp[,6], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = vp[,7], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))
gsea(signature = vp[,8], geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))

clinical_dat=recurrence_annotation[colnames(data_normalized_scaled),]
clinical_dat$TumorMacrophageGSEA=c(-4.64,-2.24,-3.35,-3.49,5.17,-4.27,1.61,2.73)
library(survival)
library(survminer)
clinical_dat$SurvObj=with(clinical_dat,Surv(Time_to_recurrence_months,recurrence))
res.cox1 <- coxph(SurvObj ~ ., data =  clinical_dat[,c(10,11)])
summary(res.cox1)
ggforest(res.cox1, data = clinical_dat,fontsize = 1,main = "Tumor Macrophage GSEA Hazard Ratio")
clinical_dat$TumorMacrophage=clinical_dat$TumorMacrophageGSEA
clinical_dat$TumorMacrophage[which(clinical_dat$TumorMacrophageGSEA<0)]="low"
clinical_dat$TumorMacrophage[which(clinical_dat$TumorMacrophageGSEA>0)]="high"
clinical_dat$TumorMacrophage=as.factor(clinical_dat$TumorMacrophage)
km.by.resp <- survfit(SurvObj ~ TumorMacrophage, data = clinical_dat, conf.type = "log-log")
km.by.resp
ggsurvplot(km.by.resp, data = clinical_dat, pval = T,ylab="Fraction Without Recurrence",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(90, 0.9))

library("pheatmap")
select <- names(ledge)
select=unique(c("LILRB5","TREM2","APOE",names(ledge)))
df <- as.data.frame(colData(dds)[,c("condition")])
colnames(df)="condition"
pheatmap(t(apply(vp[select,],1,function(x){(x-mean(x))/sd(x)}))*3, cluster_rows=FALSE, show_rownames=T,cluster_cols=FALSE, annotation_col=df,fontsize_row = 8)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(t(apply(vp[select,],1,function(x){(x-mean(x))/sd(x)}))*3, n = 30)
pheatmap(t(apply(vp[select,],1,function(x){(x-mean(x))/sd(x)}))*3, cluster_rows=FALSE, show_rownames=T,cluster_cols=FALSE, annotation_col=df,fontsize_row = 8,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'yellow', 'red'))(length(mat_breaks)))

