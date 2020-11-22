#kaplan-meier and fisher test for enrichment of tumor mac signature (GSEA) with recurrence/ time-to-recurrence
#same test for IHC density
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

Mac_signature_genexp=markers.cd45pos$gene[which(markers.cd45pos$cluster=="Macrophage")]
Mac_signature_vp=markers.cd45pos.viper$gene[which(markers.cd45pos.viper$cluster=="Macrophage")]


## load Brian Rini dataset
rini_metadata=read.csv("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/RCC_allsamples_fullpipeline/ccRCC_organized_data_and_code_v3/validation_dataset_metadata.csv")
rini_data=read.table("validation_dataset_bulkRNASeq_counts.gct",sep="\t",header=T)
rini_data=rini_data[which(rini_data[,2] %in% names(which(table(rini_data[,2])==1))),]
rownames(rini_data)=rini_data[,2]
rini_data=rini_data[,3:ncol(rini_data)]
colnames(rini_data)=unlist(lapply(strsplit(colnames(rini_data),"\\."),function(x){x[3]}))
rini_metadata=rini_metadata[which(rini_metadata$rna.seq != "#N/A" & rini_metadata$rna.seq != "0"),]
rini_metadata$rna.seq=unlist(lapply(strsplit(as.character(rini_metadata$rna.seq),"-"),function(x){x[3]}))
rini_metadata=rini_metadata[!grepl("DUPLICATE",rini_metadata$Batch.ID),]
rownames(rini_metadata)=rini_metadata$rna.seq
rini_data=rini_data[,rownames(rini_metadata)]
rini_metadata$Date.of.recurrence=as.character(rini_metadata$Date.of.recurrence)
rini_metadata$Date.of.recurrence[which(rini_metadata$Date.of.recurrence=="")]=as.character(rini_metadata$Date.of.last.contact[which(rini_metadata$Date.of.recurrence=="")])
rini_metadata$Date.of.recurrence=as.Date(rini_metadata$Date.of.recurrence, "%m/%d/%y")
rini_metadata$Date.of.resection.or.CR=as.Date(rini_metadata$Date.of.resection.or.CR, "%m/%d/%y")
rini_metadata$Time_without_recurrence=rini_metadata$Date.of.recurrence-rini_metadata$Date.of.resection.or.CR
rini_metadata$Time_without_recurrence=rini_metadata$Time_without_recurrence/365*12
rini_metadata=rini_metadata[which(!is.na(rini_metadata$Time_without_recurrence)),]
rini_metadata$Time_without_recurrence=as.numeric(rini_metadata$Time_without_recurrence)
#filter out rows
rini_metadata=rini_metadata[which(rini_metadata$Metastasectomy.or..met..Biopsy=="No"),]
rini_data=rini_data[,rownames(rini_metadata)]
rini_metadata$Recurrence


#DESeq2 analysis
data=rini_data
data=data[which(rowSums(data)>0),]
group=rini_metadata$Recurrence
data_normalized=log10(apply(data,2,function(x){x/sum(x)*1000000})+1)
data_normalized_scaled=data_normalized[which(apply(data_normalized,1,sd)>0),]
data_normalized_scaled=t(apply(data_normalized_scaled,1,function(x){(x-mean(x))/sd(x)}))
pca=prcomp(t(data_normalized))
plot(pca$x[,1],pca$x[,2],pch=19,cex=1,col=as.numeric(group)+1,xlab="PC1",ylab="PC2")
text(pca$x[,1], pca$x[,2]-.5, colnames(data_normalized),cex=0.5)

rini_metadata=rini_metadata[setdiff(rownames(rini_metadata),c("1","62","104","28")),]
rini_data=rini_data[,rownames(rini_metadata)]
data=rini_data
data=data[which(rowSums(data)>0),]
group=rini_metadata$Recurrence
data_normalized=log10(apply(data,2,function(x){x/sum(x)*1000000})+1)
data_normalized_scaled=data_normalized[which(apply(data_normalized,1,sd)>0),]
data_normalized_scaled=t(apply(data_normalized_scaled,1,function(x){(x-mean(x))/sd(x)}))
pca=prcomp(t(data_normalized))
plot(pca$x[,1],pca$x[,2],pch=19,cex=1,col=as.numeric(group)+1,xlab="PC1",ylab="PC2")
text(pca$x[,1], pca$x[,2]-.5, colnames(data_normalized),cex=0.5)


library(DESeq2)
colData <- matrix(NA,nrow=ncol(data),ncol=2)
rownames(colData) <- colnames(data)
colnames(colData) <- c("condition","type")
colData[,"condition"]=as.character(rini_metadata$Recurrence)
colData[,"type"] <- "paired-end"
colData <- as.data.frame(colData)

dds <- DESeqDataSetFromMatrix(countData = data,colData = colData,design= ~ condition)
dds <- DESeq(dds)
res <- results(dds, name="condition_Yes_vs_No")
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered),file="res.csv")
res <- read.csv("res.csv", header=TRUE, colClasses=c("character",rep("numeric",6)))

##TUMOR MACROPHAGE SIGNATURE
library(atools)
set.seed(1234)
res=res[which(!is.na(res$pvalue)),]
x=res$log2FoldChange
#x=res$pvalue
names(x)=res$X
x=x[which(!is.na(x))]
g=rep(100,length(Mac_signature_genexp))
names(g)=Mac_signature_genexp
set.seed(1234)
ledge=gsea(signature = x, geneset =g, twoTails = F, pout = TRUE, per = 1000,colSig = c(.15,0.45, 0.3, 1), colHit = c(.05,0.58, 0.1, 2))$ledge

means=apply(data_normalized,1,mean)
sds=apply(data_normalized,1,sd)
data_normalized_scaled=(data_normalized-means)/sds
data_normalized_scaled=data_normalized_scaled[which(sds>0),]
set.seed(1234)
gseaBySample=apply(data_normalized_scaled,2,function(x){gsea(signature = x, geneset =g, twoTails = F, pout = T, per = 500)$nes})

clinical_dat=rini_metadata
clinical_dat$TumorMacrophageGSEA=gseaBySample
p=ggplot(clinical_dat, aes(x=Recurrence, y=TumorMacrophageGSEA,fill=Recurrence)) +
  geom_boxplot() + 
  #geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,fill="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Tumor-Macrophage Enrichment: Recurrence vs No Recurrence")+ylab("NES")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=10, face="bold"),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),legend.position = "none")
library(survival)
library(survminer)
clinical_dat$Recurrence=as.character(clinical_dat$Recurrence)
clinical_dat$Recurrence[which(clinical_dat$Recurrence=="Yes")]=T
clinical_dat$Recurrence[which(clinical_dat$Recurrence=="No")]=F
clinical_dat$Recurrence=as.logical(clinical_dat$Recurrence)
clinical_dat$SurvObj=with(clinical_dat,Surv(Time_without_recurrence,Recurrence))
res.cox1 <- coxph(SurvObj ~ ., data =  clinical_dat[,c(63,64)])
summary(res.cox1)
ggforest(res.cox1, data = clinical_dat,fontsize = 1,main = "Tumor Macrophage GSEA Hazard Ratio")
clinical_dat$TumorMacrophage=clinical_dat$TumorMacrophageGSEA
clinical_dat$TumorMacrophage[which(clinical_dat$TumorMacrophageGSEA<0)]="low"
clinical_dat$TumorMacrophage[which(clinical_dat$TumorMacrophageGSEA>0)]="high"
clinical_dat$TumorMacrophage=as.factor(clinical_dat$TumorMacrophage)
km.by.resp <- survfit(SurvObj ~ TumorMacrophage, data = clinical_dat, conf.type = "log-log")
km.by.resp
ggsurvplot(km.by.resp, data = clinical_dat, pval = T,ylab="Fraction Without Recurrence",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))
#optimal cut-point by maximally selected rank statistic
res.cut <- surv_cutpoint(clinical_dat, time = "Time_without_recurrence", event = "Recurrence",variables = c("TumorMacrophageGSEA"))
summary(res.cut)
plot(res.cut, "TumorMacrophageGSEA", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(Time_without_recurrence, Recurrence) ~TumorMacrophageGSEA, data = res.cat)
ggsurvplot(fit, conf.int=F,pval = T,ylab="Fraction Without Recurrence",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(0, 0))

library("pheatmap")
select <- names(ledge)
select=unique(c("TREM2","APOE","C1QA","C1QB","C1QC",names(ledge)))
df <- as.data.frame(colData(dds)[,c("condition")])
colnames(df)="recurrence"
pheatmap(data_normalized_scaled[select,order(df$recurrence)], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df,fontsize_row = 8)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(data_normalized_scaled[select,order(df$recurrence)], n = 30)
pheatmap(data_normalized_scaled[select,order(df$recurrence)], cluster_rows=FALSE, show_rownames=T,show_colnames = F,cluster_cols=FALSE, annotation_col=df,fontsize_row = 7,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)))


#radomforest and AUC
clinical_dat$Recurrence=as.factor(clinical_dat$Recurrence)
res.cat$Recurrence=as.factor(res.cat$Recurrence)
res.cat$TumorMacrophageGSEA=as.factor(res.cat$TumorMacrophageGSEA)
library(randomForest)
library(pROC)
#dat=clinical_dat[which((clinical_dat$Time_without_recurrence>60 & clinical_dat$Recurrence==FALSE)|(clinical_dat$Time_without_recurrence<60 & clinical_dat$Recurrence==TRUE)),]
dat=clinical_dat
plot(dat$TumorMacrophageGSEA,dat$Time_without_recurrence,pch=16,col=dat$Recurrence)
model=randomForest(Recurrence~TumorMacrophageGSEA+Fuhrman.grade+Formal.staging+Time_without_recurrence,data=dat,ntree=10000,importance=T)
#model=randomForest(Recurrence~Fuhrman.grade+Formal.staging,data=dat,ntree=10000,importance=T)
rf.roc<-roc(dat$Recurrence,model$votes[,2],ci=T,ci.alpha=0.95,stratified=F,plot=T,print.auc=T)
plot.roc(rf.roc,print.auc=T,asp=NA)
