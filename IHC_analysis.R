library(phenoptr)
ihc_data=data.frame(patient=c("PatientA","PatientB","PatientC","Patient1","Patient2","Patient3","Patient4","Patient5","Patient6","Patient7","Patient8"),
                    stage=c("pT3a","pT1b","pT1b","pT1b","pT3a","pT3a","pT1b","pT3a","pT1b","pT3a","pT1b"),grade=c(3,2,2,2,2,1,1,4,2,3,2))

C1Q.macrophage.oddsratio=c()
TREM2.macrophage.oddsratio=c()
APOE.macrophage.oddsratio=c()

tumor_freqs=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other")
stroma_freqs=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other")
normal_freqs=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other")
normal_freqs=rbind(normal_freqs,rep(NaN,length(normal_freqs)))

#Patient A: stage 3, grade 3
#IHC_data_allsamples/Batch_Sp18_10879
#APOE Stroma 2.72 Tumor 2.02
#C1Q Stroma 3.02 Tumor 2.5
#CA9 Stroma 4.12 Tumor
#CD3 Stroma 0.28 Tumor
#CD68+CD163+ Stroma 3.2 Tumor 2.2
#TREM2 Stroma 0.42 Tumor 0.28
files=list_cell_seg_files("IHC_data_allsamples/Batch_Sp18_10879")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell C1q (Opal 520) Mean`>2.5),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell ApoE (Opal 620) Mean`>2.02),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell C1q (Opal 520) Mean`>2.5))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell ApoE (Opal 620) Mean`>2.02))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>4.12)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.28)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.5)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.02)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.5 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.02 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.5 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.5 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.02 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.02 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.5 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.02 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=csd
####
####

#Patient B: stage 1, grade 2
#IHC_data_allsamples/Batch_SP18-18677
#APOE Stroma 3.95 Tumor 4.4
#C1Q Stroma 3.44 Tumor 2.82
#CA9 Stroma Tumor
#CD3 Stroma Tumor
#CD68+CD163+ Stroma 3.24 Tumor 2.43
#TREM2 Stroma 0.605 Tumor 0.58
files=list_cell_seg_files("IHC_data_allsamples/Batch_SP18-18677")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell C1q (Opal 520) Mean`>2.82),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.58),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell ApoE (Opal 620) Mean`>4.4),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell C1q (Opal 520) Mean`>2.82))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.58))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell ApoE (Opal 620) Mean`>4.4))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>4.12)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.28)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.82)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.58)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>4.4)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.82 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.58 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>4.4 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.82 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.58 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.82 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.4 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.58 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.4 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.82 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.58 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.4 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=rbind(csd_distance_to_tumor,csd)
####
####

#Patient C: stage 1, grade 2
#IHC_data_allsamples/Batch_19282
#APOE Stroma 2.98 Tumor 3.11
#C1Q Stroma 2.73 Tumor 2.24
#CA9 Stroma Tumor
#CD3 Stroma Tumor
#CD68+CD163+ Stroma 1.8 Tumor 1.33
#TREM2 Stroma 0.3 Tumor 0.24
files=list_cell_seg_files("IHC_data_allsamples/Batch_19282")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.33,csd$`Entire Cell C1q (Opal 520) Mean`>2.73),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.33,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.24),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.33,csd$`Entire Cell ApoE (Opal 620) Mean`>3.11),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.33,csd$`Entire Cell C1q (Opal 520) Mean`>2.73))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.33,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.24))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.33,csd$`Entire Cell ApoE (Opal 620) Mean`>3.11))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>4.12)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.28)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.33)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.73)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.24)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>3.11)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.73 & csd$`Entire Cell Macs (Opal 570) Mean`>1.33)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.24 & csd$`Entire Cell Macs (Opal 570) Mean`>1.33)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>3.11 & csd$`Entire Cell Macs (Opal 570) Mean`>1.33)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.73 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.24 & csd$`Entire Cell Macs (Opal 570) Mean`>1.33)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.73 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.11 & csd$`Entire Cell Macs (Opal 570) Mean`>1.33)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.24 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.11 & csd$`Entire Cell Macs (Opal 570) Mean`>1.33)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.73 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.24 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.11 & csd$`Entire Cell Macs (Opal 570) Mean`>1.33)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent Normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=rbind(csd_distance_to_tumor,csd)
####
####

#Patient 1: stage 1, grade 2
#IHC_data_allsamples/Batch_SP19-994
#APOE Stroma 1.22 Tumor 1.53
#C1Q Stroma 2.74 Tumor 2.2
#CA9 Stroma 1.19 Tumor
#CD3 Stroma 0.16 Tumor
#CD68+CD163+ Stroma 1.02 Tumor 1.89
#TREM2 Stroma 0.392 Tumor 0.378
files=list_cell_seg_files("IHC_data_allsamples/Batch_SP19-994")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.89,csd$`Entire Cell C1q (Opal 520) Mean`>2.2),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.89,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.378),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.89,csd$`Entire Cell ApoE (Opal 620) Mean`>1.53),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.89,csd$`Entire Cell C1q (Opal 520) Mean`>2.2))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.89,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.378))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.89,csd$`Entire Cell ApoE (Opal 620) Mean`>1.53))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>1.19)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.16)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.89)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.2)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.378)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>1.53)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.2 & csd$`Entire Cell Macs (Opal 570) Mean`>1.89)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.378 & csd$`Entire Cell Macs (Opal 570) Mean`>1.89)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>1.53 & csd$`Entire Cell Macs (Opal 570) Mean`>1.89)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.2 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.378 & csd$`Entire Cell Macs (Opal 570) Mean`>1.89)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.2 & csd$`Entire Cell ApoE (Opal 620) Mean`>1.53 & csd$`Entire Cell Macs (Opal 570) Mean`>1.89)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.378 & csd$`Entire Cell ApoE (Opal 620) Mean`>1.53 & csd$`Entire Cell Macs (Opal 570) Mean`>1.89)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.2 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.378 & csd$`Entire Cell ApoE (Opal 620) Mean`>1.53 & csd$`Entire Cell Macs (Opal 570) Mean`>1.89)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=rbind(csd_distance_to_tumor,csd)
####
####

#Patient 2: stage 3, grade 2
#IHC_data_allsamples/Batch_SP19-1727
#APOE Stroma 4.09 Tumor 4.02
#C1Q Stroma 3.67 Tumor 2.53
#CA9 Stroma 2.1 Tumor
#CD3 Stroma 0.245 Tumor
#CD68+CD163+ Stroma 3.11 Tumor 1.47
#TREM2 Stroma 0.53 Tumor 0.29
files=list_cell_seg_files("IHC_data_allsamples/Batch_SP19-1727")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.47,csd$`Entire Cell C1q (Opal 520) Mean`>3.67),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.47,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.53),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.47,csd$`Entire Cell ApoE (Opal 620) Mean`>4.02),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.47,csd$`Entire Cell C1q (Opal 520) Mean`>3.67))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.47,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.53))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.47,csd$`Entire Cell ApoE (Opal 620) Mean`>4.02))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>2.1)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.245)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.47)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.67)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.53)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>4.02)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.67 & csd$`Entire Cell Macs (Opal 570) Mean`>1.47)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.53 & csd$`Entire Cell Macs (Opal 570) Mean`>1.47)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>4.02 & csd$`Entire Cell Macs (Opal 570) Mean`>1.47)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.67 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.53 & csd$`Entire Cell Macs (Opal 570) Mean`>1.47)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.67 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.02 & csd$`Entire Cell Macs (Opal 570) Mean`>1.47)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.53 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.02 & csd$`Entire Cell Macs (Opal 570) Mean`>1.47)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.67 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.53 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.02 & csd$`Entire Cell Macs (Opal 570) Mean`>1.47)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent Normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=rbind(csd_distance_to_tumor,csd)
####
####

#Patient 3: stage 3, grade 1
#IHC_data_allsamples/Batch_SP19-1994
#APOE Stroma 2.81 Tumor 2.71
#C1Q Stroma 2.62 Tumor 2.12
#CA9 Stroma 1.61 Tumor
#CD3 Stroma 0.3 Tumor
#CD68+CD163+ Stroma 2.57 Tumor 2.43
#TREM2 Stroma 0.42 Tumor 0.35
files=list_cell_seg_files("IHC_data_allsamples/Batch_SP19-1994")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43, csd$`Entire Cell C1q (Opal 520) Mean`>2.62),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.35),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell ApoE (Opal 620) Mean`>2.71),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell C1q (Opal 520) Mean`>2.62))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.35))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.43,csd$`Entire Cell ApoE (Opal 620) Mean`>2.71))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>1.61)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.3)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.62)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.35)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.71)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.62 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.35 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.71 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.62 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.35 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.62 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.71 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.35 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.71 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.62 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.35 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.71 & csd$`Entire Cell Macs (Opal 570) Mean`>2.43)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent Normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=rbind(csd_distance_to_tumor,csd)
####
####

#Patient 4: stage 1, grade 1
#IHC_data_allsamples/Batch_SP19-3695
#APOE Stroma 4.21 Tumor 5.41
#C1Q Stroma 3.89 Tumor 2.85
#CA9 Stroma 2.19 Tumor
#CD3 Stroma 0.21 Tumor
#CD68+CD163+ Stroma 1.13 Tumor 1.37
#TREM2 Stroma 0.585 Tumor 0.485
files=list_cell_seg_files("IHC_data_allsamples/Batch_SP19-3695")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.37, csd$`Entire Cell C1q (Opal 520) Mean`>3.89),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.37,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.485),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.37,csd$`Entire Cell ApoE (Opal 620) Mean`>5.41),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.37,csd$`Entire Cell C1q (Opal 520) Mean`>3.89))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.37,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.485))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.37,csd$`Entire Cell ApoE (Opal 620) Mean`>5.41))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>2.19)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.21)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.37)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.89)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.485)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>5.41)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.89 & csd$`Entire Cell Macs (Opal 570) Mean`>1.37)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.485 & csd$`Entire Cell Macs (Opal 570) Mean`>1.37)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>5.41 & csd$`Entire Cell Macs (Opal 570) Mean`>1.37)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.89 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.485 & csd$`Entire Cell Macs (Opal 570) Mean`>1.37)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.89 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.41 & csd$`Entire Cell Macs (Opal 570) Mean`>1.37)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.485 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.41 & csd$`Entire Cell Macs (Opal 570) Mean`>1.37)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.89 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.485 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.41 & csd$`Entire Cell Macs (Opal 570) Mean`>1.37)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent Normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=rbind(csd_distance_to_tumor,csd)
####
####

#Patient 5: stage 3, grade 4
#IHC_data_allsamples/Batch_SP19-9557
#APOE Stroma 3.47 Tumor 3.07
#C1Q Stroma 2.12 Tumor 1.77
#CA9 Stroma 1.44 Tumor
#CD3 Stroma 0.186 Tumor
#CD68+CD163+ Stroma 3.25 Tumor 2.61
#TREM2 Stroma 0.348 Tumor 0.316
files=list_cell_seg_files("IHC_data_allsamples/Batch_SP19-9557")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.61, csd$`Entire Cell C1q (Opal 520) Mean`>1.77),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.61,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.316),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.61,csd$`Entire Cell ApoE (Opal 620) Mean`>3.07),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.61,csd$`Entire Cell C1q (Opal 520) Mean`>1.77))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.61,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.316))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.61,csd$`Entire Cell ApoE (Opal 620) Mean`>3.07))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>1.44)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.186)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>2.61)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>1.77)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.316)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>3.07)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>1.77 & csd$`Entire Cell Macs (Opal 570) Mean`>2.61)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.316 & csd$`Entire Cell Macs (Opal 570) Mean`>2.61)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>3.07 & csd$`Entire Cell Macs (Opal 570) Mean`>2.61)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>1.77 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.316 & csd$`Entire Cell Macs (Opal 570) Mean`>2.61)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>1.77 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.07 & csd$`Entire Cell Macs (Opal 570) Mean`>2.61)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.316 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.07 & csd$`Entire Cell Macs (Opal 570) Mean`>2.61)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>1.77 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.316 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.07 & csd$`Entire Cell Macs (Opal 570) Mean`>2.61)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=rbind(csd_distance_to_tumor,csd)
####
####

#Patient 6: stage 1, grade 2
#IHC_data_allsamples/Batch_SP19-12625
#APOE Stroma 2.29 Tumor 2.24
#C1Q Stroma 2.48 Tumor 1.82
#CA9 Stroma 1.69 Tumor
#CD3 Stroma 0.11 Tumor
#CD68+CD163+ Stroma 1.96 Tumor 1.96
#TREM2 Stroma 0.83 Tumor 0.86
files=list_cell_seg_files("IHC_data_allsamples/Batch_SP19-12625")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.96, csd$`Entire Cell C1q (Opal 520) Mean`>1.82),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.96,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.86),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.96,csd$`Entire Cell ApoE (Opal 620) Mean`>2.24),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.96,csd$`Entire Cell C1q (Opal 520) Mean`>1.82))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.96,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.86))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.96,csd$`Entire Cell ApoE (Opal 620) Mean`>2.24))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>1.69)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.11)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.96)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>1.82)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.86)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.24)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>1.82 & csd$`Entire Cell Macs (Opal 570) Mean`>1.96)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.86 & csd$`Entire Cell Macs (Opal 570) Mean`>1.96)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.24 & csd$`Entire Cell Macs (Opal 570) Mean`>1.96)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>1.82 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.86 & csd$`Entire Cell Macs (Opal 570) Mean`>1.96)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>1.82 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.24 & csd$`Entire Cell Macs (Opal 570) Mean`>1.96)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.86 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.24 & csd$`Entire Cell Macs (Opal 570) Mean`>1.96)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>1.82 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.86 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.24 & csd$`Entire Cell Macs (Opal 570) Mean`>1.96)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=rbind(csd_distance_to_tumor,csd)
####
####

#Patient 7: stage 3, grade 3
#IHC_data_allsamples/Batch_SP19-13364
#APOE Stroma 2.87 Tumor 2.32
#C1Q Stroma 2.2 Tumor 2.08
#CA9 Stroma 1.39 Tumor
#CD3 Stroma 0.0755 Tumor
#CD68+CD163+ Stroma 1.39 Tumor 1.24
#TREM2 Stroma 0.1065 Tumor 0.0965
files=list_cell_seg_files("IHC_data_allsamples/Batch_SP19-13364")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.24, csd$`Entire Cell C1q (Opal 520) Mean`>2.2),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.24,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.0965),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.24,csd$`Entire Cell ApoE (Opal 620) Mean`>2.32),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.24,csd$`Entire Cell C1q (Opal 520) Mean`>2.2))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.24,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.0965))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.24,csd$`Entire Cell ApoE (Opal 620) Mean`>2.32))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>1.39)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.0755)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.24)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.2)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.0965)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.32)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.2 & csd$`Entire Cell Macs (Opal 570) Mean`>1.24)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.0965 & csd$`Entire Cell Macs (Opal 570) Mean`>1.24)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.32 & csd$`Entire Cell Macs (Opal 570) Mean`>1.24)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.2 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.0965 & csd$`Entire Cell Macs (Opal 570) Mean`>1.24)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.2 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.32 & csd$`Entire Cell Macs (Opal 570) Mean`>1.24)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.0965 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.32 & csd$`Entire Cell Macs (Opal 570) Mean`>1.24)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.2 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.0965 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.32 & csd$`Entire Cell Macs (Opal 570) Mean`>1.24)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent Normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=rbind(csd_distance_to_tumor,csd)
####
####

#Patient 8: stage 1, grade 2
#IHC_data_allsamples/Batch_SP19-14528
#APOE Stroma 4.98 Tumor 4.98
#C1Q Stroma 4.96 Tumor 4.56
#CA9 Stroma 2.77 Tumor
#CD3 Stroma 0.57 Tumor
#CD68+CD163+ Stroma 7.97 Tumor 7.35
#TREM2 Stroma 0.723 Tumor 0.623
files=list_cell_seg_files("IHC_data_allsamples/Batch_SP19-14528")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>7.35, csd$`Entire Cell C1q (Opal 520) Mean`>4.96),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>7.35,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.623),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>7.35,csd$`Entire Cell ApoE (Opal 620) Mean`>4.98),sep="\t",quote=F)
C1Q.macrophage.oddsratio=c(C1Q.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>7.35,csd$`Entire Cell C1q (Opal 520) Mean`>4.96))$estimate)
TREM2.macrophage.oddsratio=c(TREM2.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>7.35,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.623))$estimate)
APOE.macrophage.oddsratio=c(APOE.macrophage.oddsratio,fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>7.35,csd$`Entire Cell ApoE (Opal 620) Mean`>4.98))$estimate)
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>2.77)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.57)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>7.35)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>4.96)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.623)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>4.98)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>4.96 & csd$`Entire Cell Macs (Opal 570) Mean`>7.35)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.623 & csd$`Entire Cell Macs (Opal 570) Mean`>7.35)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>4.98 & csd$`Entire Cell Macs (Opal 570) Mean`>7.35)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>4.96 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.623 & csd$`Entire Cell Macs (Opal 570) Mean`>7.35)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>4.96 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.98 & csd$`Entire Cell Macs (Opal 570) Mean`>7.35)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.623 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.98 & csd$`Entire Cell Macs (Opal 570) Mean`>7.35)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>4.96 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.623 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.98 & csd$`Entire Cell Macs (Opal 570) Mean`>7.35)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}

####
####average distance between APOE+Macs to tumor cells vs APOE-Macs to tumor cells
csd$Phenotype=csd$celltype
#csd=csd[which(csd$`Tissue Category` == "Stroma"),]
distances=lapply(unique(csd$`Sample Name`),function(i){find_nearest_distance(csd[which(csd$`Sample Name`==i),],phenotypes="CA9+")})
d=distances[[1]]
for(i in 2:length(distances)){
  d=rbind(d,distances[[i]])
}
csd=cbind(csd,d)
csd$MacrophageType=as.character(csd$Phenotype)
#csd=csd[which(csd$MacrophageType %in% c("APOE+CD68/CD163+","C1Q+APOE+CD68/CD163+","C1Q+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","CD68/CD163+","TREM2+APOE+CD68/CD163+","TREM2+CD68/CD163+")),]
csd$MacrophageType[which(csd$MacrophageType %in% c("CD68/CD163+"))]="CD68/CD163+"
csd$MacrophageType[which(csd$MacrophageType %in% c("C1Q+TREM2+CD68/CD163+"))]="CD68/CD163+APOE+"
csd=csd[which(csd$MacrophageType %in% c("CD68/CD163+","CD68/CD163+APOE+")),]
ggplot(csd, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+")],csd$`Distance to CA9+`[which(csd$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor=rbind(csd_distance_to_tumor,csd)
ggplot(csd_distance_to_tumor, aes(`Distance to CA9+`, color=MacrophageType)) +
  geom_density(size=1)
ggplot(csd_distance_to_tumor, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_boxplot()
mean(na.exclude(csd_distance_to_tumor$`Distance to CA9+`[which(csd_distance_to_tumor$MacrophageType=="CD68/CD163+")]))
mean(na.exclude(csd_distance_to_tumor$`Distance to CA9+`[which(csd_distance_to_tumor$MacrophageType=="CD68/CD163+APOE+")]))
t.test(csd_distance_to_tumor$`Distance to CA9+`[which(csd_distance_to_tumor$MacrophageType=="CD68/CD163+")],csd_distance_to_tumor$`Distance to CA9+`[which(csd_distance_to_tumor$MacrophageType=="CD68/CD163+APOE+")])
csd_distance_to_tumor$`Distance to CA9+`=log10(csd_distance_to_tumor$`Distance to CA9+`)
ggplot(csd_distance_to_tumor, aes(MacrophageType,`Distance to CA9+`,color=MacrophageType)) +
  geom_violin()
####
####

###ggplot visualizations
#odds ratios of enrichment for each marker in CD68/CD163+ macrophages
#frequencies of each cell population in tumor, stromal, adjacent normal, grouped by stage or grade, and in individual patients
ihc_data$C1Q=C1Q.macrophage.oddsratio
ihc_data$TREM2=TREM2.macrophage.oddsratio
ihc_data$APOE=APOE.macrophage.oddsratio
colnames(stroma_freqs)=paste("Stroma",stroma_freqs[1,],sep="_")
colnames(tumor_freqs)=paste("Tumor",tumor_freqs[1,],sep="_")
colnames(normal_freqs)=paste("Normal",normal_freqs[1,],sep="_")
ihc_data=cbind(ihc_data,stroma_freqs[2:nrow(stroma_freqs),],tumor_freqs[2:nrow(tumor_freqs),],normal_freqs[3:nrow(normal_freqs),])
ihc_data$grade=as.character(ihc_data$grade)

theme_update(plot.title = element_text(hjust = 0.5))

x=ihc_data
p=ggplot(x, aes(x=stage,y=C1Q,fill=stage)) +
  geom_boxplot() + ylim(1,NA) + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("C1Q+ Odds Ratio in Macrophages vs Non-Macrophages")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())+ geom_hline(yintercept=1, linetype="dashed", color = "red")

p=ggplot(x, aes(x=stage,y=TREM2,fill=stage)) +
  geom_boxplot() + ylim(1,NA) + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("TREM2+ Odds Ratio in Macrophages vs Non-Macrophages")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())+ geom_hline(yintercept=1, linetype="dashed", color = "red")

p=ggplot(x, aes(x=stage,y=APOE,fill=stage)) +
  geom_boxplot() + ylim(1,NA) + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("APOE+ Odds Ratio in Macrophages vs Non-Macrophages")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())+ geom_hline(yintercept=1, linetype="dashed", color = "red")

x=ihc_data[,c(2,4,5,6)]
library(reshape)
x <- melt(x, id=c("stage"))
colnames(x)=c("stage","marker","odds_ratio")
p=ggplot(x, aes(x=marker,y=odds_ratio,fill=marker)) +
  geom_boxplot() + ylim(1,NA) + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Novel Marker Odds Ratios in Macrophages vs Non-Macrophages")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())+ geom_hline(yintercept=1, linetype="dashed", color = "red")


x=ihc_data[,c(2,9:16,23:30,37:44)]   #rows c(2,3,4,7,8,10,11)
x <- melt(x, id=c("stage"))
x$tissue="Stroma"
x$tissue[(nrow(x)/3+1):(2*nrow(x)/3)]="Tumor"
x$tissue[(2*nrow(x)/3+1):(nrow(x))]="Adjacent Normal"
x=x[which(x$tissue!="Tumor"),]
x$tissue[which(x$tissue=="Stroma")]="Tumor Stroma"
x$phenotype=unlist(lapply(strsplit(as.character(x$variable),"_"),function(y){y[2]}))
colnames(x)=c("stage","label","frequency","tissue","phenotype")
x$frequency[which(x$frequency=="NaN")]=NA
x$frequency=as.numeric(as.character(x$frequency))
x$phenotype=factor(x$phenotype,levels=c("CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+"))
p=ggplot(x, aes(x=phenotype,y=frequency,fill=tissue)) +
  geom_boxplot() + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Population Frequencies in Tumor vs Normal")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank()) #+ scale_y_continuous(trans='log10')
y=x[which(x$phenotype=="C1Q+CD68/CD163+"),]
p=ggplot(y, aes(x=stage,y=frequency,fill=tissue)) +
  geom_boxplot() + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("C1Q+CD68/CD163+ Frequencies in Tumor vs Normal")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())
y=x[which(x$phenotype=="TREM2+CD68/CD163+"),]
p=ggplot(y, aes(x=stage,y=frequency,fill=tissue)) +
  geom_boxplot() + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("TREM2+CD68/CD163+ Frequencies in Tumor vs Normal")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())
y=x[which(x$phenotype=="APOE+CD68/CD163+"),]
p=ggplot(y, aes(x=stage,y=frequency,fill=tissue)) +
  geom_boxplot() + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("APOE+CD68/CD163+ Frequencies in Tumor vs Normal")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())
y=x[which(x$phenotype=="C1Q+TREM2+CD68/CD163+"),]
p=ggplot(y, aes(x=stage,y=frequency,fill=tissue)) +
  geom_boxplot() + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("C1Q+TREM2+CD68/CD163+ Frequencies in Tumor vs Normal")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())
y=x[which(x$phenotype=="C1Q+APOE+CD68/CD163+"),]
p=ggplot(y, aes(x=stage,y=frequency,fill=tissue)) +
  geom_boxplot() + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("C1Q+APOE+CD68/CD163+ Frequencies in Tumor vs Normal")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())
y=x[which(x$phenotype=="TREM2+APOE+CD68/CD163+"),]
p=ggplot(y, aes(x=stage,y=frequency,fill=tissue)) +
  geom_boxplot() + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("TREM2+APOE+CD68/CD163+ Frequencies in Tumor vs Normal")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())
y=x[which(x$phenotype=="C1Q+TREM2+APOE+CD68/CD163+"),]
p=ggplot(y, aes(x=stage,y=frequency,fill=tissue)) +
  geom_boxplot() + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("C1Q+TREM2+APOE+CD68/CD163+ Frequencies in Tumor vs Normal")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank())

x=ihc_data[c(2,3,5,6,7,10,11),c(9:16,37:44)]
x=apply(x,2,function(y){as.numeric(as.character(y))})
write.table(apply(x[,1:8]-x[,9:16],2,mean),sep="\t",quote=F)
wilcox.test(x[,1],x[,9],paired=T)
wilcox.test(x[,2],x[,10],paired=T)
wilcox.test(x[,3],x[,11],paired=T)
wilcox.test(x[,4],x[,12],paired=T)
wilcox.test(x[,5],x[,13],paired=T)
wilcox.test(x[,6],x[,14],paired=T)
wilcox.test(x[,7],x[,15],paired=T)
wilcox.test(x[,8],x[,16],paired=T)






####
###Vinson COHORT####
####
recurrence_annotation_ihc=data.frame(sample=c("VC001","VC002","VC003","VC004","VC006","VC009","VC010","VC011"),
                                     recurrence=c(F,F,F,F,T,T,T,T),
                                     Gender=c("F","M","M","F","F","M","M","F"),
                                     Tumor_size_cm=c(8.3,9,4,5.2,8,10,15,10),
                                     Grade=c("1","2","2/3","2","3","2/3","3","3"),
                                     Margins=rep("Negative",8),
                                     Stage=c("pT2","pT3b","pT3b","pT3b","pT2","pT3a","pT3a","pT3a"),
                                     Time_to_recurrence_months=c(110,86,113,35,12,82,8,12))
tumor_freqs=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other")
stroma_freqs=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other")
normal_freqs=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other")
normal_freqs=rbind(normal_freqs,rep(NaN,length(normal_freqs)))

#Patient RCC2, VC001, no recurrence, 110months
#IHC_data_allsamples/Vinson_cohort/Batch_RCC2
#APOE Stroma 3.72 Tumor 2.72  
#C1Q Stroma 3.79 Tumor 1.89
#CA9 Stroma 2.19 Tumor 
#CD3 Stroma 0.06 Tumor 
#CD68+CD163+ Stroma 1 Tumor 1.427
#TREM2 Stroma 0.26 Tumor 0.2
files=list_cell_seg_files("IHC_data_allsamples/Vinson_cohort/Batch_RCC2")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.427, csd$`Entire Cell C1q (Opal 520) Mean`>3.79),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.427,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.26),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.427,csd$`Entire Cell ApoE (Opal 620) Mean`>2.72),sep="\t",quote=F)
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.427,csd$`Entire Cell C1q (Opal 520) Mean`>3.79))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.427,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.26))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.427,csd$`Entire Cell ApoE (Opal 620) Mean`>2.72))
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>2.19)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.06)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.427)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.79)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.26)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.72)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.79 & csd$`Entire Cell Macs (Opal 570) Mean`>1.427)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.26 & csd$`Entire Cell Macs (Opal 570) Mean`>1.427)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.72 & csd$`Entire Cell Macs (Opal 570) Mean`>1.427)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.79 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.26 & csd$`Entire Cell Macs (Opal 570) Mean`>1.427)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.79 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.72 & csd$`Entire Cell Macs (Opal 570) Mean`>1.427)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.26 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.72 & csd$`Entire Cell Macs (Opal 570) Mean`>1.427)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.79 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.26 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.72 & csd$`Entire Cell Macs (Opal 570) Mean`>1.427)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}


#Patient RCC3, VC002, no recurrence, 86 months
#IHC_data_allsamples/Vinson_cohort/Batch_RCC3
#APOE Stroma 4.88 Tumor 5.88  
#C1Q Stroma 3.61 Tumor 3.07
#CA9 Stroma 2.05 Tumor 
#CD3 Stroma 0.107 Tumor 
#CD68+CD163+ Stroma 1.36 Tumor 1.86
#TREM2 Stroma 0.278 Tumor 0.28
files=list_cell_seg_files("IHC_data_allsamples/Vinson_cohort/Batch_RCC3")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.86, csd$`Entire Cell C1q (Opal 520) Mean`>3.61),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.86,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.86,csd$`Entire Cell ApoE (Opal 620) Mean`>5.88),sep="\t",quote=F)
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.86,csd$`Entire Cell C1q (Opal 520) Mean`>3.61))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.86,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.86,csd$`Entire Cell ApoE (Opal 620) Mean`>5.88))
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>2.05)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.107)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.86)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.61)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.623)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>5.88)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.61 & csd$`Entire Cell Macs (Opal 570) Mean`>1.86)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell Macs (Opal 570) Mean`>1.86)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>5.88 & csd$`Entire Cell Macs (Opal 570) Mean`>1.86)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.61 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell Macs (Opal 570) Mean`>1.86)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.61 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.88 & csd$`Entire Cell Macs (Opal 570) Mean`>1.86)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.88 & csd$`Entire Cell Macs (Opal 570) Mean`>1.86)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.61 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.88 & csd$`Entire Cell Macs (Opal 570) Mean`>1.86)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}


#Patient RCC4, VC003, no recurrence, 113 months
#IHC_data_allsamples/Vinson_cohort/Batch_RCC4
#APOE Stroma 5.52 Tumor 5.3
#C1Q Stroma 2.95 Tumor 2.35
#CA9 Stroma 1.29 Tumor 
#CD3 Stroma 0.18 Tumor 
#CD68+CD163+ Stroma 3.23 Tumor 2.36
#TREM2 Stroma 0.87 Tumor 0.56
files=list_cell_seg_files("IHC_data_allsamples/Vinson_cohort/Batch_RCC4")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.36, csd$`Entire Cell C1q (Opal 520) Mean`>2.95),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.36,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.87),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.36,csd$`Entire Cell ApoE (Opal 620) Mean`>5.3),sep="\t",quote=F)
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.36,csd$`Entire Cell C1q (Opal 520) Mean`>2.95))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.36,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.87))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.36,csd$`Entire Cell ApoE (Opal 620) Mean`>5.3))
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>1.29)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.18)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>2.36)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.95)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.87)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>5.3)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.95 & csd$`Entire Cell Macs (Opal 570) Mean`>2.36)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.87 & csd$`Entire Cell Macs (Opal 570) Mean`>2.36)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>5.3 & csd$`Entire Cell Macs (Opal 570) Mean`>2.36)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.95 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.87 & csd$`Entire Cell Macs (Opal 570) Mean`>2.36)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.95 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.3 & csd$`Entire Cell Macs (Opal 570) Mean`>2.36)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.87 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.3 & csd$`Entire Cell Macs (Opal 570) Mean`>2.36)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.95 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.87 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.3 & csd$`Entire Cell Macs (Opal 570) Mean`>2.36)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}


#Patient RCC5, VC004, no recurrence, 35 months
#IHC_data_allsamples/Vinson_cohort/Batch_RCC5
#APOE Stroma 2.2 Tumor 3.59 
#C1Q Stroma 2.464 Tumor 3.51
#CA9 Stroma 1.85 Tumor 
#CD3 Stroma 0.12 Tumor 
#CD68+CD163+ Stroma 1.95 Tumor 2.2
#TREM2 Stroma 0.25 Tumor 0.2
files=list_cell_seg_files("IHC_data_allsamples/Vinson_cohort/Batch_RCC5")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2, csd$`Entire Cell C1q (Opal 520) Mean`>3.51),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.25),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell ApoE (Opal 620) Mean`>3.59),sep="\t",quote=F)
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell C1q (Opal 520) Mean`>3.51))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.25))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.2,csd$`Entire Cell ApoE (Opal 620) Mean`>3.59))
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>1.85)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.12)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.51)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.25)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>3.59)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.51 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.25 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>3.59 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.51 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.25 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.51 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.59 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.25 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.59 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.51 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.25 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.59 & csd$`Entire Cell Macs (Opal 570) Mean`>2.2)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}


#Patient RCC7, VC006, recurrence, 12 months
#IHC_data_allsamples/Vinson_cohort/Batch_RCC7
#APOE Stroma 5.95 Tumor 5.15
#C1Q Stroma 3.48 Tumor 2.48
#CA9 Stroma 2.09 Tumor 
#CD3 Stroma 0.265 Tumor 
#CD68+CD163+ Stroma 1.8 Tumor 1.5
#TREM2 Stroma 0.49 Tumor 0.42
files=list_cell_seg_files("IHC_data_allsamples/Vinson_cohort/Batch_RCC7")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.5, csd$`Entire Cell C1q (Opal 520) Mean`>2.48),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.5,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.42),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.5,csd$`Entire Cell ApoE (Opal 620) Mean`>5.15),sep="\t",quote=F)
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.5,csd$`Entire Cell C1q (Opal 520) Mean`>2.48))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.5,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.42))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.5,csd$`Entire Cell ApoE (Opal 620) Mean`>5.15))
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>2.09)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.265)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.5)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.48)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.42)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>5.15)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.48 & csd$`Entire Cell Macs (Opal 570) Mean`>1.5)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.42 & csd$`Entire Cell Macs (Opal 570) Mean`>1.5)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>5.15 & csd$`Entire Cell Macs (Opal 570) Mean`>1.5)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.48 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.42 & csd$`Entire Cell Macs (Opal 570) Mean`>1.5)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.48 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.15 & csd$`Entire Cell Macs (Opal 570) Mean`>1.5)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.42 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.15 & csd$`Entire Cell Macs (Opal 570) Mean`>1.5)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.48 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.42 & csd$`Entire Cell ApoE (Opal 620) Mean`>5.15 & csd$`Entire Cell Macs (Opal 570) Mean`>1.5)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}


#Patient RCC10, VC009, recurrence, 82 months
#IHC_data_allsamples/Vinson_cohort/Batch_RCC10
#APOE Stroma 4.12 Tumor 2.98
#C1Q Stroma 4.35 Tumor 3.85
#CA9 Stroma 2.43 Tumor 
#CD3 Stroma 0.12 Tumor 
#CD68+CD163+ Stroma 2.25 Tumor 2.25
#TREM2 Stroma 0.37 Tumor 0.28
files=list_cell_seg_files("IHC_data_allsamples/Vinson_cohort/Batch_RCC10")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.25, csd$`Entire Cell C1q (Opal 520) Mean`>3.85),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.25,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.25,csd$`Entire Cell ApoE (Opal 620) Mean`>2.98),sep="\t",quote=F)
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.25,csd$`Entire Cell C1q (Opal 520) Mean`>3.85))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.25,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>2.25,csd$`Entire Cell ApoE (Opal 620) Mean`>2.98))
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>2.43)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.12)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>2.25)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.85)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.98)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.85 & csd$`Entire Cell Macs (Opal 570) Mean`>2.25)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell Macs (Opal 570) Mean`>2.25)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>2.98 & csd$`Entire Cell Macs (Opal 570) Mean`>2.25)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.85 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell Macs (Opal 570) Mean`>2.25)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.85 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.98 & csd$`Entire Cell Macs (Opal 570) Mean`>2.25)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.98 & csd$`Entire Cell Macs (Opal 570) Mean`>2.25)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>3.85 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.28 & csd$`Entire Cell ApoE (Opal 620) Mean`>2.98 & csd$`Entire Cell Macs (Opal 570) Mean`>2.25)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}


#Patient RCC11, VC010, recurrence, 8 months
#IHC_data_allsamples/Vinson_cohort/Batch_RCC11
#APOE Stroma 7.15 Tumor 4.86
#C1Q Stroma 3.96 Tumor 2.96
#CA9 Stroma 2.29 Tumor 
#CD3 Stroma 0.15 Tumor 
#CD68+CD163+ Stroma 2.19 Tumor 1.71
#TREM2 Stroma 0.71 Tumor 0.63
files=list_cell_seg_files("IHC_data_allsamples/Vinson_cohort/Batch_RCC11")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.71, csd$`Entire Cell C1q (Opal 520) Mean`>2.96),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.71,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.63),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.71,csd$`Entire Cell ApoE (Opal 620) Mean`>4.86),sep="\t",quote=F)
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.71,csd$`Entire Cell C1q (Opal 520) Mean`>2.96))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.71,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.63))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.71,csd$`Entire Cell ApoE (Opal 620) Mean`>4.86))
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>2.29)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.15)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.71)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.96)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.63)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>4.86)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.96 & csd$`Entire Cell Macs (Opal 570) Mean`>1.71)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.63 & csd$`Entire Cell Macs (Opal 570) Mean`>1.71)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>4.86 & csd$`Entire Cell Macs (Opal 570) Mean`>1.71)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.96 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.63 & csd$`Entire Cell Macs (Opal 570) Mean`>1.71)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.96 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.86 & csd$`Entire Cell Macs (Opal 570) Mean`>1.71)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.63 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.86 & csd$`Entire Cell Macs (Opal 570) Mean`>1.71)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.96 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.63 & csd$`Entire Cell ApoE (Opal 620) Mean`>4.86 & csd$`Entire Cell Macs (Opal 570) Mean`>1.71)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}


#Patient RCC12, VC011, recurrence, 12 months
#IHC_data_allsamples/Vinson_cohort/Batch_RCC12
#APOE Stroma 5.498 Tumor 3.498
#C1Q Stroma 3.79 Tumor 2.92
#CA9 Stroma 1.52 Tumor 
#CD3 Stroma 0.119 Tumor 
#CD68+CD163+ Stroma 2.19 Tumor 1.4
#TREM2 Stroma 0.42 Tumor 0.41
files=list_cell_seg_files("IHC_data_allsamples/Vinson_cohort/Batch_RCC12")
csd=read_cell_seg_data(files[1])
for(i in 2:length(files)){
  t=read_cell_seg_data(files[i])
  csd=rbind(csd,t)
}
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.4, csd$`Entire Cell C1q (Opal 520) Mean`>2.92),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.4,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.41),sep="\t",quote=F)
write.table(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.4,csd$`Entire Cell ApoE (Opal 620) Mean`>3.498),sep="\t",quote=F)
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.4,csd$`Entire Cell C1q (Opal 520) Mean`>2.92))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.4,csd$`Entire Cell TREM2 (Opal 650) Mean`>0.41))
fisher.test(table(csd$`Entire Cell Macs (Opal 570) Mean`>1.4,csd$`Entire Cell ApoE (Opal 620) Mean`>3.498))
csd$celltype="Other"
csd$celltype[which(csd$`Entire Cell CA9 (Opal 540) Mean`>1.52)]="CA9+"
csd$celltype[which(csd$`Entire Cell CD3 (Opal 690) Mean`>0.119)]="CD3+"
csd$celltype[which(csd$`Entire Cell Macs (Opal 570) Mean`>1.4)]="CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.92)]="C1Q+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.41)]="TREM2+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>3.498)]="APOE+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.92 & csd$`Entire Cell Macs (Opal 570) Mean`>1.4)]="C1Q+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.41 & csd$`Entire Cell Macs (Opal 570) Mean`>1.4)]="TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell ApoE (Opal 620) Mean`>3.498 & csd$`Entire Cell Macs (Opal 570) Mean`>1.4)]="APOE+CD68/CD163+"
#
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.92 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.41 & csd$`Entire Cell Macs (Opal 570) Mean`>1.4)]="C1Q+TREM2+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.92 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.498 & csd$`Entire Cell Macs (Opal 570) Mean`>1.4)]="C1Q+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell TREM2 (Opal 650) Mean`>0.41 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.498 & csd$`Entire Cell Macs (Opal 570) Mean`>1.4)]="TREM2+APOE+CD68/CD163+"
csd$celltype[which(csd$`Entire Cell C1q (Opal 520) Mean`>2.92 & csd$`Entire Cell TREM2 (Opal 650) Mean`>0.41 & csd$`Entire Cell ApoE (Opal 620) Mean`>3.498 & csd$`Entire Cell Macs (Opal 570) Mean`>1.4)]="C1Q+TREM2+APOE+CD68/CD163+"
csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
#csd$celltype=factor(csd$celltype,levels=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other"))
table(csd$celltype,csd$`Tissue Category`)
write.table(table(csd$celltype,csd$`Tissue Category`),sep="\t",quote=F)
x=as.matrix(table(csd$celltype,csd$`Tissue Category`))
x=as.data.frame(apply(x,2,function(y){y/sum(y)}))
write.table(x,sep="\t",quote=F)
stroma_freqs=rbind(stroma_freqs,x$Stroma)
tumor_freqs=rbind(tumor_freqs,x$Tumor)
if(ncol(x)==3){normal_freqs=rbind(normal_freqs,x$`Adjacent normal`)} else{normal_freqs=rbind(normal_freqs,rep(NaN,ncol(normal_freqs)))}


###ggplot visualizations
#frequencies of each cell population in tumor, stromal, adjacent normal, grouped by stage or grade, and in individual patients
colnames(stroma_freqs)=paste("Stroma",stroma_freqs[1,],sep="_")
colnames(tumor_freqs)=paste("Tumor",tumor_freqs[1,],sep="_")
recurrence_annotation_ihc=cbind(recurrence_annotation_ihc,stroma_freqs[2:nrow(stroma_freqs),])
colnames(recurrence_annotation_ihc)[9:22]=c("CA9+","CD3+","CD68/CD163+","C1Q+CD68/CD163+","TREM2+CD68/CD163+","APOE+CD68/CD163+","C1Q+TREM2+CD68/CD163+","C1Q+APOE+CD68/CD163+","TREM2+APOE+CD68/CD163+","C1Q+TREM2+APOE+CD68/CD163+","C1Q+","TREM2+","APOE+","Other")

x=recurrence_annotation_ihc[,c(2,11:21)]
x <- melt(x, id=c("recurrence"))
colnames(x)=c("recurrence","phenotype","frequency")
x$frequency=as.numeric(as.character(x$frequency))
p=ggplot(x, aes(x=phenotype,y=frequency,fill=recurrence)) +
  geom_boxplot() + #geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Tumor Frequencies in Recurrence vs Non-Recurrence")
p + theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=10, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          axis.title.x = element_blank()) #+ scale_y_continuous(trans='log10')
apply(recurrence_annotation_ihc[,11:21],2,function(x){wilcox.test(as.numeric(as.character(x[1:4])),as.numeric(as.character(x[5:8])))$p.value})

#survival analysis
library(survival)
library(survminer)
recurrence_annotation_ihc[,9:22]=apply(recurrence_annotation_ihc[,9:22],2,function(x){as.numeric(as.character(x))})
recurrence_annotation_ihc$SurvObj=with(recurrence_annotation_ihc,Surv(Time_to_recurrence_months,recurrence))
res.cox1 <- coxph(SurvObj ~ ., data =  recurrence_annotation_ihc[,c(19,23)])
summary(res.cox1)
ggforest(res.cox1, data = recurrence_annotation_ihc,fontsize = 1,main = "C1Q+ Hazard Ratio")
res.cox1 <- coxph(SurvObj ~ ., data =  recurrence_annotation_ihc[,c(12,23)])
summary(res.cox1)
ggforest(res.cox1, data = recurrence_annotation_ihc,fontsize = 1,main = "C1Q+CD68/CD163+ Hazard Ratio")
#optimal cut-point by maximally selected rank statistic
res.cut <- surv_cutpoint(recurrence_annotation_ihc, time = "Time_to_recurrence_months", event = "recurrence",variables = c("C1Q+"))
summary(res.cut)
plot(res.cut, "C1Q+", palette = "npg")
recurrence_annotation_ihc$C1Q_fraction="low"
recurrence_annotation_ihc$C1Q_fraction[which(recurrence_annotation_ihc$`C1Q+`>res.cut$cutpoint$cutpoint)]="high"
fit <- survfit(SurvObj ~ C1Q_fraction, data = recurrence_annotation_ihc, conf.type = "log-log")
ggsurvplot(fit, data = recurrence_annotation_ihc, pval = T,ylab="Fraction Without Recurrence",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(90, 0.9))
#optimal cut-point by maximally selected rank statistic
res.cut <- surv_cutpoint(recurrence_annotation_ihc, time = "Time_to_recurrence_months", event = "recurrence",variables = c("C1Q+CD68/CD163+"))
summary(res.cut)
plot(res.cut, "C1Q+CD68/CD163+", palette = "npg")
recurrence_annotation_ihc$C1Q_CD68CD163_fraction="low"
recurrence_annotation_ihc$C1Q_CD68CD163_fraction[which(recurrence_annotation_ihc$`C1Q+CD68/CD163+`>res.cut$cutpoint$cutpoint)]="high"
fit <- survfit(SurvObj ~ C1Q_CD68CD163_fraction, data = recurrence_annotation_ihc, conf.type = "log-log")
ggsurvplot(fit, data = recurrence_annotation_ihc, pval = T,ylab="Fraction Without Recurrence",xlab="Time (months)",linetype=c("dashed","solid"),palette="jco",risk.table=T,pval.coord = c(90, 0.9))
