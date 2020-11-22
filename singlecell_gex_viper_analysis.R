#Stage annotations
#"pT3a"-- PatientA, Patient2, Patient3, Patient5(metastatic), Patient7
#"pT1b"-- PatientB, PatientC, Patient1, Patient4, Patient6, Patient8

#Grade annotations
#"1"-- Patient3, Patient4, 
#"2"-- PatientB, PatientC, Patient1, Patient2, Patient6, Patient8
#"3"-- PatientA, Patient7
#"4"-- Patient5
#######

set.seed(1234)
library(reticulate)
library(SingleR)
library(Seurat)
library(cowplot)
library(dplyr)
library(iterClust)
library(cluster)
library(umap)
library(reshape)
library(pheatmap)
library(viper)
library("org.Hs.eg.db")
library(clustree)
library(factoextra)
library(leiden)
library(MAST)
library(Hmisc)
library(ggplot2)
library(scales)
library(flowCore)
library(ggcyto)
library(infercnv)
library(ggrepel)
library(plyr)

##bug-fixed version of CreateBigSingleRObject from singleR pipeline. 
#in singleR version the internal function would default to fine.tune=T regardless of input parameter settings
CreateBigSingleRObjectv2=function (counts, annot = NULL, project.name, xy, clusters, N = 10000, 
                                   min.genes = 200, technology = "10X", species = "Human", citation = "", 
                                   ref.list = list(), normalize.gene.length = F, variable.genes = "de", 
                                   fine.tune = T, reduce.file.size = T, do.signatures = F, do.main.types = T, 
                                   temp.dir = getwd(), numCores = SingleR.numCores){
  n = ncol(counts)
  s = seq(1, n, by = N)
  dir.create(paste0(temp.dir, "/singler.temp/"), showWarnings = FALSE)
  for (i in s) {
    print(i)
    A = seq(i, min(i + N - 1, n))
    singler = CreateSinglerObject(counts[, A], annot = annot[A], 
                                  project.name = project.name, min.genes = min.genes, 
                                  technology = technology, species = species, citation = citation, 
                                  do.signatures = do.signatures, clusters = NULL, numCores = numCores,fine.tune=fine.tune)
    save(singler, file = paste0(temp.dir, "/singler.temp/", project.name, ".", i, ".RData"))
  }
  singler.objects.file <- list.files(paste0(temp.dir, "/singler.temp/"), pattern = "RData", full.names = T)
  singler.objects = list()
  for (i in 1:length(singler.objects.file)) {
    load(singler.objects.file[[i]])
    singler.objects[[i]] = singler}
  singler = SingleR.Combine(singler.objects, order = colnames(counts), clusters = clusters, xy = xy)
  singler
}

#takes a data matrix mat and clust as input, where each column of clust is a unique louvain clustering
#subsamples 1000 cells 100 times to compute silhouette score. 
#outputs list of mean silhouette scores and standard deviations of silhouette scores for each clustering. 
sil_subsample_v2=function(mat,clust){
  out=as.data.frame(matrix(rep(NA,100*ncol(clust)),nrow=100))
  for(x in 1:100){
    i=sample(1:ncol(mat),min(1000,ncol(mat)))
    d=as.dist(1 - cor(mat[,i], method = "pearson"))
    for(j in 1:ncol(clust)){
      if(length(table(clust[i,j]))==1){out[x,j]=0}
      else{
        sil=silhouette(as.numeric(clust[i,j]),d)
        out[x,j]=mean(sil[, "sil_width"])}}
  }
  means=apply(out,2,mean)
  sd=apply(out,2,sd)
  return(list(means,sd))
}

#takes a data matrix mat and clust as input, where each column of clust is a unique louvain clustering
#subsamples 1000 cells 100 times to compute silhouette score. balanced subsampling to select proportional numbers of cells from each cluster
#outputs list of mean silhouette scores and standard deviations of silhouette scores for each clustering. 
sil_subsample_v3=function(mat,clust){
  out=as.data.frame(matrix(rep(NA,100*ncol(clust)),nrow=100))
  for(x in 1:100){
    for(j in 1:ncol(clust)){
      i=c()
      a=clust[,j]
      for(lab in unique(a)){
        i=c(i,sample((1:ncol(mat))[which(a==lab)],length(which(a==lab))/length(a)*min(1000,ncol(mat))))
      }
      d=as.dist(1 - cor(mat[,i], method = "pearson"))
      if(length(table(clust[i,j]))==1){out[x,j]=0}
      else{
        sil=silhouette(as.numeric(clust[i,j]),d)
        out[x,j]=mean(sil[, "sil_width"])}}
  }
  means=apply(out,2,median)
  sd=apply(out,2,mad)
  return(list(means,sd))
}

#Function to draw heatmap of genes in seurat object dat, with cluster labels clust
#seurat object must have labels "tissue" and "patient" 
geneHeatmap=function(dat,clust,genes,genes_by_cluster=T,n_top_genes_per_cluster=5,viper=F,color_palette=NA,scaled=F){
  identities <- levels(clust)
  if(is.na(color_palette)){my_color_palette <- hue_pal()(length(identities))}
  else{my_color_palette=color_palette}
  features=genes
  i=sample(1:ncol(dat),min(10000,ncol(dat)),replace = F)
  if(viper==F){x=dat[["SCT"]]@data[features,i]}
  if(viper==T){x=dat[["RNA"]]@data[features,i]}
  #df <- data.frame(clust[i],dat$tissue[i],dat$patient[i])
  df <- data.frame(clust[i],dat$tissue[i])
  rownames(df)=colnames(x)
  #colnames(df)=c("cluster","tissue","patient")
  colnames(df)=c("cluster","tissue")
  anno_colors <- list(cluster = my_color_palette,tissue=c("cornflowerblue","coral3"))
  names(anno_colors$cluster) <- levels(df$cluster)
  names(anno_colors$tissue)<-c("Normal","Tumor")
  #o=order(df$cluster,df$tissue,df$patient)
  o=order(df$cluster,df$tissue)
  x=x[,o]
  df=df[o,]
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  #mat_breaks <- quantile_breaks(as.matrix(apply(x,1,function(x){(x-mean(x))/sd(x)})), n = 30)
  if(scaled==F){t=as.matrix(apply(x,1,function(x){(x-mean(x))/sd(x)}))}
  if(scaled==T){t=as.matrix(x)}
  mat_breaks <- c(quantile_breaks(t[which(t<0)], n = 10),0,quantile_breaks(t[which(t>0)], n = 10))
  mat_breaks=mat_breaks[2:(length(mat_breaks)-1)] #restrict range of data to quantiles 5%-95%, extreme values excluded
  if(genes_by_cluster){
    anno_colors$group=anno_colors$cluster
    anno_row=data.frame(group=unlist(lapply(unique(df$cluster),function(x){rep(x,n_top_genes_per_cluster)})))
    gene_names=rownames(x)
    rownames(x)=1:nrow(x)
    if(!scaled){pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_row = anno_row,annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 10,show_colnames = F,annotation_colors = anno_colors,scale="row",gaps_row=(2:length(unique(clust))-1)*n_top_genes_per_cluster,annotation_names_row = F,labels_row=gene_names,row_annotation_legend=F)}
    if(scaled){pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_row = anno_row,annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 10,show_colnames = F,annotation_colors = anno_colors,gaps_row=(2:length(unique(clust))-1)*n_top_genes_per_cluster,annotation_names_row = F,labels_row=gene_names,row_annotation_legend=F)}
  }
  else{
    if(!scaled){pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 8,show_colnames = F,annotation_colors = anno_colors,scale="row")}
    if(scaled){pheatmap(x, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 8,show_colnames = F,annotation_colors = anno_colors)}
  }
}

#' Identifies MRs for given data using stouffer integration.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param cluster Vector of cluster lables. If not included, integrates the entire matrix.
#' @param weights A named vector of sample weights. If included, stouffer integration is weighted.
#' @return Returns the stouffer integrated scores for each protien.
StoufferMRs <- function(dat.mat, cluster, weights) {
  # generate dummy weights if missing
  if (missing(weights)) {
    weights <- as.numeric(rep(1, ncol(dat.mat))); names(weights) <- colnames(dat.mat)
  }
  # perform integration across full matrix if cluster was missing
  if (missing(cluster)) {
    sInt <- rowSums(t(t(dat.mat) * weights))
    sInt <- rowSums(t(t(dat.mat) * weights)) / sqrt(sum(weights ** 2))
    return(sInt)
  }
  # separate cluster specific matrices
  k <- length(table(cluster))
  mrs <- list()
  for (i in 1:k) { # for each cluster
    clust.cells <- names(cluster)[which(cluster == i)]
    clust.mat <- dat.mat[, clust.cells]
    clust.weights <- weights[clust.cells]
    clust.mrs <- StoufferMRs(clust.mat, weights = clust.weights)
    mrs[[paste('c', i, sep = '')]] <- sort(clust.mrs, decreasing = TRUE)
  }
  return(mrs)
}

#' Identifies MRs based on ANOVA analysis for a given clustering.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param clustering Clustering vector
#' @return A named vector of p-values for each protein
AnovaMRs <- function(dat.mat, clustering) {
  pVals <- c()
  group.vec <- clustering[colnames(dat.mat)]
  # perform an anova for each protein, storing pValues in a vector
  for (i in 1:nrow(dat.mat)) {
    aov.df <- data.frame('weights' = dat.mat[i,], 'group' = group.vec)
    #print(aov.df)
    aov.test <- aov(weights ~ group, aov.df)
    pVal <- summary(aov.test)[[1]][1,5]
    pVals <- c(pVals, pVal)
  }
  # name and return the vector
  names(pVals) <- rownames(dat.mat)
  return(pVals)
}

#' Performs a bootstrap t-test between two sample vectors x and y. Returns a log p-value.
#' 
#' @param x Vector of test values.
#' @param y Vector of reference values.
#' @param bootstrap.num Number of bootstraps to use. Default of 100.
#' @return A signed log p-value.
LogBootstrapTTest <- function(x, y, bootstrap.num = 100) {
  x.n <- length(x); y.n <- length(y)
  log.pValue <- c()
  ## perform test for each bootstrap
  for (i in 1:bootstrap.num) {
    # create bootstraps
    x.boot <- sample(1:x.n, size = x.n, replace = TRUE)
    x.boot <- x[x.boot]
    y.boot <- sample(1:y.n, size = y.n, replace = TRUE)
    y.boot <- y[y.boot]
    # perform t.test
    test.res <- t.test(x = x.boot, y = y.boot, alternative = "two.sided")
    # generate log p-value
    log.p <- 2*pt(q = abs(test.res$statistic), df = floor(test.res$parameter), log.p = TRUE, lower.tail = FALSE)*(-sign(test.res$statistic))
    log.pValue <- c(log.pValue, log.p)
  }
  # return mean log p-value
  return(mean(log.pValue))
}

#' Performs a t-test between two sample vectors x and y. Returns a log p-value.
#' @param x Vector of test values.
#' @param y Vector of reference values.
#' @param bootstrap.num Number of bootstraps to use. Default of 100.
#' @return A signed log p-value.
LogTTest <- function(x, y) {
  test.res <- t.test(x, y, alternative = 'two.sided')
  log.p <- 2*pt(q = abs(test.res$statistic), df = floor(test.res$parameter), log.p = TRUE, lower.tail = FALSE)*(-sign(test.res$statistic))
  return(log.p)
}

#' Identifies MRs based on a bootstraped Ttest between clusters.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param clustering Vector of cluster labels.
#' @param bootstrap.num Number of bootstraps to use. Default of 10 
#' @return Returns a list of lists; each list is a vector of sorted log p-values for each cluster.
BTTestMRs <- function(dat.mat, clustering, bootstrap.num = 100) {
  # set initial variables
  clustering <- clustering
  k <- length(table(clustering))
  mrs <- list()
  # identify MRs for each cluster
  for (i in 1:k) {
    print(paste('Identifying MRs for cluster ', i, '...', sep = ''))
    mrs.mat <- matrix(0L, nrow = nrow(dat.mat), ncol = bootstrap.num)
    rownames(mrs.mat) <- rownames(dat.mat)
    # split test and ref matrices
    clust <- names(table(clustering))[i]
    clust.vect <- which(clustering == clust)
    test.mat <- dat.mat[, clust.vect]; ref.mat <- dat.mat[, -clust.vect]
    t.n <- ncol(test.mat); r.n <- ncol(ref.mat)
    # for each bootstrap
    for (b in 1:bootstrap.num) {
      test.boot <- test.mat[, sample(colnames(test.mat), size = t.n, replace = TRUE)]
      ref.boot <- ref.mat[, sample(colnames(ref.mat), size = t.n, replace = TRUE)]
      # for each gene
      for (g in rownames(dat.mat)) {
        mrs.mat[g, b] <- LogTTest(test.boot[g,], ref.boot[g,])
      }
    }
    # sort and add to list
    mList <- sort(rowMeans(mrs.mat), decreasing = TRUE)
    mrs[[clust]] <- mList
  }
  # return 
  return(mrs)
}

#' Returns the master regulators for the given data.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param method 'Stouffer' or 'ANOVA'
#' @param clustering Optional argument for a vector of cluster labels.
#' @param numMRs Number of MRs to return per cluster. Default of 50.
#' @param bottom Switch to return downregulated proteins in MR list. Default FALSE>
#' @param weights Optional argument for weights, which can be used in the Stouffer method.
#' @return Returns a list of master regulators, or a list of lists if a clustring is specified.
GetMRs <- function(dat.mat, clustering, method, numMRs = 50, bottom = FALSE, weights, ...) {
  if (method == 'ANOVA') {
    mr.vals <- AnovaMRs(dat.mat, clustering)
  } else if (method == 'Stouffer') {
    # generate dummy weights if not specified
    if (missing(weights)) {
      weights <- rep(1, ncol(dat.mat))
      names(weights) <- colnames(dat.mat)
    }
    # recursive calls for each cluster
    if (missing(clustering)) { # no clustering specified
      mr.vals <- StoufferMRs(dat.mat, weights)
    } else {
      k <- length(table(clustering))
      mrs <- list()
      for (i in 1:k) {
        # get cluster specific matrix and weights
        clust.cells <- names(which(clustering == i))
        clust.mat <- dat.mat[, clust.cells]
        print(dim(clust.mat))
        clust.weights <- weights[clust.cells]
        # find mrs and add to list
        clust.mrs <- GetMRs(clust.mat, method = method, weights = clust.weights, numMRs = numMRs, bottom = bottom)
        print(head(clust.mrs))
        mrs[[paste('c', i, sep = '')]] <- clust.mrs
      }
      return(mrs)
    }
  } else {
    print('Invalid method: must be "Stouffer" or "ANOVA".')
  }
  # return appropriate portion of MR list
  mr.vals <- sort(mr.vals, decreasing = TRUE)
  if (bottom) {
    return(c(mr.vals[1:numMRs], tail(mr.vals, numMRs)))
  } else {
    return(mr.vals[1:numMRs])
  }
}

#' Identifies MRs on a cell-by-cell basis and returns a merged, unique list of all such MRs.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param numMRs Default number of MRs to identify in each cell. Default of 25.
#' @return Returns a list of master regulators, the unique, merged set from all cells.
CBCMRs <- function(dat.mat, numMRs = 25) {
  # identify MRs
  cbc.mrs <- apply(dat.mat, 2, function(x) { names(sort(x, decreasing = TRUE))[1:numMRs] })
  cbc.mrs <- unique(unlist(as.list(cbc.mrs)))
  # return
  return(cbc.mrs)
}

#' Unwraps a nested MR list: previous functions return cluster specific master regulators as a list of lists. This funciton will unwrap that object into one, unique list.
#'
#' @param MRs List of lists, with MR names as sub-list names and MR activity as sub-list entries.
#' @param top If specified, will subset the top X regulators from each set.
#' @return Returns a de-duplicated list of MRs.
MR_UnWrap <- function(MRs, top) {
  if (missing(top)) {
    return( unique(unlist(lapply(MRs, names), use.names = FALSE)) )
  } else {
    mr.unwrap <- lapply(MRs, function(x) {
      names(sort(x, decreasing = TRUE))[ 1:min(top, length(x)) ]
    })
    return( unique(unlist(mr.unwrap, use.names = FALSE)) )
  }
}

#' Performs a rank transformation on a given matrix.
#' 
#' @param dat.mat Matrix of data, usually gene expression (genes X samples).
#' @return Rank transformed matrix.
RankTransform <- function(dat.mat) {
  rank.mat <- apply(dat.mat, 2, rank)
  median <- apply(rank.mat, 1, median)
  mad <- apply(rank.mat, 1, mad)
  rank.mat <- (rank.mat - median) / mad
  return(rank.mat)
}

#' Make Cluster Metacells for ARACNe. Will take a clustering and produce saved meta cell matrices.
#'
#' @param dat.mat Matrix of raw gene expression (genes X samples).
#' @param dist.mat Distance matrix to be used for neighbor calculation. We recommend using a viper similarity matrix.
#' @param numNeighbors Number of neighbors to use for each meta cell. Default of 5.
#' @param clustering Vector of cluster labels. 
#' @param subSize Size to subset the data too. Since 200 cells is adequate for ARACNe runs, this allows for speedup. Default of 200.
#' @param out.dir Directory for sub matrices to be saved in.
#' @param out.name Optional argument for preface of file names. 
MakeCMfA <- function(dat.mat, numNeighbors = 10, clustering, subSize = 200, out.dir, out.name = '',sizeThresh=100) {
  # generate cluster matrices
  clust.mats <- ClusterMatrices(dat.mat, clustering,sizeThresh = sizeThresh)
  clust.mats=clust.mats[which(!unlist(lapply(clust.mats,is.null)))]
  # produce metaCell matrix and save for each cluster matrix
  k <- length(clust.mats)
  meta.mats <- list()
  for (i in 1:k) {
    mat <- clust.mats[[i]]
    meta.mat <- MetaCells(mat, numNeighbors)
    file.name <- paste(out.dir, out.name, '_clust-', i, '-metaCells_all', sep = '')
    ARACNeTable(meta.mat, file.name, subset = FALSE)
    meta.mats[[i]] <- meta.mat
    if(subSize < ncol(meta.mat)){meta.mat=meta.mat[, sample(colnames(meta.mat), subSize)]}
    meta.mat <- CPMTransform(meta.mat)
    file.name <- paste(out.dir, out.name, '_clust-', i, '-metaCells', sep = '')
    ARACNeTable(meta.mat, file.name, subset = FALSE)
  }
  return(meta.mats)
}

#' Generates a meta cell matrix for given data.
#' 
#' @param dat.mat Raw gene expression matrix (genes X samples).
#' @param dist.mat Distance matrix to be used for neighbor inference.
#' @param numNeighbors Number of neighbors to use for each meta cell. Default of 10.
#' @param subSize If specified, number of metaCells to be subset from the final matrix. No subsetting occurs if not incldued.
#' @return A matrix of meta cells (genes X samples).
MetaCells <- function(dat.mat, numNeighbors = 10, subSize) {
  # prune distance matrix if necessary
  #dist.mat <- as.matrix(dist.mat)
  #dist.mat <- dist.mat[colnames(dat.mat), colnames(dat.mat)]
  #dist.mat <- as.dist(dist.mat)
  dist.mat=as.dist(1 - cor(dat.mat, method = "pearson"))
  # KNN function
  KNN <- function(dist.mat, k){
    dist.mat <- as.matrix(dist.mat)
    n <- nrow(dist.mat)
    neighbor.mat <- matrix(0L, nrow = n, ncol = k)
    for (i in 1:n) {
      neighbor.mat[i,] <- order(dist.mat[i,])[2:(k + 1)]
    }
    return(neighbor.mat)
  }
  knn.neighbors <- KNN(dist.mat, numNeighbors)
  # create imputed matrix
  imp.mat <- matrix(0, nrow = nrow(dat.mat), ncol = ncol(dat.mat))
  rownames(imp.mat) <- rownames(dat.mat); colnames(imp.mat) <- colnames(dat.mat)
  for (i in 1:ncol(dat.mat)) {
    neighbor.mat <- dat.mat[,c(i, knn.neighbors[i,])]
    imp.mat[,i] <- rowSums(neighbor.mat)
  }
  # subset if requested and return
  if (missing(subSize)) {
    return(imp.mat)
  } else if (subSize > ncol(imp.mat)) {
    return(imp.mat)
  } else {
    return(imp.mat[, sample(colnames(imp.mat), subSize) ])
  }
}


#' Generates cluster-specific matrices for given data based on a clustering object.
#' 
#' @param dat.mat Data matrix to be split (features X samples).
#' @param clust Clustering object.
#' @param savePath If specified, matrices will be saved. Otherwise, a list of matrices will be returned.
#' @param savePref Preface for file names, if saving.
#' @param sizeThresh Smallest size cluster for which a matrix will be created. Default 300.
#' @return If files are NOT saved, returnes a list of matrices, one for each cluster. Otherwise, returns nothing.
ClusterMatrices <- function(dat.mat, clust, savePath, savePref, sizeThresh = 100) {
  ## set savePath if it is specified
  if(!missing(savePath)) {
    if (!missing(savePref)) {
      savePath <- paste(savePath, savePref, sep = '')
    }
  } else {
    clust.mats <- list()
  }
  ## generate matrices
  clust.table <- table(clust)
  for (i in 1:length(clust.table)) {
    if (clust.table[i] > sizeThresh) {
      clust.cells <- which(clust == names(clust.table)[i])
      clust.mat <- dat.mat[, clust.cells]
      clust.mat <- clust.mat[ rowSums(clust.mat) >= 1 ,]
      if (missing(savePath)) {
        clust.mats[[i]] <- clust.mat
      } else {
        saveRDS(clust.mat, file = paste(savePath, '_', names(clust.table)[i], '.rds', sep = ''))
      }
    }
  }
  ## return if not saving
  if (missing(savePath)) {
    return(clust.mats)
  }
}

#' Performs a CPM normalization on the given data. 
#' 
#' @param dat.mat Matrix of gene expression data (genes X samples).
#' @param l2 Optional log2 normalization switch. Default of False.
#' @return Returns CPM normalized matrix
CPMTransform <- function(dat.mat, l2 = FALSE) {
  cpm.mat <- t(t(dat.mat) / (colSums(dat.mat) / 1e6))
  if (l2) {
    cpm.mat <- log2(cpm.mat + 1)
  }
  return(cpm.mat)
}

#' Saves a matrix in a format for input to ARACNe
#'
#' @param dat.mat Matrix of data (genes X samples).
#' @param out.file Output file where matrix will be saved.
#' @param subset Switch for subsetting the matrix to 500 samples. Default TRUE.
ARACNeTable <- function(dat.mat, out.file, subset = TRUE) {
  dat.mat <- dat.mat[!duplicated(rownames(dat.mat)), ]
  saveRDS(dat.mat, file = paste(out.file,".rds",sep=""))
  if (subset) {
    dat.mat <- dat.mat[, sample(colnames(dat.mat), min(ncol(dat.mat), 500)) ]
  }
  sample.names <- colnames(dat.mat)
  gene.ids <- rownames(dat.mat)
  m <- dat.mat
  mm <- rbind( c("gene", sample.names), cbind(gene.ids, m))
  write.table( x = mm , file = paste(out.file,".tsv",sep=""), 
               sep="\t", quote = F , row.names = F , col.names = F )
}

#' Processes ARACNe results into a regulon object compatible with VIPER.
#'
#' @param a.file ARACNe final network .tsv.
#' @param exp.mat Matrix of expression from which the network was generated (genes X samples).
#' @param out.dir Output directory for networks to be saved to.
#' @param out.name Optional argument for prefix of the file name.
RegProcess <- function(a.file, exp.mat, out.dir, out.name = '.') {
  require(viper)
  processed.reg <- aracne2regulon(afile = a.file, eset = exp.mat, format = '3col')
  saveRDS(processed.reg, file = paste(out.dir, out.name, 'unPruned.rds', sep = ''))
  pruned.reg <- pruneRegulon(processed.reg, 50, adaptive = FALSE, eliminate = TRUE)
  saveRDS(pruned.reg, file = paste(out.dir, out.name, 'pruned.rds', sep = ''))
}



#batch1
#Patient A - pT3a, Grade III
#Patient B - pT1b, Grade II
#Patient C - pT1b, Grade II

###Sample IDs for RCC CyTEK. 
#TvsAdjNormal: 1==Tumor, 2==AdjNormal
#PatientID= Patient number 1-10 as in excel
#patientdate MonthDayYear

#batch2
#1	ccRCC. pT1a, grade 2, negative margins CN1-4
#2	ccRCC. pT3a, grade 2, negative margins CN5-8
#3	ccRCC. pT3a, grade 1, negative margins CN9-12
#4	ccRCC. pT1a, grade 1, negative margins CN13-16
#6	ccRCC. pT3aN0M1, grade 4, positive margins CN21-24
#7	Oncocytoma   CN26-28
#8	clear cell RCC, pT1a, grade 2, neg margins CN29-32
#9	clear cell RCC, pT2a, grade 3, neg margins CN33-36
#10	clear cell RCC, pT1a, grade 2, neg margins CN37-40

###
###LOAD AND ANNOTATE ALL DATA
###
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180426_DRAEK_NIVI_2_HUMAN_10X/CN004/outs/filtered_gene_bc_matrices/GRCh38")
pbmca1 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmca1$patient="PatientA"
pbmca1$tissue="Normal"
pbmca1$cd45="CD45+"
pbmca1$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180426_DRAEK_NIVI_2_HUMAN_10X/CN005/outs/filtered_gene_bc_matrices/GRCh38")
pbmca2 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmca2$patient="PatientA"
pbmca2$tissue="Tumor"
pbmca2$cd45="CD45+"
pbmca2$batch="batch1"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180730_CHARLES_NIVI_3_HUMAN_10X/CN009/outs/filtered_gene_bc_matrices/GRCh38")
pbmcb1 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcb1$patient="PatientB"
pbmcb1$tissue="Normal"
pbmcb1$cd45="CD45+"
pbmcb1$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180730_CHARLES_NIVI_3_HUMAN_10X/CN010/outs/filtered_gene_bc_matrices/GRCh38")
pbmcb2 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcb2$patient="PatientB"
pbmcb2$tissue="Tumor"
pbmcb2$cd45="CD45+"
pbmcb2$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180730_CHARLES_NIVI_3_HUMAN_10X/CN011/outs/filtered_gene_bc_matrices/GRCh38")
pbmcb3 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcb3$patient="PatientB"
pbmcb3$tissue="Tumor"
pbmcb3$cd45="CD45+"
pbmcb3$batch="batch1"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180807_CHARLES_NIVI_3_HUMAN_10X/CN012/outs/filtered_gene_bc_matrices/GRCh38")
pbmcc1 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcc1$patient="PatientC"
pbmcc1$tissue="Normal"
pbmcc1$cd45="CD45+"
pbmcc1$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180807_CHARLES_NIVI_3_HUMAN_10X/CN013/outs/filtered_gene_bc_matrices/GRCh38")
pbmcc2 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcc2$patient="PatientC"
pbmcc2$tissue="Tumor"
pbmcc2$cd45="CD45+"
pbmcc2$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180807_CHARLES_NIVI_3_HUMAN_10X/CN014/outs/filtered_gene_bc_matrices/GRCh38")
pbmcc3 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcc3$patient="PatientC"
pbmcc3$tissue="Tumor"
pbmcc3$cd45="CD45+"
pbmcc3$batch="batch1"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN1/outs/filtered_feature_bc_matrix")
pbmc1 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc1$patient="Patient1"
pbmc1$tissue="Tumor"
pbmc1$cd45="CD45+"
pbmc1$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN2/outs/filtered_feature_bc_matrix")
pbmc2 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc2$patient="Patient1"
pbmc2$tissue="Tumor"
pbmc2$cd45="CD45-"
pbmc2$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN3/outs/filtered_feature_bc_matrix")
pbmc3 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc3$patient="Patient1"
pbmc3$tissue="Normal"
pbmc3$cd45="CD45+"
pbmc3$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN4/outs/filtered_feature_bc_matrix")
pbmc4 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc4$patient="Patient1"
pbmc4$tissue="Normal"
pbmc4$cd45="CD45-"
pbmc4$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN5/outs/filtered_feature_bc_matrix")
pbmc5 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc5$patient="Patient2"
pbmc5$tissue="Tumor"
pbmc5$cd45="CD45+"
pbmc5$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN6/outs/filtered_feature_bc_matrix")
pbmc6 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc6$patient="Patient2"
pbmc6$tissue="Tumor"
pbmc6$cd45="CD45-"
pbmc6$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN7/outs/filtered_feature_bc_matrix")
pbmc7 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc7$patient="Patient2"
pbmc7$tissue="Normal"
pbmc7$cd45="CD45+"
pbmc7$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN8/outs/filtered_feature_bc_matrix")
pbmc8 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc8$patient="Patient2"
pbmc8$tissue="Normal"
pbmc8$cd45="CD45-"
pbmc8$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN9/outs/filtered_feature_bc_matrix")
pbmc9 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc9$patient="Patient3"
pbmc9$tissue="Tumor"
pbmc9$cd45="CD45+"
pbmc9$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN10/outs/filtered_feature_bc_matrix")
pbmc10 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc10$patient="Patient3"
pbmc10$tissue="Tumor"
pbmc10$cd45="CD45-"
pbmc10$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN11/outs/filtered_feature_bc_matrix")
pbmc11 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc11$patient="Patient3"
pbmc11$tissue="Normal"
pbmc11$cd45="CD45+"
pbmc11$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN12/outs/filtered_feature_bc_matrix")
pbmc12 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc12$patient="Patient3"
pbmc12$tissue="Normal"
pbmc12$cd45="CD45-"
pbmc12$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190603_CHARLES_NIVI_9_HUMAN_10X/CN13/outs/filtered_feature_bc_matrix")
pbmc13 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc13$patient="Patient4"
pbmc13$tissue="Tumor"
pbmc13$cd45="CD45+"
pbmc13$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190603_CHARLES_NIVI_9_HUMAN_10X/CN14/outs/filtered_feature_bc_matrix")
pbmc14 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc14$patient="Patient4"
pbmc14$tissue="Tumor"
pbmc14$cd45="CD45-"
pbmc14$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190603_CHARLES_NIVI_9_HUMAN_10X/CN15/outs/filtered_feature_bc_matrix")
pbmc15 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc15$patient="Patient4"
pbmc15$tissue="Normal"
pbmc15$cd45="CD45+"
pbmc15$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190603_CHARLES_NIVI_9_HUMAN_10X/CN16/outs/filtered_feature_bc_matrix")
pbmc16 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc16$patient="Patient4"
pbmc16$tissue="Normal"
pbmc16$cd45="CD45-"
pbmc16$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN21b/outs/filtered_feature_bc_matrix")
pbmc_data=pbmc_data[,sample(colnames(pbmc_data),10000)]
pbmc21 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc21$patient="Patient5"
pbmc21$tissue="Tumor"
pbmc21$cd45="CD45+"
pbmc21$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN22b/outs/filtered_feature_bc_matrix")
pbmc_data=pbmc_data[,sample(colnames(pbmc_data),10000)]
pbmc22 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc22$patient="Patient5"
pbmc22$tissue="Tumor"
pbmc22$cd45="CD45-"
pbmc22$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN23/outs/filtered_feature_bc_matrix")
pbmc_data=pbmc_data[,sample(colnames(pbmc_data),10000)]
pbmc23 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc23$patient="Patient5"
pbmc23$tissue="Normal"
pbmc23$cd45="CD45+"
pbmc23$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN24/outs/filtered_feature_bc_matrix")
pbmc_data=pbmc_data[,sample(colnames(pbmc_data),10000)]
pbmc24 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc24$patient="Patient5"
pbmc24$tissue="Normal"
pbmc24$cd45="CD45-"
pbmc24$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN29/outs/filtered_feature_bc_matrix")
pbmc29 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc29$patient="Patient6"
pbmc29$tissue="Tumor"
pbmc29$cd45="CD45+"
pbmc29$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN30/outs/filtered_feature_bc_matrix")
pbmc30 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc30$patient="Patient6"
pbmc30$tissue="Tumor"
pbmc30$cd45="CD45-"
pbmc30$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN31/outs/filtered_feature_bc_matrix")
pbmc31 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc31$patient="Patient6"
pbmc31$tissue="Normal"
pbmc31$cd45="CD45+"
pbmc31$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN32/outs/filtered_feature_bc_matrix")
pbmc32 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc32$patient="Patient6"
pbmc32$tissue="Normal"
pbmc32$cd45="CD45-"
pbmc32$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN33/outs/filtered_feature_bc_matrix")
pbmc33 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc33$patient="Patient7"
pbmc33$tissue="Tumor"
pbmc33$cd45="CD45+"
pbmc33$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN34/outs/filtered_feature_bc_matrix")
pbmc34 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc34$patient="Patient7"
pbmc34$tissue="Tumor"
pbmc34$cd45="CD45-"
pbmc34$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN35/outs/filtered_feature_bc_matrix")
pbmc35 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc35$patient="Patient7"
pbmc35$tissue="Normal"
pbmc35$cd45="CD45+"
pbmc35$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN36/outs/filtered_feature_bc_matrix")
pbmc36 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc36$patient="Patient7"
pbmc36$tissue="Normal"
pbmc36$cd45="CD45-"
pbmc36$batch="batch2"
rm(pbmc_data)

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190825_VICTOR_VICTOR_1_HUMAN_10X/CN37/outs/filtered_feature_bc_matrix")
pbmc37 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc37$patient="Patient8"
pbmc37$tissue="Tumor"
pbmc37$cd45="CD45+"
pbmc37$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190825_VICTOR_VICTOR_1_HUMAN_10X/CN38/outs/filtered_feature_bc_matrix")
pbmc38 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc38$patient="Patient8"
pbmc38$tissue="Tumor"
pbmc38$cd45="CD45-"
pbmc38$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190825_VICTOR_VICTOR_1_HUMAN_10X/CN39/outs/filtered_feature_bc_matrix")
pbmc39 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc39$patient="Patient8"
pbmc39$tissue="Normal"
pbmc39$cd45="CD45+"
pbmc3$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190825_VICTOR_VICTOR_1_HUMAN_10X/CN40/outs/filtered_feature_bc_matrix")
pbmc40 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc40$patient="Patient8"
pbmc40$tissue="Normal"
pbmc40$cd45="CD45-"
pbmc40$batch="batch2"
rm(pbmc_data)


##
##combine all samples and save in one file
##
pbmc.big=merge(pbmc1,y=c(pbmc2,pbmc3,pbmc4,pbmc5,pbmc6,pbmc7,pbmc8,pbmc9,pbmc10),project = "RCC_SC")
pbmc.big=merge(pbmc.big,y=c(pbmc11,pbmc12,pbmc13,pbmc14,pbmc15,pbmc16,pbmc21,pbmc22,pbmc23,pbmc24),project="RCC_SC")
pbmc.big=merge(pbmc.big,y=c(pbmc29,pbmc30,pbmc31,pbmc32,pbmc33,pbmc34,pbmc35,pbmc36,pbmc37,pbmc38,pbmc39,pbmc40),project="RCC_SC")
pbmc.big<- merge(pbmc.big, y = c(pbmca1,pbmca2,pbmcb1,pbmcb2,pbmcb3,pbmcc1,pbmcc2,pbmcc3), project = "RCC_SC")
rm(pbmca1,pbmca2,pbmcb1,pbmcb2,pbmcb3,pbmcc1,pbmcc2,pbmcc3)
rm(pbmc1,pbmc2,pbmc3,pbmc4,pbmc5,pbmc6,pbmc7,pbmc8,pbmc9,pbmc10,pbmc11,pbmc12,pbmc13,pbmc14,pbmc15,pbmc16,pbmc21,pbmc22,pbmc23,pbmc24,pbmc29,pbmc30,pbmc31,pbmc32,pbmc33,pbmc34,pbmc35,pbmc36,pbmc37,pbmc38,pbmc39,pbmc40)
big_list=SplitObject(pbmc.big, split.by = "cd45")
rm(pbmc.big)
pbmc.big=big_list[[1]]
pbmc.big <- PercentageFeatureSet(pbmc.big, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(pbmc.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0,group.by="tissue")
print(ncol(pbmc.big))
pbmc.big <- subset(pbmc.big, subset = percent.mt < 10 & nCount_RNA > 1500 & nCount_RNA < 15000)
VlnPlot(pbmc.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0,group.by="tissue")
print(ncol(pbmc.big))
pbmc.big <- SCTransform(pbmc.big,return.only.var.genes = F,verbose = T,conserve.memory = T,ncells=10000)
pbmc.big <- RunPCA(pbmc.big, features = VariableFeatures(object = pbmc.big),npcs = 30)
ElbowPlot(pbmc.big)
pbmc.big <- RunUMAP(pbmc.big, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
pbmc.big.singler=CreateBigSingleRObjectv2(counts = pbmc.big[["SCT"]]@counts,annot=NULL,project.name="RCC1",N=10000,
                                          min.genes=0,technology='10X',
                                          species='Human',citation='',normalize.gene.length=F,
                                          variable.genes='de',fine.tune=F,
                                          reduce.file.size=T,do.signatures=F,
                                          do.main.types=T,
                                          temp.dir="RCC_singler/singler_rcc1",xy = pbmc.big$umap@cell.embeddings,clusters = NULL)
pbmc.big$hpca_labels=pbmc.big.singler$singler[[1]][[1]][[2]]
pbmc.big$hpca_main_labels=pbmc.big.singler$singler[[1]][[4]][[2]]
pbmc.big$blueprint_labels=pbmc.big.singler$singler[[2]][[1]][[2]]
pbmc.big$blueprint_main_labels=pbmc.big.singler$singler[[2]][[4]][[2]]
pbmc.big$hpca_pvals=pbmc.big.singler$singler[[1]][[1]][[3]]
pbmc.big$hpca_main_pvals=pbmc.big.singler$singler[[1]][[4]][[3]]
pbmc.big$blueprint_pvals=pbmc.big.singler$singler[[2]][[1]][[3]]
pbmc.big$blueprint_main_pvals=pbmc.big.singler$singler[[2]][[4]][[3]]
rm(pbmc.big.singler)
big_list[[1]]=pbmc.big
rm(pbmc.big)

pbmc.big=big_list[[2]]
pbmc.big <- PercentageFeatureSet(pbmc.big, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(pbmc.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0,group.by="tissue")
print(ncol(pbmc.big))
pbmc.big <- subset(pbmc.big, subset = percent.mt < 10 & nCount_RNA > 1500 & nCount_RNA < 15000)
VlnPlot(pbmc.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0,group.by="tissue")
print(ncol(pbmc.big))
pbmc.big <- SCTransform(pbmc.big,return.only.var.genes = F,verbose = T,conserve.memory = T,ncells=10000)
pbmc.big <- RunPCA(pbmc.big, features = VariableFeatures(object = pbmc.big),npcs = 30)
ElbowPlot(pbmc.big)
pbmc.big <- RunUMAP(pbmc.big, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
pbmc.big.singler=CreateBigSingleRObjectv2(counts = pbmc.big[["SCT"]]@counts,annot=NULL,project.name="RCC1",N=10000,
                                          min.genes=0,technology='10X',
                                          species='Human',citation='',normalize.gene.length=F,
                                          variable.genes='de',fine.tune=F,
                                          reduce.file.size=T,do.signatures=F,
                                          do.main.types=T,
                                          temp.dir="RCC_singler/singler_rcc1",xy = pbmc.big$umap@cell.embeddings,clusters = NULL)
pbmc.big$hpca_labels=pbmc.big.singler$singler[[1]][[1]][[2]]
pbmc.big$hpca_main_labels=pbmc.big.singler$singler[[1]][[4]][[2]]
pbmc.big$blueprint_labels=pbmc.big.singler$singler[[2]][[1]][[2]]
pbmc.big$blueprint_main_labels=pbmc.big.singler$singler[[2]][[4]][[2]]
pbmc.big$hpca_pvals=pbmc.big.singler$singler[[1]][[1]][[3]]
pbmc.big$hpca_main_pvals=pbmc.big.singler$singler[[1]][[4]][[3]]
pbmc.big$blueprint_pvals=pbmc.big.singler$singler[[2]][[1]][[3]]
pbmc.big$blueprint_main_pvals=pbmc.big.singler$singler[[2]][[4]][[3]]
rm(pbmc.big.singler)
big_list[[2]]=pbmc.big
rm(pbmc.big)
saveRDS(big_list, file = "seurat_human_rcc_allpatients.rds")
rm(big_list)


#####
#####Do patient-level gene expression analysis for all patients
#####
big_list=readRDS("seurat_human_rcc_allpatients.rds")
seurat_list_cd45pos=SplitObject(big_list[[1]],split.by="patient")
seurat_list_cd45neg=SplitObject(big_list[[2]],split.by="patient")
rm(big_list)

for(iter in 1:length(seurat_list_cd45pos)){
  print(iter)
  s=seurat_list_cd45pos[[iter]]
  patientnumber=unique(s$patient)
  
  s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")
  tiff(paste(patientnumber,"QCPlot_CD45pos.tiff",sep = "_"), width = 8, height = 6, units = 'in', res = 600)
  plot(VlnPlot(s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0,group.by="tissue"))
  dev.off()
  s <- subset(s, subset = percent.mt < 10 & nCount_RNA > 1000 & nCount_RNA < 15000)
  s <- SCTransform(s,return.only.var.genes = F,verbose = T,conserve.memory = T)
  s <- RunPCA(s, features = VariableFeatures(object = s))
  s <- RunUMAP(s, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
  s <- FindNeighbors(s, dims = 1:30, verbose = FALSE)
  s <- FindClusters(s, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
  clust=s@meta.data[,which(grepl("SCT_snn_res.",colnames(s@meta.data)))]
  mat=as.data.frame(t(s$pca@cell.embeddings))
  out=sil_subsample_v3(mat,clust)
  means=out[[1]]
  sd=out[[2]]
  x=seq(0.01,1,by=0.01)
  tiff(paste(patientnumber,"louvain_resolution_CD45pos.tiff",sep = "_"), width = 6, height = 6, units = 'in', res = 600)
  errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
  lines(x,means)
  best=tail(x[which(means==max(means))],n=1)
  legend("topright",paste("Best",best,sep = " = "))
  dev.off()
  s$seurat_clusters=s@meta.data[,which(colnames(s@meta.data)==paste("SCT_snn_res.",best,sep=""))]
  Idents(s) <- "seurat_clusters"
  tiff(paste(patientnumber,"louvain_split_umap_CD45pos.tiff",sep = "_"), width = 8, height = 6, units = 'in', res = 600)
  plot(DimPlot(s,reduction="umap",group.by="seurat_clusters",split.by="tissue"))
  dev.off()
  tiff(paste(patientnumber,"louvain_umap_CD45pos.tiff",sep = "_"), width = 6, height = 6, units = 'in', res = 600)
  plot(DimPlot(s, reduction = "umap",label = TRUE) + NoLegend())
  dev.off()
  markers <- FindAllMarkers(s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
  tiff(paste(patientnumber,"louvain_heatmap_CD45pos.tiff",sep = "_"), width = 8, height = 8, units = 'in', res = 600)
  geneHeatmap(s,s$seurat_clusters,top10$gene)
  dev.off()
  l=s$blueprint_labels
  l[which(s$blueprint_pvals>0.1)]=NA
  l[which(l %in% names(which(table(l)<20)))]=NA
  s$l=l
  Idents(s) <- "l"
  tiff(paste(patientnumber,"singler_umap_CD45pos.tiff",sep = "_"), width = 8, height = 6, units = 'in', res = 600)
  plot(DimPlot(s, reduction = "umap",label=TRUE))
  dev.off()
  saveRDS(s, file = paste(patientnumber,"CD45pos.rds",sep="_"))
  
  s_meta=MakeCMfA(dat.mat=as.matrix(s[["SCT"]]@counts),clustering=s$seurat_clusters,out.dir="metacell_aracne_inputs/",out.name=paste(patientnumber,"CD45pos",sep="_"))
  meta=s_meta[[1]]
  for(i in 2:length(s_meta)){
    meta=merge(meta,s_meta[[i]],by=0,all=T)
    rownames(meta)=meta[,1]
    meta=meta[,2:ncol(meta)]
  }
  meta[is.na(meta)] <- 0
  saveRDS(meta, file = paste(patientnumber,"CD45pos_metacells.rds",sep="_"))
}

for(iter in 1:length(seurat_list_cd45neg)){
  print(iter)
  s2=seurat_list_cd45neg[[iter]]
  patientnumber=unique(s2$patient)
  
  s2 <- PercentageFeatureSet(s2, pattern = "^MT-", col.name = "percent.mt")
  tiff(paste(patientnumber,"QCPlot_CD45neg.tiff",sep = "_"), width = 8, height = 6, units = 'in', res = 600)
  plot(VlnPlot(s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0,group.by="tissue"))
  dev.off()
  s2 <- subset(s2, subset = percent.mt < 10 & nCount_RNA > 1000 & nCount_RNA < 15000)
  s2 <- SCTransform(s2,return.only.var.genes = F,verbose = T,conserve.memory = T)
  s2 <- RunPCA(s2, features = VariableFeatures(object = s2))
  s2 <- RunUMAP(s2, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
  s2 <- FindNeighbors(s2, dims = 1:30, verbose = FALSE)
  s2 <- FindClusters(s2, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
  clust=s2@meta.data[,which(grepl("SCT_snn_res.",colnames(s2@meta.data)))]
  mat=as.data.frame(t(s2$pca@cell.embeddings))
  out=sil_subsample_v3(mat,clust)
  means=out[[1]]
  sd=out[[2]]
  x=seq(0.01,1,by=0.01)
  tiff(paste(patientnumber,"louvain_resolution_CD45neg.tiff",sep = "_"), width = 6, height = 6, units = 'in', res = 600)
  errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
  lines(x,means)
  best=tail(x[which(means==max(means))],n=1)
  legend("topright",paste("Best",best,sep = " = "))
  dev.off()
  s2$seurat_clusters=s2@meta.data[,which(colnames(s2@meta.data)==paste("SCT_snn_res.",best,sep=""))]
  Idents(s2) <- "seurat_clusters"
  tiff(paste(patientnumber,"louvain_split_umap_CD45neg.tiff",sep = "_"), width = 8, height = 6, units = 'in', res = 600)
  plot(DimPlot(s2,reduction="umap",group.by="seurat_clusters",split.by="tissue"))
  dev.off()
  tiff(paste(patientnumber,"louvain_umap_CD45neg.tiff",sep = "_"), width = 6, height = 6, units = 'in', res = 600)
  plot(DimPlot(s2, reduction = "umap",label = TRUE) + NoLegend())
  dev.off()
  markers <- FindAllMarkers(s2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
  tiff(paste(patientnumber,"louvain_heatmap_CD45neg.tiff",sep = "_"), width = 8, height = 8, units = 'in', res = 600)
  geneHeatmap(s2,s2$seurat_clusters,top10$gene)
  dev.off()
  l=s2$blueprint_labels
  l[which(s2$blueprint_pvals>0.1)]=NA
  l[which(l %in% names(which(table(l)<20)))]=NA
  s2$l=l
  Idents(s2) <- "l"
  tiff(paste(patientnumber,"singler_umap_CD45neg.tiff",sep = "_"), width = 8, height = 6, units = 'in', res = 600)
  plot(DimPlot(s2, reduction = "umap",label=TRUE))
  dev.off()
  saveRDS(s2, file = paste(patientnumber,"CD45neg.rds",sep="_"))
  
  s2_meta=MakeCMfA(dat.mat=as.matrix(s2[["SCT"]]@counts),clustering=s2$seurat_clusters,out.dir="metacell_aracne_inputs/",out.name=paste(patientnumber,"CD45neg",sep="_"))
  meta=s2_meta[[1]]
  for(i in 2:length(s2_meta)){
    meta=merge(meta,s2_meta[[i]],by=0,all=T)
    rownames(meta)=meta[,1]
    meta=meta[,2:ncol(meta)]
  }
  meta[is.na(meta)] <- 0
  saveRDS(meta, file = paste(patientnumber,"CD45neg_metacells.rds",sep="_"))
}
rm(seurat_list_cd45pos,seurat_list_cd45neg)



###
#####MERGE + BATCH-CORRECT SINGLE-CELLS USING SEURAT ALGORITHM
###
list_cd45pos=list()
list_cd45neg=list()
patients=c("Patient1","Patient2","Patient3","Patient4","Patient5","Patient6","Patient7","Patient8","PatientA","PatientB","PatientC")
for(pt in patients){
  s_allcells=readRDS(paste(pt,"_CD45pos.rds",sep=""))
  list_cd45pos=c(list_cd45pos,s_allcells)
  rm(s_allcells)
  if(!(pt %in% c("PatientA","PatientB","PatientC"))){
    s2_allcells=readRDS(paste(pt,"_CD45neg.rds",sep=""))
    list_cd45neg=c(list_cd45neg,s2_allcells)
    rm(s2_allcells)}
}

names(list_cd45pos)=patients
cd45pos.features <- SelectIntegrationFeatures(object.list = list_cd45pos, nfeatures = 2000)
cd45pos.list <- PrepSCTIntegration(object.list = list_cd45pos, anchor.features = cd45pos.features, verbose = T)
rm(list_cd45pos)
cd45pos.anchors <- FindIntegrationAnchors(object.list = cd45pos.list, normalization.method = "SCT", anchor.features = cd45pos.features, verbose = T,reference = 3)
rm(cd45pos.features,cd45pos.list)
cd45pos.integrated <- IntegrateData(anchorset = cd45pos.anchors, normalization.method = "SCT", verbose = T)
rm(cd45pos.anchors)

names(list_cd45neg)=patients[1:8]
cd45neg.features <- SelectIntegrationFeatures(object.list = list_cd45neg, nfeatures = 2000)
cd45neg.list <- PrepSCTIntegration(object.list = list_cd45neg, anchor.features = cd45neg.features, verbose = T)
rm(list_cd45neg)
cd45neg.anchors <- FindIntegrationAnchors(object.list = cd45neg.list, normalization.method = "SCT", anchor.features = cd45neg.features, verbose = T,reference = 3)
rm(cd45neg.features,cd45neg.list)
cd45neg.integrated <- IntegrateData(anchorset = cd45neg.anchors, normalization.method = "SCT", verbose = T)
rm(cd45neg.anchors)

cd45pos.integrated <- RunPCA(cd45pos.integrated, features = VariableFeatures(object = cd45pos.integrated))
cd45pos.integrated <- RunUMAP(cd45pos.integrated, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
cd45pos.integrated <- FindNeighbors(cd45pos.integrated, dims = 1:50, verbose = FALSE)
cd45pos.integrated <- FindClusters(cd45pos.integrated, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=cd45pos.integrated@meta.data[,which(grepl("integrated_snn_res.",colnames(cd45pos.integrated@meta.data)))]
mat=as.data.frame(t(cd45pos.integrated$pca@cell.embeddings))
out=sil_subsample_v2(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
tiff("cd45pos_merged_louvain_resolution_seurat.tiff", width = 6, height = 6, units = 'in', res = 600)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
dev.off()
cd45pos.integrated$seurat_clusters=cd45pos.integrated@meta.data[,which(colnames(cd45pos.integrated@meta.data)==paste("integrated_snn_res.",best,sep=""))]
Idents(cd45pos.integrated) <- "seurat_clusters"
#
# iterative louvain sub-clustering
cd45pos.integrated$original_seurat_clusters=cd45pos.integrated$seurat_clusters
cluster_output=as.character(cd45pos.integrated$seurat_clusters)
for(i in unique(as.character(cd45pos.integrated$seurat_clusters))){
  print(i)
  s <- subset(cd45pos.integrated, subset = seurat_clusters == i)
  if(ncol(s)/ncol(cd45pos.integrated)>0.05){
    s <- RunPCA(s, features = VariableFeatures(object = s))
    s <- RunUMAP(s, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
    s <- FindNeighbors(s, dims = 1:50, verbose = FALSE)
    s <- FindClusters(s, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
    clust=s@meta.data[,which(grepl("integrated_snn_res.",colnames(s@meta.data)))]
    mat=as.data.frame(t(s$pca@cell.embeddings))
    out=sil_subsample_v2(mat,clust)
    means=out[[1]]
    sd=out[[2]]
    x=seq(0.01,1,by=0.01)
    best=tail(x[which(means==max(means))],n=1)
    s$seurat_clusters=s@meta.data[,which(colnames(s@meta.data)==paste("integrated_snn_res.",best,sep=""))]
    if(max(means)>0.25 & min(table(s$seurat_clusters))/sum(table(s$seurat_clusters))>0.05){cluster_output[which(cluster_output==i)]=paste(i,s$seurat_clusters)}
  }
}
cd45pos.integrated$seurat_clusters=as.factor(as.numeric(as.factor(cluster_output)))
Idents(cd45pos.integrated) <- "seurat_clusters"
#end iterative louvain sub-clustering
#
cd45pos.integrated$seurat_clusters=mapvalues(cd45pos.integrated$seurat_clusters, from = 1:length(unique(cd45pos.integrated$seurat_clusters)), to = c("CD4 T-cell","Treg","CD8 T-cell", "NK cell 1", "NK cell 2","Macrophage","Monocyte","B cell","Misc.","Mast cell","Plasma cell"))
Idents(cd45pos.integrated) <- "seurat_clusters"
#
tiff("cd45pos_merged_louvain_split_umap_seurat.tiff", width = 10, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos.integrated,reduction="umap",group.by="seurat_clusters",split.by="tissue",cols = c("darkgoldenrod2","darkorange3","mediumorchid2","forestgreen","chartreuse3","cyan4","cyan2","lightsteelblue3","cornsilk3","deeppink3","lightsteelblue4")) + NoLegend())
dev.off()
tiff("cd45pos_merged_louvain_split_patient_umap_seurat.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos.integrated,reduction="umap",group.by="seurat_clusters",split.by="patient",ncol=3,cols = c("darkgoldenrod2","darkorange3","mediumorchid2","forestgreen","chartreuse3","cyan4","cyan2","lightsteelblue3","cornsilk3","deeppink3","lightsteelblue4")))
dev.off()
tiff("cd45pos_merged_louvain_umap_seurat.tiff", width = 6, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos.integrated, reduction = "umap",label = TRUE,label.size=7,repel=T,cols = c("darkgoldenrod2","darkorange3","mediumorchid2","forestgreen","chartreuse3","cyan4","cyan2","lightsteelblue3","cornsilk3","deeppink3","lightsteelblue4")) + NoLegend())
dev.off()
markers <- FindAllMarkers(cd45pos.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tiff("cd45pos_merged_louvain_heatmap_seurat.tiff", width = 8, height = 8, units = 'in', res = 600)
#geneHeatmap(cd45pos.integrated,cd45pos.integrated$seurat_clusters,top10$gene)
geneHeatmap(cd45pos.integrated,cd45pos.integrated$seurat_clusters,top10$gene,color_palette=c("darkgoldenrod2","darkorange3","mediumorchid2","forestgreen","chartreuse3","cyan4","cyan2","lightsteelblue3","cornsilk3","deeppink3","lightsteelblue4"))
dev.off()
l=cd45pos.integrated$blueprint_labels
l[which(cd45pos.integrated$blueprint_pvals>0.1)]=NA
l[which(l %in% names(which(table(l)<150)))]=NA
cd45pos.integrated$l=l
Idents(cd45pos.integrated) <- "l"
tiff("cd45pos_merged_singler_umap_seurat.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos.integrated, reduction = "umap",label=TRUE,repel=T,label.size=5)+NoLegend())
dev.off()
saveRDS(cd45pos.integrated, file = "cd45pos_merged_seurat.rds")
rm(list_cd45pos,cd45pos.integrated)

cd45neg.integrated <- RunPCA(cd45neg.integrated, features = VariableFeatures(object = cd45neg.integrated))
cd45neg.integrated <- RunUMAP(cd45neg.integrated, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
cd45neg.integrated <- FindNeighbors(cd45neg.integrated, dims = 1:50, verbose = FALSE)
cd45neg.integrated <- FindClusters(cd45neg.integrated, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=cd45neg.integrated@meta.data[,which(grepl("integrated_snn_res.",colnames(cd45neg.integrated@meta.data)))]
mat=as.data.frame(t(cd45neg.integrated$pca@cell.embeddings))
out=sil_subsample_v2(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
tiff("cd45neg_merged_louvain_resolution_seurat.tiff", width = 6, height = 6, units = 'in', res = 600)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
dev.off()
cd45neg.integrated$seurat_clusters=cd45neg.integrated@meta.data[,which(colnames(cd45neg.integrated@meta.data)==paste("integrated_snn_res.",best,sep=""))]
Idents(cd45neg.integrated) <- "seurat_clusters"
#
cd45neg.integrated$seurat_clusters=mapvalues(cd45neg.integrated$seurat_clusters, from = 1:length(unique(cd45neg.integrated$seurat_clusters))-1, to = c("Epithelial","Endothelial","Fibroblast", "M2 Macrophage", "Adipocyte"))
Idents(cd45neg.integrated) <- "seurat_clusters"
#
tiff("cd45neg_merged_louvain_split_umap_seurat.tiff", width = 10, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg.integrated,reduction="umap",group.by="seurat_clusters",split.by="tissue",cols = c("chocolate2","chartreuse3","gold2","lemonchiffon3","pink1")) + NoLegend())
dev.off()
tiff("cd45neg_merged_louvain_split_patient_umap_seurat.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg.integrated,reduction="umap",group.by="seurat_clusters",split.by="patient",ncol=3,cols = c("chocolate2","chartreuse3","gold2","lemonchiffon3","pink1")))
dev.off()
tiff("cd45neg_merged_louvain_umap_seurat.tiff", width = 6, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg.integrated, reduction = "umap",label = TRUE,label.size=7,repel=T,cols = c("chocolate2","chartreuse3","gold2","lemonchiffon3","pink1")) + NoLegend())
dev.off()
markers <- FindAllMarkers(cd45neg.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tiff("cd45neg_merged_louvain_heatmap_seurat.tiff", width = 8, height = 8, units = 'in', res = 600)
geneHeatmap(cd45neg.integrated,cd45neg.integrated$seurat_clusters,top10$gene,color_palette = c("chocolate2","chartreuse3","gold2","lemonchiffon3","pink1"))
dev.off()
l=cd45neg.integrated$blueprint_labels
l[which(cd45neg.integrated$blueprint_pvals>0.1)]=NA
l[which(l %in% names(which(table(l)<1400)))]=NA
cd45neg.integrated$l=l
Idents(cd45neg.integrated) <- "l"
tiff("cd45neg_merged_singler_umap_seurat.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg.integrated, reduction = "umap",label=TRUE,repel=T,label.size=5)+NoLegend())
dev.off()
saveRDS(cd45neg.integrated, file = "cd45neg_merged_seurat.rds")
rm(list_cd45neg,cd45neg.integrated)


###
#####MERGE + BATCH-CORRECT META-CELLS USING SEURAT ALGORITHM
###
list_cd45pos=list()
list_cd45neg=list()
patients=c("Patient1","Patient2","Patient3","Patient4","Patient5","Patient6","Patient7","Patient8","PatientA","PatientB","PatientC")
for(pt in patients){
  print(pt)
  s_allcells=readRDS(paste(pt,"_CD45pos.rds",sep=""))
  s_meta=readRDS(paste(pt,"_CD45pos_metacells.rds",sep=""))
  i=sample(colnames(s_meta),ncol(s_meta)/3)
  s_meta=s_meta[,i]
  temp=CreateSeuratObject(counts = s_meta)
  temp=SCTransform(temp,return.only.var.genes = F,verbose = T,conserve.memory = T)
  temp$l=s_allcells$l[colnames(temp)]
  temp$blueprint_labels=s_allcells$blueprint_labels[colnames(temp)]
  temp$blueprint_pvals=s_allcells$blueprint_pvals[colnames(temp)]
  temp$tissue=s_allcells$tissue[colnames(temp)]
  temp$cd45=s_allcells$cd45[colnames(temp)]
  temp$patient=s_allcells$patient[colnames(temp)]
  list_cd45pos=c(list_cd45pos,temp)
  rm(s_allcells,s_meta,temp)
  if(!(pt %in% c("PatientA","PatientB","PatientC"))){
    s2_allcells=readRDS(paste(pt,"_CD45neg.rds",sep=""))
    s2_meta=readRDS(paste(pt,"_CD45neg_metacells.rds",sep=""))
    i=sample(colnames(s2_meta),ncol(s2_meta)/3)
    s2_meta=s2_meta[,i]
    temp=CreateSeuratObject(counts = s2_meta)
    temp=SCTransform(temp,return.only.var.genes = F,verbose = T,conserve.memory = T)
    temp$l=s2_allcells$l[colnames(temp)]
    temp$blueprint_labels=s2_allcells$blueprint_labels[colnames(temp)]
    temp$blueprint_pvals=s2_allcells$blueprint_pvals[colnames(temp)]
    temp$tissue=s2_allcells$tissue[colnames(temp)]
    temp$cd45=s2_allcells$cd45[colnames(temp)]
    temp$patient=s2_allcells$patient[colnames(temp)]
    list_cd45neg=c(list_cd45neg,temp)
    rm(s2_allcells,s2_meta,temp)
  }
}
saveRDS(list_cd45pos,"list_cd45pos_meta.rds")
saveRDS(list_cd45neg,"list_cd45neg_meta.rds")
rm(list_cd45pos,list_cd45neg)

list_cd45pos=readRDS("list_cd45pos_meta.rds")
names(list_cd45pos)=patients
cd45pos.features <- SelectIntegrationFeatures(object.list = list_cd45pos, nfeatures = 2500)
cd45pos.list <- PrepSCTIntegration(object.list = list_cd45pos, anchor.features = cd45pos.features, verbose = T)
rm(list_cd45pos)
cd45pos.anchors <- FindIntegrationAnchors(object.list = cd45pos.list, normalization.method = "SCT", anchor.features = cd45pos.features, verbose = T,reference = 3)
rm(cd45pos.features,cd45pos.list)
cd45pos.integrated <- IntegrateData(anchorset = cd45pos.anchors, normalization.method = "SCT", verbose = T)
rm(cd45pos.anchors)

cd45pos.integrated <- RunPCA(cd45pos.integrated, features = VariableFeatures(object = cd45pos.integrated))
cd45pos.integrated <- RunUMAP(cd45pos.integrated, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
cd45pos.integrated <- FindNeighbors(cd45pos.integrated, dims = 1:30, verbose = FALSE)
cd45pos.integrated <- FindClusters(cd45pos.integrated, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=cd45pos.integrated@meta.data[,which(grepl("integrated_snn_res.",colnames(cd45pos.integrated@meta.data)))]
mat=as.data.frame(t(cd45pos.integrated$pca@cell.embeddings))
out=sil_subsample_v2(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
tiff("cd45pos_merged_louvain_resolution_seurat_meta.tiff", width = 6, height = 6, units = 'in', res = 600)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
dev.off()
cd45pos.integrated$seurat_clusters=cd45pos.integrated@meta.data[,which(colnames(cd45pos.integrated@meta.data)==paste("integrated_snn_res.",best,sep=""))]
Idents(cd45pos.integrated) <- "seurat_clusters"
tiff("cd45pos_merged_louvain_split_umap_seurat_meta.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos.integrated,reduction="umap",group.by="seurat_clusters",split.by="tissue"))
dev.off()
tiff("cd45pos_merged_louvain_split_patient_umap_seurat_meta.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos.integrated,reduction="umap",group.by="seurat_clusters",split.by="patient",ncol=3))
dev.off()
tiff("cd45pos_merged_louvain_umap_seurat_meta.tiff", width = 6, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos.integrated, reduction = "umap",label = TRUE) + NoLegend())
dev.off()
markers <- FindAllMarkers(cd45pos.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tiff("cd45pos_merged_louvain_heatmap_seurat_meta.tiff", width = 8, height = 8, units = 'in', res = 600)
geneHeatmap(cd45pos.integrated,cd45pos.integrated$seurat_clusters,top10$gene)
dev.off()
l=cd45pos.integrated$blueprint_labels
l[which(cd45pos.integrated$blueprint_pvals>0.1)]=NA
l[which(l %in% names(which(table(l)<150)))]=NA
cd45pos.integrated$l=l
Idents(cd45pos.integrated) <- "l"
tiff("cd45pos_merged_singler_umap_seurat_meta.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos.integrated, reduction = "umap",label=TRUE))
dev.off()
saveRDS(cd45pos.integrated, file = "cd45pos_merged_seurat_meta.rds")
rm(list_cd45pos,cd45pos.integrated)

list_cd45neg=readRDS("list_cd45neg_meta.rds")
names(list_cd45neg)=patients[1:8]
cd45neg.features <- SelectIntegrationFeatures(object.list = list_cd45neg, nfeatures = 2500)
cd45neg.list <- PrepSCTIntegration(object.list = list_cd45neg, anchor.features = cd45neg.features, verbose = T)
rm(list_cd45neg)
cd45neg.anchors <- FindIntegrationAnchors(object.list = cd45neg.list, normalization.method = "SCT", anchor.features = cd45neg.features, verbose = T,reference = 3)
rm(cd45neg.features,cd45neg.list)
cd45neg.integrated <- IntegrateData(anchorset = cd45neg.anchors, normalization.method = "SCT", verbose = T)
rm(cd45neg.anchors)

cd45neg.integrated <- RunPCA(cd45neg.integrated, features = VariableFeatures(object = cd45neg.integrated))
cd45neg.integrated <- RunUMAP(cd45neg.integrated, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
cd45neg.integrated <- FindNeighbors(cd45neg.integrated, dims = 1:30, verbose = FALSE)
cd45neg.integrated <- FindClusters(cd45neg.integrated, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=cd45neg.integrated@meta.data[,which(grepl("integrated_snn_res.",colnames(cd45neg.integrated@meta.data)))]
mat=as.data.frame(t(cd45neg.integrated$pca@cell.embeddings))
out=sil_subsample_v2(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
tiff("cd45neg_merged_louvain_resolution_seurat_meta.tiff", width = 6, height = 6, units = 'in', res = 600)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
dev.off()
cd45neg.integrated$seurat_clusters=cd45neg.integrated@meta.data[,which(colnames(cd45neg.integrated@meta.data)==paste("integrated_snn_res.",best,sep=""))]
Idents(cd45neg.integrated) <- "seurat_clusters"
tiff("cd45neg_merged_louvain_split_umap_seurat_meta.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg.integrated,reduction="umap",group.by="seurat_clusters",split.by="tissue"))
dev.off()
tiff("cd45neg_merged_louvain_split_patient_umap_seurat_meta.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg.integrated,reduction="umap",group.by="seurat_clusters",split.by="patient",ncol=3))
dev.off()
tiff("cd45neg_merged_louvain_umap_seurat_meta.tiff", width = 6, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg.integrated, reduction = "umap",label = TRUE) + NoLegend())
dev.off()
markers <- FindAllMarkers(cd45neg.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tiff("cd45neg_merged_louvain_heatmap_seurat_meta.tiff", width = 8, height = 8, units = 'in', res = 600)
geneHeatmap(cd45neg.integrated,cd45neg.integrated$seurat_clusters,top10$gene)
dev.off()
l=cd45neg.integrated$blueprint_labels
l[which(cd45neg.integrated$blueprint_pvals>0.1)]=NA
l[which(l %in% names(which(table(l)<1400)))]=NA
cd45neg.integrated$l=l
Idents(cd45neg.integrated) <- "l"
tiff("cd45neg_merged_singler_umap_seurat_meta.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg.integrated, reduction = "umap",label=TRUE))
dev.off()
saveRDS(cd45neg.integrated, file = "cd45neg_merged_seurat_meta.rds")
rm(list_cd45neg,cd45neg.integrated)

###load viper gene lists
surfacemarkers=read.table("single-cell-pipeline-master/ARACNe/human_hugo/surface-hugo.txt",sep="\t")
surfacemarkers=surfacemarkers[,1]
tfs=read.table("single-cell-pipeline-master/ARACNe/human_hugo/tfs-hugo.txt",sep="\t")
tfs=tfs[,1]
cotfs=read.table("single-cell-pipeline-master/ARACNe/human_hugo/cotfs-hugo.txt",sep="\t")
cotfs=cotfs[,1]
sig=read.table("single-cell-pipeline-master/ARACNe/human_hugo/sig-hugo.txt",sep="\t")
sig=sig[,1]


#####
#####Do patient-level VIPER analysis for all patients
#####
#run VIPER on batch-corrected meta-cells
cd45neg.integrated=readRDS("cd45neg_merged_seurat_meta.rds")
Idents(cd45neg.integrated) <- "seurat_clusters"
cd45neg_patients=unique(cd45neg.integrated$patient)
for(pt in cd45neg_patients){
  print(paste(pt,"CD45neg VIPER",sep=" "))
  filenames <- list.files("sc_nets/", pattern="Patient*", full.names=TRUE)
  filenames=filenames[grepl("neg",filenames)]
  nets=lapply(filenames,readRDS)
  dat=cd45neg.integrated@assays$integrated@scale.data[,which(cd45neg.integrated$patient==pt)]
  ncols=length(which(cd45neg.integrated$patient==pt))
  if(ncols>8000){
    s_vp_meta_1<-viper(dat[,1:4000], nets, method = 'none')
    dat=dat[,4001:ncol(dat)]
    s_vp_meta_2<-viper(dat[,1:4000], nets, method = 'none')
    dat=dat[,4001:ncol(dat)]
    s_vp_meta_3<-viper(dat, nets, method = 'none')
    s_vp_meta=cbind(s_vp_meta_1,s_vp_meta_2,s_vp_meta_3)
    rm(s_vp_meta_1,s_vp_meta_2,s_vp_meta_3)
  }
  if(ncols<=8000){
    s_vp_meta <- viper(dat, nets, method = 'none')
  }
  rm(dat)
  saveRDS(s_vp_meta, file = paste(pt,"_cd45neg_vp.rds",sep=""))
  rm(s_vp_meta)
}
rm(cd45neg.integrated,cd45neg_meta)

cd45pos.integrated=readRDS("cd45pos_merged_seurat_meta.rds")
Idents(cd45pos.integrated) <- "seurat_clusters"
cd45pos_patients=unique(cd45pos.integrated$patient)
for(pt in cd45pos_patients){
  print(paste(pt,"CD45pos VIPER",sep=" "))
  filenames <- list.files("sc_nets/", pattern="Patient*", full.names=TRUE)
  filenames=filenames[grepl("pos",filenames)]
  nets=lapply(filenames,readRDS)
  dat=cd45pos.integrated@assays$integrated@scale.data[,which(cd45pos.integrated$patient==pt)]
  ncols=length(which(cd45pos.integrated$patient==pt))
  if(ncols>8000){
    s1_vp_meta_1<-viper(dat[,1:4000], nets, method = 'none')
    dat=dat[,4001:ncol(dat)]
    s1_vp_meta_2<-viper(dat[,1:4000], nets, method = 'none')
    dat=dat[,4001:ncol(dat)]
    s1_vp_meta_3<-viper(dat, nets, method = 'none')
    s1_vp_meta=cbind(s1_vp_meta_1,s1_vp_meta_2,s1_vp_meta_3)
    rm(s1_vp_meta_1,s1_vp_meta_2,s1_vp_meta_3)
  }
  if(ncols<=8000){
    s1_vp_meta <- viper(dat, nets, method = 'none')
  }
  rm(dat)
  saveRDS(s1_vp_meta, file = paste(pt,"_cd45pos_vp.rds",sep=""))
  rm(s1_vp_meta)
}
rm(cd45pos.integrated,cd45pos_meta)


for(pt in cd45pos_patients){
  print(paste(pt,"CD45pos VIPER",sep=" "))
  s1_vp_meta=readRDS(paste(pt,"_cd45pos_vp.rds",sep=""))
  s1=readRDS(paste(pt,"_CD45pos.rds",sep=""))
  #cbcMRs <- CBCMRs(s1_vp_meta)
  #s1_vp_meta <- s1_vp_meta[cbcMRs,]
  s1_vp_meta_seurat <- CreateSeuratObject(counts = s1_vp_meta)
  s1_vp_meta_seurat@assays$RNA@scale.data=as.matrix(s1_vp_meta_seurat@assays$RNA@data)
  s1_vp_meta_seurat <- RunPCA(s1_vp_meta_seurat,features=rownames(s1_vp_meta_seurat))
  s1_vp_meta_seurat <- RunUMAP(s1_vp_meta_seurat, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
  s1_vp_meta_seurat <- FindNeighbors(s1_vp_meta_seurat, dims = 1:30, verbose = FALSE)
  s1_vp_meta_seurat <- FindClusters(s1_vp_meta_seurat, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
  clust=s1_vp_meta_seurat@meta.data[,which(grepl("RNA_snn_res.",colnames(s1_vp_meta_seurat@meta.data)))]
  mat=s1_vp_meta_seurat@assays$RNA@scale.data
  out=sil_subsample_v2_viper(mat,clust)
  means=out[[1]]
  sd=out[[2]]
  x=seq(0.01,1,by=0.01)
  tiff(paste(pt,"cd45pos_viper_louvain_resolution_anchor.tiff",sep="_"), width = 8, height = 6, units = 'in', res = 600)
  errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
  lines(x,means)
  best=tail(x[which(means==max(means))],n=1)
  legend("topright",paste("Best",best,sep = " = "))
  dev.off()
  s1_vp_meta_seurat$seurat_clusters=s1_vp_meta_seurat@meta.data[,which(colnames(s1_vp_meta_seurat@meta.data)==paste("RNA_snn_res.",best,sep=""))]
  Idents(s1_vp_meta_seurat) <- "seurat_clusters"
  tiff(paste(pt,"cd45pos_viper_louvain_umap_anchor.tiff",sep="_"), width = 8, height = 8, units = 'in', res = 600)
  plot(DimPlot(s1_vp_meta_seurat, reduction = "umap",label = TRUE) + NoLegend())
  dev.off()
  MRs <- BTTestMRs(s1_vp_meta_seurat@assays$RNA@scale.data, s1_vp_meta_seurat$seurat_clusters)
  top10=MR_UnWrap(MRs, top = 5)
  s1_vp_meta_seurat$tissue=s1$tissue[names(s1_vp_meta_seurat$seurat_clusters)]
  s1_vp_meta_seurat$patient=s1$patient[names(s1_vp_meta_seurat$seurat_clusters)]
  s1_vp_meta_seurat$gene_clustering=s1$seurat_clusters[names(s1_vp_meta_seurat$seurat_clusters)]
  s1_vp_meta_seurat$l=s1$l[names(s1_vp_meta_seurat$seurat_clusters)]
  tiff(paste(pt,"cd45pos_viper_louvain_heatmap_anchor.tiff",sep="_"), width = 8, height = 8, units = 'in', res = 600)
  geneHeatmap_viper(s1_vp_meta_seurat,s1_vp_meta_seurat$seurat_clusters,top10)
  dev.off()
  saveRDS(s1_vp_meta_seurat, file = paste(pt,"cd45pos_viper_meta_seurat.rds",sep=""))
  rm(s1_vp_meta_seurat,s1_vp_meta)
  
  if(pt %in% cd45neg_patients){
    print(paste(pt,"CD45neg VIPER",sep=" "))
    s2_vp_meta=readRDS(paste(pt,"_cd45neg_vp.rds",sep=""))
    s2=readRDS(paste(pt,"_CD45neg.rds",sep=""))
    #cbcMRs <- CBCMRs(s2_vp_meta)
    #s2_vp_meta <- s2_vp_meta[cbcMRs ,]
    s2_vp_meta_seurat <- CreateSeuratObject(counts = s2_vp_meta)
    s2_vp_meta_seurat@assays$RNA@scale.data=as.matrix(s2_vp_meta_seurat@assays$RNA@data)
    s2_vp_meta_seurat <- RunPCA(s2_vp_meta_seurat,features=rownames(s2_vp_meta_seurat))
    s2_vp_meta_seurat <- RunUMAP(s2_vp_meta_seurat, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
    s2_vp_meta_seurat <- FindNeighbors(s2_vp_meta_seurat, dims = 1:30, verbose = FALSE)
    s2_vp_meta_seurat <- FindClusters(s2_vp_meta_seurat, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
    clust=s2_vp_meta_seurat@meta.data[,which(grepl("RNA_snn_res.",colnames(s2_vp_meta_seurat@meta.data)))]
    mat=s2_vp_meta_seurat@assays$RNA@scale.data
    out=sil_subsample_v2_viper(mat,clust)
    means=out[[1]]
    sd=out[[2]]
    x=seq(0.01,1,by=0.01)
    tiff(paste(pt,"cd45neg_viper_louvain_resolution_anchor.tiff",sep="_"), width = 8, height = 6, units = 'in', res = 600)
    errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
    lines(x,means)
    best=tail(x[which(means==max(means))],n=1)
    legend("topright",paste("Best",best,sep = " = "))
    dev.off()
    s2_vp_meta_seurat$seurat_clusters=s2_vp_meta_seurat@meta.data[,which(colnames(s2_vp_meta_seurat@meta.data)==paste("RNA_snn_res.",best,sep=""))]
    Idents(s2_vp_meta_seurat) <- "seurat_clusters"
    tiff(paste(pt,"cd45neg_viper_louvain_umap_anchor.tiff",sep="_"), width = 8, height = 8, units = 'in', res = 600)
    plot(DimPlot(s2_vp_meta_seurat, reduction = "umap",label = TRUE) + NoLegend())
    dev.off()
    MRs <- BTTestMRs(s2_vp_meta_seurat@assays$RNA@scale.data, s2_vp_meta_seurat$seurat_clusters)
    top10=MR_UnWrap(MRs, top = 5)
    s2_vp_meta_seurat$tissue=s2$tissue[names(s2_vp_meta_seurat$seurat_clusters)]
    s2_vp_meta_seurat$patient=s2$patient[names(s2_vp_meta_seurat$seurat_clusters)]
    s2_vp_meta_seurat$gene_clustering=s2$seurat_clusters[names(s2_vp_meta_seurat$seurat_clusters)]
    s2_vp_meta_seurat$l=s2$l[names(s2_vp_meta_seurat$seurat_clusters)]
    tiff(paste(pt,"cd45neg_viper_louvain_heatmap_anchor.tiff",sep="_"), width = 8, height = 8, units = 'in', res = 600)
    geneHeatmap_viper(s2_vp_meta_seurat,s2_vp_meta_seurat$seurat_clusters,top10)
    dev.off()
    saveRDS(s2_vp_meta_seurat, file = paste(pt,"cd45neg_viper_meta_seurat.rds",sep=""))
    rm(s2_vp_meta_seurat,s2_vp_meta)}
}


#####
#####merged VIPER analysis
#####
list_cd45pos=list()
list_cd45neg=list()
patients=c("Patient1","Patient2","Patient3","Patient4","Patient5","Patient6","Patient7","Patient8","PatientA","PatientB","PatientC")
for(pt in patients){
  s_allcells=readRDS(paste(pt,"cd45pos_viper_meta_seurat.rds",sep=""))
  list_cd45pos=c(list_cd45pos,s_allcells)
  rm(s_allcells)
  if(!(pt %in% c("PatientA","PatientB","PatientC"))){
    s2_allcells=readRDS(paste(pt,"cd45neg_viper_meta_seurat.rds",sep=""))
    list_cd45neg=c(list_cd45neg,s2_allcells)
    rm(s2_allcells)}
}
pos_genes=lapply(list_cd45pos,rownames)
pos_genes_intersect=pos_genes[[1]]
for(i in 2:length(pos_genes)){
  pos_genes_intersect=intersect(pos_genes_intersect,pos_genes[[i]])
}
neg_genes=lapply(list_cd45neg,rownames)
neg_genes_intersect=neg_genes[[1]]
for(i in 2:length(neg_genes)){
  neg_genes_intersect=intersect(neg_genes_intersect,neg_genes[[i]])
}
list_cd45pos=lapply(list_cd45pos,function(x){x[pos_genes_intersect,]})
list_cd45neg=lapply(list_cd45neg,function(x){x[neg_genes_intersect,]})


cd45pos_gene=readRDS("cd45pos_merged_seurat.rds")
cd45pos=merge(list_cd45pos[[1]],y=list_cd45pos[2:length(list_cd45pos)],project = "RCC_SC_VIPER")
rm(list_cd45pos)
cbcMRs <- CBCMRs(cd45pos@assays$RNA@data)
cd45pos <- cd45pos@assays$RNA@data[cbcMRs,]
cd45pos <- CreateSeuratObject(counts = cd45pos)
cd45pos@assays$RNA@scale.data=as.matrix(cd45pos@assays$RNA@data)
cd45pos <- RunPCA(cd45pos,features=rownames(cd45pos))
cd45pos <- RunUMAP(cd45pos, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
cd45pos <- FindNeighbors(cd45pos, dims = 1:50, verbose = FALSE)
cd45pos <- FindClusters(cd45pos, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=cd45pos@meta.data[,which(grepl("RNA_snn_res.",colnames(cd45pos@meta.data)))]
mat=cd45pos@assays$RNA@scale.data
out=sil_subsample_v2_viper(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
tiff("cd45pos_merged_louvain_resolution_seurat_viper_meta.tiff", width = 6, height = 6, units = 'in', res = 600)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
dev.off()
cd45pos$tissue=cd45pos_gene$tissue[names(cd45pos$seurat_clusters)]
cd45pos$patient=cd45pos_gene$patient[names(cd45pos$seurat_clusters)]
cd45pos$gene_clustering=cd45pos_gene$seurat_clusters[names(cd45pos$seurat_clusters)]
cd45pos$l=cd45pos_gene$l[names(cd45pos$seurat_clusters)]
cd45pos$seurat_clusters=cd45pos@meta.data[,which(colnames(cd45pos@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(cd45pos) <- "seurat_clusters"
#
# iterative louvain sub-clustering
cd45pos$original_seurat_clusters=cd45pos$seurat_clusters
cluster_output=as.character(cd45pos$seurat_clusters)
for(i in unique(as.character(cd45pos$seurat_clusters))){
  print(i)
  s <- subset(cd45pos, subset = seurat_clusters == i)
  if(ncol(s)/ncol(cd45pos)>0.05){
    s <- RunPCA(s, features = rownames(s))
    s <- RunUMAP(s, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
    s <- FindNeighbors(s, dims = 1:50, verbose = FALSE)
    s <- FindClusters(s, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
    clust=s@meta.data[,which(grepl("RNA_snn_res.",colnames(s@meta.data)))]
    mat=s@assays$RNA@scale.data
    out=sil_subsample_v2_viper(mat,clust)
    means=out[[1]]
    sd=out[[2]]
    x=seq(0.01,1,by=0.01)
    #best=tail(x[5:length(x)][which(means[5:length(means)]==max(means[5:length(means)]))],n=1)
    best=tail(x[which(means==max(means))],n=1)
    s$seurat_clusters=s@meta.data[,which(colnames(s@meta.data)==paste("RNA_snn_res.",best,sep=""))]
    if(max(means)>0.25 & min(table(s$seurat_clusters))/sum(table(s$seurat_clusters))>0.05){cluster_output[which(cluster_output==i)]=paste(i,s$seurat_clusters)}
  }
}
cd45pos$seurat_clusters=as.factor(as.numeric(as.factor(cluster_output)))
Idents(cd45pos) <- "seurat_clusters"
#end iterative louvain sub-clustering
#
cd45pos$seurat_clusters=mapvalues(cd45pos$seurat_clusters, from = 1:length(unique(cd45pos$seurat_clusters)), to = c("CD4 T-cell","Treg", "NK cell 1", "NK cell 2","Macrophage","Monocyte 1", "Monocyte 2","Monocyte 3","CD8 T-cell 1","CD8 T-cell 2", "B cell","Mast cell"))
Idents(cd45pos) <- "seurat_clusters"
#
tiff("cd45pos_merged_louvain_split_umap_seurat_viper_meta.tiff", width = 10, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos,reduction="umap",group.by="seurat_clusters",split.by="tissue",cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"))+NoLegend())
dev.off()
tiff("cd45pos_merged_louvain_split_patient_umap_seurat_viper_meta.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos,reduction="umap",group.by="seurat_clusters",split.by="patient",ncol=3,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3")))
dev.off()
tiff("cd45pos_merged_louvain_umap_seurat_viper_meta.tiff", width = 6, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos, reduction = "umap",label = TRUE,label.size=7,repel=T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3")) + NoLegend())
dev.off()
MRs <- BTTestMRs(cd45pos@assays$RNA@scale.data, cd45pos$seurat_clusters)
top10=MR_UnWrap(MRs, top = 5)
tiff("cd45pos_merged_louvain_heatmap_seurat_viper_meta.tiff", width = 8, height = 8, units = 'in', res = 600)
geneHeatmap(cd45pos,cd45pos$seurat_clusters,top10,viper = T,scaled = T,color_palette = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"))
dev.off()
saveRDS(cd45pos, file = "cd45pos_merged_seurat_viper.rds")
Idents(cd45pos) <- "l"
tiff("cd45pos_merged_singler_umap_seurat_viper_meta.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos, reduction = "umap",label=TRUE,repel=T,label.size=5)+NoLegend())
dev.off()
#surfaceMRs <- BTTestMRs(cd45pos@assays$RNA@scale.data[intersect(surfacemarkers,rownames(cd45pos)),], cd45pos$seurat_clusters)
#top10=MR_UnWrap(surfaceMRs, top = 10)
#tiff("cd45pos_merged_louvain_surfaceheatmap_seurat_viper_meta.tiff", width = 8, height = 10, units = 'in', res = 600)
#geneHeatmap_viper(cd45pos,cd45pos$seurat_clusters,top10)
#dev.off()
rm(list_cd45pos,cd45pos,cd45pos_gene)


cd45neg_gene=readRDS("cd45neg_merged_seurat.rds")
cd45neg=merge(list_cd45neg[[1]],y=list_cd45neg[2:length(list_cd45neg)],project = "RCC_SC_VIPER")
rm(list_cd45neg)
cbcMRs <- CBCMRs(cd45neg@assays$RNA@data)
cd45neg <- cd45neg@assays$RNA@data[cbcMRs,]
cd45neg <- CreateSeuratObject(counts = cd45neg)
cd45neg@assays$RNA@scale.data=as.matrix(cd45neg@assays$RNA@data)
cd45neg <- RunPCA(cd45neg,features=rownames(cd45neg))
cd45neg <- RunUMAP(cd45neg, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
cd45neg <- FindNeighbors(cd45neg, dims = 1:50, verbose = FALSE)
cd45neg <- FindClusters(cd45neg, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=cd45neg@meta.data[,which(grepl("RNA_snn_res.",colnames(cd45neg@meta.data)))]
mat=cd45neg@assays$RNA@scale.data
out=sil_subsample_v2_viper(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
tiff("cd45neg_merged_louvain_resolution_seurat_viper_meta.tiff", width = 6, height = 6, units = 'in', res = 600)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
dev.off()
cd45neg$tissue=cd45neg_gene$tissue[names(cd45neg$seurat_clusters)]
cd45neg$patient=cd45neg_gene$patient[names(cd45neg$seurat_clusters)]
cd45neg$gene_clustering=cd45neg_gene$seurat_clusters[names(cd45neg$seurat_clusters)]
cd45neg$l=cd45neg_gene$l[names(cd45neg$seurat_clusters)]
cd45neg$seurat_clusters=cd45neg@meta.data[,which(colnames(cd45neg@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(cd45neg) <- "seurat_clusters"
#
cd45neg$seurat_clusters=mapvalues(cd45neg$seurat_clusters, from = 1:length(unique(cd45neg$seurat_clusters))-1, to = c("Epithelial 1","Fibroblast", "Endothelial", "M2 Macrophage","Epithelial 2","Epithelial 3", "Epithelial 4"))
Idents(cd45neg) <- "seurat_clusters"
#
tiff("cd45neg_merged_louvain_split_umap_seurat_viper_meta.tiff", width = 10, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg,reduction="umap",group.by="seurat_clusters",split.by="tissue",cols = c("chocolate2","gold2","chartreuse3","lemonchiffon3","darkslategray3","coral3","darksalmon"))+NoLegend())
dev.off()
tiff("cd45neg_merged_louvain_split_patient_umap_seurat_viper_meta.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg,reduction="umap",group.by="seurat_clusters",split.by="patient",ncol=3,cols = c("chocolate2","gold2","chartreuse3","lemonchiffon3","darkslategray3","coral3","darksalmon")))
dev.off()
tiff("cd45neg_merged_louvain_umap_seurat_viper_meta.tiff", width = 6, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg, reduction = "umap",label = TRUE,label.size=7,repel=T,cols = c("chocolate2","gold2","chartreuse3","lemonchiffon3","darkslategray3","coral3","darksalmon")) + NoLegend())
dev.off()
MRs <- BTTestMRs(cd45neg@assays$RNA@scale.data, cd45neg$seurat_clusters)
top10=MR_UnWrap(MRs, top = 5)
tiff("cd45neg_merged_louvain_heatmap_seurat_viper_meta.tiff", width = 8, height = 8, units = 'in', res = 600)
geneHeatmap(cd45neg,cd45neg$seurat_clusters,top10,viper = T,scaled = T,color_palette = c("chocolate2","gold2","chartreuse3","lemonchiffon3","darkslategray3","coral3","darksalmon"))
dev.off()
saveRDS(cd45neg, file = "cd45neg_merged_seurat_viper.rds")
Idents(cd45neg) <- "l"
tiff("cd45neg_merged_singler_umap_seurat_viper_meta.tiff", width = 8, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg, reduction = "umap",label=TRUE,repel=T,label.size=5)+NoLegend())
dev.off()
#surfaceMRs <- BTTestMRs(cd45neg@assays$RNA@scale.data[intersect(surfacemarkers,rownames(cd45neg)),], cd45neg$seurat_clusters)
#top10=MR_UnWrap(surfaceMRs, top = 10)
#tiff("cd45neg_merged_louvain_surfaceheatmap_seurat_viper_meta.tiff", width = 8, height = 10, units = 'in', res = 600)
#geneHeatmap_viper(cd45neg,cd45neg$seurat_clusters,top10)
#dev.off()
rm(list_cd45neg,cd45neg,cd45neg_gene)


###Violinplots
Idents(cd45pos)=cd45pos$seurat_clusters
Idents(cd45neg)=cd45neg$seurat_clusters
Idents(cd45pos.integrated)=cd45pos.integrated$seurat_clusters
Idents(cd45neg.integrated)=cd45neg.integrated$seurat_clusters
tiff("cd45pos_vp_markergene_violinplots.tiff", width = 12, height = 9, units = 'in', res = 600)
#VlnPlot(cd45pos,c("LILRB5","APOE","FOXP3","LAG3","TOX2","CTLA4","PDCD1","CD8A","CD8B"),ncol = 3,same.y.lims = T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"),pt.size = 0,combine=T)+geom_hline(yintercept = 0)
p1=VlnPlot(cd45pos,"LILRB5",same.y.lims = T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
p2=VlnPlot(cd45pos,"APOE",same.y.lims = T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
p3=VlnPlot(cd45pos,"FOXP3",same.y.lims = T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
p4=VlnPlot(cd45pos,"LAG3",same.y.lims = T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
p5=VlnPlot(cd45pos,"TOX2",same.y.lims = T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
p6=VlnPlot(cd45pos,"CTLA4",same.y.lims = T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
p7=VlnPlot(cd45pos,"PDCD1",same.y.lims = T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
p8=VlnPlot(cd45pos,"CD8A",same.y.lims = T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
p9=VlnPlot(cd45pos,"CD8B",same.y.lims = T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3)
dev.off()
tiff("cd45neg_vp_markergene_violinplots.tiff", width = 12, height = 3, units = 'in', res = 600)
#VlnPlot(cd45neg,c("PAX8","PAX2","CA9"),ncol = 3,same.y.lims=T,cols = c("chocolate2","gold2","chartreuse3","lemonchiffon3","darkslategray3","coral3","darksalmon"),pt.size = 0,combine=T)+geom_hline(yintercept = 0)
p1=VlnPlot(cd45neg,"PAX8",same.y.lims = T,cols=c("chocolate2","gold2","chartreuse3","lemonchiffon3","darkslategray3","coral3","darksalmon"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
p2=VlnPlot(cd45neg,"PAX2",same.y.lims = T,cols=c("chocolate2","gold2","chartreuse3","lemonchiffon3","darkslategray3","coral3","darksalmon"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
cd45neg.integrated=cd45neg.integrated[,colnames(cd45neg)]
cd45neg.integrated$seurat_clusters=cd45neg$seurat_clusters
Idents(cd45neg.integrated)=cd45neg.integrated$seurat_clusters
p3=VlnPlot(cd45neg.integrated,"CA9",same.y.lims = T,cols=c("chocolate2","gold2","chartreuse3","lemonchiffon3","darkslategray3","coral3","darksalmon"),pt.size = 0,combine=T)+NoLegend()+geom_hline(yintercept = 0)
plot_grid(p1,p2,p3,ncol=3)
dev.off()
tiff("cd45pos_geneexp_markergene_violinplots.tiff", width = 12, height = 6, units = 'in', res = 600)
VlnPlot(cd45pos.integrated,c("C1QA","C1QB","C1QC","APOE","TREM2"),cols=c("cornflowerblue","coral3"),ncol = 3,same.y.lims=T,split.by="tissue",pt.size = 0) #cols = c("darkgoldenrod2","darkorange3","mediumorchid2","forestgreen","chartreuse3","cyan4","cyan2","lightsteelblue3","cornsilk3","deeppink3","lightsteelblue4")
dev.off()


##UMAP plots colored by alternative clustering
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
cd45pos_vp$gene_clustering=cd45pos$seurat_clusters
cd45neg_vp$gene_clustering=cd45neg$seurat_clusters
cd45pos$seurat_clusters=cd45pos_vp$seurat_clusters
cd45neg$seurat_clusters=cd45neg_vp$seurat_clusters
Idents(cd45pos) <- "seurat_clusters"
Idents(cd45neg) <- "seurat_clusters"
Idents(cd45pos_vp)="gene_clustering"
Idents(cd45neg_vp)="gene_clustering"
tiff("cd45pos_merged_louvain_umap_seurat_vpplot_ColoredByGexCluster.tiff", width = 6, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos_vp, reduction = "umap",label = TRUE,label.size=7,repel=T,cols = c("darkgoldenrod2","darkorange3","mediumorchid2","forestgreen","chartreuse3","cyan4","cyan2","lightsteelblue3","cornsilk3","deeppink3","lightsteelblue4")) + NoLegend())
dev.off()
tiff("cd45neg_merged_louvain_umap_seurat_vpplot_ColoredByGexCluster.tiff", width = 6, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg_vp, reduction = "umap",label = TRUE,label.size=7,repel=T,cols = c("chocolate2","chartreuse3","gold2","lemonchiffon3","pink1")) + NoLegend())
dev.off()
tiff("cd45neg_merged_louvain_umap_seurat_geneexpplot_ColoredByVPCluster.tiff", width = 6, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45neg, reduction = "umap",label = TRUE,label.size=7,repel=T,cols = c("chocolate2","gold2","chartreuse3","lemonchiffon3","darkslategray3","coral3","darksalmon")) + NoLegend())
dev.off()
tiff("cd45pos_merged_louvain_umap_seurat_geneexpplot_ColoredByVPCluster.tiff", width = 6, height = 6, units = 'in', res = 600)
plot(DimPlot(cd45pos, reduction = "umap",label = TRUE,label.size=7,repel=T,cols = c("darkgoldenrod2","darkorange3","forestgreen","chartreuse3","cyan4","cyan2","deepskyblue","dodgerblue3","mediumorchid2","maroon2","lightsteelblue3","deeppink3")) + NoLegend())
dev.off()
rm(cd45pos,cd45neg,cd45pos_vp,cd45neg_vp)



###INFERCNV:: all patients -- gene expression clusters
cd45pos=readRDS("cd45pos_merged_seurat.rds")
Idents(cd45pos)="seurat_clusters"
cd45pos=cd45pos[,sample(colnames(cd45pos),ncol(cd45pos)/10)]
cd45pos$seurat_clusters=paste("cd45pos",cd45pos$seurat_clusters,sep="_")
cd45neg=readRDS("cd45neg_merged_seurat.rds")
Idents(cd45neg)="seurat_clusters"
cd45neg=cd45neg[,sample(colnames(cd45neg),ncol(cd45neg)/10)]
cd45neg$seurat_clusters=paste("cd45neg",cd45neg$seurat_clusters,sep="_")
cd45pos2=CreateSeuratObject(counts=cd45pos@assays$SCT@counts)
cd45pos2$seurat_clusters=cd45pos$seurat_clusters
cd45pos2$cd45=cd45pos$cd45
cd45neg2=CreateSeuratObject(counts=cd45neg@assays$SCT@counts)
cd45neg2$seurat_clusters=cd45neg$seurat_clusters
cd45neg2$cd45=cd45neg$cd45
combined=merge(cd45pos2,cd45neg2)
write.table(combined@assays$RNA@counts,file="raw_counts_matrix.txt",sep="\t",quote=F)
write.table(combined$seurat_clusters,file="annotations.txt",sep="\t",quote=F,col.names = F)
infercnv_combined = CreateInfercnvObject(raw_counts_matrix= 'raw_counts_matrix.txt',
                                         annotations_file='annotations.txt',
                                         gene_order_file='chromosome_locations_noDupGenename_armloc.txt',
                                         delim = '\t',
                                         ref_group_names=unique(combined$seurat_clusters[which(combined$cd45=="CD45+")]))
infercnv_combined = infercnv::run(infercnv_combined,
                                  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                  out_dir='infercnv_results_combined', 
                                  cluster_by_groups=TRUE, 
                                  denoise=TRUE,
                                  noise_logistic = TRUE, 
                                  scale_data = FALSE,
                                  HMM=F, 
                                  up_to_step = 50, 
                                  k_obs_groups = 1, 
                                  HMM_type = 'i6',
                                  analysis_mode = 'subclusters', 
                                  no_plot = F, 
                                  plot_steps = T,
                                  tumor_subcluster_partition_method = 'qnorm', 
                                  num_threads = 11, 
                                  debug = F)
saveRDS(infercnv_combined,"infercnv_combined.rds")
head(infercnv_combined@gene_order)
infercnv_combined@expr.data[1:5,1:5]

###INFERCNV:: all patients -- viper clusters
cd45pos=readRDS("cd45pos_merged_seurat.rds")
cd45pos_vp=readRDS("cd45pos_merged_seurat_viper.rds")
Idents(cd45pos)="seurat_clusters"
cd45pos$genexp_clusters=paste("cd45pos",cd45pos$seurat_clusters,sep="_")
Idents(cd45pos_vp)="seurat_clusters"
cd45pos=cd45pos[,colnames(cd45pos_vp)]
cd45pos=cd45pos[,sample(colnames(cd45pos),ncol(cd45pos)/10)]
cd45pos_vp=cd45pos_vp[,colnames(cd45pos)]
cd45pos$seurat_clusters=paste("cd45pos",cd45pos_vp$seurat_clusters,sep="_")
cd45neg=readRDS("cd45neg_merged_seurat.rds")
cd45neg_vp=readRDS("cd45neg_merged_seurat_viper.rds")
Idents(cd45neg)="seurat_clusters"
cd45neg$genexp_clusters=paste("cd45neg",cd45neg$seurat_clusters,sep="_")
Idents(cd45neg_vp)="seurat_clusters"
cd45neg=cd45neg[,colnames(cd45neg_vp)]
cd45neg=cd45neg[,sample(colnames(cd45neg),ncol(cd45neg)/10)]
cd45neg_vp=cd45neg_vp[,colnames(cd45neg)]
cd45neg$seurat_clusters=paste("cd45neg",cd45neg_vp$seurat_clusters,sep="_")
cd45pos2=CreateSeuratObject(counts=cd45pos@assays$SCT@counts)
cd45pos2$seurat_clusters=cd45pos$seurat_clusters
cd45pos2$genexp_clusters=cd45pos$genexp_clusters
cd45pos2$cd45=cd45pos$cd45
cd45neg2=CreateSeuratObject(counts=cd45neg@assays$SCT@counts)
cd45neg2$seurat_clusters=cd45neg$seurat_clusters
cd45neg2$genexp_clusters=cd45neg$genexp_clusters
cd45neg2$cd45=cd45neg$cd45
combined=merge(cd45pos2,cd45neg2)
combined=combined[,order(combined$seurat_clusters)]
write.table(combined@assays$RNA@counts[,order(combined$seurat_clusters)],file="raw_counts_matrix.txt",sep="\t",quote=F)
write.table(combined$seurat_clusters[order(combined$seurat_clusters)],file="annotations.txt",sep="\t",quote=F,col.names = F)
infercnv_combined = CreateInfercnvObject(raw_counts_matrix= 'raw_counts_matrix.txt',
                                         annotations_file='annotations.txt',
                                         gene_order_file='/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/chromosome_locations_noDupGenename_armloc.txt',
                                         delim = '\t',
                                         ref_group_names=unique(combined$seurat_clusters[which(combined$cd45=="CD45+")]))
infercnv_combined = infercnv::run(infercnv_combined,
                                  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                  out_dir='infercnv_results_combined_viper_sameorderasgeneexp', 
                                  cluster_by_groups=TRUE, 
                                  denoise=TRUE,
                                  noise_logistic = TRUE, 
                                  scale_data = FALSE,
                                  HMM=F, 
                                  up_to_step = 50, 
                                  k_obs_groups = 1, 
                                  HMM_type = 'i6',
                                  analysis_mode = 'samples', 
                                  no_plot = F, 
                                  plot_steps = T,
                                  tumor_subcluster_partition_method = 'qnorm', 
                                  num_threads = 11, 
                                  debug = F,
                                  cluster_references=F)
write.table(combined@assays$RNA@counts[,order(combined$seurat_clusters)],file="raw_counts_matrix.txt",sep="\t",quote=F)
write.table(combined$genexp_clusters[order(combined$seurat_clusters)],file="annotations.txt",sep="\t",quote=F,col.names = F)
infercnv_combined = CreateInfercnvObject(raw_counts_matrix= 'raw_counts_matrix.txt',
                                         annotations_file='annotations.txt',
                                         gene_order_file='/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/chromosome_locations_noDupGenename_armloc.txt',
                                         delim = '\t',
                                         ref_group_names=unique(combined$genexp_clusters[which(combined$cd45=="CD45+")]))
infercnv_combined = infercnv::run(infercnv_combined,
                                  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                  out_dir='infercnv_results_combined_geneexp_sameorderasviper', 
                                  cluster_by_groups=TRUE, 
                                  denoise=TRUE,
                                  noise_logistic = TRUE, 
                                  scale_data = FALSE,
                                  HMM=F, 
                                  up_to_step = 50, 
                                  k_obs_groups = 1, 
                                  HMM_type = 'i6',
                                  analysis_mode = 'samples', 
                                  no_plot = F, 
                                  plot_steps = T,
                                  tumor_subcluster_partition_method = 'qnorm', 
                                  num_threads = 11, 
                                  debug = F,
                                  cluster_references=F)



####Identify known tumor-immune receptor-ligand co-upregulation -- Gene Exp Ligands to VIPER receptors
#database of known receptor-ligand interactions http://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/
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
receptor_ligand_pairs=read.csv("PairsLigRec.csv")
pairs=receptor_ligand_pairs[,c(2,4)]
pairs[,1]=as.character(pairs[,1])
pairs[,2]=as.character(pairs[,2])
pairs=rbind(pairs,c("PDCD1","CD274"))
markers.cd45pos.viper <- FindAllMarkers(cd45pos_vp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "t")
markers.cd45neg.viper <- FindAllMarkers(cd45neg_vp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "t")
markers.cd45pos <- FindAllMarkers(cd45pos, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
markers.cd45neg <- FindAllMarkers(cd45neg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
markers.cd45neg=markers.cd45neg[markers.cd45neg$cluster!="M2 Macrophage",]
markers.cd45neg.viper=markers.cd45neg.viper[markers.cd45neg.viper$cluster!="M2 Macrophage",]
pairs=pairs[which(pairs[,1] %in% union(markers.cd45pos$gene,markers.cd45neg$gene) & pairs[,2] %in% union(markers.cd45pos.viper$gene,markers.cd45neg.viper$gene)),]
pos_clusters=unique(cd45pos$seurat_clusters)
neg_clusters=unique(cd45neg$seurat_clusters)
neg_clusters=setdiff(neg_clusters,"M2 Macrophage")
cluster_labels=c(as.character(pos_clusters),as.character(neg_clusters))
receptor_source=c()
ligand_source=c()
for(i in 1:nrow(pairs)){
  ligand=as.character(pairs[i,1])
  receptor=as.character(pairs[i,2])
  if(ligand %in% rownames(cd45pos)){pos_clusters_ligand_median=lapply(pos_clusters,function(clust){median(cd45pos@assays$integrated@scale.data[ligand,which(cd45pos$seurat_clusters==clust)])})}
  if(!(ligand %in% rownames(cd45pos))){pos_clusters_ligand_median=lapply(pos_clusters,function(clust){-Inf})}
  if(ligand %in% rownames(cd45neg)){neg_clusters_ligand_median=lapply(neg_clusters,function(clust){median(cd45neg@assays$integrated@scale.data[ligand,which(cd45neg$seurat_clusters==clust)])})}
  if(!(ligand %in% rownames(cd45neg))){neg_clusters_ligand_median=lapply(neg_clusters,function(clust){-Inf})}
  if(receptor %in% rownames(cd45pos_vp)){pos_clusters_receptor_median=lapply(pos_clusters,function(clust){median(cd45pos_vp@assays$RNA@scale.data[receptor,which(cd45pos_vp$seurat_clusters==clust)])})}
  if(!(receptor %in% rownames(cd45pos_vp))){pos_clusters_receptor_median=lapply(pos_clusters,function(clust){-Inf})}
  if(receptor %in% rownames(cd45neg_vp)){neg_clusters_receptor_median=lapply(neg_clusters,function(clust){median(cd45neg_vp@assays$RNA@scale.data[receptor,which(cd45neg_vp$seurat_clusters==clust)])})}
  if(!(receptor %in% rownames(cd45neg_vp))){neg_clusters_receptor_median=lapply(neg_clusters,function(clust){-Inf})}
  ligand_medians=c(unlist(pos_clusters_ligand_median),unlist(neg_clusters_ligand_median))
  receptor_medians=c(unlist(pos_clusters_receptor_median),unlist(neg_clusters_receptor_median))
  receptor_source=c(receptor_source,cluster_labels[which(receptor_medians==max(receptor_medians))])
  ligand_source=c(ligand_source,cluster_labels[which(ligand_medians==max(ligand_medians))])
}
pairs$ligand_cluster=ligand_source
pairs$receptor_cluster=receptor_source
write.csv(pairs,file="receptor_ligand_pairs_geneexpLigand_viperReceptor.csv")
pairs[which((pairs[,3] %in% c("Epithelial 1", "Epithelial 3", "Epithelial 4") & pairs[,4] %in% c("Treg","CD4 T-cell","CD8 T-cell 1", "CD8 T-cell 2")) | (pairs[,4] %in% c("Epithelial 1", "Epithelial 3", "Epithelial 4") & pairs[,3] %in% c("Treg","CD4 T-cell","CD8 T-cell 1", "CD8 T-cell 2"))),]
pairs[which((pairs[,3] %in% c("Macrophage") & pairs[,4] %in% c("Treg","CD4 T-cell","CD8 T-cell 1", "CD8 T-cell 2")) | (pairs[,4] %in% c("Macrophage") & pairs[,3] %in% c("Treg","CD4 T-cell","CD8 T-cell 1", "CD8 T-cell 2"))),]
pairs[which((pairs[,3] %in% c("Macrophage") & pairs[,4] %in% c("Epithelial 1", "Epithelial 3", "Epithelial 4")) | (pairs[,4] %in% c("Macrophage") & pairs[,3] %in% c("Epithelial 1", "Epithelial 3", "Epithelial 4"))),]


#######
####GENE EXPRESSION plots of cluster frequencies in each patient, separated by tumor/normal and labelled by stage/grade
#Grade annotations
#"1"-- Patient3, Patient4, 
#"2"-- PatientB, PatientC, Patient1, Patient2, Patient6, Patient8
#"3"-- PatientA, Patient7
#"4"-- Patient5
#Stage annotations
#"pT3a"-- PatientA, Patient2, Patient3, Patient5(metastatic), Patient7
#"pT1b"-- PatientB, PatientC, Patient1, Patient4, Patient6, Patient8
cd45pos=readRDS("cd45pos_merged_seurat.rds")
cd45neg=readRDS("cd45neg_merged_seurat.rds")
dat_norm=as.data.frame.matrix(table(cd45pos$seurat_clusters[which(cd45pos$tissue=="Normal")],cd45pos$patient[which(cd45pos$tissue=="Normal")]))
dat_tumor=as.data.frame.matrix(table(cd45pos$seurat_clusters[which(cd45pos$tissue=="Tumor")],cd45pos$patient[which(cd45pos$tissue=="Tumor")]))
dat2_norm=as.data.frame.matrix(table(cd45neg$seurat_clusters[which(cd45neg$tissue=="Normal")],cd45neg$patient[which(cd45neg$tissue=="Normal")]))
dat2_tumor=as.data.frame.matrix(table(cd45neg$seurat_clusters[which(cd45neg$tissue=="Tumor")],cd45neg$patient[which(cd45neg$tissue=="Tumor")]))
table(cd45pos$seurat_clusters,cd45pos$l)
table(cd45neg$seurat_clusters,cd45neg$l)
Idents(cd45pos) <- "seurat_clusters"
Idents(cd45neg) <- "seurat_clusters"

dat=cbind(dat_norm,dat_tumor)
stage=colnames(dat)
stage[which(stage %in% c("Patient1","Patient4","Patient6","Patient8","PatientB","PatientC"))]="pT1a"
stage[which(stage %in% c("Patient2","Patient3","Patient5","Patient7","PatientA"))]="pT3a"
grade=colnames(dat)
grade[which(grade %in% c("Patient3","Patient4"))]="grade1"
grade[which(grade %in% c("PatientB", "PatientC", "Patient1", "Patient2", "Patient6", "Patient8"))]="grade2"
grade[which(grade %in% c("PatientA", "Patient7"))]="grade3"
grade[which(grade %in% c("Patient5"))]="grade4"
df=data.frame(patient=colnames(dat),tissue=c(rep("Normal",ncol(dat)/2),rep("Tumor",ncol(dat)/2)),stage=stage,grade=grade)
dat=apply(dat,2,function(x){x/sum(x)})
colnames(dat)=make.unique(colnames(dat))
rownames(df)=colnames(dat)
o=order(df$tissue,df$stage,df$grade,df$patient)
dat=dat[,o]
df=df[o,]
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(dat, n = 20)
tiff("clusterfreqs_cd45pos_gene.tiff", width = 8, height = 6, units = 'in', res = 600)
pheatmap(dat, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 6,show_colnames = F,display_numbers = T)
dev.off()
##tumor vs normal
t_vs_n_pvals_gene_cd45pos=apply(dat,1,function(x){t.test(x[which(df$tissue=="Tumor")],x[which(df$tissue=="Normal")],paired=T)$p.value})
t_vs_n_pvals_gene_cd45pos=p.adjust(t_vs_n_pvals_gene_cd45pos,"BH")
t_vs_n_difference_gene_cd45pos=apply(dat,1,function(x){t.test(x[which(df$tissue=="Tumor")],x[which(df$tissue=="Normal")],paired=T)$estimate})
#stage1 vs stage3
stage_pvals_gene_cd45pos=apply(dat,1,function(x){t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$p.value})
stage_pvals_gene_cd45pos=p.adjust(stage_pvals_gene_cd45pos,"BH")
stage_difference_gene_cd45pos=apply(dat,1,function(x){t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$estimate[1]-t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$estimate[2]})
######### boxplots
dat_barplot=as.data.frame(t(dat[,(ncol(dat)/2+1):ncol(dat)]-dat[,1:(ncol(dat)/2)]))
dat_barplot$stage=df[rownames(dat_barplot),3]
library(reshape)
x <- melt(dat_barplot, id=c("stage"))
colnames(x)=c("stage","celltype","TumorMinusNormalFreq")
tiff("clusterfreqs_cd45pos_gene_stage_boxplot.tiff", width = 8, height = 4, units = 'in', res = 600)
ggplot(x, aes(x=celltype, y=TumorMinusNormalFreq,fill=stage)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
dat_barplot=as.data.frame(t(dat[,(ncol(dat)/2+1):ncol(dat)]-dat[,1:(ncol(dat)/2)]))
dat_barplot$grade=df[rownames(dat_barplot),4]
library(reshape)
x <- melt(dat_barplot, id=c("grade"))
colnames(x)=c("grade","celltype","TumorMinusNormalFreq")
tiff("clusterfreqs_cd45pos_gene_grade_boxplot.tiff", width = 8, height = 4, units = 'in', res = 600)
ggplot(x, aes(x=celltype, y=TumorMinusNormalFreq,fill=grade)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
##########

dat=cbind(dat2_norm,dat2_tumor)
stage=colnames(dat)
stage[which(stage %in% c("Patient1","Patient4","Patient6","Patient8","PatientB","PatientC"))]="pT1a"
stage[which(stage %in% c("Patient2","Patient3","Patient5","Patient7","PatientA"))]="pT3a"
grade=colnames(dat)
grade[which(grade %in% c("Patient3","Patient4"))]="grade1"
grade[which(grade %in% c("PatientB", "PatientC", "Patient1", "Patient2", "Patient6", "Patient8"))]="grade2"
grade[which(grade %in% c("PatientA", "Patient7"))]="grade3"
grade[which(grade %in% c("Patient5"))]="grade4"
df=data.frame(patient=colnames(dat),tissue=c(rep("Normal",ncol(dat)/2),rep("Tumor",ncol(dat)/2)),stage=stage,grade=grade)
dat=apply(dat,2,function(x){x/sum(x)})
colnames(dat)=make.unique(colnames(dat))
rownames(df)=colnames(dat)
o=order(df$tissue,df$stage,df$grade,df$patient)
dat=dat[,o]
df=df[o,]
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(dat, n = 20)
#rownames(dat)=c("epithelial/mesangial","endothelial","fibroblasts","M2 macs","adipocytes")
tiff("clusterfreqs_cd45neg_gene.tiff", width = 8, height = 6, units = 'in', res = 600)
pheatmap(dat, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 6,show_colnames = F,display_numbers = T)
dev.off()
##tumor vs normal
t_vs_n_pvals_gene_cd45neg=apply(dat,1,function(x){t.test(x[which(df$tissue=="Tumor")],x[which(df$tissue=="Normal")],paired=T)$p.value})
t_vs_n_pvals_gene_cd45neg=p.adjust(t_vs_n_pvals_gene_cd45neg,"BH")
t_vs_n_difference_gene_cd45neg=apply(dat,1,function(x){t.test(x[which(df$tissue=="Tumor")],x[which(df$tissue=="Normal")],paired=T)$estimate})
#stage1 vs stage3
stage_pvals_gene_cd45neg=apply(dat,1,function(x){t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$p.value})
stage_pvals_gene_cd45neg=p.adjust(stage_pvals_gene_cd45neg,"BH")
stage_difference_gene_cd45neg=apply(dat,1,function(x){t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$estimate[1]-t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$estimate[2]})
######### boxplots
dat_barplot=as.data.frame(t(dat[,(ncol(dat)/2+1):ncol(dat)]-dat[,1:(ncol(dat)/2)]))
dat_barplot$stage=df[rownames(dat_barplot),3]
library(reshape)
x <- melt(dat_barplot, id=c("stage"))
colnames(x)=c("stage","celltype","TumorMinusNormalFreq")
tiff("clusterfreqs_cd45neg_gene_stage_boxplot.tiff", width = 8, height = 4, units = 'in', res = 600)
ggplot(x, aes(x=celltype, y=TumorMinusNormalFreq,fill=stage)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
dat_barplot=as.data.frame(t(dat[,(ncol(dat)/2+1):ncol(dat)]-dat[,1:(ncol(dat)/2)]))
dat_barplot$grade=df[rownames(dat_barplot),4]
library(reshape)
x <- melt(dat_barplot, id=c("grade"))
colnames(x)=c("grade","celltype","TumorMinusNormalFreq")
tiff("clusterfreqs_cd45neg_gene_grade_boxplot.tiff", width = 8, height = 4, units = 'in', res = 600)
ggplot(x, aes(x=celltype, y=TumorMinusNormalFreq,fill=grade)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
##########

####VIPER plots of cluster frequencies in each patient, separated by tumor/normal and labelled by stage/grade
cd45pos_vp=readRDS("cd45pos_merged_seurat_viper.rds")
cd45neg_vp=readRDS("cd45neg_merged_seurat_viper.rds")
dat_norm=as.data.frame.matrix(table(cd45pos_vp$seurat_clusters[which(cd45pos_vp$tissue=="Normal")],cd45pos_vp$patient[which(cd45pos_vp$tissue=="Normal")]))
dat_tumor=as.data.frame.matrix(table(cd45pos_vp$seurat_clusters[which(cd45pos_vp$tissue=="Tumor")],cd45pos_vp$patient[which(cd45pos_vp$tissue=="Tumor")]))
dat2_norm=as.data.frame.matrix(table(cd45neg_vp$seurat_clusters[which(cd45neg_vp$tissue=="Normal")],cd45neg_vp$patient[which(cd45neg_vp$tissue=="Normal")]))
dat2_tumor=as.data.frame.matrix(table(cd45neg_vp$seurat_clusters[which(cd45neg_vp$tissue=="Tumor")],cd45neg_vp$patient[which(cd45neg_vp$tissue=="Tumor")]))
table(cd45pos_vp$seurat_clusters,cd45pos_vp$l)
table(cd45neg_vp$seurat_clusters,cd45neg_vp$l)
Idents(cd45pos_vp) <- "seurat_clusters"
Idents(cd45neg_vp) <- "seurat_clusters"

dat=cbind(dat_norm,dat_tumor)
stage=colnames(dat)
stage[which(stage %in% c("Patient1","Patient4","Patient6","Patient8","PatientB","PatientC"))]="pT1a"
stage[which(stage %in% c("Patient2","Patient3","Patient5","Patient7","PatientA"))]="pT3a"
grade=colnames(dat)
grade[which(grade %in% c("Patient3","Patient4"))]="grade1"
grade[which(grade %in% c("PatientB", "PatientC", "Patient1", "Patient2", "Patient6", "Patient8"))]="grade2"
grade[which(grade %in% c("PatientA", "Patient7"))]="grade3"
grade[which(grade %in% c("Patient5"))]="grade4"
df=data.frame(patient=colnames(dat),tissue=c(rep("Normal",ncol(dat)/2),rep("Tumor",ncol(dat)/2)),stage=stage,grade=grade)
dat=apply(dat,2,function(x){x/sum(x)})
colnames(dat)=make.unique(colnames(dat))
rownames(df)=colnames(dat)
o=order(df$tissue,df$stage,df$grade,df$patient)
dat=dat[,o]
df=df[o,]
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(dat, n = 20)
#rownames(dat)=c("CD4 T-cells","Tregs","NK cells 1","NK cells 2","Macrophages","Monocytes1","Monocytes2","Monocytes3","CD8 T-cells 1","CD8 T-Cells 2","B Cells","Mast Cells")
tiff("clusterfreqs_cd45pos_viper.tiff", width = 8, height = 6, units = 'in', res = 600)
pheatmap(dat, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 6,show_colnames = F,display_numbers = T)
dev.off()
##tumor vs normal
t_vs_n_pvals_viper_cd45pos=apply(dat,1,function(x){t.test(x[which(df$tissue=="Tumor")],x[which(df$tissue=="Normal")],paired=T)$p.value})
t_vs_n_pvals_viper_cd45pos=p.adjust(t_vs_n_pvals_viper_cd45pos,"BH")
t_vs_n_difference_viper_cd45pos=apply(dat,1,function(x){t.test(x[which(df$tissue=="Tumor")],x[which(df$tissue=="Normal")],paired=T)$estimate})
#stage1 vs stage3
stage_pvals_viper_cd45pos=apply(dat,1,function(x){t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$p.value})
stage_pvals_viper_cd45pos=p.adjust(stage_pvals_viper_cd45pos,"BH")
stage_difference_viper_cd45pos=apply(dat,1,function(x){t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$estimate[1]-t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$estimate[2]})
######### boxplots
dat_barplot=as.data.frame(t(dat[,(ncol(dat)/2+1):ncol(dat)]-dat[,1:(ncol(dat)/2)]))
dat_barplot$stage=df[rownames(dat_barplot),3]
library(reshape)
x <- melt(dat_barplot, id=c("stage"))
colnames(x)=c("stage","celltype","TumorMinusNormalFreq")
tiff("clusterfreqs_cd45pos_vp_stage_boxplot.tiff", width = 8, height = 4, units = 'in', res = 600)
ggplot(x, aes(x=celltype, y=TumorMinusNormalFreq,fill=stage)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
dat_barplot=as.data.frame(t(dat[,(ncol(dat)/2+1):ncol(dat)]-dat[,1:(ncol(dat)/2)]))
dat_barplot$grade=df[rownames(dat_barplot),4]
library(reshape)
x <- melt(dat_barplot, id=c("grade"))
colnames(x)=c("grade","celltype","TumorMinusNormalFreq")
tiff("clusterfreqs_cd45pos_vp_grade_boxplot.tiff", width = 8, height = 4, units = 'in', res = 600)
ggplot(x, aes(x=celltype, y=TumorMinusNormalFreq,fill=grade)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
##########

dat=cbind(dat2_norm,dat2_tumor)
stage=colnames(dat)
stage[which(stage %in% c("Patient1","Patient4","Patient6","Patient8","PatientB","PatientC"))]="pT1a"
stage[which(stage %in% c("Patient2","Patient3","Patient5","Patient7","PatientA"))]="pT3a"
grade=colnames(dat)
grade[which(grade %in% c("Patient3","Patient4"))]="grade1"
grade[which(grade %in% c("PatientB", "PatientC", "Patient1", "Patient2", "Patient6", "Patient8"))]="grade2"
grade[which(grade %in% c("PatientA", "Patient7"))]="grade3"
grade[which(grade %in% c("Patient5"))]="grade4"
df=data.frame(patient=colnames(dat),tissue=c(rep("Normal",ncol(dat)/2),rep("Tumor",ncol(dat)/2)),stage=stage,grade=grade)
dat=apply(dat,2,function(x){x/sum(x)})
colnames(dat)=make.unique(colnames(dat))
rownames(df)=colnames(dat)
o=order(df$tissue,df$stage,df$grade,df$patient)
dat=dat[,o]
df=df[o,]
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(dat, n = 20)
#rownames(dat)=c("epithelial/mesangial 1","fibroblast","endothelial","M2 macs","epithelial/mesangial 2","epithelial/mesangial 3","epithelial/mesangial 4")
tiff("clusterfreqs_cd45neg_viper.tiff", width = 8, height = 6, units = 'in', res = 600)
pheatmap(dat, cluster_rows=FALSE,show_rownames=T,cluster_cols=FALSE, annotation_col=df,breaks=mat_breaks,color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),fontsize_row = 6,show_colnames = F,display_numbers = T)
dev.off()
##tumor vs normal
t_vs_n_pvals_viper_cd45neg=apply(dat,1,function(x){t.test(x[which(df$tissue=="Tumor")],x[which(df$tissue=="Normal")],paired=T)$p.value})
t_vs_n_pvals_viper_cd45neg=p.adjust(t_vs_n_pvals_viper_cd45neg,"BH")
t_vs_n_difference_viper_cd45neg=apply(dat,1,function(x){t.test(x[which(df$tissue=="Tumor")],x[which(df$tissue=="Normal")],paired=T)$estimate})
#stage1 vs stage3
stage_pvals_viper_cd45neg=apply(dat,1,function(x){t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$p.value})
stage_pvals_viper_cd45neg=p.adjust(stage_pvals_viper_cd45neg,"BH")
stage_difference_viper_cd45neg=apply(dat,1,function(x){t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$estimate[1]-t.test(x[which(df$stage=="pT3a" & df$tissue=="Tumor")],x[which(df$stage=="pT1a" & df$tissue=="Tumor")])$estimate[2]})
######### boxplots
dat_barplot=as.data.frame(t(dat[,(ncol(dat)/2+1):ncol(dat)]-dat[,1:(ncol(dat)/2)]))
dat_barplot$stage=df[rownames(dat_barplot),3]
library(reshape)
x <- melt(dat_barplot, id=c("stage"))
colnames(x)=c("stage","celltype","TumorMinusNormalFreq")
tiff("clusterfreqs_cd45neg_vp_stage_boxplot.tiff", width = 8, height = 4, units = 'in', res = 600)
ggplot(x, aes(x=celltype, y=TumorMinusNormalFreq,fill=stage)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
dat_barplot=as.data.frame(t(dat[,(ncol(dat)/2+1):ncol(dat)]-dat[,1:(ncol(dat)/2)]))
dat_barplot$grade=df[rownames(dat_barplot),4]
library(reshape)
x <- melt(dat_barplot, id=c("grade"))
colnames(x)=c("grade","celltype","TumorMinusNormalFreq")
tiff("clusterfreqs_cd45neg_vp_grade_boxplot.tiff", width = 8, height = 4, units = 'in', res = 600)
ggplot(x, aes(x=celltype, y=TumorMinusNormalFreq,fill=grade)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
##########

FeaturePlot(cd45pos,c("APOE"),cols = c("white","red"))
VlnPlot(cd45pos,"APOE",split.by="tissue",pt.size=0)


#Run infercnv on only tumor cells with normal cluster reference, 
###INFERCNV:: all patients -- viper clusters
cd45neg=readRDS("cd45neg_merged_seurat.rds")
cd45neg_vp=readRDS("cd45neg_merged_seurat_viper.rds")
Idents(cd45neg)="seurat_clusters"
Idents(cd45neg_vp)="seurat_clusters"
cd45neg_vp=cd45neg_vp[,which(cd45neg_vp$seurat_clusters %in% c("Epithelial 1", "Epithelial 2", "Epithelial 3", "Epithelial 4"))]
cd45neg=cd45neg[,colnames(cd45neg_vp)]
cd45neg=cd45neg[,sample(colnames(cd45neg),ncol(cd45neg)/10)]
cd45neg_vp=cd45neg_vp[,colnames(cd45neg)]
cd45neg$seurat_clusters=cd45neg_vp$seurat_clusters
write.table(cd45neg@assays$SCT@counts,file="raw_counts_matrix.txt",sep="\t",quote=F)
write.table(cd45neg$seurat_clusters,file="annotations.txt",sep="\t",quote=F,col.names = F)
infercnv_combined = CreateInfercnvObject(raw_counts_matrix= 'raw_counts_matrix.txt',
                                         annotations_file='annotations.txt',
                                         gene_order_file='chromosome_locations_noDupGenename_armloc.txt',
                                         delim = '\t',
                                         ref_group_names="Epithelial 2")
infercnv_combined = infercnv::run(infercnv_combined,
                                  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                  out_dir='infercnv_results_combined_viper_normalEpithelialAsReference', 
                                  cluster_by_groups=TRUE, 
                                  denoise=TRUE,
                                  noise_logistic = TRUE, 
                                  scale_data = FALSE,
                                  HMM=F, 
                                  up_to_step = 50, 
                                  k_obs_groups = 1, 
                                  HMM_type = 'i6',
                                  analysis_mode = 'subclusters', 
                                  no_plot = F, 
                                  plot_steps = T,
                                  tumor_subcluster_partition_method = 'qnorm', 
                                  num_threads = 11, 
                                  debug = F)


