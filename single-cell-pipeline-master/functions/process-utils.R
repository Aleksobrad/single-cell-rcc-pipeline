#' Filters out genes with no expression and low quality cells.
#' 
#' @param raw.mat Matrix of raw gene expression data (genes X samples).
#' @param minCount Minimum number of reads in a cell. Default of 1000.
#' @param maxCount Maximum number of reads in a cell. Default of 100000.
#' @param minGeneReads Minimum number of reads for a gene to be kept. Default of 1 (any gene with no reads will be removed).
#' @return Quality controlled matrix.
QCTransform <- function(raw.mat, minCount = 1000, maxCount = 100000, minGeneReads = 1) {
  filt.mat <- raw.mat[, colSums(raw.mat) > minCount & colSums(raw.mat) < maxCount]
  filt.mat <- filt.mat[ rowSums(raw.mat) >= minGeneReads ,]
  rem.genes <- nrow(raw.mat) - nrow(filt.mat); rem.cells <- ncol(raw.mat) - ncol(filt.mat)
  print(paste('Removed ', rem.genes, ' genes and ', rem.cells, ' cells.', sep =''))
  return(filt.mat)
}

#' Generates basic quality control plots from raw gene expression data.
#' 
#' @param raw.mat Matrix of raw gene expression data (genes X samples).
#' @param mt.genes Path to .csv file with ENSG and Hugo names for mitochondrial genes.
#' @param plot.path Optional argumetn of save path for plot.'
QCPlots <- function(raw.mat, mt.genes, plot.path) {
  # packages
  require(ggplot2)
  require(ggpubr)
  ## sequencing depth plot
  p1.dat <- data.frame('Depth' = colSums(raw.mat), 'Sample' = as.factor(rep('raw', ncol(raw.mat))))
  p1 <- ggplot(p1.dat, aes(x=Sample, y=Depth)) + geom_violin(color = '#F8766D', fill = '#F8766D') + 
    theme_bw()
  ## detected gene plot
  p2.dat <- data.frame('dgenes' = colSums(raw.mat > 0), 'Sample' = as.factor(rep('raw', ncol(raw.mat))))
  p2 <- ggplot(p2.dat, aes(x=Sample, y=dgenes)) + geom_violin(color = '#00BA38', fill = '#00BA38') + 
    ylab('Datected Genes') + theme_bw()
  ## mt percentage plot
  mt.perc <- MTPercent(raw.mat, mt.genes)
  p3.dat <- data.frame('mt' = mt.perc, 'Sample' = as.factor(rep('raw', length(mt.perc))))
  p3 <- ggplot(p3.dat, aes(x=Sample, y=mt)) + geom_violin(color = '#619CFF', fill = '#619CFF') +
    ylab('MT%') + theme_bw()
  ## arrange and plot
  if (!missing(plot.path)) {
    ggarrange(plotlist = list(p1, p2, p3), ncol = 3) %>% ggexport(filename = plot.path, height = 700, width = 1000)
  } else {
    ggarrange(plotlist = list(p1, p2, p3), ncol = 3)
  }
}

#' Returns vector of mitochondrial percentages for the given samples.
#'
#' @param dat.mat Matrix of raw gene expression data (genes X samples).
#' @param mt.genes List of mitochondrial genes
#' @retun Returns named vector of mitochondrial gene percentage.
MTPercent <- function(dat.mat, mt.genes) {
  mt.count <- colSums(dat.mat[ intersect(rownames(dat.mat), mt.genes) ,])
  total.count <- colSums(dat.mat)
  mt.percent <- mt.count / total.count
  head(mt.percent)
  return( mt.percent )
}

#' Filters data based on percentage of mitochondrial gens.
#'
#' @param raw.mat Matrix of raw gene expression data (genes X samples)
#' @param mt.genes Path to .csv file with ENSG and Hugo names for mitochondrial genes.
#' @param mt.thresh Threshold above which cells will be removed. Default of 0.15
MTFilter <- function(dat.mat, mt.genes, mt.thresh = 0.1) {
  ## find mt percentages
  mt.perc <- MTPercent(raw.mat, mt.genes)
  ## filter matrix
  thresh.cells <- names(mt.perc)[which(mt.perc < mt.thresh)]
  rem.cells <- ncol(dat.mat) - length(thresh.cells)
  print(paste('Removed', rem.cells, 'cells with too many MT reads', sep = ' '))
  dat.mat <- dat.mat[, thresh.cells ]
  return(dat.mat)
}

#' Performas a CPM normalization on the given data. 
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

#' Transforms gene names from Ensembl to hgnc.
#' 
#' @param dat.mat Matrix of data with ENSEMBL names (genes X samples).
#' @return Data with HGNC names. Some data will likely be lost in the conversion.
Ensemble2GeneName<-function(dat.mat) {
  # packages
  require(biomaRt)
  # get the mart
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://useast.ensembl.org")
  name.map <- getBM(attributes=c('hgnc_symbol','hgnc_id','ensembl_gene_id'), 
                       filters = 'ensembl_gene_id', values = (rownames(dat.mat)), mart = ensembl)
  # remove rows with no match for gene name, then match and replace
  name.map <- name.map[ which(!is.na(name.map$hgnc_symbol)) , ]
  name.map <- name.map[ which(name.map$hgnc_symbol != '') , ]
  convert.dat <- dat.mat[ name.map$ensembl_gene_id , ]
  rownames(convert.dat) <- name.map$hgnc_symbol
  return(convert.dat)
}

#' Transforms gene names from Ensembl to Entrez
#' 
#' @param dat.mat Matrix of data with ENSEMBL names (genes X samples).
#' @return Data with Entrez names. Some data will likely be lost in the conversion.
Ensemble2Entrez <- function(dat.mat) {
  # packages
  require(biomaRt)
  # get the mart
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://useast.ensembl.org")
  name.map <- getBM(attributes=c('entrezgene_id','ensembl_gene_id'), 
                    filters = 'ensembl_gene_id', values = (rownames(dat.mat)), mart = ensembl)
  # remove rows with no match for entrez, then match and replace
  name.map <- name.map[ which(!is.na(name.map$entrezgene)) , ]
  name.map <- name.map[ which(name.map$entrezgene != '') , ]
  convert.dat <- dat.mat[ name.map$ensembl_gene_id , ]
  rownames(convert.dat) <- name.map$entrezgene
  return(convert.dat)
}

#' Transforms gene names from Entrez to Ensembl
#' 
#' @param dat.mat Matrix of data with Entrez names (genes X samples).
#' @return Data with Ensembl names. Some data will likely be lost in the conversion.
Entrez2Ensemble <- function(dat.mat) {
  # packages
  require(biomaRt)
  # get the mart
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://useast.ensembl.org")
  name.map <- getBM(attributes=c('entrezgene_id','ensembl_gene_id'), 
                    filters = 'entrezgene_id', values = (rownames(dat.mat)), mart = ensembl)
  # remove rows with no match for entrez, then match and replace
  name.map <- name.map[ which(!is.na(name.map$ensembl_gene_id)) , ]
  name.map <- name.map[ which(name.map$ensembl_gene_id != '') , ]
  convert.dat <- dat.mat[ as.character( name.map$entrezgene ) , ]
  rownames(convert.dat) <- name.map$ensembl_gene_id
  return(convert.dat)
}

#' Saves a matrix in a format for input to ARACNe
#'
#' @param dat.mat Matrix of data (genes X samples).
#' @param out.file Output file where matrix will be saved.
#' @param subset Switch for subsetting the matrix to 500 samples. Default TRUE.
ARACNeTable <- function(dat.mat, out.file, subset = TRUE) {
  dat.mat <- dat.mat[!duplicated(rownames(dat.mat)), ]
  if (subset) {
    dat.mat <- dat.mat[, sample(colnames(dat.mat), min(ncol(dat.mat), 500)) ]
  }
  sample.names <- colnames(dat.mat)
  gene.ids <- rownames(dat.mat)
  m <- dat.mat
  mm <- rbind( c("gene", sample.names), cbind(gene.ids, m))
  write.table( x = mm , file = out.file , 
               sep="\t", quote = F , row.names = F , col.names = F )
}

#' Builds a UMAP transformation based on the shared, cell-by-cell master regulators. Utilizes custom UMAP parameters ('pearson', 25 neighbors).
#'
#' @param dat.mat Matrix of protein activity (protein X samples).
#' @param num.mrs Number of top master regulators to take from each cell. Default of 50.
#' @return Returns a UMAP based on the cell-by-cell master regulators.
cbcMR_UMAP <- function(dat.mat, num.mrs = 50) {
  require(umap)
  # set UMAP parameters
  umap_custom <- umap.defaults
  umap_custom$n_neighbors <- 25
  umap_custom$metric <- 'pearson'
  # identify MRs
  cbc.mrs <- apply(dat.mat, 2, function(x) { names(sort(x, decreasing = TRUE))[1:num.mrs] })
  cbc.mrs <- unique(unlist(as.list(cbc.mrs)))
  # generate and return UMAP
  cbc.umap <- umap(t(dat.mat[ match(cbc.mrs, rownames(dat.mat)) , ]), config = umap_custom, init = "random")
  return(cbc.umap)
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

#' Generates a UMAP based on a set of proteins.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param mrs List of proteins to use as master regulators.
#' @return UMAP object.
CustomUMAP <- function(dat.mat) {
  require(umap)
  # set UMAP parameters
  set.seed(1)
  umap_custom <- umap.defaults
  umap_custom$n_neighbors <- 25
  umap_custom$metric <- 'pearson'
  # compute umap
  c.umap <- umap(t(dat.mat), config = umap_custom, init = 'random')
  return(c.umap)
}

#' Read in data in 10x format
#'
#' @param dat.path Path to directory containing 'matrix.mtx', 'genes.tsv', and 'barcodes.tsv'.
#' @return Matrix of raw gene expression (genes X samples).
read10X <- function(dat.path) {
  require(Matrix)
  # read in data
  raw.mat <- as.matrix(readMM( paste(dat.path, 'matrix.mtx', sep = '') ))
  genes <- read.table( paste(dat.path, 'genes.tsv', sep = ''), sep = '\t')
  barcodes <- read.table( paste(dat.path, 'barcodes.tsv', sep = ''), sep = '\t')
  # set names and return
  colnames(raw.mat) <- barcodes[,1]
  rownames(raw.mat) <- genes[,1]
  return(raw.mat)
}