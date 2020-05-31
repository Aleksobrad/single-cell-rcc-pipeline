#' PAM clustering over a range of k values.
#' 
#' @param dist.mat A distance matrix for the data to be clustered.
#' @param kmin The minimum k to be used (default of 2).
#' @param kmax The maximum k to be used (default of 5).
#' @param verbose Switch to print status updates
#' @return A list of clustering objects for each value of k.
PamKRange <- function(dist.mat, kmin = 2, kmax = 5, verbose = TRUE) {
  ## packages
  require('cluster')  
  ## generate clusterings for each k
  clusterings <- list()
  for (k in kmin:kmax) {
    if (verbose) { print(paste('Clustering with k=', k, '...', sep = ''))}
    clustering <- pam(dist.mat, k)
    clusterings[[paste('k', k, sep='')]] <- clustering
  }
  return(clusterings)
}

#' Generation of silhouette score for a list of clusters.
#'
#' @param clusterings List of clustering objects.
#' @param dist.mat Distance matrix for the data clustered in clusterings.
#' @param plotPath If specified, will save a plot of the silhouette scores.
#' @return List of silhouette scores.
SilScoreEval <- function(clusterings, dist.mat, plotPath) {
  ## packages
  require('cluster')
  ## find silhouette scores for each cluster
  L <- length(clusterings)
  sil.scores <- c()
  k.vals <- c()
  for (i in 1:L) {
    # identify k for this clustering
    k <- length(table(clusterings[[i]]$clustering))
    k.vals <- c(k.vals, k)
    # find silhouette score
    sil <- silhouette(clusterings[[i]], dist.mat)
    sil.scores <- c(sil.scores, mean(sil[,3]))
  }
  ## plot if requested
  if (!missing(plotPath)) {
    require('ggplot2')
    plot.dat <- data.frame('k' = k.vals, 'Silhouette.Scores' = sil.scores)
    ggplot(plot.dat, aes(x=k, y=Silhouette.Scores)) + geom_point() + geom_line() +
      ggsave(plotPath, height=2, width=3)
  }
  ## return scores
  return(sil.scores)
}

#' Generates cluster-specific matrices for given data based on a clustering object.
#' 
#' @param dat.mat Data matrix to be split (features X samples).
#' @param clust Clustering object.
#' @param savePath If specified, matrices will be saved. Otherwise, a list of matrices will be returned.
#' @param savePref Preface for file names, if saving.
#' @param sizeThresh Smallest size cluster for which a matrix will be created. Default 300.
#' @return If files are NOT saved, returnes a list of matrices, one for each cluster. Otherwise, returns nothing.
ClusterMatrices <- function(dat.mat, clust, savePath, savePref, sizeThresh = 300) {
  ## set savePath if it is specified
  if(!missing(savePath)) {
    if (!missing(savePref)) {
      savePath <- paste(savePath, savePref, sep = '')
    }
  } else {
    clust.mats <- list()
  }
  ## generate matrices
  clust.table <- table(clust$clustering)
  for (i in 1:length(clust.table)) {
    if (clust.table[i] > sizeThresh) {
      clust.cells <- names(which(clust$clustering == names(clust.table)[i]))
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

#' Identifies a vector of colors for a given number of clusters. Internal function.
#'
#' @param k Number of clusters.
#' @param offset Optional argument to shift colors along color wheel. 
#' @return A vector of hues
ClusterColors <- function(k, offset = 0) {
  hues = seq(15, 375, length = k + 1) + offset
  return(hcl(h = hues, l = 65, c = 100)[1:k])
}

#' Genreates breaks for a color scale based on quantiles.
#'
#' @param dat.mat Data matrix (features X samples).
#' @param n Number of breaks to generate. Default of 10.
#' @return Numeric vector of break values.
QuantileBreaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(from = 0, to = 1, length.out = n))
  return(unique(breaks))
}

#' Generates a pheatmap from the given data and clustering object.
#' 
#' @param dat.mat Data matrix to be used (features X samples).
#' @param clust Vector of cluster labels.
#' @param plotTitle Title of the plot.
#' @param plotPath Optional argument for savign the plot, rather than displaying it.
#' @param quantile If TRUE, colors using quantile breaks; otherwise, uses linear breaks. Default of TRUE.
#' @return NULL
ClusterHeatmap <- function(dat.mat, clust, plotTitle, plotPath, quantile = TRUE) {
  require(pheatmap)
  require(RColorBrewer)
  ## generate sorting, plot data, mat breaks, and annotation color
  sorted.cells <- sort(clust)
  pheatmap.mat <- dat.mat[, names(sorted.cells)]
  num.breaks <- 100
  if (quantile) {
    mat.breaks <- QuantileBreaks(pheatmap.mat, num.breaks)
  } else {
    mat.breaks <- seq(min(dat.mat), max(dat.mat), (max(dat.mat) - min(dat.mat)) / num.breaks)
  }
  annot_color <- ClusterColors(length(unique(clust))); names(annot_color) <- sort(unique(clust))
  ## save if specified
  if (!missing(plotPath)) {
    jpeg(filename = plotPath) 
  }
  ## specify
  pheatmap(pheatmap.mat, annotation_col = data.frame('Cluster' = as.factor(clust)), 
           annotation_colors = list('Cluster' = annot_color),
           main = plotTitle, width = 6, height = 8, scale = 'row', fontsize_row = 4,
           cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
           color = colorRampPalette(rev(brewer.pal(11, 'RdBu')))(num.breaks), breaks = mat.breaks)
  if (!missing(plotPath)) {
    dev.off()
  }
}

#' Generates a ggplot2 scatter plot of the clusters in UMAP space.
#' 
#' @param umap.obj UMAP object, generated from the CustomUMAP function (or manually using the umap package).
#' @param clust Vector of cluster labels.
#' @param plotTitle Titel of the plot.
#' @return NULL
ClusterScatter <- function(umap.obj, clust, plotTitle) {
  require(ggplot2)
  plot.dat <- data.frame('UMAP1' = umap.obj$layout[,1], 'UMAP2' = umap.obj$layout[,2], 'Clusters' = clust)
  ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = Clusters), alpha=0.5, size = 2) +
    theme_bw() + theme(legend.justification = c("right", "top")) + ggtitle(plotTitle) +
    scale_color_manual(values = ClusterColors(length(unique(clust))))
}

#' Iterative clustering using PAM
#' 
#' @param dat.mat Matrix of data (features X samples).
#' @param dist.func Functiont to compute a distance matrix from the data.
#' @param iter.max Maximum number of iterations to use
#' @param sil.thresh Minimum cluster quality threshold (default of 0.25).
#' @return A data frame with samples in rows and iterations in columns, containing cluster labels for each iteration.
IterPAM <- function(dat.mat, dist.func, iter.max = 5, sil.thresh = 0.25) {
  require(cluster)
  require(psych)
  ## helper function for clustering
  IterClustHelper <- function(dist.mat) {
    clust <- PamKRange(dist.mat, kmin = 2, kmax = 5, verbose = FALSE)
    sil.scores <- unlist(lapply(clust, function(x) {x$silinfo$avg.width} ))
    best.sil <- which.max(sil.scores)
    opt.clust <- clust[[names(best.sil)]]
    return(opt.clust)
  }

  ## set up result matrices
  label.mat <- matrix(0L, ncol = 1, nrow = ncol(dat.mat))
  sil.mat <- matrix(0L, ncol = 1, nrow = ncol(dat.mat))
  rownames(label.mat) <- colnames(dat.mat); rownames(sil.mat) <- colnames(dat.mat)
  ## perform initial clustering
  print('Iteration 1...')
  dist.mat <- dist.func(dat.mat)
  i1.clust <- IterClustHelper(dist.mat)
  if (i1.clust$silinfo$avg.width < sil.thresh) {
    print('Data is homogeneous with given silhouette score threshold.')
    return()
  }
  label.mat[,1] <- i1.clust$clustering
  sil.mat[,1] <- i1.clust$silinfo$widths[,3]
  iter.num <- 1
  ## iterative clustering
  while (iter.num < iter.max) {
    print(paste('Iteration ', iter.num + 1, '...', sep = ''))
    ## prepare tracking variables
    k <- length(table(label.mat[,iter.num]))
    spacer <- 0
    homo.clusts <- 0
    ## prepare new vector for labels and silhouette scores
    label.vect <- rep(0, nrow(label.mat))
    sil.vect <- rep(0, nrow(sil.mat))
    ## sub cluster each cluster from previous iteration
    for (i in 1:k) {
      ## subset the data matrix and cluster
      clust.cells <- which(label.mat[, iter.num] == i)
      clust.mat <- dat.mat[, clust.cells]
      dist.mat <- dist.func(clust.mat)
      subClust <- IterClustHelper(dist.mat)
      ## check for homogeneity
      if (subClust$silinfo$avg.width < sil.thresh) {
        ## if there was no robust clustering, transfer previous labels down, with spacer
        label.vect[clust.cells] <- rep(spacer + 1, length(clust.cells))
        sil.vect[clust.cells] <- sil.mat[clust.cells, iter.num]
        spacer <- spacer + 1
        homo.clusts <- homo.clusts + 1
      } else {
        label.vect[clust.cells] <- subClust$clustering + spacer
        sil.vect[clust.cells] <- subClust$silinfo$widths[,3]
        spacer <- spacer + length(table(subClust$clustering))
      }
    }
    ## if all the clusterings in this iteration were homogeneous, return
    if (homo.clusts == k) {
      cNames <- paste0('iter', 1:iter.num)
      colnames(label.mat) <- cNames; colnames(sil.mat) <- cNames
      label.mat <- dfOrder(label.mat, 1:iter.num)
      sil.mat <- sil.mat[rownames(label.mat), ]
      return(list('clustering' = label.mat, 'silhouette' = sil.mat))
    } else { ## if there were new clusterings, update the matrices and iterate
      label.mat <- cbind(label.mat, label.vect)
      sil.mat <- cbind(sil.mat, sil.vect)
      iter.num <- iter.num + 1
    }
  }
  ## if max iterations were reached, return
  cNames <- paste0('iter', 1:iter.num)
  colnames(label.mat) <- cNames; colnames(sil.mat) <- cNames
  label.mat <- dfOrder(label.mat, 1:iter.num)
  sil.mat <- sil.mat[rownames(label.mat), ]
  return(list('clustering' = label.mat, 'silhouette' = sil.mat))
}

#' Louvain clustering using the MUDAN package.
#' 
#' @param dat.mat Matrix of data to be clustered (features X samples).
#' @param k Number of neighbors used in the community detection. Default of 100.
#' @return Clustering object.
LouvainClust <- function(dat.mat, k = 100) {
  require(MUDAN)
  set.seed(1)
  res <- MUDAN::getComMembership(t(dat.mat), k = k, method = igraph::cluster_infomap, verbose = FALSE) 
  return(res)
}

#' Creates a grid of plots for the given markers that exist in the given matrix.
#'
#' @param umap UMAP object generated from the 'UMAP' package
#' @param cluster Vector of cluster labels for each sample
#' @param pAct Matrix of protein activity (features X samples).
#' @param markers List of markers of interest. 
#' @param plotTitle Title for the plot.
#' @return NULL
MarkerGrid <- function(umap, cluster, pAct, markers, plotTitle) {
  require(ggplot2)
  require(ggpubr)
  # create plot data frame
  plot.dat <- data.frame('UMAP1' = umap$layout[,1], 'UMAP2' = umap$layout[,2],
                       'clust' = cluster)
  marker.mat <- t(pAct[ intersect(markers, rownames(pAct)) ,])
  plot.dat <- cbind(plot.dat, as.data.frame(marker.mat))
  # create plots
  clust.plot <- ggplot(plot.dat, aes(x = UMAP1, y = UMAP2, color = clust)) + geom_point() + ggtitle('Clusters')
  plot.list <- list('Clusters' = clust.plot)
  for (i in 4:ncol(plot.dat)) { # creat a plot for each marker and add to list
    m <- colnames(plot.dat)[i]
    m.plot <- ggplot(plot.dat, aes_string(x = 'UMAP1', y = 'UMAP2', color = m)) + geom_point() + ggtitle(m) + 
      scale_colour_gradientn(colours = c('blue', 'white', 'red'))
    plot.list[[m]] <- m.plot
  }
  # determine dimensions
  nCol <- min(3, length(plot.list))
  nRow <- ceiling(length(plot.list) / nCol)
  # arrange plots
  marker.plot <- ggarrange(plotlist = plot.list, ncol = nCol, nrow = nRow)
  print(annotate_figure(marker.plot, top = text_grob(plotTitle, size = 24)))
}
