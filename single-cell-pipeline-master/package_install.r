## install packages from CRAN
install.packages('devtools')
install.packages('BiocManager')
install.packages('dplyr')
install.packages('Matrix')
install.packages('gdata')
install.packages('scater')
install.packages('knitr')
install.packages('kableExtra')
install.packages('reshape2')
install.packages('ggpubr')
install.packages('RcolorBrewer')
install.packages('umap')
install.packages('ggplot2')
install.packages('optparse')

## install packages from BioConductor
BiocManager::install("DropletUtils")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("viper")

## install packages from Github
devtools::install_github("JEFworks/MUDAN")