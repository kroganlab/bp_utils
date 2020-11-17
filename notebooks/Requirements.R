required.packages <- c('circlize', 'data.table', 'ggplot2', 'pbapply')
new.packages <- setdiff(required.packages, installed.packages())
if(length(new.packages) > 0)
  install.packages(new.packages)

required.bioconductor.packages <- c("ComplexHeatmap", "artMS", "MSstats")
new.packages <- setdiff(required.bioconductor.packages, installed.packages())
if(length(new.packages) > 0){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install(new.packages)
}

