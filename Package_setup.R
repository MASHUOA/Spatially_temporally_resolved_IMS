#t-SNE package folder

tsnefolder<-"~/FIt-SNE/"
source(paste0(tsnefolder,"fast_tsne.R"))
FAST_TSNE_SCRIPT_DIR <- tsnefolder

if (!require(Cardinal)){
  remotes::install_github("kuwisdelu/Cardinal",force=T)
}

if (!require(HiTMaP)){
  install.packages("remotes")
  remotes::install_github("MASHUOA/HiTMaP",force=T) 
}
#closeAllConnections()
register(SnowParam())
options(Cardinal.verbose=FALSE)
options(Cardinal.progress=T)
RNGkind("Mersenne-Twister")
parallel=try(detectCores()/2)
if (parallel<1 | is.null(parallel)){parallel=1}
BPPARAM=SnowParam()
bpprogressbar(BPPARAM)=TRUE
setCardinalBPPARAM(BPPARAM)

#Cardinal
install.packages("BiocManager")
BiocManager::install("Cardinal")

#MetaboAnalystR
metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}
metanr_packages()
install.packages("devtools")
library(devtools)
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)

#HiMAP
install.packages("remotes")
#library(devtools)
library(remotes)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github("MASHUOA/HiTMaP",force=T)
remotes::install_github("kuwisdelu/Cardinal",force=T)
3
no
#Update all dependencies
BiocManager::install(ask = F)
yes
library(HiTMaP)

install.packages("dplyr")
install.packages("stringr")
install.packages("readr")
install.packages("hexbin")
install.packages("magick")
install.packages("FactoMineR")
install.packages("M3C")
install.packages("plotly")
install.packages("rsvd")
install.packages("egg")
install.packages("reshape2")
install.packages("flextable")
install.packages("grid")
install.packages("wesanderson")
install.packages("ggpubr")