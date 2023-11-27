Annotation and visualisation of spatially and temporally resolved,
isotopically-labelled imaging mass spectrometry metabolomics data
================
Dingchang Shi, George Guo

- [1. Initialization settings and
  functions](#1-initialization-settings-and-functions)
  - [1.1 Raw data, processed spectrum, and
    result](#11-raw-data-processed-spectrum-and-result)
  - [1.2 Global functions (shown in the
    RMD)](#12-global-functions-shown-in-the-rmd)
  - [1.3 Package Setup](#13-package-setup)
- [2. Working directory setup](#2-working-directory-setup)
- [3. Lens coordinate rescaling and
  alignment](#3-lens-coordinate-rescaling-and-alignment)
- [4. Lens-to-lens hexagonal binning](#4-lens-to-lens-hexagonal-binning)
- [5. Time-scale transformation of the lens
  data](#5-time-scale-transformation-of-the-lens-data)
- [6. Dimensionality reduction, top loading feature selection and
  K-means
  segmentation](#6-dimensionality-reduction-top-loading-feature-selection-and-k-means-segmentation)
  - [6.1 Principle conponent analysis and top loading feature
    selection](#61-principle-conponent-analysis-and-top-loading-feature-selection)
  - [6.2 UMAP (all bio rep)](#62-umap-all-bio-rep)
  - [6.3 t-SNE (all bio rep)](#63-t-sne-all-bio-rep)
- [7. Metabolite library construction
  (function)](#7-metabolite-library-construction-function)
- [8. 1st round metabolomic annotation and pathway enrichment for all
  detected m/z
  features](#8-1st-round-metabolomic-annotation-and-pathway-enrichment-for-all-detected-mz-features)
  - [8.1 Metabolite annotation and pathway
    enrichment](#81-metabolite-annotation-and-pathway-enrichment)
  - [8.2 pathway summary](#82-pathway-summary)
- [9. Feature scoring and FDR controlled
  filtering](#9-feature-scoring-and-fdr-controlled-filtering)
- [10. 2nd round metabolomic annotation and pathway enrichment within
  FDR filtered
  features](#10-2nd-round-metabolomic-annotation-and-pathway-enrichment-within-fdr-filtered-features)
- [11. Feature visualization](#11-feature-visualization)
  - [11.1 Selected bio rep](#111-selected-bio-rep)
  - [11.2 Visualization](#112-visualization)
  - [11.3 Metabolite annotation in pathway
    format](#113-metabolite-annotation-in-pathway-format)
- [12. Shrunk data](#12-shrunk-data)
  - [12.1 PCA visualisation](#121-pca-visualisation)
  - [12.2 PCA K means](#122-pca-k-means)
  - [12.3 UMAP K means](#123-umap-k-means)
  - [12.4 t-SNE K means](#124-t-sne-k-means)

# 1. Initialization settings and functions

## 1.1 Raw data, processed spectrum, and result

Raw data link for steps 3 and 4:

[*Project on
metaspace*](https://metaspace2020.eu/api_auth/review?prj=54e09ade-8cb1-11ee-adab-831686c45448&token=pwemUgFOuNiB)

Result link for step 5 and later can be retrieved from the GitHub repo
by using the R code below:

Raw data and Result files for all steps can accessed temporarily via
dropbox for reviewing process: [Dropbox
transfer](https://www.dropbox.com/t/5zOCpQDGbzaevoZy)

## 1.2 Global functions (shown in the RMD)

``` r
post_mod <-
  function (fig, point_size) {
    for (n in 1:length(fig[["x"]][["data"]])) {
      fig[["x"]][["data"]][[n]][["marker"]]$size = point_size
    }
    return(fig)
  }

# m/z match with tolerance
`%~%` <- function(x, y) {
  x = as.numeric(x)
  y = as.numeric(y)
  sapply(x, function(.x, tol = 0.000012) {
    any(sapply(y, function(.y)
      isTRUE(all.equal(
        .x, .y, tolerance = tol
      ))))
  })
}

knitr::opts_chunk$set(
  echo = TRUE,
  fig.height = 12,
  fig.width = 12,
  message = FALSE,
  warning = FALSE
)

#Adopted from metaboanalyst to make the local rsever available

init.global.vars <- function(anal.type) {
  # other global variables
  msg.vec <<- ""
  err.vec <<- ""
  # for network analysis
  module.count <<- 0
  # counter for naming different json file (pathway viewer)
  smpdbpw.count <<- 0
  mdata.all <<- list()
  mdata.siggenes <<- vector("list")
  meta.selected <<- TRUE
  anal.type <<- anal.type
  if (.on.public.web) {
    # disable parallel prcessing for public server
    library(BiocParallel)
    register(SerialParam())
  } else {
    if ("stat" %in% anal.type |
        "msetqea" %in% anal.type |
        "pathqea" %in% anal.type |
        "roc" %in% anal.type)
      # start Rserve engine for Rpackage
      load_Rserve()
 }
  
  # plotting required by all
  Cairo::CairoFonts(
    regular = "Arial:style=Regular",
    bold = "Arial:style=Bold",
    italic = "Arial:style=Italic",
    bolditalic = "Arial:style=Bold Italic",
    symbol = "Symbol"
  )
  
  # sqlite db path for gene annotation
  if (file.exists("/home/glassfish/sqlite/")) {
    #.on.public.web
    url.pre <<- "/home/glassfish/sqlite/"
    
  } else if (file.exists("/home/jasmine/Downloads/sqlite/")) {
    #jasmine's local
    url.pre <<- "/home/jasmine/Downloads/sqlite/"
    
    #api.base <<- "localhost:8686"
  } else if (file.exists("/Users/soufanom/Documents/Projects/gene-id-mapping/")) {
    # soufan laptop
    url.pre <<-
      "/Users/soufanom/Documents/Projects/gene-id-mapping/"
    
  } else if (file.exists("~/Documents/Projects/gene-id-mapping/")) {
    url.pre <<- "~/Documents/Projects/gene-id-mapping/"
  } else if (file.exists("/Users/xia/Dropbox/sqlite/")) {
    # xia local
    url.pre <<- "/Users/xia/Dropbox/sqlite/"
    
  } else if (file.exists("/media/zzggyy/disk/sqlite/")) {
    url.pre <<- "/media/zzggyy/disk/sqlite/"
    #zgy local)
  } else if (file.exists("/home/zgy/sqlite/")) {
    url.pre <<- "/home/zgy/sqlite/"
    #zgy local)
  } else if (file.exists("/home/le/sqlite/")) {
    # le local
    url.pre <<- "/home/le/sqlite/"
    
  } else if (file.exists("/home/qiang/Music/")) {
    # qiang local
    url.pre <<- "/home/qiang/sqlite/"
    
  } else{
    url.pre <<-
      paste0(dirname(
        system.file(
          "database",
          "sqlite/GeneID_25Species_JE/ath_genes.sqlite",
          package = "MetaboAnalystR"
        )
      ), "/")
  }
  
  api.base <<- "http://api.xialab.ca"
  #api.base <<- "132.216.38.6:8987"
  
}

load_Rserve <- function() {
  installed <- c("Rserve") %in% rownames(installed.packages())
  
  if (installed) {
    # first need to start up an Rserve instance
    suppressMessages(library(Rserve))
    Rserve::Rserve(args = "--no-save")
  } else{
    print("Please install Rserve R package!")
  }
}

.on.public.web <<- F

# load pixel level data for all runs
load_pixel_label <-
  function(combinedimdata,
           datafile_base,
           workdir,
           coordata_file = "coordata.csv",
           pixel_idx_col = base::row.names,
           label_col = "pattern",
           ...) {
    library(Cardinal)
    library(stringr)
    library(HiTMaP)
    library(stringr)
    library(dplyr)
    datafile <- str_remove(datafile_base, "\\.imzML$")
    if (length(workdir) != length(datafile)) {
      workdir = rep(workdir[1], length(datafile))
    }
    
    datafile_imzML = datafile
    coordata_file_tb <- NULL
    for (z in 1:length(datafile)) {
      name <- basename(datafile[z])
      name <- gsub(".imzML$", "", name)
      name <- gsub("/$", "", name)
      setwd(workdir[z])
      folder <- base::dirname(datafile[z])
      #imdata <- Cardinal::readImzML(datafile[z],preprocessing = F,attach.only = T,resolution = 200,rotate = rotate[z],as="MSImageSet",BPPARAM = BPPARAM)
      if (!str_detect(datafile[z], ".imzML$")) {
        datafile_imzML[z] <- paste0(datafile[z], ".imzML")
      }
      
      coordata_file_tb[[datafile[z]]] <-
        read.csv(paste0(workdir[z], "/", datafile[z], " ID/", coordata_file))
      coordata_file_tb[[datafile[z]]]$run <- datafile[z]
      coordata_file_tb[[datafile[z]]]$pixel_idx <-
        pixel_idx_col(coordata_file_tb[[datafile[z]]])
    }
    
    coordata_file_bind <- do.call(rbind, coordata_file_tb)
    coordata_file_bind$run <-
      as.factor(tolower(coordata_file_bind$run))
    coordata_file_bind$pixel_idx <-
      as.numeric(coordata_file_bind$pixel_idx)
    Pixel_run <- run(combinedimdata)
    Pixel_run <- data.frame(run = Pixel_run)
    Pixel_run <-
      Pixel_run %>% group_by(run) %>% summarise(pixel_idx = 1:length(run))
    
    label_run <-
      merge(
        Pixel_run,
        coordata_file_bind,
        by = c("run", "pixel_idx"),
        all.x = T,
        sort = F
      )
    
    Pixel_label <- label_run[[label_col]]
    return(Pixel_label)
  }

match_mz = function(x,
                    ref,
                    tol = 1e-10,
                    unit = c("Da", "ppm")) {
  sref = sort(ref)
  
  i = findInterval(x, sref, all.inside = TRUE)
  
  dif1 = abs(x - sref[i])
  
  dif2 = abs(x - sref[i + 1])
  
  dif = dif1 > dif2
  
  dif1[dif] = dif2[dif]
  
  i[dif] = i[dif] + 1
  
  res_indx <- rep(NA, length(x))
  
  if (unit[1] == "ppm") {
    tol_ppm <- tol
    tol <- tol_ppm * x / 1000000
  }
  res_indx[dif1 <= tol] <- sref[i[dif1 <= tol]]
  res_indx
}

# 13C metabolites library modification and adoption
mummichog.lib.mod <-
  function(lib = "bta_kegg",
           lib.new = "bta_kegg_13C",
           C13_number = c(3, 6),
           C13_deltamass = 1.003355,
           wd = getwd(),
           force_update_lib = F,
           adducts_list = "all",
           method = c("new_cpd", "new_adduct")) {
    filenm <- paste(wd, "/", lib, ".qs", sep = "")
    mum.url <-
      paste(
        "https://www.metaboanalyst.ca/resources/libs/mummichog/",
        paste(lib, ".qs", sep = ""),
        sep = ""
      )
    if (`|`(force_update_lib, !file.exists(filenm))) {
      download.file(mum.url,
                    destfile = filenm,
                    method = "libcurl",
                    mode = "wb")
    }
    
    mummichog.lib <- qs::qread(filenm)
    cmpd.map <- MetaboAnalystR:::.get.my.lib("compound_db.qs")
    hit.inx <-
      match(tolower(mummichog.lib$cpd.lib$id), tolower(cmpd.map$kegg))
    nms <- cmpd.map[hit.inx, "smiles"]
    c_count <- stringr::str_count(nms, "C|c")
    nms_pub <- cmpd.map[hit.inx, "pubchem_id"]
    
    df <-
      data.frame(
        smile = nms,
        c_count,
        mw_c_ratio = mummichog.lib$cpd.lib$mw / c_count,
        id = mummichog.lib$cpd.lib$id,
        mw = mummichog.lib$cpd.lib$mw,
        name = mummichog.lib$cpd.lib$name
      )
    
    mummichog.lib.mod <- NULL
    mummichog.lib.new <- mummichog.lib
    
    if (method[1] == "new_cpd") {
      for (cnum in C13_number) {
        (df$c_count >= cnum) -> mod.item
        (df$mw >= 26.09864 * as.numeric(cnum)) -> mod.item2
        mod.item[is.na(mod.item)] <- F
        mod.item2[is.na(mod.item2)] <- F
        mod.item_final <-
          as.numeric(ifelse(mod.item, mod.item, mod.item2))
        mummichog.lib.mod[[cnum]] <- mummichog.lib
        lapply(mummichog.lib.mod[[cnum]][["pathways"]][["cpds"]], function(x, cnum) {
          paste0(x, "_13C", cnum)
        }, cnum) -> mummichog.lib.mod[[cnum]][["pathways"]][["cpds"]]
        mummichog.lib.mod[[cnum]][["cpd.lib"]][["id"]] <-
          paste0(mummichog.lib.mod[[cnum]][["cpd.lib"]][["id"]], "_13C", cnum)
        mummichog.lib.mod[[cnum]][["cpd.lib"]][["name"]] <-
          paste0(mummichog.lib.mod[[cnum]][["cpd.lib"]][["name"]], "_13C", cnum)
        mummichog.lib.mod[[cnum]][["cpd.lib"]][["mw"]] <-
          mummichog.lib.mod[[cnum]][["cpd.lib"]][["mw"]] + C13_deltamass * as.numeric(cnum) *
          mod.item_final
        mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["dpj_positive"]] <-
          mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["dpj_positive"]] + C13_deltamass *
          as.numeric(cnum) * mod.item_final
        mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["positive"]] <-
          mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["positive"]] + C13_deltamass *
          as.numeric(cnum) * mod.item_final
        mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["negative"]] <-
          mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["negative"]] + C13_deltamass *
          as.numeric(cnum) * mod.item_final
        
      }
      for (cnum in C13_number) {
        lapply(1:length(mummichog.lib.new$pathways$cpds), function(x) {
          c(mummichog.lib.new$pathways$cpds[[x]],
            mummichog.lib.mod[[cnum]]$pathways$cpds[[x]])
        }) -> mummichog.lib.new$pathways$cpds
        
        mummichog.lib.new$cpd.lib$id <-
          c(mummichog.lib.new$cpd.lib$id,
            mummichog.lib.mod[[cnum]]$cpd.lib$id)
        mummichog.lib.new$cpd.lib$name <-
          c(mummichog.lib.new$cpd.lib$name,
            mummichog.lib.mod[[cnum]]$cpd.lib$name)
        mummichog.lib.new$cpd.lib$mw <-
          c(mummichog.lib.new$cpd.lib$mw,
            mummichog.lib.mod[[cnum]]$cpd.lib$mw)
        mummichog.lib.new[["cpd.lib"]][["adducts"]]$dpj_positive <-
          rbind(mummichog.lib.new[["cpd.lib"]][["adducts"]]$dpj_positive,
                mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]]$dpj_positive)
        mummichog.lib.new[["cpd.lib"]][["adducts"]]$positive <-
          rbind(mummichog.lib.new[["cpd.lib"]][["adducts"]]$positive,
                mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]]$positive)
        mummichog.lib.new[["cpd.lib"]][["adducts"]]$negative <-
          rbind(mummichog.lib.new[["cpd.lib"]][["adducts"]]$negative,
                mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]]$negative)
        
      }
    }
    
    cpd.lib <- mummichog.lib.new$cpd.lib
    
    ms_modes <- c('dpj_positive', 'positive', 'negative')
    
    adducts <- list()
    
    for (ms_mode in ms_modes) {
      if (adducts_list[1] == "all") {
        adducts[[ms_mode]] <-
          MetaboAnalystR:::Compound_function_mzlist(ms_mode, cpd.lib$mw)
      } else{
        resdf <-
          MetaboAnalystR:::Compound_function_mzlist(ms_mode, cpd.lib$mw)
        if (sum(adducts_list %in% colnames(resdf)) != 0) {
          adducts[[ms_mode]] <-
            resdf[, adducts_list[adducts_list %in% colnames(resdf)]]
        } else{
          adducts[[ms_mode]] <- resdf
        }
      }
    }
    
    if (method[1] == "new_adduct") {
      ms_modes <- c('dpj_positive', 'positive', 'negative')
      for (ms_mode in ms_modes) {
        new_adduct <- NULL
        for (cnum in (C13_number)) {
          new_adduct[[cnum]] <-
            adducts[[ms_mode]] + C13_deltamass * as.numeric(cnum)
          colnames(new_adduct[[cnum]]) <-
            stringr::str_replace(colnames(new_adduct[[cnum]]),
                                 "^M",
                                 paste0("[M_13C", cnum, "]"))
        }
        new_adduct <- do.call(cbind, new_adduct)
        adducts[[ms_mode]] <- cbind(adducts[[ms_mode]], new_adduct)
      }
    }
    cpd.lib$adducts <- adducts
    
    cpd.tree <- list()
    
    for (ms_mode in ms_modes) {
      l2 <- list()
      
      l2[[49]] <- ""
      
      l2[[2001]] <- ""
      
      mz.mat <- cpd.lib$adducts[[ms_mode]]
      
      floor.mzs <- floor(mz.mat)
      
      for (i in 1:nrow(floor.mzs)) {
        neighbourhood <- floor.mzs[i, ]
        
        for (n in neighbourhood) {
          if ((n > 50) & (n < 2000)) {
            l2[[n]] <- append(l2[[n]], i)
            
          }
        }
      }
      cpd.tree[[ms_mode]] <- lapply(l2, unique)
      
    }
    # set up the variables
    mummichog.lib.new <- list(
      pathways = mummichog.lib.new$pathways,
      cpd.tree = cpd.tree,
      cpd.lib = cpd.lib
    )
    qs::qsave(mummichog.lib.new, paste0(wd, "/", lib.new, ".qs"))
    return(mummichog.lib.new)
  }


# annotate and summarize the pathway mapping result
annotate_PSEA <- function(mSet, mummichog.lib.new) {
  lapply(mSet[["path.hits"]], function(x) {
    data.frame(
      cpd_id = paste0(x, collapse = ";"),
      cpd_names = paste0(mummichog.lib.new$cpd.lib$name[match(x, mummichog.lib.new$cpd.lib$id)], collapse = ";")
    )
  }) -> res
  res <- do.call(rbind, res)
  cbind(mSet[["mummi.resmat"]], res) -> resdf
  return(resdf)
}
```

## 1.3 Package Setup

The t-SNE on large scale data requires
<https://github.com/KlugerLab/FIt-SNE> to be installed. FFTW
(<http://www.fftw.org/>) is also required.

Please follow the instruction from the repo to get “fftRtsne” function
available and working properly.

For detailed installation guide of MetaboAnalystR please visit:
<https://github.com/xia-lab/MetaboAnalystR>

``` r
tsnefolder<-"~/FIt-SNE/"
source(paste0(tsnefolder,"fast_tsne.R"))

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
```

# 2. Working directory setup

``` r
library(MetaboAnalystR)
library(Cardinal)

#change to your own project directory here
project_dir <- "~/Manuscript/test/"

large_files_dir <- "/large_files/"
spatial_alignment_dir <- "/spatial_alignment/"
dim_reduction_dir <-  "/dim_reduction/"
PCA_dir <- "/pca_anal_new/"
UMAP_dir <- "/UMAP/"
tSNE_dir <- "/tSNE/"

FDR_dir <- "/SQRTP_0.1/"
visualization_dir <- "/visualization/"
shrink_dir<- "/data_shrink/"
heximg_dir <- "/heximg/"

dir.create(paste0(project_dir, spatial_alignment_dir))
dir.create(paste0(project_dir, FDR_dir))
dir.create(paste0(project_dir, dim_reduction_dir))
dir.create(paste0(project_dir,visualization_dir))

#data infromation 
fileinfo_v2 <- readr::read_csv(paste0(project_dir,"/fileinfo_v2.csv"))
```

# 3. Lens coordinate rescaling and alignment

``` r
wd <- project_dir
outwd<- paste0(project_dir, spatial_alignment_dir)

#load IMS data of full m/z features in all runs
combinedimdata <- readRDS(paste0(project_dir, large_files_dir,"/combinedimdata_norm_full.rds"))
example_imdata<-combinedimdata[,run(combinedimdata) == "4hrsilglucose_2nd_20201109_trim"]
library(Cardinal)
library(readr)
library(hexbin)
library(dplyr)

fileinfo_v2 <- readr::read_csv(paste0(wd,"/fileinfo_v2.csv"))
fileinfo_v2$filenames<-tolower(fileinfo_v2$filenames)
fileinfo_v2<-fileinfo_v2[fileinfo_v2$filenames %in% levels(run(combinedimdata)),]
coordata<-combinedimdata@elementMetadata@coord@listData
coordata$run<-run(combinedimdata)
coordata<-as.data.frame(coordata)

pixeldf<-coordata
pixeldf$x->pixeldf$x_scale
pixeldf$y->pixeldf$y_scale

#lens rescaling 
for (imsrun in unique(pixeldf$run)){
  pixeldf$x[pixeldf$run == imsrun]->x
  pixeldf$y[pixeldf$run == imsrun]->y
  x_scale<-(x-(max(x)-min(x))/2)/(max(y)-min(y))
  y_scale<-(y-(max(y)-min(y))/2)/(max(y)-min(y))
  pixeldf$x_scale[pixeldf$run == imsrun]<-x_scale 
  pixeldf$y_scale[pixeldf$run == imsrun]<-y_scale
}

pixeldf %>% group_by(run) %>% summarize(xmin=range(x_scale)[1],ymin=range(y_scale)[1])->pixeldf_sum


bin<-hexbin(pixeldf$x_scale,pixeldf$y_scale, xbins=round(sqrt(nrow(pixeldf)/7/length(unique(pixeldf$run)))), IDs= TRUE)
bin@cID->pixeldf$new_label_scale
# write.csv(pixeldf,paste0(wd,"all_coordata_label_sel_arc_bind_scale.csv"))
```

# 4. Lens-to-lens hexagonal binning

``` r
wd <- project_dir
outwd<- paste0(project_dir, spatial_alignment_dir)

combinedimdata->combinedimdata_full
final_selection_df<-read.csv(paste0(wd,"/final_selection.csv"))

fileinfo_v2$filenames<-tolower(fileinfo_v2$filenames)
fileinfo_v2<-fileinfo_v2[fileinfo_v2$filenames %in% levels(run(combinedimdata_full)),]


coordata<-combinedimdata_full@elementMetadata@coord@listData
coordata$run<-run(combinedimdata_full)
coordata<-as.data.frame(coordata)
#hexbin all features

library(hexbin)
label_sel="all"
merged_data_df<-NULL
coordata_label_sel_arc<-NULL

    setCardinalBPPARAM(BPPARAM=SerialParam())
    setCardinalVerbose(verbose=F)
    for (imsrun in fileinfo_v2$filenames){
      message(imsrun)
      time_cls<-fileinfo_v2$class[match(imsrun,fileinfo_v2$filenames)]
      combinedimdata_full@elementMetadata@coord@listData->coordata
      coordata$run<-run(combinedimdata_full)
      coordata<-as.data.frame(coordata)
      coordata_run<-coordata[coordata$run==imsrun,]
      coordata_run<-pixeldf[pixeldf$run==imsrun,]
      rownames(coordata_run)<-1:nrow(coordata_run)
      library(hexbin)
      combinedimdata_full[,run(combinedimdata_full)==imsrun]->combinedimdata_full_sample_bin
    
    bincell_df<-NULL
    for (bincell in unique(coordata_run$new_label_scale)){
      combinedimdata_full_sample_bin[,which(coordata_run$new_label_scale==bincell)] %>% summarizeFeatures( FUN = "mean") %>% process()->meanspec
      meanspec@featureData@listData[["mean"]]->bincell_df[[bincell]]
    }
    
    bincell_df_bind<-do.call(cbind,bincell_df)
    colnames(bincell_df_bind)<-paste(imsrun,time_cls,unique(coordata_run$new_label_scale),sep = "@")
    merged_data_df[[imsrun]]<-bincell_df_bind
    coordata_label_sel_arc[[imsrun]]<-coordata_run
    }

coordata_label_sel_arc_bind<-do.call(rbind,coordata_label_sel_arc)
merged_data_df_bind<-do.call(cbind,merged_data_df)
rownames(merged_data_df_bind)<-combinedimdata_full@featureData@mz


#Generate hexbin data
write.csv(merged_data_df_bind,paste0(outwd,"/",label_sel,"_merged_data_df_bind2",".csv"))

#Generate Master Table 
write.csv(coordata_label_sel_arc_bind,paste0(outwd,"/",label_sel,"_coordata_label_sel_arc_bind2",".csv"))


message(label_sel," DONE")
```

# 5. Time-scale transformation of the lens data

``` r
wd<-paste0(project_dir, spatial_alignment_dir)
setwd(wd)
library(readr)
library(stringr)

all_data_original<-read_csv(paste0(wd,"/all_merged_data_df_bind2.csv"))
all_data_original<-as.data.frame(all_data_original)
rownames(all_data_original)<-all_data_original[,1]
all_data_original[,-1]->all_data_original
col.nms.split<-str_split_fixed(colnames(all_data_original),"@",3)
all_data_original_t<-t(all_data_original)

group_sum <- rowsum(all_data_original_t, group = interaction (col.nms.split[,2],col.nms.split[,3]))
group_sum_t<-as.data.frame(t(group_sum))

a = paste(col.nms.split[,2], col.nms.split[,3], sep=".")
for (col in colnames(group_sum_t)) { 
  group_sum_t [,col]=  group_sum_t [,col] / sum(a==col)
  # message(col,sum(a==col))
}

grouped_bio_rep_all_pixel_mean<-group_sum_t
grouped_bio_rep_all_pixel_mean<-as.data.frame(grouped_bio_rep_all_pixel_mean)


col.nms.split_2<-str_split_fixed(colnames(grouped_bio_rep_all_pixel_mean),"\\.",2)
Label<-col.nms.split_2[,1]
grouped_bio_rep_all_pixel_mean<-rbind(Label, grouped_bio_rep_all_pixel_mean)
row.names(grouped_bio_rep_all_pixel_mean)[1]<-"Label"
grouped_bio_rep_all_pixel_mean<-grouped_bio_rep_all_pixel_mean[,order(grouped_bio_rep_all_pixel_mean[1,])]


cls_merge<-list()
for (cls in unique(Label)){
    group_sum_t[,col.nms.split_2[,1]==cls]->temp
    colnames(temp)<-str_remove(colnames(temp),paste0(cls,"\\."))
    row.names(temp)<-paste0(row.names(temp),"@",cls)
  cls_merge[[cls]]<-(temp)
  
}

library(plyr)
cls_merge_bind<-rbind.fill(cls_merge)
rownms<-unlist(lapply(cls_merge,rownames))
rownames(cls_merge_bind)<-rownms

write.csv(cls_merge_bind, paste0(wd, "grouped_bio_rep_all_pixel_mean_label_time_merge.csv"))
```

# 6. Dimensionality reduction, top loading feature selection and K-means segmentation

check if the file needed exists in the right folder; otherwise download
from the repository

``` r
file.exists(paste0(project_dir, spatial_alignment_dir, "grouped_bio_rep_all_pixel_mean_label_time_merge.csv"))
```

## 6.1 Principle conponent analysis and top loading feature selection

``` r
library(dplyr)
library(stringr)
library(FactoMineR)

#read reshaped data set
cls_merge_bind<-readr::read_csv(paste0(project_dir, spatial_alignment_dir, "grouped_bio_rep_all_pixel_mean_label_time_merge.csv"))
cls_merge_bind<-as.data.frame(cls_merge_bind)
rownames(cls_merge_bind)<-cls_merge_bind[,1]
cls_merge_bind[,1]<-NULL
cls_merge_mt<-t(as.matrix(cls_merge_bind))
cls_merge_mt[is.na(cls_merge_mt)]<-0
cls_merge_mt_rowsum<-rowsum(cls_merge_mt,rep(1,nrow(cls_merge_mt)))
cls_merge_mt_non_zero<-cls_merge_mt[,cls_merge_mt_rowsum!=0]


PCA<-FactoMineR::PCA(cls_merge_mt_non_zero, graph = TRUE, ncp = 5)

PCA.contri<-data.frame(PCA[["var"]][["contrib"]])
PCA.contri.feature <-stringr::str_split_fixed(rownames(PCA.contri),"@", 2) 
colnames(PCA.contri.feature )<- c("mz","class")
PCA.contri<-cbind(PCA.contri.feature, PCA.contri)

colnames(PCA.contri)[!(colnames(PCA.contri) %in% c("mz", "class"))]->dimname
dir.create(paste0(wd,PCA_dir))

#Top loading feature selection 

Dim.feature.contri.list<- list()
Dim.feature.contri.summmary.list<- list()
for (dim in dimname){
  PCA.contri[,c("mz", "class",dim)]->dim_df
  colnames(dim_df)<-c("mz", "class","Contribution")
  dim_df$logcontrib<-log(dim_df$Contribution)
  dim_df$mz<-as.numeric(dim_df$mz)
  dim_df  %>% mutate(new_bin = ntile(Contribution, n=30))->dim_df
  dim_df$Dim=dim
    write.csv(dim_df[dim_df$new_bin %in% c(30),],paste0(project_dir, dim_reduction_dir, PCA_dir,dim,"_anno_freq_top1bin",".csv"),row.names = F)
  
  feature_slt<- unique(dim_df[dim_df$new_bin %in% c(30),"mz"])
  Dim.feature.contri.list[[dim]]<-feature_slt
  Dim.feature.contri.summmary.list[[dim]]<-dim_df
}
Dim.feature.contri.summmary.list.bind<-do.call(rbind,Dim.feature.contri.summmary.list)

Dim.feature.contri.summmary.list.bind %>% group_by(mz, Dim) %>% mutate(Rank = rank(Contribution)) -> Dim.feature.contri.summmary.list.bind

save(list = c("PCA","PCA.contri","PCA.desc","Dim.feature.contri.list","cls_merge_mt","Dim.feature.contri.summmary.list.bind"), file = paste0(project_dir, dim_reduction_dir, PCA_dir,"Time_merged_PCA_all_bio.rda"))
```

### 6.1.1 PCA Single dimension K means

``` r
library(dplyr)
library(Cardinal)

load(paste0(project_dir, dim_reduction_dir, PCA_dir,"Time_merged_PCA_all_bio.rda"))
dims.pca <- PCA[["ind"]][["coord"]]
dims.pca<- as.data.frame(dims.pca)

all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir,"all_coordata_label_sel_arc_bind2.csv"),row.names = 1)

example_imdata<-combinedimdata[,run(combinedimdata) == "4hrsilglucose_2nd_20201109_trim"]

dir.create(paste0(project_dir, dim_reduction_dir, PCA_dir,"/kmean_seg_single_demo/"))

  
for (a in 1:5) {

  dims.pca<-dims.pca[as.numeric(rownames(dims.pca)) %in% all_coordata_label_sel_arc_bind$new_label_scale[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim"],]

dims.pca$Rank<-dense_rank(dims.pca[,a])




png(paste0(project_dir, dim_reduction_dir, PCA_dir,"/kmean_seg_single_demo/","Dim.",a,"_pca_single_dim_lens_projections.png"),height = 5,width = 6,units = "in",res = 600)



match(all_coordata_label_sel_arc_bind$new_label_scale, as.numeric(rownames(dims.pca)))-> rank_index

all_coordata_label_sel_arc_bind$Rank<- dims.pca$Rank[rank_index]

factor(all_coordata_label_sel_arc_bind$Rank[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim"])->final_selection


fig.arbitrary_seg<-image(example_imdata, factor(final_selection) ~ x * y,superpose=F, key=F, xlab= NULL, col = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(max(dims.pca$Rank)) )



fig.arbitrary_seg[["par"]][["ann"]]=F
fig.arbitrary_seg[["par"]][["bty"]]="n"
fig.arbitrary_seg[["par"]][["pty"]]="s"
fig.arbitrary_seg[["par"]][["xaxt"]]="n"
fig.arbitrary_seg[["par"]][["yaxt"]]="n"
fig.arbitrary_seg[["par"]][["fg"]]="white"
fig.arbitrary_seg[["par"]][["oma"]]=c(0, 0, 0, 0)

print(fig.arbitrary_seg)
dev.off()
}
```

## 6.2 UMAP (all bio rep)

``` r
library(readr)
library(Cardinal)
library(rsvd)
library(dplyr)
library(plotly)
library(M3C)

dir.create(paste0(project_dir, dim_reduction_dir, UMAP_dir))


all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir,"all_coordata_label_sel_arc_bind2.csv"),row.names = 1)

cls_merge_bind<-readr::read_csv(paste0(project_dir, spatial_alignment_dir,"/grouped_bio_rep_all_pixel_mean_label_time_merge.csv"))
cls_merge_bind<-as.data.frame(cls_merge_bind)
rownames(cls_merge_bind)<-cls_merge_bind[,1]
cls_merge_bind[,1]<-NULL
cls_merge_mt<-(as.matrix(cls_merge_bind))
cls_merge_mt[is.na(cls_merge_mt)]<-0
cls_merge_mt==0->cls_merge_mt_blk_indx
rowSums(cls_merge_mt_blk_indx)->fblk_count
fblk_count_thres<-fblk_count>=max(fblk_count)*0.90
cls_merge_mt_test<-t(cls_merge_mt[!fblk_count_thres,])


##Time series umap
features.umap = M3C::umap((cls_merge_mt)) 
Time_merged_umap<-data.frame(features.umap$data$X1,features.umap$data$X2,rownames(cls_merge_mt_test))
colnames(Time_merged_umap)<-c("x","y","bin")
saveRDS(Time_merged_umap, paste0(project_dir, dim_reduction_dir, UMAP_dir, "/Time_merged_umap_all_bio.rds"))
```

### 6.2.1 UMAP Single dimension lens projection

``` r
library(dplyr)
library(Cardinal)

Time_merged_umap <- readRDS(paste0(project_dir, dim_reduction_dir, UMAP_dir,"Time_merged_umap_all_bio.rds"))


all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir,"all_coordata_label_sel_arc_bind2.csv"),row.names = 1)

example_imdata<-combinedimdata[,run(combinedimdata) == "4hrsilglucose_2nd_20201109_trim"]

dir.create(paste0(project_dir, dim_reduction_dir, UMAP_dir,"/kmean_seg_single_demo/"))

  
for (a in c("x", "y")) {

 
Time_merged_umap <-Time_merged_umap[as.numeric(Time_merged_umap$bin) %in% all_coordata_label_sel_arc_bind$new_label_scale[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim"],]

Time_merged_umap$Rank<-dense_rank(Time_merged_umap[,a])




png(paste0(project_dir, dim_reduction_dir, UMAP_dir,"/kmean_seg_single_demo/","umap_",a,"_single_dim_lens_projections.png"),height = 5,width = 6,units = "in",res = 600)



match(all_coordata_label_sel_arc_bind$new_label_scale, as.numeric(Time_merged_umap$bin))-> rank_index

all_coordata_label_sel_arc_bind$Rank<- Time_merged_umap$Rank[rank_index]

factor(all_coordata_label_sel_arc_bind$Rank[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim"])->final_selection


fig.arbitrary_seg<-image(example_imdata, factor(final_selection) ~ x * y,superpose=F, key=F, xlab= NULL, col = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(max(Time_merged_umap$Rank)) )



fig.arbitrary_seg[["par"]][["ann"]]=F
fig.arbitrary_seg[["par"]][["bty"]]="n"
fig.arbitrary_seg[["par"]][["pty"]]="s"
fig.arbitrary_seg[["par"]][["xaxt"]]="n"
fig.arbitrary_seg[["par"]][["yaxt"]]="n"
fig.arbitrary_seg[["par"]][["fg"]]="white"
fig.arbitrary_seg[["par"]][["oma"]]=c(0, 0, 0, 0)

print(fig.arbitrary_seg)
dev.off()
}
```

## 6.3 t-SNE (all bio rep)

``` r
dir.create(paste0(project_dir, dim_reduction_dir, tSNE_dir))


all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir,"all_coordata_label_sel_arc_bind2.csv"),row.names = 1)

cls_merge_bind<-readr::read_csv(paste0(project_dir, spatial_alignment_dir,"/grouped_bio_rep_all_pixel_mean_label_time_merge.csv"))
cls_merge_bind<-as.data.frame(cls_merge_bind)
rownames(cls_merge_bind)<-cls_merge_bind[,1]
cls_merge_bind[,1]<-NULL
cls_merge_mt<-(as.matrix(cls_merge_bind))
cls_merge_mt[is.na(cls_merge_mt)]<-0
cls_merge_mt==0->cls_merge_mt_blk_indx
rowSums(cls_merge_mt_blk_indx)->fblk_count
fblk_count_thres<-fblk_count>=max(fblk_count)*0.90
cls_merge_mt_test<-t(cls_merge_mt[!fblk_count_thres,])


##Time series tSNE
features.tsne= fftRtsne (cls_merge_mt_test, dims = 2, max_iter = 1000)
features.tsne<-as.data.frame(features.tsne)
features.tsne$Bin<-rownames(cls_merge_mt_test)
Time_merged_tsne<-data.frame(features.tsne$V1,features.tsne$V2, features.tsne$Bin)
saveRDS(Time_merged_tsne, paste0(wd, "/Time_merged_tsne_all_bio.rds"))
```

### 6.3.1 tSNE Single dimension K means

``` r
library(dplyr)
library(Cardinal)

tSNE_dir<-"V:/Bioinformatics/FIt-SNE/"
source(paste0(tSNE_dir,"fast_tsne.R"))
FAST_TSNE_SCRIPT_DIR <- tSNE_dir

Time_merged_tsne <- readRDS(paste0(project_dir, dim_reduction_dir, tSNE_dir,"Time_merged_tsne_all_bio.rds"))
colnames(Time_merged_tsne)<- c("x", "y", "bin")

all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir,"all_coordata_label_sel_arc_bind2.csv"),row.names = 1)

example_imdata<-combinedimdata[,run(combinedimdata) == "4hrsilglucose_2nd_20201109_trim"]

dir.create(paste0(project_dir, dim_reduction_dir, tSNE_dir,"/kmean_seg_single_demo/"))

  
for (a in c("x", "y")) {

 
Time_merged_tsne  <-Time_merged_tsne [as.numeric(Time_merged_tsne $bin) %in% all_coordata_label_sel_arc_bind$new_label_scale[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim"],]

Time_merged_tsne $Rank<-dense_rank(Time_merged_tsne[,a])




png(paste0(project_dir, dim_reduction_dir, tSNE_dir,"/kmean_seg_single_demo/","tsne_",a,"_single_dim_lens_projections.png"),height = 5,width = 6,units = "in",res = 600)



match(all_coordata_label_sel_arc_bind$new_label_scale, as.numeric(Time_merged_tsne $bin))-> rank_index

all_coordata_label_sel_arc_bind$Rank<- Time_merged_tsne$Rank[rank_index]

factor(all_coordata_label_sel_arc_bind$Rank[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim"])->final_selection


fig.arbitrary_seg<-image(example_imdata, factor(final_selection) ~ x * y,superpose=F, key=F, xlab= NULL, col = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(max(Time_merged_tsne$Rank)) )



fig.arbitrary_seg[["par"]][["ann"]]=F
fig.arbitrary_seg[["par"]][["bty"]]="n"
fig.arbitrary_seg[["par"]][["pty"]]="s"
fig.arbitrary_seg[["par"]][["xaxt"]]="n"
fig.arbitrary_seg[["par"]][["yaxt"]]="n"
fig.arbitrary_seg[["par"]][["fg"]]="white"
fig.arbitrary_seg[["par"]][["oma"]]=c(0, 0, 0, 0)

print(fig.arbitrary_seg)
dev.off()
}
```

# 7. Metabolite library construction (function)

``` r
mummichog.lib.mod<-function(lib="bta_kegg",lib.new="bta_kegg_13C",C13_number=c(3,6),C13_deltamass=1.003355,wd=getwd(),force_update_lib=F,adducts_list="all",method=c("new_cpd","new_adduct")){

filenm <- paste(wd,"/",lib, ".qs", sep="")
mum.url <- paste("https://www.metaboanalyst.ca/resources/libs/mummichog/", paste(lib, ".qs", sep=""), sep="")
if(`|`(force_update_lib,!file.exists(filenm))){
download.file(mum.url, destfile = filenm, method="libcurl", mode = "wb")
}

mummichog.lib <- qs::qread(filenm)
cmpd.map <- MetaboAnalystR:::.get.my.lib("compound_db.qs")
    hit.inx <- match(tolower(mummichog.lib$cpd.lib$id), tolower(cmpd.map$kegg))
    nms <- cmpd.map[hit.inx, "smiles"]
    c_count <-stringr::str_count(nms,"C|c")
    nms_pub <- cmpd.map[hit.inx, "pubchem_id"]
    
df<-data.frame(smile=nms,c_count,mw_c_ratio=mummichog.lib$cpd.lib$mw/c_count,id=mummichog.lib$cpd.lib$id,mw=mummichog.lib$cpd.lib$mw,name=mummichog.lib$cpd.lib$name)
    
mummichog.lib.mod<-NULL
mummichog.lib.new<-mummichog.lib

if (method[1]=="new_cpd"){
for (cnum in C13_number){
  (df$c_count>=cnum)->mod.item
  (df$mw>=26.09864*as.numeric(cnum))->mod.item2
  mod.item[is.na(mod.item)]<-F
  mod.item2[is.na(mod.item2)]<-F
  mod.item_final<-as.numeric(ifelse(mod.item,mod.item,mod.item2))
  mummichog.lib.mod[[cnum]]<-mummichog.lib
  lapply(mummichog.lib.mod[[cnum]][["pathways"]][["cpds"]],function(x,cnum){
    paste0(x,"_13C",cnum)
  },cnum)->mummichog.lib.mod[[cnum]][["pathways"]][["cpds"]]
  mummichog.lib.mod[[cnum]][["cpd.lib"]][["id"]]<-paste0(mummichog.lib.mod[[cnum]][["cpd.lib"]][["id"]],"_13C",cnum)
    mummichog.lib.mod[[cnum]][["cpd.lib"]][["name"]]<-paste0(mummichog.lib.mod[[cnum]][["cpd.lib"]][["name"]],"_13C",cnum)
  mummichog.lib.mod[[cnum]][["cpd.lib"]][["mw"]]<-mummichog.lib.mod[[cnum]][["cpd.lib"]][["mw"]] + C13_deltamass*as.numeric(cnum)*mod.item_final
  mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["dpj_positive"]]<-mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["dpj_positive"]]+ C13_deltamass*as.numeric(cnum)*mod.item_final
  mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["positive"]]<-mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["positive"]]+ C13_deltamass*as.numeric(cnum)*mod.item_final
  mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["negative"]]<-mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["negative"]]+ C13_deltamass*as.numeric(cnum)*mod.item_final
  
}
for (cnum in C13_number){
  
lapply(1:length(mummichog.lib.new$pathways$cpds),function(x){
  c(mummichog.lib.new$pathways$cpds[[x]],mummichog.lib.mod[[cnum]]$pathways$cpds[[x]])
})->mummichog.lib.new$pathways$cpds
  
mummichog.lib.new$cpd.lib$id<-c(mummichog.lib.new$cpd.lib$id,mummichog.lib.mod[[cnum]]$cpd.lib$id)
mummichog.lib.new$cpd.lib$name<-c(mummichog.lib.new$cpd.lib$name,mummichog.lib.mod[[cnum]]$cpd.lib$name)
mummichog.lib.new$cpd.lib$mw<-c(mummichog.lib.new$cpd.lib$mw,mummichog.lib.mod[[cnum]]$cpd.lib$mw)
mummichog.lib.new[["cpd.lib"]][["adducts"]]$dpj_positive<-rbind(mummichog.lib.new[["cpd.lib"]][["adducts"]]$dpj_positive,mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]]$dpj_positive)
mummichog.lib.new[["cpd.lib"]][["adducts"]]$positive<-rbind(mummichog.lib.new[["cpd.lib"]][["adducts"]]$positive,mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]]$positive)
mummichog.lib.new[["cpd.lib"]][["adducts"]]$negative<-rbind(mummichog.lib.new[["cpd.lib"]][["adducts"]]$negative,mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]]$negative)
  
}
}





#MetaboAnalystR:::CreateLibFromKEGG(mummichog.lib.new$cpd.lib, mummichog.lib.new$pathways, lib.new)
  cpd.lib <- mummichog.lib.new$cpd.lib;
  ms_modes <- c('dpj_positive', 'positive', 'negative');
  adducts <- list();
  for (ms_mode in ms_modes){
    if(adducts_list[1]=="all"){
      adducts[[ms_mode]] <- MetaboAnalystR:::Compound_function_mzlist(ms_mode, cpd.lib$mw)
    }else{
      resdf <- MetaboAnalystR:::Compound_function_mzlist(ms_mode, cpd.lib$mw)
      if(sum(adducts_list %in% colnames(resdf))!=0){
      adducts[[ms_mode]]<-resdf[,adducts_list[adducts_list %in% colnames(resdf)]]
      }else{
      adducts[[ms_mode]]<-resdf
      }
    }
  }

if (method[1]=="new_adduct"){
  
    ms_modes <- c('dpj_positive', 'positive', 'negative')
  for (ms_mode in ms_modes){
    new_adduct<-NULL
  for (cnum in (C13_number)){
    new_adduct[[cnum]]<-adducts[[ms_mode]]+C13_deltamass*as.numeric(cnum)
    colnames(new_adduct[[cnum]])<-stringr::str_replace(colnames(new_adduct[[cnum]]),"^M",paste0("[M_13C",cnum,"]"))
  }
    new_adduct<-do.call(cbind,new_adduct)
    adducts[[ms_mode]]<-cbind(adducts[[ms_mode]],new_adduct)
  }
}
  cpd.lib$adducts <- adducts;
  
  # create a dictionary for look up in the range of 50-2000
  # now need to create ladder (tree) for each new mz
  # key is the mass 50 to 2000, values are the compounds (if any of their modified mw gives the value)
  # now create cpd tree for each mass pos
  # note, this can be slow, but this can be created before hand
  # for each species and for each mode
  # note l2 only stores the index of the cpd.lib
  
  cpd.tree <- list();
  for (ms_mode in ms_modes){
    l2 <- list();
    l2[[49]] <- "";
    l2[[2001]] <- "";
    mz.mat <- cpd.lib$adducts[[ms_mode]];
    floor.mzs <- floor(mz.mat);
    for(i in 1:nrow(floor.mzs)){
      neighbourhood <- floor.mzs[i,];
      for(n in neighbourhood){
        if((n>50) & (n<2000)){
          l2[[n]] <- append(l2[[n]], i);
        }
      }
    }
    cpd.tree[[ms_mode]] <- lapply(l2, unique);
  }
  
  # set up the variables
  mummichog.lib.new <- list(
    pathways = mummichog.lib.new$pathways,
    cpd.tree = cpd.tree,
    cpd.lib = cpd.lib
  )
qs::qsave(mummichog.lib.new,paste0(wd,"/",lib.new,".qs"))
return(mummichog.lib.new)
}
```

# 8. 1st round metabolomic annotation and pathway enrichment for all detected m/z features

## 8.1 Metabolite annotation and pathway enrichment

``` r
outwd = paste0(project_dir, dim_reduction_dir,PCA_dir)

for  (dim in paste0("Dim.", 1:5)){
  
dim.sig<-readr::read_csv(paste0(outwd, dim, "_anno_freq_top1bin.csv"))
dir.create(paste0(outwd,"/Metabo_",dim))
setwd(paste0(outwd,"/Metabo_",dim))
message(getwd())

mSet<-NULL

mSet<-InitDataObjects("mass_table", "mummichog", FALSE)
anal.type<<-mSet$analSet$type
mSet<-SetPeakFormat(mSet, "pvalue")
mSet<-UpdateInstrumentParameters(mSet, 10.0, "negative", "yes");
mSet<-SetRTincluded(mSet, "no")
mSet<-Read.TextData(mSet, paste0(project_dir, "/bin_data_time_sumspec_tb.csv"), "colu", "disc")
mSet<-SanityCheckMummichogData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckMummichogData(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "none", "none", "none", ratio=FALSE, ratioNum=20)



mSet<-SetPeakEnrichMethod(mSet, "mum", "v2")
mSet<-PreparePeakTable4PSEA(mSet)

feature_sig <-  names(mSet[["analSet"]][["tt"]][["t.score"]]) %in% dim.sig$mz
mSet[["analSet"]][["tt"]][["p.value"]][feature_sig] <-0.01
mSet[["analSet"]][["tt"]][["p.value"]][!feature_sig] <-1
nData<-mSet$dataSet$mummi.proc
feature_sig <-  nData$m.z %in% dim.sig$mz
match_mz(nData$m.z,dim.sig$mz,10,"ppm")->mapped_mz
feature_sig <- !is.na(mapped_mz)
nData[feature_sig,1] <-0.000001
nData[!feature_sig,1] <-1
nData->mSet$dataSet$mummi.proc

mSet<-SetMummichogPval(mSet, 0.05)
mummichog.lib.new<-mummichog.lib.mod(lib="bta_kegg",lib.new="bta_kegg_13C_adduct",C13_number=c(6),C13_deltamass=1.003355,wd=paste0(outwd,"/Metabo_",dim),force_update_lib=F,adducts_list=c("M+Cl[-]","M-H[-]"),method = "new_adduct")  
mSet$paramSet$ContainsMS2<-F
mSet<-PerformPSEA(mSet, "bta_kegg_13C_adduct", "current", 3 , 100)
mSet<-SaveTransformedData(mSet)
mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_0_", "png", 300, width=NA)
write.csv(annotate_PSEA(mSet = mSet,mummichog.lib.new),"mSet_mummi_resmat_anno.csv")
  
}

df_C13<-NULL
for (p_thres in paste0("Metabo_Dim.",1:5)){
df_C13[[as.character(p_thres)]]<-read.csv(paste0(outwd,"/",p_thres,"/","/","mSet_mummi_resmat_anno.csv"))
df_C13[[as.character(p_thres)]]$source<-p_thres
}
```

## 8.2 pathway summary

``` r
#read KEGG library

bta_kegg<- qs::qread( paste0(outwd, "/Metabo_Dim.1/bta_kegg.qs") )
cpd.info<- data.frame(cpd.name = bta_kegg[["cpd.lib"]][["name"]], id = bta_kegg[["cpd.lib"]][["id"]])
```

``` r
library(reshape2)
library(dplyr)
library(stringr)

df_C13<-NULL
for (p_thres in paste0("Metabo_Dim.",1:5)){
df_C13[[as.character(p_thres)]]<-read.csv(paste0(outwd,"/",p_thres,"/","mSet_mummi_resmat_anno.csv"))
df_C13[[as.character(p_thres)]]$source<-p_thres
}
df_C13<-do.call(rbind,df_C13)
df_C13$Gamma_log<--log(df_C13$Gamma,10)
df_C13$Pathway<-as.factor(df_C13$X)
df_C13$Pval_quartile<-df_C13$source
df_C13$PCA_dim<-as.numeric(str_remove(df_C13$source,"Metabo_Dim."))
df_C13_selected<-df_C13[df_C13$X %in% df_C13$X[df_C13$FET<=0.05],]
df_C13_selected$sig_level <- gtools::stars.pval(df_C13_selected$FET)



pathway_summary<-df_C13_selected[,c("X", "PCA_dim", "sig_level","Hits.sig")]
colnames(pathway_summary)<-c("Pathway", "Dim", "Sig_level","Hits.sig")
pathway_summary$Dim<-paste0("Dim",".",pathway_summary$Dim)
pathway_summary_sig<-dcast(pathway_summary,Pathway ~ Dim, value.var="Sig_level")
pathway_summary_sig[is.na(pathway_summary_sig)]<-""


sig.cpd.all <- list()
for (dim in paste0("Dim.", 1:5)){

matched.cpd<-read.csv(paste0(outwd, "/Metabo_", dim, "/mummichog_matched_compound_all.csv"))

sig.cpd.dim <- read.csv(paste0(outwd, "/", dim, "_anno_freq_top1bin.csv"))
sig.cpd.id <- match_mz(matched.cpd$Query.Mass, sig.cpd.dim$mz, 10,"ppm")
sig.cpd.dim <- matched.cpd[!is.na(sig.cpd.id),]
sig.cpd.dim$Dim <- paste0(dim)
sig.cpd.all[[dim]]  <- sig.cpd.dim
}
sig.cpd.all<- do.call(rbind, sig.cpd.all)

pathway.sig <- unique(df_C13_selected$X)
pathway.uni.id<- match ( pathway.sig, df_C13_selected$X)
pathway.sig<-  df_C13_selected[pathway.uni.id,]

pathway.cpd<- list()
for (pathway in  unique(df_C13_selected$X)){
  pathway.cpd[[pathway]] <- unlist ( str_split(pathway.sig[which (pathway.sig$X == pathway), which(colnames(pathway.sig)=="cpd_id")],";") )

}

lapply(names(pathway.cpd), function(x){
  sig.cpd.all[sig.cpd.all$Matched.Compound %in% pathway.cpd[[x]],]->df
  df$Pathway<-(x)
  df
  })->Pathway_matched_sig
Pathway_matched_sig<-do.call(rbind,Pathway_matched_sig)

Pathway_matched_sig$CPD_name<- cpd.info$cpd.name[match(Pathway_matched_sig$Matched.Compound, cpd.info$id)]
Pathway_matched_sig$Isotype=NA
Pathway_matched_sig$Isotype[str_detect(Pathway_matched_sig$Matched.Form, "13C")] <- "SIL_C6"
Pathway_matched_sig$Isotype[!str_detect(Pathway_matched_sig$Matched.Form, "13C")] <- "Normal"

df<- Pathway_matched_sig %>% group_by(Pathway, Dim) %>% dplyr::summarise(MZ_Normal = length(unique(Query.Mass[!str_detect(Matched.Form,"M_13C")])),MZ_SIL = length(unique(Query.Mass[str_detect(Matched.Form,"M_13C")]))) 


df_cast_SIL<- dcast(df,Pathway ~ Dim, value.var="MZ_SIL")
df_cast_SIL$Category<- "m/z SIL"
df_cast_SIL[is.na(df_cast_SIL)]<-0
df_cast_SIL<- df_cast_SIL%>% relocate(Category, .after = Pathway)


df_cast_Normal<- dcast(df,Pathway ~ Dim, value.var="MZ_Normal")
df_cast_Normal$Category<- "m/z Normal"
df_cast_Normal[is.na(df_cast_Normal)]<-0
df_cast_Normal<- df_cast_Normal%>% relocate(Category, .after = Pathway)

df_output <- rbind(df_cast_SIL,df_cast_Normal)

pathway_summary_sig$Category<- "Sig Level"
pathway_summary_sig<- pathway_summary_sig%>% relocate(Category, .after = Pathway)

df_output <- list(df_cast_Normal,df_cast_SIL, pathway_summary_sig )
df_output<-do.call(rbind,df_output)
# df_output<- df_output %>% relocate(Dim.10, .after = Dim.9)

write.csv(df_output, paste0(outwd, "pathway_final_output_1_to_5.csv"), row.names = F)
# 
# Pathway_matched_sig <- Pathway_matched_sig %>% group_by(Matched.Compound, Matched.Form) %>%  dplyr::summarise(mean = mean (Query.Mass), ) 
```

``` r
Pathway_matched_sig_summary<-Pathway_matched_sig
Pathway_matched_sig_summary<-Pathway_matched_sig_summary[!duplicated(Pathway_matched_sig_summary[,c( "Matched.Compound", "Matched.Form","Dim" )]),]
```

# 9. Feature scoring and FDR controlled filtering

``` r
library(HiTMaP)
library(stringr)
  wd="V:/Bioinformatics/test_scoring_cpd_mummichog/"
  raw_wd<-"G:/imzML" 
  #spectrum_to_match<- read.csv(paste0(wd,"raw_data_time_sumspec.csv"))
    cpd_list<-read.csv(paste0(wd,"pca_anal_thesis_data_test_tol_10/Metabo_Dim.1/mummichog_matched_compound_all.csv"))
    cmpd.map<-read.csv(paste0(wd,"cmpd.map_anno.csv"))
    fileinfo_v2 <- read.csv(paste0(raw_wd,"/fileinfo_v2.csv"))
    #SMPdb<-read.csv(paste0(wd,"SMPdb.csv"),row.names = 1)
    #spectrum_to_match<-read.csv(paste0(wd,"bin_data_time_sumspec.csv"),row.names=1)
    locs_file<-fileinfo_v2$filenames
    locs<-paste0("G:/imzML/",locs_file," ID/")

    cpd_list$formula<-cmpd.map$formula[match(cpd_list$Matched.Compound,cmpd.map$kegg_id)]
    cpd_list$Isotype<-""  
    cpd_list<-cpd_list[!is.na(cpd_list$formula),]
    cpd_list->cpd_list_org
    for (score_method in c("SQRTP")){
    for (timepoint in unique(locs_file)){
    #for (timepoint in "20hrSILGlucose_03_20211123"){
    #"16hrSILGlucose_01_20211123_trim"->timepoint
    for (FDR_cutoff in c(0.1))  {
      dir.create(paste0(wd,"/",score_method,"_",FDR_cutoff))
    peaklist<-as.data.frame(readr::read_csv(paste0(locs<-paste0("G:/imzML/",timepoint," ID/"),"overall_Spectrum.csv")))
library(HiTMaP)  
library(BiocParallel)
   HiTMaP:::Cpd_spectrum_match_rescore(cpd_list,peaklist,
                                       wd=paste0(wd,"/",score_method,"_",FDR_cutoff),
                                       adjust_score=T,
                                       SPECTRUM_batch=timepoint,
                                       BPPARAM=SerialParam(),
                                     atom_isotope_sub=NULL,
                                     ppm=10,
                                     score_method=score_method,
                                     FDR_cutoff=FDR_cutoff,Decoy_search=T)
    }
      }
        }
    
    
    peak_res_sum<-NULL
    peak_res_sum_2nd<-NULL
    for (score_method in c("SQRTP")){
      dir.create(paste0(wd,"/",score_method))
    for (timepoint in unique(locs_file)){
    for (FDR_cutoff in c(0.1))  {
      dir.create(paste0(wd,"/",score_method,"_",FDR_cutoff))
    peak_res<-read.csv(paste0(wd,"/",score_method,"_",FDR_cutoff,"/",timepoint,"/","CPD_1st_ID_score_rank_SQRTP.csv"))
    peak_res2nd<-read.csv(paste0(wd,"/",score_method,"_",FDR_cutoff,"/",timepoint,"/","CPD_2nd_ID_score_rank_SQRTP.csv"))
    if(nrow(peak_res)>=1) peak_res$source<-paste0(score_method,"_",FDR_cutoff,"_",timepoint)
    peak_res_sum[[paste0(score_method,"_",FDR_cutoff,"_",timepoint)]]<-peak_res
    if(nrow(peak_res2nd)>=1) peak_res2nd$source<-paste0(score_method,"_",FDR_cutoff,"_",timepoint)
    peak_res_sum_2nd[[paste0(score_method,"_",FDR_cutoff,"_",timepoint)]]<-peak_res2nd
    }
    }
    }
    peak_res_sum_bind<-do.call(rbind,peak_res_sum)
    peak_res_sum_2nd_bind<-do.call(rbind,peak_res_sum_2nd)
    peak_res_sum_bind$CPD_file<-paste0(peak_res_sum_bind$formula,"_",
                                       peak_res_sum_bind$isdecoy,"_",
                                       peak_res_sum_bind$Isotype,"_",
                                       peak_res_sum_bind$source)
    peak_res_sum_2nd_bind$CPD_file<-paste0(peak_res_sum_2nd_bind$formula,"_",
                                           peak_res_sum_2nd_bind$isdecoy,"_",
                                           peak_res_sum_2nd_bind$Isotype,"_",
                                           peak_res_sum_2nd_bind$source)

    
    target_fm<-c("C6H12O6","C12H22O11","C15H24N2O17P2","C6H14O6","C7H15O10P","C6H13O9P1","C6H14O12P2")
    target_nm<-c("Glucose","Sucrose","UDP_Glucose","Sorbitol","S7P","G6P","F16P")
    peak_res_sum_bind_tgt<-peak_res_sum_bind[peak_res_sum_bind$formula_mono %in% target_fm,]
    peak_res_sum_bind_tgt_meta<-match(str_split_fixed(peak_res_sum_bind_tgt$source,"_",3)[,3],fileinfo_v2$filenames)
    peak_res_sum_bind_tgt<-cbind(peak_res_sum_bind_tgt,fileinfo_v2[peak_res_sum_bind_tgt_meta,])
    peak_res_sum_bind_tgt$FDR_pass<-peak_res_sum_bind_tgt$CPD_file %in% unique(peak_res_sum_2nd_bind$CPD_file)
    peak_res_sum_bind_tgt<-unique(peak_res_sum_bind_tgt[,c("Score" ,"adduct", "formula","formula_mono","Isotype","isdecoy","mins","bioreplicate","class" ,"FDR_pass" )])
    peak_res_sum_bind_tgt$cpd<-target_nm[match(peak_res_sum_bind_tgt$formula_mono,target_fm)]
    
        data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
    }
library(ggplot2)
       df2 <- data_summary(peak_res_sum_bind_tgt, varname="Score", 
                    groupnames=c("cpd","Isotype","isdecoy","class","adduct"))
    df2 <- merge(df2,peak_res_sum_bind_tgt[,c("cpd","Isotype","isdecoy","class","adduct","FDR_pass")],by=c("cpd","Isotype","isdecoy","class","adduct"))
    p<-ggplot(df2) +geom_point(aes(x=class,y=Score,colour=FDR_pass,size=FDR_pass)) + geom_line(aes(x=class,y=Score,color=interaction(isdecoy))) + 
  geom_errorbar(aes(x=class,ymin=Score-sem, ymax=Score+sem), width=.2,
                 position=position_dodge(.9)) + 
facet_grid(rows = vars(Isotype,adduct),cols = vars(cpd))
    png(paste0(wd,"/","/target_fm_summary_",FDR_cutoff,".png"),height = 2400*2,width = 1500*length(target_nm),res = 300)

print(p)
    dev.off()    
    write.csv(peak_res_sum_bind,paste0(wd,"/",score_method,"_","peak_res_sum_bind_",FDR_cutoff,".csv"),row.names = F)
    write.csv(peak_res_sum_2nd_bind,paste0(wd,"/",score_method,"_","peak_res_sum_2nd_bind_",FDR_cutoff,".csv"),row.names = F)
#     
#     for (score_method in c("Equal-intensity-SQRT",
#                            "Mix-SQRT","Equal-SQRT",
#                            "balanced-SQRT",
#                            "SQRT","SQRTP")){
#       dir.create(paste0(wd,"/",score_method))
#     for (timepoint in colnames(spectrum_to_match)){
#     peaklist<-data.frame(m.z=as.numeric(rownames(spectrum_to_match)),intensities=spectrum_to_match[,timepoint])
#     
# library(HiTMaP)  
# library(BiocParallel)
#    HiTMaP:::Cpd_spectrum_match_rescore(cpd_list,peaklist,wd=paste0(wd,"/",score_method),
#                                        SPECTRUM_batch=timepoint,
#                                        BPPARAM=SerialParam(),
#                                      atom_isotope_sub=NULL,
#                                      ppm=10,
#                                      score_method=score_method,
#                                      adjust_score=F,
#                                      FDR_cutoff=0.15,Decoy_search=T)
#    }
#     }
```

# 10. 2nd round metabolomic annotation and pathway enrichment within FDR filtered features

``` r
library(stringr)
library(readr)

FDR_wd<- paste0(project_dir, FDR_dir)
main_dir <- FDR_wd
dir_list <- list.dirs(main_dir,recursive = FALSE)

all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir, "all_coordata_label_sel_arc_bind2.csv"),row.names = 1)

match.features <- read_csv(paste0(outwd,"/Metabo_Dim.1/mummichog_matched_compound_all.csv" ))
match.features$adduct <- NA
match.features$adduct[str_detect(match.features$Matched.Form,"Cl")]<- "M+Cl"
match.features$adduct[!str_detect(match.features$Matched.Form,"Cl")]<- "M-H"

match.features$Isotype<- NA
match.features$Isotype[str_detect(match.features$Matched.Form,"13C6")]<- "SIL"
match.features$Isotype[!str_detect(match.features$Matched.Form,"13C6")]<- "Normal"

folder_name<- list.files(FDR_wd)
folder.slt <- c() 

for (x in 1:length(folder_name)){

folder.exist = str_detect(unique(all_coordata_label_sel_arc_bind$run),tolower(folder_name[x]))

if (sum(folder.exist) > 0){
  folder.slt<- append(folder.slt,folder_name[x])
}
  else {
    message("Data in this run was not included.")
  }
}

cpd.list<- list() 
  
for (dir in folder.slt){
  cpd.df <- read_csv(paste0(FDR_wd,dir,"/CPD_2nd_ID_score_rank_SQRTP.csv"))
  run <- str_remove(dir, paste0(main_dir,"/"))
  
  cpd.df$run<-run
  cpd.list[[run]] <-  cpd.df
}

cpd.list.all<- do.call(rbind, cpd.list)



file.info<-fileinfo_v2

cpd_filter <- unique(cpd.list.all[,c("Isotype", "Metabolite.Name", "adduct")])

cpd_filter <- unique(cpd.list.all[,c("mz","Isotype", "Metabolite.Name", "adduct")])

cpd_filter$match <- paste(cpd_filter$Metabolite.Name, cpd_filter$Isotype, cpd_filter$adduct, sep = "_")



match.features$Isotype<- NA
match.features$Isotype[str_detect(match.features$Matched.Form,"13C")] <- "SIL"
match.features$Isotype[!str_detect(match.features$Matched.Form,"13C")] <- "Normal"

match.features$adduct<- NA
match.features$adduct[str_detect(match.features$Matched.Form, "Cl")]  <- "M+Cl"
match.features$adduct[!str_detect(match.features$Matched.Form, "Cl")]  <- "M-H"

match.features$match<- paste(match.features$Matched.Compound, match.features$Isotype, match.features$adduct, sep = "_")


match.features_HQ<-match.features[match.features$match %in% cpd_filter$match,]

cpd_filter <- unique(match.features_HQ$Query.Mass)
file.path <- outwd
for (dim in paste0("Dim.", 1:5)){
  top.features <- read_csv(paste0(file.path ,dim,"_anno_freq_top1bin.csv"))
  
top.features.flt <- top.features[which(top.features$mz %in%  cpd_filter),]

write.csv(top.features.flt, paste0(file.path, dim,"_anno_freq_top1bin_FDR_filter_0.1.csv"), row.names = FALSE)
}

df <- read.csv(paste0(project_dir,"/bin_data_time_sumspec.csv") , row.names = 1)
df_filter <- df[as.numeric(rownames(df)) %in% cpd_filter,] 

label<- c(rep("Early", 5), rep("Late", 5) )
df_filter_label <-rbind(label, df_filter )

write.csv(df_filter_label, paste0(project_dir, "/bin_data_time_sumspec_FDR_0.1_tb.csv") )
write.csv(df_filter,paste0(project_dir,  "/bin_data_time_sumspec_FDR_0.1.csv") )


for  (dim in paste0("Dim.", 1:5)){
  
dim.sig<-readr::read_csv(paste0(outwd, dim, "_anno_freq_top1bin_FDR_filter_0.1.csv"))
dir.create(paste0(outwd,"/Metabo_",dim,"FDR_0.1"))
setwd(paste0(outwd,"/Metabo_",dim,"FDR_0.1"))
message(getwd())

mSet<-NULL
library(MetaboAnalystR)
mSet<-InitDataObjects("mass_table", "mummichog", FALSE)
anal.type<<-mSet$analSet$type
mSet<-SetPeakFormat(mSet, "pvalue")
mSet<-UpdateInstrumentParameters(mSet, 10.0, "negative", "yes");
mSet<-SetRTincluded(mSet, "no")
mSet<-Read.TextData(mSet, paste0(project_dir, "/bin_data_time_sumspec_FDR_0.1_tb.csv"), "colu", "disc")
mSet<-SanityCheckMummichogData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckMummichogData(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "none", "none", "none", ratio=FALSE, ratioNum=20)



mSet<-SetPeakEnrichMethod(mSet, "mum", "v2")
mSet<-PreparePeakTable4PSEA(mSet)

feature_sig <-  names(mSet[["analSet"]][["tt"]][["t.score"]]) %in% dim.sig$mz
mSet[["analSet"]][["tt"]][["p.value"]][feature_sig] <-0.01
mSet[["analSet"]][["tt"]][["p.value"]][!feature_sig] <-1
nData<-mSet$dataSet$mummi.proc
feature_sig <-  nData$m.z %in% dim.sig$mz
match_mz(nData$m.z,dim.sig$mz,10,"ppm")->mapped_mz
feature_sig <- !is.na(mapped_mz)
nData[feature_sig,1] <-0.000001
nData[!feature_sig,1] <-1
nData->mSet$dataSet$mummi.proc

mSet<-SetMummichogPval(mSet, 0.05)
mummichog.lib.new<-mummichog.lib.mod(lib="bta_kegg",lib.new="bta_kegg_13C_adduct",C13_number=c(6),C13_deltamass=1.003355,wd=paste0(outwd,"/Metabo_",dim, "FDR_0.1"),force_update_lib=F,adducts_list=c("M+Cl[-]","M-H[-]"),method = "new_adduct")  
mSet$paramSet$ContainsMS2<-F
mSet<-PerformPSEA(mSet, "bta_kegg_13C_adduct", "current", 3 , 100)
mSet<-SaveTransformedData(mSet)
mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_0_", "png", 300, width=NA)
write.csv(annotate_PSEA(mSet = mSet,mummichog.lib.new),"mSet_mummi_resmat_anno.csv")
  
}

df_C13<-NULL
for (p_thres in paste0("Metabo_Dim.",1:5,"FDR_0.1")){
df_C13[[as.character(p_thres)]]<-read.csv(paste0(outwd,"/",p_thres,"/","/","mSet_mummi_resmat_anno.csv"))
df_C13[[as.character(p_thres)]]$source<-p_thres
}

df_C13<-do.call(rbind,df_C13)


library(reshape2)
library(dplyr)


df_C13<-NULL
for (p_thres in paste0("Metabo_Dim.",1:5)){
df_C13[[as.character(p_thres)]]<-read.csv(paste0(outwd,"/",p_thres,"FDR_0.1","/","mSet_mummi_resmat_anno.csv"))
df_C13[[as.character(p_thres)]]$source<-p_thres
}
df_C13<-do.call(rbind,df_C13)
df_C13$Gamma_log<--log(df_C13$Gamma,10)
df_C13$Pathway<-as.factor(df_C13$X)
df_C13$Pval_quartile<-df_C13$source
df_C13$PCA_dim<-as.numeric(str_remove(df_C13$source,"Metabo_Dim."))
df_C13_selected<-df_C13[df_C13$X %in% df_C13$X[df_C13$FET<=0.05],]
df_C13_selected$sig_level <- gtools::stars.pval(df_C13_selected$FET)



pathway_summary<-df_C13_selected[,c("X", "PCA_dim", "sig_level","Hits.sig")]
colnames(pathway_summary)<-c("Pathway", "Dim", "Sig_level","Hits.sig")
pathway_summary$Dim<-paste0("Dim",".",pathway_summary$Dim)
pathway_summary_sig<-dcast(pathway_summary,Pathway ~ Dim, value.var="Sig_level")
pathway_summary_sig[is.na(pathway_summary_sig)]<-""


sig.cpd.all <- list()
for (dim in paste0("Dim.", 1:5)){

matched.cpd<-read.csv(paste0(outwd, "/Metabo_", dim,"FDR_0.1", "/mummichog_matched_compound_all.csv"))

sig.cpd.dim <- read.csv(paste0(outwd, "/", dim, "_anno_freq_top1bin_FDR_filter_0.1.csv"))



sig.cpd.id <- match_mz(matched.cpd$Query.Mass, sig.cpd.dim$mz, 10,"ppm")
  
   




sig.cpd.dim <- matched.cpd[!is.na(sig.cpd.id),]
sig.cpd.dim$Dim <- paste0(dim)
sig.cpd.all[[dim]]  <- sig.cpd.dim
}
sig.cpd.all<- do.call(rbind, sig.cpd.all)

pathway.sig <- unique(df_C13_selected$X)
pathway.uni.id<- match ( pathway.sig, df_C13_selected$X)
pathway.sig<-  df_C13_selected[pathway.uni.id,]

pathway.cpd<- list()
for (pathway in  unique(df_C13_selected$X)){
  pathway.cpd[[pathway]] <- unlist ( str_split(pathway.sig[which (pathway.sig$X == pathway), which(colnames(pathway.sig)=="cpd_id")],";") )

}

lapply(names(pathway.cpd), function(x){
  sig.cpd.all[sig.cpd.all$Matched.Compound %in% pathway.cpd[[x]],]->df
  df$Pathway<-(x)
  df
  })->Pathway_matched_sig
Pathway_matched_sig<-do.call(rbind,Pathway_matched_sig)

df<- Pathway_matched_sig %>% group_by(Pathway, Dim) %>% dplyr::summarise(MZ_Normal = length(unique(Query.Mass[!str_detect(Matched.Form,"M_13C")])),MZ_SIL = length(unique(Query.Mass[str_detect(Matched.Form,"M_13C")]))) 


df_cast_SIL<- dcast(df,Pathway ~ Dim, value.var="MZ_SIL")
df_cast_SIL$Category<- "m/z SIL"
df_cast_SIL[is.na(df_cast_SIL)]<-0
df_cast_SIL<- df_cast_SIL%>% relocate(Category, .after = Pathway)


df_cast_Normal<- dcast(df,Pathway ~ Dim, value.var="MZ_Normal")
df_cast_Normal$Category<- "m/z Normal"
df_cast_Normal[is.na(df_cast_Normal)]<-0
df_cast_Normal<- df_cast_Normal%>% relocate(Category, .after = Pathway)

df_output <- rbind(df_cast_SIL,df_cast_Normal)

pathway_summary_sig$Category<- "Sig Level"
pathway_summary_sig<- pathway_summary_sig%>% relocate(Category, .after = Pathway)

df_output <- list(df_cast_Normal,df_cast_SIL, pathway_summary_sig )
df_output<-do.call(rbind,df_output)
# df_output<- df_output %>% relocate(Dim.10, .after = Dim.9)

write.csv(df_output, paste0(outwd, "pathway_final_output_1_to_5_FDR_0.1_remapping.csv"), row.names = F)
```

# 11. Feature visualization

## 11.1 Selected bio rep

``` r
library(Cardinal)
library(stringr)
library(dplyr)

heximg_dir <- "/heximg/"
dir.create(paste0(project_dir, visualization_dir, heximg_dir))

run_demo_slt<-c("ff_trim_02_08042021","5minsilglucose_2nd_20201116_trim","15minsilglucose_2nd_20201116_trim","30minsilglucose_2nd_20201116_trim","1hrsilglucose_201130_trim","2hrsilglucose_2nd_20201109_trim","4hr_01_20201109_trim","8hr_01_20201109_trim","16hrsilglucose_03_20211123_trim","20hrsilglucose_02_20211123_trim")

combinedimdata_demo<-combinedimdata[,run(combinedimdata) %in% run_demo_slt]

all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir,  spatial_alignment_dir,"all_coordata_label_sel_arc_bind2.csv"),row.names = 1)
pixeldf<-all_coordata_label_sel_arc_bind
```

``` r
library(Cardinal)

combinedimdata_full<- combinedimdata

fileinfo_v2 <- readr::read_csv(paste0(project_dir,"/fileinfo_v2.csv"))

# final_selection_df<-read.csv(paste0(wd,"/final_selection.csv"))

fileinfo_v2$filenames<-tolower(fileinfo_v2$filenames)
fileinfo_v2<-fileinfo_v2[fileinfo_v2$filenames %in% levels(run(combinedimdata_demo)),]
unique(fileinfo_v2$filenames)

coordata<-combinedimdata_demo@elementMetadata@coord@listData
coordata$run<-run(combinedimdata_demo)
coordata<-as.data.frame(coordata)
unique(coordata$run)
#hexbin all features

combinedimdata_full<-combinedimdata_demo


library(hexbin)
label_sel="slt"
merged_data_df<-NULL
coordata_label_sel_arc<-NULL

    setCardinalBPPARAM(BPPARAM=SerialParam())
    setCardinalVerbose(verbose=F)
    for (imsrun in fileinfo_v2$filenames){
      message(imsrun)
      time_cls<-fileinfo_v2$class[match(imsrun,fileinfo_v2$filenames)]
      combinedimdata_full@elementMetadata@coord@listData->coordata
      coordata$run<-run(combinedimdata_full)
      coordata<-as.data.frame(coordata)
      coordata_run<-coordata[coordata$run==imsrun,]
      coordata_run<-pixeldf[pixeldf$run==imsrun,]
      rownames(coordata_run)<-1:nrow(coordata_run)
      library(hexbin)
      combinedimdata_full[,run(combinedimdata_full)==imsrun]->combinedimdata_full_sample_bin
    
    bincell_df<-NULL
    for (bincell in unique(coordata_run$new_label_scale)){
      combinedimdata_full_sample_bin[,which(coordata_run$new_label_scale==bincell)] %>% summarizeFeatures( FUN = "mean") %>% process()->meanspec
      meanspec@featureData@listData[["mean"]]->bincell_df[[bincell]]
    }
    
    bincell_df_bind<-do.call(cbind,bincell_df)
    colnames(bincell_df_bind)<-paste(imsrun,time_cls,unique(coordata_run$new_label_scale),sep = "@")
    merged_data_df[[imsrun]]<-bincell_df_bind
    coordata_label_sel_arc[[imsrun]]<-coordata_run
    }

coordata_label_sel_arc_bind<-do.call(rbind,coordata_label_sel_arc)
merged_data_df_bind<-do.call(cbind,merged_data_df)
rownames(merged_data_df_bind)<-combinedimdata_full@featureData@mz


#Generate hexbin data
write.csv(merged_data_df_bind,paste0(project_dir, visualization_dir , heximg_dir,label_sel,"_merged_data_df_bind",".csv"))

#Generate Master Table 
write.csv(coordata_label_sel_arc_bind,paste0(project_dir, visualization_dir , heximg_dir,label_sel,"_coordata_label_sel_arc_bind",".csv"))


message(label_sel," DONE")
```

``` r
setwd(wd)
library(readr)
all_data_original<-read_csv(paste0(project_dir, visualization_dir , heximg_dir,"/slt_merged_data_df_bind.csv"))
all_data_original<-as.data.frame(all_data_original)
rownames(all_data_original)<-all_data_original[,1]
all_data_original[,-1]->all_data_original
library(stringr)
col.nms.split<-str_split_fixed(colnames(all_data_original),"@",3)

all_data_original_t<-t(all_data_original)

group_sum <- rowsum(all_data_original_t, group = interaction (col.nms.split[,2],col.nms.split[,3]))

group_sum_t<-as.data.frame(t(group_sum))

a = paste(col.nms.split[,2], col.nms.split[,3], sep=".")

for (col in colnames(group_sum_t)) { 
  group_sum_t [,col]=  group_sum_t [,col] / sum(a==col)
  # message(col,sum(a==col))
}

grouped_bio_rep_all_pixel_mean<-group_sum_t

grouped_bio_rep_all_pixel_mean<-as.data.frame(grouped_bio_rep_all_pixel_mean)

library(stringr)
col.nms.split_2<-str_split_fixed(colnames(grouped_bio_rep_all_pixel_mean),"\\.",2)

Label<-col.nms.split_2[,1]
grouped_bio_rep_all_pixel_mean<-rbind(Label, grouped_bio_rep_all_pixel_mean)
row.names(grouped_bio_rep_all_pixel_mean)[1]<-"Label"
# grouped_bio_rep_all_pixel_mean<-grouped_bio_rep_all_pixel_mean[,order(grouped_bio_rep_all_pixel_mean[1,])]


cls_merge<-list()
for (cls in unique(Label)){
    group_sum_t[,col.nms.split_2[,1]==cls]->temp
    colnames(temp)<-str_remove(colnames(temp),paste0(cls,"\\."))
    row.names(temp)<-paste0(row.names(temp),"@",cls)
  cls_merge[[cls]]<-(temp)
  
}

library(plyr)
cls_merge_bind<-rbind.fill(cls_merge)
rownms<-unlist(lapply(cls_merge,rownames))
rownames(cls_merge_bind)<-rownms

write.csv(cls_merge_bind, paste0(project_dir, visualization_dir , heximg_dir, "grouped_bio_rep_slt_pixel_mean_label_time_merge_new.csv"))
```

## 11.2 Visualization

``` r
mum.url <- paste("https://www.metaboanalyst.ca/resources/libs/mummichog/bta_kegg.qs")

download.file(mum.url, paste0(project_dir, visualization_dir, heximg_dir, "bta_kegg.qs"), 
        mode = "wb")

bta_kegg<-qs::qread(paste0(project_dir, visualization_dir, heximg_dir, "bta_kegg.qs"))
cpd_list<- data.frame(id =bta_kegg[["cpd.lib"]][["id"]], name = bta_kegg[["cpd.lib"]][["name"]])

cpd.match <- read_csv(paste0(project_dir, dim_reduction_dir, PCA_dir, "/Metabo_Dim.1/mummichog_matched_compound_all.csv"))


cpd.match$name <- cpd_list$name[match(cpd.match$Matched.Compound, cpd_list$id)]

target_cpd<- c("C00031","C00668","C05378","C05382", "C00029", "C00794")
cpd.target<- cpd.match[which(cpd.match$Matched.Compound %in% target_cpd),] 
cpd.target.sil<- cpd.target[str_detect(cpd.target$Matched.Form,"13C6"),]


library(RColorBrewer)
fun_color_range <- colorRampPalette(brewer.pal(9, "Spectral"))
my_colors <- fun_color_range(10)



library(Cardinal)
library(stringr)
library(dplyr)

all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir, "all_coordata_label_sel_arc_bind2.csv"),row.names = 1)

cls_merge_bind<-readr::read_csv(paste0(project_dir, visualization_dir, heximg_dir,"grouped_bio_rep_slt_pixel_mean_label_time_merge_new.csv"))

cls_merge_bind<-as.data.frame(cls_merge_bind)
rownames(cls_merge_bind)<-cls_merge_bind[,1]
cls_merge_bind[,1]<-NULL


cpd.target<-cpd.target[!duplicated(paste(cpd.target$Matched.Compound, cpd.target$Matched.Form, sep = ";")),] 



for (mass in cpd.target$Query.Mass) {
glucose_intensity <- cls_merge_bind[stringr::str_detect(rownames(cls_merge_bind), paste0(mass)),]
glucose_intensity[is.na(glucose_intensity)]<-0

glucose_intensity<-as.data.frame (t (glucose_intensity))
colnames(glucose_intensity)<- paste("Time",str_remove(colnames(glucose_intensity),paste0(mass,"@"))) 
glucose_intensity<-glucose_intensity[,order(colnames(glucose_intensity))]

for (a in c(0, 0.05,0.1 , 0.12,0.15)){
  for (col in 1:10){
  
  target_class <- data.frame(intensity=glucose_intensity[,colnames(glucose_intensity)[col]], row.names = rownames(glucose_intensity) )

all_coordata_label_sel_arc_bind$intensity<-target_class$intensity[match(all_coordata_label_sel_arc_bind$new_label_scale, as.numeric(rownames(target_class)))]
  
all_coordata_label_sel_arc_bind$intensity[is.na(all_coordata_label_sel_arc_bind$intensity)]<-0

all_coordata_label_sel_arc_bind$signal_exist<- NA
all_coordata_label_sel_arc_bind$signal_exist[which(all_coordata_label_sel_arc_bind$intensity <= a)]<- 0
all_coordata_label_sel_arc_bind$signal_exist[is.na(all_coordata_label_sel_arc_bind$signal_exist)]<-1

factor(all_coordata_label_sel_arc_bind$signal_exist[all_coordata_label_sel_arc_bind$run %in% unique(run(example_imdata))])->final_selection

fig.arbitrary_seg<- image(example_imdata, factor(final_selection) ~ x * y,superpose=F, key=F, col = c("#FFFFFF00", my_colors[col]))


fig.arbitrary_seg[["par"]][["ann"]]=F
fig.arbitrary_seg[["par"]][["bty"]]="n"
fig.arbitrary_seg[["par"]][["pty"]]="s"
fig.arbitrary_seg[["par"]][["xaxt"]]="n"
fig.arbitrary_seg[["par"]][["yaxt"]]="n"
fig.arbitrary_seg[["par"]][["fg"]]="white"
fig.arbitrary_seg[["par"]][["oma"]]=c(0, 0, 0, 0)

dir.create(paste0(project_dir, visualization_dir, heximg_dir,"image_time_stack_selected/"))
png(paste0(project_dir, visualization_dir, heximg_dir,"/image_time_stack_selected/", mass,"_time_",col,"_", a,".png"),height = 5 ,width = 6,units = "in",res = 300, bg="transparent")

print(fig.arbitrary_seg)
dev.off()
}


col=10:1

img_set<-magick::image_read(paste0(project_dir, visualization_dir, heximg_dir,"/image_time_stack_selected/", mass,"_time_",col,"_", a,".png"))
library(magick)
image_mosaic(img_set)-> image_stack
image_write(image_stack,paste0(project_dir, visualization_dir, heximg_dir,"/image_time_stack_selected/", mass,"_time_stack_",a,".png"))

}
}
```

## 11.3 Metabolite annotation in pathway format

### 11.3.1 ion image plotting (skip if ion image available )

``` r
library(stringr)

run_demo_slt<-c("ff_trim_02_08042021","5minsilglucose_2nd_20201116_trim","15minsilglucose_2nd_20201116_trim","30minsilglucose_2nd_20201116_trim","1hrsilglucose_201130_trim","2hrsilglucose_2nd_20201109_trim","4hr_01_20201109_trim","8hr_01_20201109_trim","16hrsilglucose_03_20211123_trim","20hrsilglucose_02_20211123_trim")

combinedimdata_demo<-combinedimdata[,run(combinedimdata) %in% run_demo_slt]

run(combinedimdata_demo)->run_vec
run_vec_replace<-str_extract(as.character(run_vec),"^.{0,2}min|^.{0,2}hr|ff")
run_vec_replace<-ifelse(str_detect(run_vec_replace,"min"),str_remove(run_vec_replace,"min"),ifelse(str_detect(run_vec_replace,"hr"),as.numeric(str_remove(run_vec_replace,"hr"))*60,"0"))
run(combinedimdata_demo)<-factor(run_vec_replace,level=as.character(sort(unique(as.numeric(run_vec_replace)))))

cpd_mapped <- read.csv(paste0(project_dir, dim_reduction_dir, PCA_dir, "/Metabo_Dim.1FDR_0.1/mummichog_matched_compound_all.csv"))


combinedimdata_demo[mz(combinedimdata_demo) %in% cpd_mapped$Query.Mass, ] -> combinedimdata_demo_cpd_match
saveRDS(combinedimdata_demo_cpd_match,paste0(project_dir, dim_reduction_dir, PCA_dir, "/combinedimdata_demo_cpd_match.rds"))
```

``` r
library(Cardinal)
readRDS(paste0(project_dir, dim_reduction_dir, PCA_dir, "/combinedimdata_demo_cpd_match.rds"))->combinedimdata_demo_cpd_match

ion_images_dir <- "/ion_images/"

dir.create(paste0(project_dir, dim_reduction_dir, PCA_dir, ion_images_dir))

outwd= paste0(project_dir, dim_reduction_dir, PCA_dir, ion_images_dir)

for (x in mz(combinedimdata_demo_cpd_match)){
png(paste0(outwd,"/",x,".png"),width = 300 * length(levels(run(combinedimdata_demo_cpd_match))),height =450,units = "px")
darkmode()
par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 0, 0),mfrow = c(1,1),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)
p<-image(combinedimdata_demo_cpd_match, mz=as.numeric(x),superpose=F, key=F,plusminus=as.numeric(x)/100000,contrast.enhance ="suppression",normalize.image = c("none"),smooth.image = c("none"), layout=c(1,length(levels(run(combinedimdata_demo_cpd_match)))))
lapply(p[["facets"]],function(x) {
  attr(x,"strip")$text<-""
  return(x)
})->p[["facets"]]
print(p)
dev.off()
}
```

### 11.3.2 Assign metabolites to pathways

``` r
df_output<- read.csv(paste0(project_dir, dim_reduction_dir, PCA_dir,  "/pathway_final_output_1_to_5_FDR_0.1_remapping.csv") )

cpd.high.quality <- read.csv(paste0(project_dir, dim_reduction_dir, PCA_dir, "/Metabo_Dim.1FDR_0.1/mummichog_matched_compound_all.csv"))

pathway.target<- unique(df_output$Pathway)
pathway.target <- pathway.db[which(pathway.info$pathway.name %in% pathway.target),]
rownames(pathway.target)<-pathway.target$pathway.name

outwd<- paste0(project_dir, dim_reduction_dir, PCA_dir)

sig.cpd.all <- list()
for (dim in paste0("Dim.", 1:5)){

sig.cpd.dim <- read.csv(paste0(outwd, "/", dim, "_anno_freq_top1bin_FDR_filter_0.1.csv"))

sig.cpd.dim <- cpd.high.quality[which(cpd.high.quality$Query.Mass %in% sig.cpd.dim$mz ),]
  sig.cpd.dim$Dim <- paste0(dim)
sig.cpd.all[[dim]]  <- sig.cpd.dim
}
sig.cpd.all<- do.call(rbind, sig.cpd.all)

pathway.sig <- unique(df_output$Pathway)
# pathway.sig<-  df_C13_selected[match ( pathway.sig, df_C13_selected$X),]

pathway.cpd<- list()
for (pathway in  pathway.target$pathway.name){
  pathway.cpd[[pathway]] <- unlist ( str_split(pathway.target[pathway,"pathway.cpd"], ";") )

}

lapply(names(pathway.cpd), function(x){
  sig.cpd.all[sig.cpd.all$Matched.Compound %in% pathway.cpd[[x]],]->df
  df$Pathway<-(x)
  df
  })->Pathway_matched_sig
Pathway_matched_sig<-do.call(rbind,Pathway_matched_sig)
```

``` r
outwd= paste0(project_dir, dim_reduction_dir, PCA_dir, ion_images_dir)
# Pathway_matched_sig_summary<-Pathway_matched_sig_summary[!duplicated(Pathway_matched_sig_summary[,c( "Matched.Compound", "Matched.Form","Dim" )]),]
library(magick)
library(flextable)
library(ggplot2)
library(grid)


Pathway_matched_sig_summary<-Pathway_matched_sig
Pathway_matched_sig_summary<-Pathway_matched_sig_summary[!duplicated(Pathway_matched_sig_summary[,c( "Matched.Compound", "Matched.Form","Dim" )]),]


Pathway_matched_sig_summary$CPD_name<-cpd.info$cpd.name[match(Pathway_matched_sig_summary$Matched.Compound,cpd.info$id)]
Pathway_matched_sig_summary$CPD_name<-str_remove_all(Pathway_matched_sig_summary$CPD_name,"^Beta-|^Alpha-|^alpha-")
Pathway_matched_sig_summary$CPD_name<-str_replace_all(Pathway_matched_sig_summary$CPD_name,"yl","yl ")
Pathway_matched_sig_summary$CPD_name_width<-str_width(Pathway_matched_sig_summary$CPD_name)
Pathway_matched_sig_summary$CPD_name<-str_wrap(Pathway_matched_sig_summary$CPD_name,width = 15,whitespace_only = T)
Pathway_matched_sig_summary$Isotype<-NA
Pathway_matched_sig_summary$Isotype[str_detect(Pathway_matched_sig_summary$Matched.Form,"13C6") ]<-"SIL_C6"
Pathway_matched_sig_summary$Isotype[!str_detect(Pathway_matched_sig_summary$Matched.Form,"13C6") ]<-"Normal"

for (pathway in unique(Pathway_matched_sig_summary$Pathway)){
  Pathway_matched_sig_summary_subset <- Pathway_matched_sig_summary[Pathway_matched_sig_summary$Pathway==pathway,]
   # Pathway_matched_sig_summary_subset_subset<-Pathway_matched_sig_summary_subset_subset %>% group_by(Query.Mass) %>% summarize(CPD_name=paste(CPD_name,collapse  = "\n\n"),Isotype=paste(Isotype,collapse  = "\n\n"))
  
  Pathway_matched_sig_summary_subset_subset <- Pathway_matched_sig_summary_subset %>% group_by(Query.Mass) %>% summarize(CPD_name = paste(unique(unlist(CPD_name)), collapse = "; "), Isotype = paste(unique(unlist(Isotype)),collapse = "\n\n"), Dim = paste(unique(Dim), collapse = ","))
  Pathway_matched_sig_summary_subset_subset$Dim<-str_remove_all(Pathway_matched_sig_summary_subset_subset$Dim,"Dim.")
 img_set<-magick::image_read(paste0(outwd,"/", Pathway_matched_sig_summary_subset_subset$Query.Mass,".png"))
 
 
 for (Cpd in 1:nrow(Pathway_matched_sig_summary_subset_subset)){
      img_info<-image_info(img_set)
      CPD_name=unique(unlist(str_split(Pathway_matched_sig_summary_subset_subset$CPD_name[Cpd],"\n\n")))
      
      Isotype=unique(unlist(str_split(Pathway_matched_sig_summary_subset_subset$Isotype[Cpd],"\n\n")))
      Dim = Pathway_matched_sig_summary_subset_subset$Dim[Cpd]
      tb<-data.frame(CPD_name=CPD_name,Isotype=Isotype, Dim = Dim)
      #colnames(tb)<-c("CPD_name","Isotype", "Dim")
      tb<-as_tibble(tb)
      set_flextable_defaults(font.size =10)
  #     ft_raster_cpd <- tb %>%
  # as.data.frame(.) %>% flextable::flextable() %>% flextable::delete_part(part = "header") %>% 
  #       bg( bg = "transparent", part = "all") %>%
  #       #flextable::delete_part(part = "header")  %>%
  #       flextable::color(color = "white") %>% border_remove %>%
  #       #flextable::delete_part(part = "footer")  %>%  
  #       set_table_properties(width = 1, layout = "autofit") 
  #     # ft_raster_cpd[["body"]][["styles"]][["cells"]][["background.color"]]->ft_raster_cpd[["body"]][["styles"]][["cells"]][["background"]]
  # 
  #    ft_raster_cpd<-ft_raster_cpd %>% 
  #       flextable::as_raster()
  #     ft_raster_iso <- tb[,"Isotype"] %>% flextable::flextable() %>%
  #       set_table_properties(width = 1, layout = "autofit")%>%
  #       #width(j=1, width = 550) %>% width(j=2, width = 100) %>%
  #       #flextable::delete_part(part = "header") %>%  
  #       flextable::color(color = "white") %>% border_remove %>%
  #       #flextable::delete_part(part = "footer") %>% 
  #       flextable::as_raster()
  #     
  #       dim_df<-str_remove_all(tb[1,"Dim"],"Dim.")
  #     ft_raster_dim <- as_tibble(dim_df) %>% flextable::flextable() %>%
  #       set_table_properties(width = 1, layout = "autofit")%>%
  #       #width(j=1, width = 550) %>% width(j=2, width = 100) %>%
  #       #flextable::delete_part(part = "header") %>%  
  #       flextable::color(color = "white") %>% border_remove %>% 
  #       #flextable::delete_part(part = "footer") %>% 
  #       flextable::as_raster()
     
      
      # legend <- image_blank(width = img_info$width[Cpd]/1.5, height=img_info$height[Cpd], color="black") 
      # legend <- image_ggplot(legend) + theme_void() + 
      #   annotation_custom(rasterGrob(ft_raster_cpd), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
      # legend1 <- image_graph(width = img_info$width[Cpd]/1.5, height=img_info$height[Cpd], bg="black")
      # print(legend)
      # dev.off()
      # 
      # legend <- image_blank(width = img_info$width[Cpd]/5, height=img_info$height[Cpd], color="black")
      # legend <- image_ggplot(legend) + theme_void() + 
      #   annotation_custom(rasterGrob(ft_raster_iso), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
      # legend2 <- image_graph(width = img_info$width[Cpd]/5, height=img_info$height[Cpd], bg="black"  )
      # print(legend)
      # dev.off()
      # 
      # 
      #  legend <- image_blank(width = img_info$width[Cpd]/5, height=img_info$height[Cpd], color="black")
      # legend <- image_ggplot(legend) + theme_void() + 
      #   annotation_custom(rasterGrob(ft_raster_dim), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
      # legend3 <- image_graph(width = img_info$width[Cpd]/5, height=img_info$height[Cpd], bg="black"  )
      # print(legend)
      # dev.off()
      
      
      #img_set[Cpd]<-image_append(c(img_set[Cpd],legend1,legend2,legend3),stack = F)
 }
    img_set_append<-image_append(img_set,stack = T)
    img_set_info<-image_info(img_set_append)
    tb<-as_tibble(Pathway_matched_sig_summary_subset_subset)
      set_flextable_defaults(font.size =10)
      ft_raster_cpd <- tb %>%
  as.data.frame(.) %>% flextable::flextable() %>% flextable::delete_part(part = "header") %>% 
        bg( bg = "transparent", part = "all") %>%
        #flextable::delete_part(part = "header")  %>%
        flextable::color(color = "white") %>% border_remove %>%
        #set_table_properties(layout = "autofit") %>%
        #flextable::delete_part(part = "footer")  %>%  
        height_all(height = 1) %>%
        hrule(rule = "exact")
      # ft_raster_cpd[["body"]][["styles"]][["cells"]][["background.color"]]->ft_raster_cpd[["body"]][["styles"]][["cells"]][["background"]]
  
     ft_raster_cpd<-ft_raster_cpd %>% 
        flextable::as_raster()
      legend <- image_blank(width = 1500, height=img_set_info$height, color="black")
      legend <- image_ggplot(legend) + theme_void() +
        annotation_custom(rasterGrob(ft_raster_cpd), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
      legend1 <- image_graph(width = 1500, height=img_set_info$height, bg="black")
      print(legend)
      dev.off()

    # header <- image_blank(width = img_set_info$width, height=img_set_info$width/10/3,color = "black")
    # header<-image_annotate(header,paste(y,z),gravity = "West",size=img_set_info$width/10/4,color = "white")
    img_set_append_anno<-image_append(c(img_set_append,legend1),stack = F)
    image_write(img_set_append_anno,paste0(outwd,pathway,"FDR_0.1_remapping.png"))
    
    
      tb<-as_tibble(Pathway_matched_sig_summary_subset_subset[,c("CPD_name","Isotype", "Dim")])
      set_flextable_defaults(font.size =10)
      ft_raster_cpd <- tb %>%
  as.data.frame(.) %>% flextable::flextable() %>% flextable::delete_part(part = "header") %>% 
        bg( bg = "transparent", part = "all") %>%
        # width(j = 1,width = 1.1) %>%
        # width(j = 2,width = 0.2) %>%
        # width(j = 3,width = 0.5) %>%
        #flextable::delete_part(part = "header")  %>%
        flextable::color(color = "white") %>% border_remove %>%
        #set_table_properties(layout = "autofit") %>%
        #flextable::delete_part(part = "footer")  %>%  
        height_all(height = 1) %>%
        hrule(rule = "exact")
      # ft_raster_cpd[["body"]][["styles"]][["cells"]][["background.color"]]->ft_raster_cpd[["body"]][["styles"]][["cells"]][["background"]]
  
     ft_raster_cpd<-ft_raster_cpd %>% 
        flextable::as_raster()
      legend <- image_blank(width = 1200, height=img_set_info$height, color="black")
      legend <- image_ggplot(legend) + theme_void() +
        annotation_custom(rasterGrob(ft_raster_cpd), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
      legend1 <- image_graph(width = 1200 , height=img_set_info$height, bg="black")
      print(legend)
      dev.off()

    # header <- image_blank(width = img_set_info$width, height=img_set_info$width/10/3,color = "black")
    # header<-image_annotate(header,paste(y,z),gravity = "West",size=img_set_info$width/10/4,color = "white")
    img_set_append_anno<-image_append(c(img_set_append,legend1),stack = F)
    image_write(img_set_append_anno,paste0(outwd,pathway,"FDR_0.1_remapping_style2.png"))
  }
```

# 12. Shrunk data

``` r
dir.create(paste0(project_dir, shrink_dir))

library(Cardinal)
all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir, "all_coordata_label_sel_arc_bind2.csv"),row.names = 1)

final_selection<- read.csv(paste0(project_dir, shrink_dir, "final_selection_shrink.csv"))


run(example_imdata)<-" "
lightmode()

fig.arbitrary_seg<- image(example_imdata, factor(final_selection$label[run(combinedimdata) == "4hrsilglucose_2nd_20201109_trim"] ) ~ x * y,
                          superpose=F,
                          key=F, xlab= NULL, col = c("#ff0027","#0092ff", "#800080",  "grey80","#20b2aa"))


fig.arbitrary_seg<-fig.arbitrary_seg

fig.arbitrary_seg[["par"]][["ann"]]=F
fig.arbitrary_seg[["par"]][["bty"]]="n"
fig.arbitrary_seg[["par"]][["pty"]]="s"
fig.arbitrary_seg[["par"]][["xaxt"]]="n"
fig.arbitrary_seg[["par"]][["yaxt"]]="n"
fig.arbitrary_seg[["par"]][["fg"]]="white"
fig.arbitrary_seg[["par"]][["oma"]]=c(0, 0, 0, 0)


png(paste0(project_dir, shrink_dir,"arbitrary_segmentation_shrink.png"),height = 3,width = 3,units = "in",res = 600)

print(fig.arbitrary_seg)
dev.off()
```

## 12.1 PCA visualisation

``` r
library(egg)
library(rsvd)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(Cardinal)

load(paste0(project_dir, dim_reduction_dir , PCA_dir,"Time_merged_PCA_all_bio.rda"))

all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir ,"all_coordata_label_sel_arc_bind2.csv"),row.names = 1)
arbitrary_seg <-read.csv( paste0(project_dir, shrink_dir,"/final_selection_shrink.csv"))

all_coordata_label_sel_arc_bind<-cbind(all_coordata_label_sel_arc_bind,arbitrary_seg= arbitrary_seg$label)


example_data <-all_coordata_label_sel_arc_bind[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim",]



dims.pca.full <- PCA[["ind"]][["coord"]]
dims.pca.full <- as.data.frame(dims.pca.full)
dims.pca = dims.pca.full

dims.pca<- dims.pca[which(as.numeric(rownames(dims.pca))%in% example_data$new_label_scale) ,]

dims.pca$region.final<-example_data$arbitrary_seg[match(as.numeric(rownames(dims.pca)), example_data$new_label_scale)] 

dims.pca.full$region.final<- dims.pca$region.final[match(rownames(dims.pca.full), rownames(dims.pca) )] 

dims.pca.full$region.final [is.na(dims.pca.full$region.final)] ="NS"


dims.pca.full$region.final<-factor(dims.pca.full$region.final, levels =c("Anterior",  "Core", "Equator","Posterior", "NS"))

dims.pca.full<-dims.pca.full[order(dims.pca.full$region.final, decreasing = T),]


png(paste0(project_dir, dim_reduction_dir, PCA_dir, "/PCA_Dim.1_vs_Dim.2.png"), width = 3, height = 3, units = "in", res=600)

fig<- ggplot(dims.pca.full, aes(x= Dim.1, y= Dim.2))+xlim(-200, 450)+ylim(-200, 100)+
      geom_point(data = dims.pca.full[dims.pca.full$region.final=="NS", ],
                  stat = "identity", position = "identity", show.legend = FALSE,
                  colour = "grey80", size = 1
                  )+
      theme_void()+ theme(legend.position = "none")



# fig<- fig + geom_point(data = dims.pca.full[!(dims.pca.full$region.final=="NS"), ],
#                   stat = "identity", position = "identity", show.legend = FALSE,
#                   size = 1, alpha = 0.6, aes(color=region.final)
#                   )+ scale_colour_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
#       theme_void()

fig<- fig + geom_point(data = dims.pca.full[!(dims.pca.full$region.final=="NS"), ],
                  stat = "identity", position = "identity", show.legend = FALSE,
                  size = 1, alpha = .6,aes(color=region.final)
                  )+ scale_colour_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
                  stat_ellipse(data = dims.pca.full[! (dims.pca.full$region.final =="NS"),], geom = "polygon",
                  type = "t", linetype = "twodash", aes(color = region.final, fill = region.final), level = 0.75, size = 0.5, alpha = 0.2) + scale_fill_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
      theme_void()+ theme(legend.position = "none")
 
print(fig)

dev.off()


png(paste0(project_dir, dim_reduction_dir, PCA_dir, "/PCA_Dim.4_vs_Dim.5.png"), width = 3, height = 3, units = "in", res=600)

fig<- ggplot(dims.pca.full, aes(x= Dim.4, y= Dim.5))+xlim(-100, 150)+ylim(-200, 200)+
      geom_point(data = dims.pca.full[dims.pca.full$region.final=="NS", ],
                  stat = "identity", position = "identity", show.legend = FALSE,
                  colour = "grey80", size = 1
                  )+
      theme_void()+ theme(legend.position = "none")



fig<- fig + geom_point(data = dims.pca.full[!(dims.pca.full$region.final=="NS"), ],
                  stat = "identity", position = "identity", show.legend = FALSE,
                  size = 1, alpha = .6,aes(color=region.final)
                  )+ scale_colour_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
                  stat_ellipse(data = dims.pca.full[! (dims.pca.full$region.final =="NS"),], geom = "polygon",
                  type = "t", linetype = "twodash", aes(color = region.final, fill = region.final), level = 0.75, size = 0.5, alpha = 0.2) +
                  scale_fill_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
      theme_void()+ theme(legend.position = "none")
print(fig)

dev.off()
```

## 12.2 PCA K means

``` r
#Dim.4 vs Dim.5

set.seed(2)
pca.kmean.res<- kmeans(dims.pca.full[,c("Dim.4", "Dim.5")], 6, nstart = 25)
Time_merged_pca<- data.frame(dims.pca.full$Dim.4, dims.pca.full$Dim.5, rownames(dims.pca.full))
colnames(Time_merged_pca)<-c("Dim.4","Dim.5","bin")

pca.kmean.res.cluster<-data.frame(pca.kmean.res$cluster)

index<- match(Time_merged_pca$bin, as.numeric(rownames(pca.kmean.res.cluster)))
Time_merged_pca$K_Label<- pca.kmean.res.cluster[,1][index]

fig<- ggplot(Time_merged_pca, aes(x= Dim.4, y= Dim.5)) + xlim(-100, 150) + ylim(-200, 200) +
        geom_point( stat = "identity", position = "identity", show.legend = FALSE,
                  size = 1, alpha=.6,aes(color = as.factor(K_Label)) )+
        scale_colour_manual(values=  c("#760F00","#4e639e","#e54616", "#dba053",  "#FF997C", "#7FBFDD"))+theme_void()

png(paste0(project_dir,dim_reduction_dir, PCA_dir, "/Scatter_K_means_Dim.4_vs_Dim.5.png"), width = 3, height = 3, units = "in", res=300)
print(fig)
dev.off()


all_coordata_label_sel_arc_bind$K_Label<-Time_merged_pca$K_Label[match(as.character(all_coordata_label_sel_arc_bind$new_label_scale),as.character(Time_merged_pca$bin))]
factor(all_coordata_label_sel_arc_bind$K_Label[which(all_coordata_label_sel_arc_bind$run =="4hrsilglucose_2nd_20201109_trim" )])->final_selection

lightmode()
fig.arbitrary_seg<- image(example_imdata, factor(final_selection) ~ x * y,superpose=F, key=F, col = c("#760F00","#4e639e","#e54616", "#dba053",  "#FF997C", "#7FBFDD"))


fig.arbitrary_seg[["par"]][["ann"]]=F
fig.arbitrary_seg[["par"]][["bty"]]="n"
fig.arbitrary_seg[["par"]][["pty"]]="s"
fig.arbitrary_seg[["par"]][["xaxt"]]="n"
fig.arbitrary_seg[["par"]][["yaxt"]]="n"
fig.arbitrary_seg[["par"]][["fg"]]="white"
fig.arbitrary_seg[["par"]][["oma"]]=c(0, 0, 0, 0)

png(paste0(project_dir,dim_reduction_dir, PCA_dir, "/Lens_K_means_Dim.4_vs_Dim.5.png"), width = 3, height = 3, units = "in", res=600)
print(fig.arbitrary_seg)
dev.off()

#Dim.1 vs Dim.2

set.seed(2)
pca.kmean.res<- kmeans(dims.pca.full[,c("Dim.1", "Dim.2")], 6, nstart = 25)
Time_merged_pca<- data.frame(dims.pca.full$Dim.1, dims.pca.full$Dim.2, rownames(dims.pca.full))
colnames(Time_merged_pca)<-c("Dim.1","Dim.2","bin")

pca.kmean.res.cluster<-data.frame(pca.kmean.res$cluster)

index<- match(Time_merged_pca$bin, as.numeric(rownames(pca.kmean.res.cluster)))
Time_merged_pca$K_Label<- pca.kmean.res.cluster[,1][index]

fig<- ggplot(Time_merged_pca, aes(x= Dim.1, y= Dim.2)) + xlim(-200, 450) + ylim(-200, 100) +
        geom_point( stat = "identity", position = "identity", show.legend = FALSE,
                  size = 1, alpha=.6, aes(color = as.factor(K_Label)) )+
        scale_colour_manual(values=  c("#760F00","#4e639e","#e54616",  "#FF997C","#dba053",  "#7FBFDD"))+theme_void()

png(paste0(project_dir,dim_reduction_dir, PCA_dir, "/Scatter_K_means_Dim.1_vs_Dim.2.png"), width = 3, height = 3, units = "in", res=300)
print(fig)
dev.off()


all_coordata_label_sel_arc_bind$K_Label<-Time_merged_pca$K_Label[match(as.character(all_coordata_label_sel_arc_bind$new_label_scale),as.character(Time_merged_pca$bin))]
factor(all_coordata_label_sel_arc_bind$K_Label[which(all_coordata_label_sel_arc_bind$run =="4hrsilglucose_2nd_20201109_trim" )])->final_selection

lightmode()
fig.arbitrary_seg<- image(example_imdata, factor(final_selection) ~ x * y,superpose=F, key=F, col = c("#760F00","#4e639e","#e54616",   "#FF997C","#dba053", "#7FBFDD"))


fig.arbitrary_seg[["par"]][["ann"]]=F
fig.arbitrary_seg[["par"]][["bty"]]="n"
fig.arbitrary_seg[["par"]][["pty"]]="s"
fig.arbitrary_seg[["par"]][["xaxt"]]="n"
fig.arbitrary_seg[["par"]][["yaxt"]]="n"
fig.arbitrary_seg[["par"]][["fg"]]="white"
fig.arbitrary_seg[["par"]][["oma"]]=c(0, 0, 0, 0)

png(paste0(project_dir,dim_reduction_dir, PCA_dir, "/Lens_K_means_Dim.1_vs_Dim.2.png"), width =3, height = 3, units = "in", res=600)
print(fig.arbitrary_seg)
dev.off()
```

## 12.3 UMAP K means

``` r
library(rsvd)
library(dplyr)
library(ggplot2)

Time_merged_umap<- readRDS(paste0(project_dir, dim_reduction_dir, UMAP_dir, "/Time_merged_umap_all_bio.rds"))
colnames(Time_merged_umap)<- c("UMAP_1","UMAP_2","bin")
Time_merged_umap_full <- Time_merged_umap

all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir ,"all_coordata_label_sel_arc_bind2.csv"),row.names = 1)
arbitrary_seg <-read.csv( paste0(project_dir, shrink_dir,"/final_selection_shrink.csv"))
all_coordata_label_sel_arc_bind<-cbind(all_coordata_label_sel_arc_bind,arbitrary_seg= arbitrary_seg$label)

example_data <-all_coordata_label_sel_arc_bind[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim",]

Time_merged_umap<-Time_merged_umap[Time_merged_umap$bin %in% all_coordata_label_sel_arc_bind$new_label_scale[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim"],]

Time_merged_umap$Region <- example_data$arbitrary_seg[match(as.numeric(Time_merged_umap$bin), as.numeric(example_data$new_label_scale)) ]

Time_merged_umap_full$Region <- Time_merged_umap$Region[match(Time_merged_umap_full$bin, Time_merged_umap$bin)]
Time_merged_umap_full$Region[is.na(Time_merged_umap_full$Region)] <- "NS"

Time_merged_umap$Region<- factor(Time_merged_umap$Region, levels =c("Anterior", "Core", "Equator", "Posterior", "NS") )

Time_merged_umap_full$Region<- factor(Time_merged_umap_full$Region, levels =c("Anterior", "Core", "Equator", "Posterior", "NS") )


png(paste0(project_dir, dim_reduction_dir, UMAP_dir,"//umap_arb_shrink.png"),height = 3,width = 3,units = "in",res = 600)

fig<- ggplot(Time_merged_umap_full, aes(x= UMAP_2, y= UMAP_1))+ 
      geom_point(data = Time_merged_umap_full[Time_merged_umap_full$Region =="NS", ],
                  stat = "identity", position = "identity", show.legend = FALSE,
                  colour = "grey80", size = 1
                  )+
      theme_void()



fig<- fig + geom_point(data = Time_merged_umap_full[! (Time_merged_umap_full$Region =="NS"),],
                  stat = "identity", position = "identity", show.legend = FALSE,
                  size = 1, alpha = .6,aes(color=Region)
                  )+ scale_colour_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
                  stat_ellipse(data = Time_merged_umap_full[! (Time_merged_umap_full$Region =="NS"),], 
                  type = "t", linetype = "twodash", geom = "polygon", aes(color = Region, fill = Region), level = 0.75, size = 0.5, alpha = 0.2) + scale_fill_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
                  theme_void()+ theme(legend.position = "none")
fig  
  #      theme_dr() +
  # theme(panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank())

# 
# fig<- fig +  
#   geom_segment(aes(x = min(Time_merged_umap$UMAP_1)-1, y = min(Time_merged_umap$UMAP_2)-2,
#                   xend = min(Time_merged_umap$UMAP_1) + 2, 
#                    yend = min(Time_merged_umap$UMAP_2)),
#                     colour = "black", linewidth = 0.6, arrow = arrow(length = unit(0.2, "cm"))) +
#   geom_segment(aes(x = min(Time_merged_umap$UMAP_1)-1, y = min(Time_merged_umap$UMAP_2)-2,
#                   xend = min(Time_merged_umap$UMAP_1) , 
#                    yend = min(Time_merged_umap$UMAP_2)+ 2),
#                     colour = "black", linewidth = 0.6, arrow = arrow(length = unit(0.2, "cm"))) + 
#     annotate("text", x = min(Time_merged_umap$UMAP_1) + 1, y = min(Time_merged_umap$UMAP_2) -1,
#            label = "UMAP_1", color = "black", size = 3, fontface = "bold") + 
#     annotate("text", x = min(Time_merged_umap$UMAP_1) -1, y = min(Time_merged_umap$UMAP_2) +1,
#            label = "UMAP_2", color = "black", size = 3, fontface = "bold")



print(fig)
dev.off()
```

``` r
set.seed(9)
umap.kmean.res<- kmeans(Time_merged_umap_full[,c("UMAP_1", "UMAP_2")], 6, nstart = 25)

Time_merged_umap<- data.frame(UMAP_1 = Time_merged_umap_full$UMAP_1, UMAP_2 = Time_merged_umap_full$UMAP_2, bin = Time_merged_umap_full$bin, K_Label=umap.kmean.res$cluster)

fig<- ggplot(Time_merged_umap, aes(x= UMAP_2, y= UMAP_1)) +
        geom_point( stat = "identity", position = "identity", show.legend = FALSE,
                  size = 1, alpha=.6,aes(color = as.factor(K_Label)) )+
        scale_colour_manual(values=  c("#dba053","#4e639e","#760F00", "#7FBFDD","#e54616",   "#FF997C"))+theme_void()

png(paste0(project_dir,dim_reduction_dir, UMAP_dir, "/Scatter_K_means_UMAP.png"), width = 3, height = 3, units = "in", res=600)
print(fig)
dev.off()


all_coordata_label_sel_arc_bind$K_Label<-Time_merged_umap$K_Label[match(as.character(all_coordata_label_sel_arc_bind$new_label_scale),as.character(Time_merged_umap$bin))]
factor(all_coordata_label_sel_arc_bind$K_Label[which(all_coordata_label_sel_arc_bind$run =="4hrsilglucose_2nd_20201109_trim" )])->final_selection

lightmode()
fig.arbitrary_seg<- image(example_imdata, factor(final_selection) ~ x * y,superpose=F, key=F, col = c("#dba053","#4e639e","#760F00", "#7FBFDD","#e54616",   "#FF997C"))


fig.arbitrary_seg[["par"]][["ann"]]=F
fig.arbitrary_seg[["par"]][["bty"]]="n"
fig.arbitrary_seg[["par"]][["pty"]]="s"
fig.arbitrary_seg[["par"]][["xaxt"]]="n"
fig.arbitrary_seg[["par"]][["yaxt"]]="n"
fig.arbitrary_seg[["par"]][["fg"]]="white"
fig.arbitrary_seg[["par"]][["oma"]]=c(0, 0, 0, 0)

png(paste0(project_dir,dim_reduction_dir, UMAP_dir, "/Lens_K_means_UMAP.png"), width = 3, height = 3, units = "in", res=600)
print(fig.arbitrary_seg)
dev.off()
```

## 12.4 t-SNE K means

``` r
library(rsvd)
library(dplyr)
library(ggplot2)
Time_merged_tsne<- readRDS(paste0(project_dir, dim_reduction_dir, tSNE_dir, "/Time_merged_tsne_all_bio.rds"))
colnames(Time_merged_tsne)<- c("tSNE_1","tSNE_2","bin")
Time_merged_tsne_full <- Time_merged_tsne

all_coordata_label_sel_arc_bind<-read.csv(paste0(project_dir, spatial_alignment_dir ,"all_coordata_label_sel_arc_bind2.csv"),row.names = 1)
arbitrary_seg <-read.csv( paste0(project_dir, shrink_dir,"/final_selection_shrink.csv"))
all_coordata_label_sel_arc_bind<-cbind(all_coordata_label_sel_arc_bind,arbitrary_seg= arbitrary_seg$label)

example_data <-all_coordata_label_sel_arc_bind[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim",]

Time_merged_tsne<-Time_merged_tsne[Time_merged_tsne$bin %in% all_coordata_label_sel_arc_bind$new_label_scale[all_coordata_label_sel_arc_bind$run == "4hrsilglucose_2nd_20201109_trim"],]

Time_merged_tsne$Region <- example_data$arbitrary_seg[match(as.numeric(Time_merged_tsne$bin), as.numeric(example_data$new_label_scale)) ]

Time_merged_tsne_full$Region <- Time_merged_tsne$Region[match(Time_merged_tsne_full$bin, Time_merged_tsne$bin)]
Time_merged_tsne_full$Region[is.na(Time_merged_tsne_full$Region)] <- "NS"

Time_merged_tsne$Region<- factor(Time_merged_tsne$Region, levels =c("Anterior", "Core", "Equator", "Posterior", "NS") )

Time_merged_tsne_full$Region<- factor(Time_merged_tsne_full$Region, levels =c("Anterior", "Core", "Equator", "Posterior", "NS") )


png(paste0(project_dir, dim_reduction_dir, tSNE_dir,"//tSNE_arb_shrink.png"),height = 3,width = 3,units = "in",res = 600)

fig<- ggplot(Time_merged_tsne_full, aes(x= tSNE_1, y= tSNE_2))+ 
      geom_point(data = Time_merged_tsne_full[Time_merged_tsne_full$Region =="NS", ],
                  stat = "identity", position = "identity", show.legend = FALSE,
                  colour = "grey80", size = 1
                  )+
      theme_void()



# fig<- fig + geom_point(data = Time_merged_tsne_full[! (Time_merged_tsne_full$Region =="NS"),],
#                   stat = "identity", position = "identity", show.legend = FALSE,
#                   size = 1, alpha = 0.6,aes(color=Region)
#                   )+ scale_colour_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
#                   stat_ellipse(data = Time_merged_tsne_full[! (Time_merged_tsne_full$Region =="NS"),], 
#                   type = "t", linetype = "twodash", aes(color = Region, fill = Region), level = 0.75, size = 0.5, alpha = 0.2) + 
#                   scale_fill_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
#                   theme_void()+ theme(legend.position = "none")



fig<- fig + geom_point(data = Time_merged_tsne_full[! (Time_merged_tsne_full$Region =="NS"),],
                  stat = "identity", position = "identity", show.legend = FALSE,
                  size = 1, alpha = .6,aes(color=Region)
                  )+ scale_colour_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
                  stat_ellipse(data = Time_merged_tsne_full[! (Time_merged_tsne_full$Region =="NS"),], 
                  type = "t", linetype = "twodash", geom = "polygon", aes(color = Region, fill = Region), level = 0.75, size = 0.5, alpha = 0.2) + scale_fill_manual(values=  c("#ff0027","#0092ff", "#800080","#20b2aa"))+
                  theme_void()+ theme(legend.position = "none")
fig


print(fig)
dev.off()
```

``` r
set.seed(9)
tsne.kmean.res<- kmeans(Time_merged_tsne_full[,c("tSNE_1", "tSNE_2")], 6, nstart = 25)

Time_merged_tsne<- data.frame(tSNE_1 = Time_merged_tsne_full$tSNE_1, tSNE_2 = Time_merged_tsne_full$tSNE_2, bin = Time_merged_tsne_full$bin, K_Label = tsne.kmean.res$cluster)

fig<- ggplot(Time_merged_tsne, aes(x= tSNE_1, y= tSNE_2)) + 
        geom_point( stat = "identity", position = "identity", show.legend = FALSE,
                  size = 1, alpha = 0.6, aes(color = as.factor(K_Label)) )+
        scale_colour_manual(values=  c("#760F00", "#dba053","#4e639e","#e54616", "#7FBFDD",  "#FF997C"))+theme_void()

png(paste0(project_dir,dim_reduction_dir, tSNE_dir, "/Scatter_K_means_tSNE.png"), width = 3, height = 3, units = "in", res=600)
print(fig)
dev.off()


all_coordata_label_sel_arc_bind$K_Label<-Time_merged_tsne$K_Label[match(as.character(all_coordata_label_sel_arc_bind$new_label_scale),as.character(Time_merged_tsne$bin))]
factor(all_coordata_label_sel_arc_bind$K_Label[which(all_coordata_label_sel_arc_bind$run =="4hrsilglucose_2nd_20201109_trim" )])->final_selection

lightmode()
fig.arbitrary_seg<- image(example_imdata, factor(final_selection) ~ x * y,superpose=F, key=F, col = c("#760F00", "#dba053","#4e639e","#e54616", "#7FBFDD",  "#FF997C"))


fig.arbitrary_seg[["par"]][["ann"]]=F
fig.arbitrary_seg[["par"]][["bty"]]="n"
fig.arbitrary_seg[["par"]][["pty"]]="s"
fig.arbitrary_seg[["par"]][["xaxt"]]="n"
fig.arbitrary_seg[["par"]][["yaxt"]]="n"
fig.arbitrary_seg[["par"]][["fg"]]="white"
fig.arbitrary_seg[["par"]][["oma"]]=c(0, 0, 0, 0)

png(paste0(project_dir,dim_reduction_dir, tSNE_dir, "/Lens_K_means_tSNE.png"), width = 3, height = 3, units = "in", res=600)
print(fig.arbitrary_seg)
dev.off()
```
