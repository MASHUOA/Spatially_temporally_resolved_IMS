
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


library(MetaboAnalystR)
library(Cardinal)

#change to your own project directory here
project_dir <- "~/data/"

Raw_data_dir <- "/raw/"

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


# test Principle conponent analysis and top loading feature selection

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



