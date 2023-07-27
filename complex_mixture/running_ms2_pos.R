rm(list = ls())
dts <- read.csv("/data/ms2_benchmark/MTBLS2207/41592_2021_1195_MOESM3_ESM_allcmpds.csv")
dt_pos <- dts[dts$Charge == 1, ]
library(OptiLCMS2ID)
#source("~/Github/OptiLCMS2ID/R/MSn_processing.R")
dt_ft <- cbind(dt_pos$m.z - 0.005, dt_pos$m.z + 0.005, dt_pos$RT.Start*60-25, dt_pos$RT.End*60+25)
mSet <- PerformMSnImport(targetFeatures = dt_ft, 
                         filesPath = "/data/ms2_benchmark/MTBLS2207/IROA/DDA/IROA_P1-6_ddMS2_pos_1Da.mzML",
                         acquisitionMode = "DDA")

mSet <- PerformDDADeconvolution(mSet, 
                                ppm1 = 10, 
                                ppm2 = 25,
                                sn = 12,
                                filtering = 2000,
                                window_size = 1,
                                intensity_thresh = 2.5e4,
                                database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
                                ncores = 8L)

mSet <- PerformSpectrumConsenus (mSet, 
                                 ppm2 = 25,
                                 concensus_fraction = 0.5,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)

mSet <- PerformDBSearchingBatch (mSet, 
                                 ppm1 = 10,
                                 ppm2 = 25,
                                 rt_tol = 5,
                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
                                 use_rt = FALSE,
                                 enableNL = TRUE,
                                 ncores = 8L)

mSet <- PerformResultsExport (mSet,
                              type = 3L,
                              topN = 10L, 
                              ncores = 8L)

save(mSet, file = "/data/ms2_benchmark/MTBLS2207/IROA/optilcms_pos_dda/mSet_res_Deco.rda")

inchik_map <- read.csv("/data/ms2_benchmark/MTBLS2207/IROA/data_pos_inchikeys.txt", header = F)
k<- mSet@MSnResults[["DecRes"]][[1]][["FeatureIdx"]] + 1
r<- c(1:length(mSet@MSnResults[["DecRes"]][[1]][["FeatureIdx"]]))
annote_res <- mSet@MSnResults[["DBAnnoteRes"]]
kc = 0;
all_score <- vector()
all_dots <- vector()
u_vec <- vector()
matched_inchis <- vector()

for(u in 1:length(r)){
  res_incks0 <- res_incks <- annote_res[[r[u]]]
  if(length(res_incks[[1]][[1]]) == 0){
    cat("u -> " , u ,"\n");
    next
  }
  res_incks <- sapply(res_incks, function(x){
    max_val <- max(x[[3]]);
    idx <- which(max_val == x[[3]])
    x[[2]][idx]
  })
  # res_dots <- sapply(res_incks0, function(x){
  #   max_val <- max(x[[3]]);
  #   idx <- which(max_val == x[[3]])
  #   x[[4]][idx]
  # })
  res_incks_val <- sapply(res_incks0, function(x){
    max_val <- max(x[[3]]);
    return(max_val)
  })
  # res_dots <- res_dots[res_incks != ""]
  res_incks_val <- (res_incks_val[res_incks != ""])
  res_incks <- (res_incks[res_incks != ""])
  
  res_incks <- (sapply(res_incks, function(x){x[1]}))
  # res_dot <- (res_dots[!is.na(res_dots)])[1]
  #res_incks_val <- (res_incks_val[res_incks != ""])[1]
  res_incks_val <- max(res_incks_val[!is.na(res_incks_val)])
  if(all(is.na(res_incks))){next;}
  map_incks <- inchik_map[which(inchik_map$V1 == dt_pos[k[u],1]),2]
  found_this <- FALSE;
  for(rs in res_incks){
    for(ms in map_incks){
      if((strsplit(rs, "-")[[1]][1] == strsplit(ms, "-")[[1]][1]) &
         (strsplit(rs, "-")[[1]][2] == strsplit(ms, "-")[[1]][2])) {
        matched_inchis <- c(matched_inchis, ms)
        kc = kc + 1;
        found_this <- TRUE;
        all_score <- c(all_score, res_incks_val)
        # all_dots <- c(all_dots, max(unlist(res_dot)))
        u_vec <- c(u_vec, u)
        break;
      }
    }
    if(found_this){break}
  }
  
  # if(!found_this){
  #   cat("NOT found_this -> ", u ,"\n")
  # }
  
}
print(kc) #160
kc/406
matched_inchis_pos_dda <- matched_inchis
save(matched_inchis_pos_dda, file = "../optilcms_pos_dda/matched_inchis_pos_dda.rda")
# idxx <- sapply(mSet@MSnResults[["DecRes"]][[1]][["Indicator"]][u_vec], function(x){any(x !=0)})
# all_score_deco <- all_score[idxx]
# all_dot_deco <- all_dots[idxx]
# 
# load("/data/ms2_benchmark/MTBLS2207/IROA/optilcms_pos_dda/mSet_res_noDeco_new.rda")
# 
# all_dot_noDeco <- sapply(u_vec, function(x) mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][[4]][1])[idxx]
# all_score_noDeco <- sapply(u_vec, function(x) mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][[3]][1])[idxx]

# save(all_score_deco, all_score_noDeco, file = "/data/ms2_benchmark/MTBLS2207/IROA/optilcms_pos_dda/all_score_raw_new.rda")

######### DDA testing -- nodeco

rm(list = ls())
dts <- read.csv("/data/ms2_benchmark/MTBLS2207/41592_2021_1195_MOESM3_ESM_allcmpds.csv")
dt_pos <- dts[dts$Charge == 1, ]
library(OptiLCMS2ID)
#source("~/Github/OptiLCMS2ID/R/MSn_processing.R")
dt_ft <- cbind(dt_pos$m.z - 0.005, dt_pos$m.z + 0.005, dt_pos$RT.Start*60-25, dt_pos$RT.End*60+25)
mSet <- PerformMSnImport(targetFeatures = dt_ft, 
                         filesPath = "/data/ms2_benchmark/MTBLS2207/IROA/DDA/IROA_P1-6_ddMS2_pos_1Da.mzML",
                         acquisitionMode = "DDA")

mSet <- PerformDDADeconvolution(mSet, 
                                ppm1 = 10, 
                                ppm2 = 25,
                                sn = 12,
                                filtering = 2000,
                                window_size = 1,
                                intensity_thresh = 2.5e4,
                                database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
                                ncores = 1L,
                                decoOn = FALSE)

mSet <- PerformSpectrumConsenus (mSet, 
                                 ppm2 = 25,
                                 concensus_fraction = 0.5,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)

mSet <- PerformDBSearchingBatch (mSet, 
                                 ppm1 = 10,
                                 ppm2 = 25,
                                 rt_tol = 5,
                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
                                 use_rt = FALSE,
                                 enableNL = TRUE,
                                 ncores = 8L)

mSet <- PerformResultsExport (mSet,
                              type = 3L,
                              topN = 10L, 
                              ncores = 8L)

#save(mSet, file = "/data/ms2_benchmark/MTBLS2207/IROA/optilcms_pos_dda/mSet_res_noDeco_new.rda")

inchik_map <- read.csv("/data/ms2_benchmark/MTBLS2207/IROA/data_pos_inchikeys.txt", header = F)
k<- mSet@MSnResults[["DecRes"]][[1]][["FeatureIdx"]] + 1
r<- c(1:length(mSet@MSnResults[["DecRes"]][[1]][["FeatureIdx"]]))
annote_res <- mSet@MSnResults[["DBAnnoteRes"]]
kc = 0;
matched_inchis <- vector()
all_score <- vector()
for(u in 1:length(r)){
  res_incks0 <- res_incks <- annote_res[[r[u]]]
  if(length(res_incks[[1]][[1]]) == 0){
    cat("u -> " , u ,"\n");
    next
  }
  res_incks <- sapply(res_incks, function(x){
    max_val <- max(x[[3]]);
    idx <- which(max_val == x[[3]])
    x[[2]][idx]
  })
  res_incks_val <- sapply(res_incks0, function(x){
    max_val <- max(x[[3]]);
    return(max_val)
  })
  res_incks <- unique(res_incks[res_incks != ""])[1]
  res_incks_val <- (res_incks_val[res_incks != ""])[1]
  if(is.na(res_incks)){next;}
  map_incks <- inchik_map[which(inchik_map$V1 == dt_pos[k[u],1]),2]
  for(rs in res_incks){
    for(ms in map_incks){
      if((strsplit(rs, "-")[[1]][1] == strsplit(ms, "-")[[1]][1]) &
         (strsplit(rs, "-")[[1]][2] == strsplit(ms, "-")[[1]][2])) {
        kc = kc + 1;
        matched_inchis <- c(matched_inchis, ms)
        all_score <- c(all_score, res_incks_val)
        break;
      }
    }
  }
}
print(kc) # find 163 features
kc/406
matched_inchis_noDeco_dda <- matched_inchis
save(matched_inchis_noDeco_dda, file = "../optilcms_pos_dda/matched_inchis_pos_noDeco_dda.rda")


########## DDA testing -- bioDB

rm(list = ls())
dts <- read.csv("/data/ms2_benchmark/MTBLS2207/41592_2021_1195_MOESM3_ESM_allcmpds.csv")
dt_pos <- dts[dts$Charge == 1, ]
library(OptiLCMS2ID)
#source("~/Github/OptiLCMS2ID/R/MSn_processing.R")
dt_ft <- cbind(dt_pos$m.z - 0.005, dt_pos$m.z + 0.005, dt_pos$RT.Start*60-25, dt_pos$RT.End*60+25)
mSet <- PerformMSnImport(targetFeatures = dt_ft, 
                         filesPath = "/data/ms2_benchmark/MTBLS2207/IROA/DDA/IROA_P1-6_ddMS2_pos_1Da.mzML",
                         acquisitionMode = "DDA")

mSet <- PerformDDADeconvolution(mSet, 
                                ppm1 = 10, 
                                ppm2 = 25,
                                sn = 12,
                                filtering = 2000,
                                window_size = 1,
                                intensity_thresh = 2.5e4,
                                database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID.sqlite",
                                ncores = 8L)

mSet <- PerformSpectrumConsenus (mSet, 
                                 ppm2 = 25,
                                 concensus_fraction = 0.5,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)

mSet <- PerformDBSearchingBatch (mSet, 
                                 ppm1 = 10,
                                 ppm2 = 25,
                                 rt_tol = 5,
                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID.sqlite",
                                 use_rt = FALSE,
                                 enableNL = TRUE,
                                 ncores = 8L)

mSet <- PerformResultsExport (mSet,
                              type = 3L,
                              topN = 10L, 
                              ncores = 8L)

save(mSet, file = "/data/ms2_benchmark/MTBLS2207/IROA/optilcms_pos_dda/mSet_res_bioDB.rda")

inchik_map <- read.csv("/data/ms2_benchmark/MTBLS2207/IROA/data_pos_inchikeys.txt", header = F)
k<- mSet@MSnResults[["DecRes"]][[1]][["FeatureIdx"]] + 1
r<- c(1:length(mSet@MSnResults[["DecRes"]][[1]][["FeatureIdx"]]))
annote_res <- mSet@MSnResults[["DBAnnoteRes"]]
kc = 0;
all_score <- vector()
all_dots <- vector()
u_vec <- vector()
matched_inchis <- vector()

for(u in 1:length(r)){
  res_incks0 <- res_incks <- annote_res[[r[u]]]
  if(length(res_incks[[1]][[1]]) == 0){
    cat("u -> " , u ,"\n");
    next
  }
  res_incks <- sapply(res_incks, function(x){
    max_val <- max(x[[3]]);
    idx <- which(max_val == x[[3]])
    x[[2]][idx]
  })
  # res_dots <- sapply(res_incks0, function(x){
  #   max_val <- max(x[[3]]);
  #   idx <- which(max_val == x[[3]])
  #   x[[4]][idx]
  # })
  res_incks_val <- sapply(res_incks0, function(x){
    max_val <- max(x[[3]]);
    return(max_val)
  })
  # res_dots <- res_dots[res_incks != ""]
  res_incks_val <- (res_incks_val[res_incks != ""])
  res_incks <- (res_incks[res_incks != ""])
  
  res_incks <- (sapply(res_incks, function(x){x[1]}))
  # res_dot <- (res_dots[!is.na(res_dots)])[1]
  #res_incks_val <- (res_incks_val[res_incks != ""])[1]
  res_incks_val <- max(res_incks_val[!is.na(res_incks_val)])
  if(all(is.na(res_incks))){next;}
  map_incks <- inchik_map[which(inchik_map$V1 == dt_pos[k[u],1]),2]
  found_this <- FALSE;
  for(rs in res_incks){
    for(ms in map_incks){
      if((strsplit(rs, "-")[[1]][1] == strsplit(ms, "-")[[1]][1]) &
         (strsplit(rs, "-")[[1]][2] == strsplit(ms, "-")[[1]][2])) {
        matched_inchis <- c(matched_inchis, ms)
        kc = kc + 1;
        found_this <- TRUE;
        all_score <- c(all_score, res_incks_val)
        # all_dots <- c(all_dots, max(unlist(res_dot)))
        u_vec <- c(u_vec, u)
        break;
      }
    }
    if(found_this){break}
  }
  
  # if(!found_this){
  #   cat("NOT found_this -> ", u ,"\n")
  # }
  
}
print(kc) #160
kc/406
matched_inchis_bioDB_dda <- matched_inchis
save(matched_inchis_bioDB_dda, file = "../optilcms_pos_dda/matched_inchis_bioDB_dda.rda")


######### DIA testing
rm(list = ls())

dts <- read.csv("/data/ms2_benchmark/MTBLS2207/41592_2021_1195_MOESM3_ESM_allcmpds.csv")
dt_pos <- dts[dts$Charge == 1, ]
library(OptiLCMS2ID)
#source("~/Github/OptiLCMS2ID/R/MSn_processing.R")
dt_ft <- cbind(dt_pos$m.z - 0.005, dt_pos$m.z + 0.005, dt_pos$RT.Start*60-25, dt_pos$RT.End*60+25)
mSet <- PerformMSnImport(targetFeatures = dt_ft, 
                         filesPath = "/data/ms2_benchmark/MTBLS2207/IROA/DIA/IROA_P1-6_DIA_test_pos1.mzML",
                         acquisitionMode = "DIA",
                         SWATH_file = "/data/ms2_benchmark/MTBLS2207/IROA/optilcms_neg_dia/DIA_SWATH_MS_experiment_file.txt")
mSet <- PerformDIADeconvolution (mSet, 
                                 min_width = 5,
                                 ppm2 = 25,
                                 sn = 12,
                                 span = 0.3,
                                 filtering = 2000, 
                                 ncores = 3L)

mSet <- PerformSpectrumConsenus (mSet, 
                                 ppm2 = 25,
                                 concensus_fraction = 0.5,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)

mSet <- PerformDBSearchingBatch (mSet, 
                                 ppm1 = 10,
                                 ppm2 = 25,
                                 rt_tol = 5,
                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
                                 use_rt = FALSE,
                                 enableNL = TRUE,
                                 ncores = 8L)

mSet <- PerformResultsExport (mSet,
                              type = 3L,
                              topN = 10L, 
                              ncores = 4L)

save(mSet, file = "/data/ms2_benchmark/MTBLS2207/IROA/optilcms_pos_dia/mSet_res_dia.rda")

inchik_map <- read.csv("/data/ms2_benchmark/MTBLS2207/IROA/data_pos_inchikeys.txt", header = F)
k<- mSet@MSnResults[["DecRes"]][[1]][["FeatureIdx"]] + 1
r<- c(1:length(mSet@MSnResults[["DecRes"]][[1]][["FeatureIdx"]]))
annote_res <- mSet@MSnResults[["DBAnnoteRes"]]
kc = 0;
all_score <- vector()
matched_inchis <- vector()
for(u in 1:length(r)){
  res_incks0 <- res_incks <- annote_res[[r[u]]]
  if(length(res_incks[[1]][[1]]) == 0){
    cat("u -> " , u ,"\n");
    next
  }
  res_incks <- sapply(res_incks, function(x){
    max_val <- max(x[[3]]);
    idx <- which(max_val == x[[3]])
    x[[2]][idx]
  })
  res_incks_val <- sapply(res_incks0, function(x){
    max_val <- max(x[[3]]);
    return(max_val)
  })
  
  res_incks <- (res_incks[res_incks != ""])
  res_incks_val <- (res_incks_val[res_incks != ""])
  res_incks <- (sapply(res_incks, function(x){x[1]}))
  #res_incks_val <- (res_incks_val[res_incks != ""])[1]
  if(all(is.na(res_incks))){next;}
  map_incks <- inchik_map[which(inchik_map$V1 == dt_pos[k[u],1]),2]
  found_this <- FALSE;
  for(rs in res_incks){
    for(ms in map_incks){
      if((strsplit(rs, "-")[[1]][1] == strsplit(ms, "-")[[1]][1]) &
         (strsplit(rs, "-")[[1]][2] == strsplit(ms, "-")[[1]][2])) {
        matched_inchis <- c(matched_inchis, ms)
        kc = kc + 1;
        found_this <- TRUE;
        break;
        #all_score <- c(all_score, res_incks_val)
      }
    }
    if(found_this){break}
  }
  
  # if(!found_this){
  #   cat("NOT found_this -> ", u ,"\n")
  # }
  
}
print(kc) # find 163 features
kc/406
matched_inchis_pos_dia <- matched_inchis
save(matched_inchis_pos_dia, file = "../optilcms_pos_dia/matched_inchis_pos_dia.rda")



######### DIA testing - biodb
rm(list = ls())

dts <- read.csv("/data/ms2_benchmark/MTBLS2207/41592_2021_1195_MOESM3_ESM_allcmpds.csv")
dt_pos <- dts[dts$Charge == 1, ]
library(OptiLCMS2ID)
#source("~/Github/OptiLCMS2ID/R/MSn_processing.R")
dt_ft <- cbind(dt_pos$m.z - 0.005, dt_pos$m.z + 0.005, dt_pos$RT.Start*60-25, dt_pos$RT.End*60+25)
mSet <- PerformMSnImport(targetFeatures = dt_ft, 
                         filesPath = "/data/ms2_benchmark/MTBLS2207/IROA/DIA/IROA_P1-6_DIA_test_pos1.mzML",
                         acquisitionMode = "DIA",
                         SWATH_file = "/data/ms2_benchmark/MTBLS2207/IROA/optilcms_neg_dia/DIA_SWATH_MS_experiment_file.txt")
mSet <- PerformDIADeconvolution (mSet, 
                                 min_width = 5,
                                 ppm2 = 25,
                                 sn = 12,
                                 span = 0.3,
                                 filtering = 2000, 
                                 ncores = 3L)

mSet <- PerformSpectrumConsenus (mSet, 
                                 ppm2 = 25,
                                 concensus_fraction = 0.5,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)

mSet <- PerformDBSearchingBatch (mSet, 
                                 ppm1 = 10,
                                 ppm2 = 25,
                                 rt_tol = 5,
                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID.sqlite",
                                 use_rt = FALSE,
                                 enableNL = TRUE,
                                 ncores = 8L)

mSet <- PerformResultsExport (mSet,
                              type = 3L,
                              topN = 10L, 
                              ncores = 4L)

save(mSet, file = "/data/ms2_benchmark/MTBLS2207/IROA/optilcms_pos_dia/mSet_res_dia_bioDB.rda")

inchik_map <- read.csv("/data/ms2_benchmark/MTBLS2207/IROA/data_pos_inchikeys.txt", header = F)
k<- mSet@MSnResults[["DecRes"]][[1]][["FeatureIdx"]] + 1
r<- c(1:length(mSet@MSnResults[["DecRes"]][[1]][["FeatureIdx"]]))
annote_res <- mSet@MSnResults[["DBAnnoteRes"]]
kc = 0;
all_score <- vector()
matched_inchis <- vector()
for(u in 1:length(r)){
  res_incks0 <- res_incks <- annote_res[[r[u]]]
  if(length(res_incks[[1]][[1]]) == 0){
    cat("u -> " , u ,"\n");
    next
  }
  res_incks <- sapply(res_incks, function(x){
    max_val <- max(x[[3]]);
    idx <- which(max_val == x[[3]])
    x[[2]][idx]
  })
  res_incks_val <- sapply(res_incks0, function(x){
    max_val <- max(x[[3]]);
    return(max_val)
  })
  
  res_incks <- (res_incks[res_incks != ""])
  res_incks_val <- (res_incks_val[res_incks != ""])
  res_incks <- (sapply(res_incks, function(x){x[1]}))
  #res_incks_val <- (res_incks_val[res_incks != ""])[1]
  if(all(is.na(res_incks))){next;}
  map_incks <- inchik_map[which(inchik_map$V1 == dt_pos[k[u],1]),2]
  found_this <- FALSE;
  for(rs in res_incks){
    for(ms in map_incks){
      if((strsplit(rs, "-")[[1]][1] == strsplit(ms, "-")[[1]][1]) &
         (strsplit(rs, "-")[[1]][2] == strsplit(ms, "-")[[1]][2])) {
        matched_inchis <- c(matched_inchis, ms)
        kc = kc + 1;
        found_this <- TRUE;
        break;
        #all_score <- c(all_score, res_incks_val)
      }
    }
    if(found_this){break}
  }
  
  # if(!found_this){
  #   cat("NOT found_this -> ", u ,"\n")
  # }
  
}
print(kc) # find 163 features
kc/406
matched_inchis_BioDB_pos_dia <- matched_inchis
save(matched_inchis_BioDB_pos_dia, file = "../optilcms_pos_dia/matched_inchis_biodb_pos_dia.rda")

