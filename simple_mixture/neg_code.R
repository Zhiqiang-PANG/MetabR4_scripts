# process IROA_dataset 2 - neg
rm(list = ls())
setwd("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/")
all_dt <- read.csv("/data/ms2_benchmark/MTBLS1861/pr0c00930_si_003.csv")

results_sum <- function(mSet, dt){
  peak_mtx <- mSet@MSnData[["peak_mtx"]]
  #peak_mtx <- mSet@MSnResults[["DecRes"]][[1]][["Peak_matrix"]]
  ft_idx <- mSet@MSnResults[["Concensus_spec"]][[1]]
  peak_mtx_ac <- peak_mtx[ft_idx+1,]
  
  kc <- 0;
  idx_vec <- r_vec <- vector()
  weight_h2o <- 18.01056
  weight_h <- 1.0080
  weight_cl <- 34.968853 
  
  for(i in 1:nrow(dt)){
    adH <- as.numeric(dt[i,7]) - weight_h
    adOH <- as.numeric(dt[i,7]) + weight_h2o - weight_h
    adH2OH <- as.numeric(dt[i,7]) - weight_h2o - weight_h
    adCl <- as.numeric(dt[i,7]) + weight_cl
    for(j in 1:nrow(peak_mtx_ac)){
      if((peak_mtx_ac[j,1]-0.005 < adH) & (peak_mtx_ac[j,2]+0.005 > adH)){
        kc <- kc+1;
        idx_vec <-c(idx_vec, i)
        r_vec <- c(r_vec, j)
        #print(i)
        #break;
      }
      if((peak_mtx_ac[j,1]-0.005 < adOH) & (peak_mtx_ac[j,2]+0.005 > adOH)){
        kc <- kc+1;
        idx_vec <-c(idx_vec, i)
        r_vec <- c(r_vec, j)
        #print(i)
        #break;
      }
      if((peak_mtx_ac[j,1]-0.005 < adH2OH) & (peak_mtx_ac[j,2]+0.005 > adH2OH)){
        kc <- kc+1;
        idx_vec <-c(idx_vec, i)
        r_vec <- c(r_vec, j)
        #print(i)
        #break;
      }
      if((peak_mtx_ac[j,1]-0.005 < adCl) & (peak_mtx_ac[j,2]+0.005 > adCl)){
        kc <- kc+1;
        idx_vec <-c(idx_vec, i)
        r_vec <- c(r_vec, j)
        #print(i)
        #break;
      }
    }
  }
  #print(kc)
  print(unique(idx_vec))
  if(length(idx_vec) == 0){return(NULL)}
  inchi_map <- read.csv("/data/ms2_benchmark/MTBLS1861/cid_inchikey_map_1_6", header = F)
  ct <- vector()
  for(k in 1:kc){
    this_cid <- dt[idx_vec[k],10]
    ms <- inchi_map[which(inchi_map$V1 == this_cid),2]
    res_inchis <- unique(mSet@MSnResults[["DBAnnoteRes"]][[r_vec[k]]][[1]][[2]])
    rs <- res_inchis[res_inchis!=""][1]
    
    if((strsplit(rs, "-")[[1]][1] == strsplit(ms, "-")[[1]][1]) &
       (strsplit(rs, "-")[[1]][2] == strsplit(ms, "-")[[1]][2])) {
      ct <- c(ct, k)
    }
  }
  counter <- 0
  counter <- unique(idx_vec[ct])
  
  return(dt[counter,])
}
## 
all_dt_1a <- all_dt[all_dt$PLATE == 1 & all_dt$NROW == "A", ];
all_dt_1b <- all_dt[all_dt$PLATE == 1 & all_dt$NROW == "B", ];
all_dt_1c <- all_dt[all_dt$PLATE == 1 & all_dt$NROW == "C", ];
all_dt_1d <- all_dt[all_dt$PLATE == 1 & all_dt$NROW == "D", ];
all_dt_1e <- all_dt[all_dt$PLATE == 1 & all_dt$NROW == "E", ];
all_dt_1f <- all_dt[all_dt$PLATE == 1 & all_dt$NROW == "F", ];
all_dt_1g <- all_dt[all_dt$PLATE == 1 & all_dt$NROW == "G", ];
all_dt_1h <- all_dt[all_dt$PLATE == 1 & all_dt$NROW == "H", ];

# write.csv(all_dt_1a, file = "all_cmps_dts/all_cmps_1_A.csv", row.names = F)
# write.csv(all_dt_1b, file = "all_cmps_dts/all_cmps_1_B.csv", row.names = F)
# write.csv(all_dt_1c, file = "all_cmps_dts/all_cmps_1_C.csv", row.names = F)
# write.csv(all_dt_1d, file = "all_cmps_dts/all_cmps_1_D.csv", row.names = F)
# write.csv(all_dt_1e, file = "all_cmps_dts/all_cmps_1_E.csv", row.names = F)
# write.csv(all_dt_1f, file = "all_cmps_dts/all_cmps_1_F.csv", row.names = F)
# write.csv(all_dt_1g, file = "all_cmps_dts/all_cmps_1_G.csv", row.names = F)
# write.csv(all_dt_1h, file = "all_cmps_dts/all_cmps_1_H.csv", row.names = F)

#
all_dt_6a <- all_dt[all_dt$PLATE == 6 & all_dt$NROW == "A", ];
all_dt_6b <- all_dt[all_dt$PLATE == 6 & all_dt$NROW == "B", ];
all_dt_6d <- all_dt[all_dt$PLATE == 6 & all_dt$NROW == "C", ];
all_dt_6e <- all_dt[all_dt$PLATE == 6 & all_dt$NROW == "E", ];
all_dt_6f <- all_dt[all_dt$PLATE == 6 & all_dt$NROW == "F", ];
all_dt_6g <- all_dt[all_dt$PLATE == 6 & all_dt$NROW == "G", ];
all_dt_6h <- all_dt[all_dt$PLATE == 6 & all_dt$NROW == "H", ];

# write.csv(all_dt_6a, file = "all_cmps_dts/all_cmps_6_A.csv", row.names = F)
# write.csv(all_dt_6b, file = "all_cmps_dts/all_cmps_6_B.csv", row.names = F)
# write.csv(all_dt_6d, file = "all_cmps_dts/all_cmps_6_D.csv", row.names = F)
# write.csv(all_dt_6e, file = "all_cmps_dts/all_cmps_6_E.csv", row.names = F)
# write.csv(all_dt_6f, file = "all_cmps_dts/all_cmps_6_F.csv", row.names = F)
# write.csv(all_dt_6g, file = "all_cmps_dts/all_cmps_6_G.csv", row.names = F)
# write.csv(all_dt_6h, file = "all_cmps_dts/all_cmps_6_H.csv", row.names = F)

# Plate 1 
#################### 1A -------------

# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE__ce_neg_1_A.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 15, mzdiff = 0.005));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt_1a <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE__ce_neg_1_A.mzML",
#                          targetFeatures = ft_dt_1a,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate1/mSets_1a.rda")
load("optilcms_neg_dda/plate1/mSets_1a.rda")
res_1a <- results_sum(mSet, all_dt_1a)
save(res_1a, file = "optilcms_neg_dda/plate1/res_1a.rda")

#################### 1B -------------
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_neg_1_B.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15,min_peakwidth = 3, max_peakwidth = 15, mzdiff = 0.005));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt_1a <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_neg_1_B.mzML",
#                          targetFeatures = ft_dt_1a,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate1/mSets_1b.rda")
load("optilcms_neg_dda/plate1/mSets_1b.rda")
res_1b <- results_sum(mSet, all_dt_1b)
save(res_1b, file = "optilcms_neg_dda/plate1/res_1b.rda")

#################### 1C -------------
# rm(mSet); rm(mSet1); rm(ft_dt_1a)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_ce_neg_1_C.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 15, mzdiff = 0.005));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_ce_neg_1_C.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate1/mSets_1c.rda")
load("optilcms_neg_dda/plate1/mSets_1c.rda")
res_1c <- results_sum(mSet, all_dt_1c)
save(res_1c, file = "optilcms_neg_dda/plate1/res_1c.rda")

#################### 1D -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_ce_neg_1_D.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 4, max_peakwidth = 15, mzdiff = 0.005));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_ce_neg_1_D.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate1/mSets_1d.rda")
load("optilcms_neg_dda/plate1/mSets_1d.rda")
res_1d <- results_sum(mSet, all_dt_1d)
save(res_1d, file = "optilcms_neg_dda/plate1/res_1d.rda")


#################### 1E -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_ce_neg_1_E.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 15, mzdiff = 0.005));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_ce_neg_1_E.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate1/mSets_1e.rda")
load("optilcms_neg_dda/plate1/mSets_1e.rda")
res_1e <- results_sum(mSet, all_dt_1e)
save(res_1e, file = "optilcms_neg_dda/plate1/res_1e.rda")


#################### 1F -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_ce_neg_1_F.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 15, mzdiff = 0.005));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_ce_neg_1_F.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate1/mSets_1f.rda")
load("optilcms_neg_dda/plate1/mSets_1f.rda")
res_1f <- results_sum(mSet, all_dt_1f)
save(res_1f, file = "optilcms_neg_dda/plate1/res_1f.rda")



#################### 1G -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_ce_neg_1_G.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 15, mzdiff = 0.005));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_ce_neg_1_G.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate1/mSets_1g.rda")
load("optilcms_neg_dda/plate1/mSets_1g.rda")
res_1g <- results_sum(mSet, all_dt_1g)
save(res_1g, file = "optilcms_neg_dda/plate1/res_1g.rda")

#################### 1H -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_neg_1_H.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 15, mzdiff = 0.005));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_neg_1_H.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate1/mSets_1h.rda")
load("optilcms_neg_dda/plate1/mSets_1h.rda")
res_1h <- results_sum(mSet, all_dt_1h)
save(res_1h, file = "optilcms_neg_dda/plate1/res_1h.rda")

# Plate 1 

##############
## Plate 6
#################### 6A -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_neg_6_A.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 12));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_neg_6_A.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate6/mSets_6a.rda")
load("optilcms_neg_dda/plate6/mSets_6a.rda")
res_6a <- results_sum(mSet, all_dt_6a)
save(res_6a, file = "optilcms_neg_dda/plate6/res_6a.rda")


#################### 6B -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_neg_6_B.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 25, min_peakwidth = 3, max_peakwidth = 12));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_neg_6_B.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 4L) ## TODO: cause error when using 8 cores
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate6/mSets_6b.rda")
load("optilcms_neg_dda/plate6/mSets_6b.rda")
res_6b <- results_sum(mSet, all_dt_6b)
save(res_6b, file = "optilcms_neg_dda/plate6/res_6b.rda")

#################### 6D -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_neg_6_D.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 12));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_neg_6_D.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate6/mSets_6d.rda")
load("optilcms_neg_dda/plate6/mSets_6d.rda")
res_6d <- results_sum(mSet, all_dt_6d)
save(res_6d, file = "optilcms_neg_dda/plate6/res_6d.rda")

#################### 6E -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_neg_6_E.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 12));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_neg_6_E.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate6/mSets_6e.rda")
load("optilcms_neg_dda/plate6/mSets_6e.rda")
res_6e <- results_sum(mSet, all_dt_6e)
save(res_6e, file = "optilcms_neg_dda/plate6/res_6e.rda")


#################### 6F -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_neg_6_F.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 12));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_neg_6_F.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate6/mSets_6f.rda")
load("optilcms_neg_dda/plate6/mSets_6f.rda")
res_6f <- results_sum(mSet, all_dt_6f)
save(res_6f, file = "optilcms_neg_dda/plate6/res_6f.rda")



#################### 6G -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_neg_6_G.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 12));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_neg_6_G.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 6L) #TODO: why give 8 cores but use only 7, and cause error
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate6/mSets_6g.rda")
load("optilcms_neg_dda/plate6/mSets_6g.rda")
res_6g <- results_sum(mSet, all_dt_6g)
save(res_6g, file = "optilcms_neg_dda/plate6/res_6g.rda")



#################### 6H -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "neg/IROA_PLATE_neg_6_H.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 3, max_peakwidth = 12));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/neg/IROA_PLATE_neg_6_H.mzML",
#                          targetFeatures = ft_dt,
#                          acquisitionMode = "DDA")
# 
# mSet <- PerformDDADeconvolution(mSet,
#                                 ppm1 = 10,
#                                 ppm2 = 25,
#                                 sn = 12,
#                                 filtering = 2000,
#                                 window_size = 1,
#                                 intensity_thresh = 1e4,
#                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                 ncores = 8L)
# 
# mSet <- PerformSpectrumConsenus (mSet,
#                                  ppm2 = 25,
#                                  concensus_fraction = 0.5,
#                                  database_path = "",
#                                  use_rt = FALSE,
#                                  user_dbCorrection = FALSE)
# 
# mSet <- PerformDBSearchingBatch (mSet,
#                                  ppm1 = 10,
#                                  ppm2 = 25,
#                                  rt_tol = 5,
#                                  database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                                  use_rt = FALSE,
#                                  enableNL = FALSE,
#                                  ncores = 8L)
# 
# mSet <- PerformResultsExport (mSet,
#                               type = 3L,
#                               topN = 10L,
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_neg_dda/plate6/mSets_6h.rda")
load("optilcms_neg_dda/plate6/mSets_6h.rda")
res_6h <- results_sum(mSet, all_dt_6h)
save(res_6h, file = "optilcms_neg_dda/plate6/res_6h.rda")





#########
#summarize
### 
# 1a-1h & 6a-6h

load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate1/res_1a.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate1/res_1b.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate1/res_1c.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate1/res_1d.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate1/res_1e.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate1/res_1f.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate1/res_1g.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate1/res_1h.rda")

load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate6/res_6a.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate6/res_6b.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate6/res_6d.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate6/res_6e.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate6/res_6f.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate6/res_6g.rda")
load("/data/ms2_benchmark/MTBLS1861/IROA_PLATE/optilcms_neg_dda/plate6/res_6h.rda")

# neg: 
# 1A - c(NAD, L-GLUTAMINE, HYPOTAURINE, INOSINE 5'-PHOSPHATE, CITRATE, L-THREONINE, PURINE, N-ACETYLNEURAMINATE, L-KYNURENINE, D-ASPARTATE, URATE)
# 1B - c(NICOTINATE, Cytidine, L-SERINE, TAURINE, D-GLUCONO-1,5-LACTONE, INOSINE, CYTOSINE, Isoleucine)
# 1C - c(GLUCONIC ACID)
# 1D - c(ADENINE, XANTHINE, 5'-METHYLTHIOADENOSINE, THYMIDINE, OROTATE, L-CYSTINE, ETHANOLAMINE PHOSPHATE, GLYCERATE, L-METHIONINE)
# 1E - c(SUCCINATE SEMIALDEHYDE, L-ALANINE, L-TRYPTOPHAN, URIDINE-5-MONOPHOSPHATE, L-PROLINE, (S)-LACTATE, URIDINE, FRUCTOSE 1,6-BIPHOSPHATE, CARNOSINE, SHIKIMATE)
# 1F - c(Guanosine, SUCCINATE, L-PHENYLALANINE, URACIL, (S)-MALATE, L-ASPARTATE, 2'-DEOXYCYTIDINE 5'-MONOPHOSPHATE, HYPOXANTHINE, CREATINE, 3,4-DIHYDROXY-L-PHENYLALANINE, L-TYROSINE )
# 1G - c((R,R)-TARTARIC ACID,  L-LYSINE, L-VALINE, L-TYROSINE, L-ASPARAGINE, HOMOSERINE, PYRIDOXINE, DAMP, Folic acid)
# 1H - c(Isocitric acid, N-ACETYL-D-TRYPTOPHAN, D-GLUCOSE-6-PHOSPHATE, GUANIDINOACETATE, CREATININE)

# 6A - c(2,5-Dihydroxybenzoic acid, DETHIOBIOTIN, ALPHA-KETOGLUTARIC ACID, N-ACETYLSEROTONIN, ITACONATE, AZELAIC ACID, 2-METHYLGLUTARIC ACID)
# 6B - c(ADIPIC ACID)
# 6D - c(2',4'-DIHYDROXYACETOPHENONE, SALICYLAMIDE, BENZYL ALCOHOL)
# 6E - c(PHENYLACETIC ACID, RESORCINOL MONOACETATE, 3-(4-HYDROXYPHENYL)LACTATE, Biotin, 5-HYDROXYINDOLEACETATE)
# 6F - c(ETHYL 3-INDOLEACETATE, L-SORBOSE, XYLITOL, Myo-Inositol, MANNOSE, RIBITOL)
# 6G - c(D-GULONIC ACID GAMA-LACTONE, SUCROSE, ALPHA-D-GLUCOSE, ALLOSE, MANNITOL, MALTOSE, D-TAGATOSE)
# 6H - c(D-PSICOSE, L-ARABITOL, D-RIBOSE, PALATINOSE)
p1 <- c(nrow(res_1a)/11, 
        nrow(res_1b)/8, 
        nrow(res_1c)/1, 
        nrow(res_1d)/9, 
        nrow(res_1e)/10, 
        nrow(res_1f)/11, 
        nrow(res_1g)/9, 
        nrow(res_1h)/5)

p6 <- c(nrow(res_6a)/7, 
        0/1, 
        nrow(res_6d)/3, 
        nrow(res_6e)/5, 
        nrow(res_6f)/6, 
        nrow(res_6g)/7, 
        nrow(res_6h)/4)
p_all <- c(p1, p6)
save(p1, p6, p_all, file = "optilcms_neg_dda/all_ratio_res.rda")
