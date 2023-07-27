# process IROA_dataset 2 - pos
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
  weight_na <- 22.9897693
  weight_nh4 <- 18.03846
  weight_h <- 1.0080
  
  for(i in 1:nrow(dt)){
    adH <- as.numeric(dt[i,7]) + weight_h
    adNa <- as.numeric(dt[i,7]) + weight_na
    adNH <- as.numeric(dt[i,7]) + weight_nh4
    for(j in 1:nrow(peak_mtx_ac)){
      if((peak_mtx_ac[j,1]-0.005 < adH) & (peak_mtx_ac[j,2]+0.005 > adH)){
        kc <- kc+1;
        idx_vec <-c(idx_vec, i)
        r_vec <- c(r_vec, j)
        #print(i)
        #break;
      }
      if((peak_mtx_ac[j,1]-0.005 < adNa) & (peak_mtx_ac[j,2]+0.005 > adNa)){
        kc <- kc+1;
        idx_vec <-c(idx_vec, i)
        r_vec <- c(r_vec, j)
        #print(i)
        #break;
      }
      if((peak_mtx_ac[j,1]-0.005 < adNH) & (peak_mtx_ac[j,2]+0.005 > adNH)){
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
all_dt_6d <- all_dt[all_dt$PLATE == 6 & all_dt$NROW == "D", ];
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
# 
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_1_A.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt_1a <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_1_A.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate1/mSets_1a.rda")
load("optilcms_pos_dda/plate1/mSets_1a.rda")
res_1a <- results_sum(mSet, all_dt_1a)
save(res_1a, file = "optilcms_pos_dda/plate1/res_1a.rda")

#################### 1B -------------
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_1_B.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt_1a <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_1_B.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate1/mSets_1b.rda")
load("optilcms_pos_dda/plate1/mSets_1b.rda")
res_1b <- results_sum(mSet, all_dt_1b)
save(res_1b, file = "optilcms_pos_dda/plate1/res_1b.rda")

#################### 1C -------------
# rm(mSet); rm(mSet1); rm(ft_dt_1a)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_1_C.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_1_C.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate1/mSets_1c.rda")
load("optilcms_pos_dda/plate1/mSets_1c.rda")
res_1c <- results_sum(mSet, all_dt_1c)
save(res_1c, file = "optilcms_pos_dda/plate1/res_1c.rda")

#################### 1D -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_1_D.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_1_D.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate1/mSets_1d.rda")
load("optilcms_pos_dda/plate1/mSets_1d.rda")
res_1d <- results_sum(mSet, all_dt_1d)
save(res_1d, file = "optilcms_pos_dda/plate1/res_1d.rda")


#################### 1E -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_1_E.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_1_E.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate1/mSets_1e.rda")
load("optilcms_pos_dda/plate1/mSets_1e.rda")
res_1e <- results_sum(mSet, all_dt_1e)
save(res_1e, file = "optilcms_pos_dda/plate1/res_1e.rda")


#################### 1F -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_1_F.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_1_F.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate1/mSets_1f.rda")
load("optilcms_pos_dda/plate1/mSets_1f.rda")
res_1f <- results_sum(mSet, all_dt_1f)
save(res_1f, file = "optilcms_pos_dda/plate1/res_1f.rda")



#################### 1G -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_1_G.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_1_G.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate1/mSets_1g.rda")
load("optilcms_pos_dda/plate1/mSets_1g.rda")
res_1g <- results_sum(mSet, all_dt_1g)
save(res_1g, file = "optilcms_pos_dda/plate1/res_1g.rda")

#################### 1H -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_1_H.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_1_H.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate1/mSets_1h.rda")
load("optilcms_pos_dda/plate1/mSets_1h.rda")
res_1h <- results_sum(mSet, all_dt_1h)
save(res_1h, file = "optilcms_pos_dda/plate1/res_1h.rda")

# Plate 1 

##############
## Plate 6
#################### 6A -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_6_A.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_6_A.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate6/mSets_6a.rda")
load("optilcms_pos_dda/plate6/mSets_6a.rda")
res_6a <- results_sum(mSet, all_dt_6a)
save(res_6a, file = "optilcms_pos_dda/plate6/res_6a.rda")


#################### 6B -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_6_B.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_6_B.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate6/mSets_6b.rda")
load("optilcms_pos_dda/plate6/mSets_6b.rda")
res_6b <- results_sum(mSet, all_dt_6b)
save(res_6b, file = "optilcms_pos_dda/plate6/res_6b.rda")

#################### 6D -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_6_D.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_6_D.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate6/mSets_6d.rda")
load("optilcms_pos_dda/plate6/mSets_6d.rda")
res_6d <- results_sum(mSet, all_dt_6d)
save(res_6d, file = "optilcms_pos_dda/plate6/res_6d.rda")

#################### 6E -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_6_E.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_6_E.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate6/mSets_6e.rda")
load("optilcms_pos_dda/plate6/mSets_6e.rda")
res_6e <- results_sum(mSet, all_dt_6e)
save(res_6e, file = "optilcms_pos_dda/plate6/res_6e.rda")


#################### 6F -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_6_F.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_6_F.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate6/mSets_6f.rda")
load("optilcms_pos_dda/plate6/mSets_6f.rda")
res_6f <- results_sum(mSet, all_dt_6f)
save(res_6f, file = "optilcms_pos_dda/plate6/res_6f.rda")



#################### 6G -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_6_G.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_6_G.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate6/mSets_6g.rda")
load("optilcms_pos_dda/plate6/mSets_6g.rda")
res_6g <- results_sum(mSet, all_dt_6g)
save(res_6g, file = "optilcms_pos_dda/plate6/res_6g.rda")



#################### 6H -------------
# rm(mSet); rm(mSet1); rm(ft_dt)
# library(OptiLCMS)
# mSet1 <- ImportRawMSData(path = "pos/IROA_PLATE_pos_6_H.mzML")
# mSet1@params <- OptiLCMS:::updateRawSpectraParam (SetPeakParam(ppm = 15, min_peakwidth = 5, max_peakwidth = 15));
# mSet1 <- PerformPeakPicking(mSet = mSet1)
# 
# ft_dt <- mSet1@peakpicking[["chromPeaks"]][,c(2,3,5,6)]
# 
# library(OptiLCMS2ID)
# 
# mSet <- PerformMSnImport(filesPath = "/data/ms2_benchmark/MTBLS1861/IROA_PLATE/pos/IROA_PLATE_pos_6_H.mzML",
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
#                                 ncores = 1L)
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
#                               database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_complete.sqlite",
#                               ncores = 8L)
# 
# save(mSet, mSet1, file = "optilcms_pos_dda/plate6/mSets_6h.rda")
load("optilcms_pos_dda/plate6/mSets_6h.rda")
res_6h <- results_sum(mSet, all_dt_6h)
save(res_6h, file = "optilcms_pos_dda/plate6/res_6h.rda")





#########
#summarize
### 
# 1a-1h
p1 <- c(nrow(res_1a)/6, 1, 1, nrow(res_1d)/3, nrow(res_1e)/5, nrow(res_1f)/8, nrow(res_1g)/7, nrow(res_1h)/7)
p6 <- c(nrow(res_6a)/3, nrow(res_6b)/2, nrow(res_6d)/1, nrow(res_6e)/2, nrow(res_6f)/1, nrow(res_6g)/2, nrow(res_6h)/1)
p_all <- c(p1, p6)
save(p1, p6, p_all, file = "optilcms_pos_dda/all_ratio_res.rda")
