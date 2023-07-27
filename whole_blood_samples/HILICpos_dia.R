setwd("/ultraData/new_wb/HILICpos//mzML/")
# HILIC POS
library(OptiLCMS)
#load("MS1/mSet.rda")
mSet <- qs::qread("../OptiLCMS_ms2/ms1/mSet1.qs")

library(OptiLCMS2ID)

mSet <- PerformMSnImport(mSet = mSet,
                         filesPath = c(list.files("MS2/DIA/", pattern = ".mzML", full.names = T)),
                         acquisitionMode = "DIA", 
                         SWATH_file = "MS2/DIA/DIA_SWATH_MS_experiment_file_optilcms.txt")


mSet <- PerformDIADeconvolution(mSet,
                                min_width = 5,
                                span = 0.3,
                                ppm2 = 25,
                                sn = 12,
                                filtering = 0,
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
                                 enableNL = FALSE,
                                 ncores = 8L)

mSet <- PerformResultsExport (mSet,
                              type = 3L,
                              topN = 10L,
                              ncores = 8L)
qs::qsave(mSet, file = "MS2/dia_res_optilcms_complete.qs")
