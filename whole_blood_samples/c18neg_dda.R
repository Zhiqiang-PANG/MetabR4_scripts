setwd("/ultraData/new_wb/C18neg/mzML/")
# C18 neg
library(OptiLCMS)

mSet1 <- ImportRawMSData(path = "MS1/")

mSet1 <- PerformPeakProfiling(mSet1, Params = SetPeakParam(ppm = 15, 
                                                           bw = 3, 
                                                           mzdiff = 0.001, 
                                                           max_peakwidth = 25, 
                                                           min_peakwidth = 5, 
                                                           noise = 500, minFraction = 0.5), 
                              ncore = 8, 
                              plotSettings = SetPlotParam(Plot = F))

annParams <- SetAnnotationParam(polarity = 'negative',
                                mz_abs_add = 0.035);

mSet1 <- PerformPeakAnnotation(mSet = mSet1,
                               annotaParam = annParams,
                               ncore =1)
mSet1 <- FormatPeakList(mSet = mSet1,
                        annParams,
                        filtIso =FALSE,
                        filtAdducts = FALSE,
                        missPercent = 1)
qs::qsave(mSet1, file = "../OptiLCMS_ms2/ms1/mSet1.qs")
ft_dt_1a <- as.matrix(data.frame(mzmin = mSet1@peakAnnotation[["camera_output"]][["mzmin"]],
                                 mzmax = mSet1@peakAnnotation[["camera_output"]][["mzmax"]],
                                 rtmin = mSet1@peakAnnotation[["camera_output"]][["rtmin"]],
                                 rtmax = mSet1@peakAnnotation[["camera_output"]][["rtmax"]]))



library(OptiLCMS2ID)

mSet <- PerformMSnImport(filesPath = c(list.files("MS2/DDA/", pattern = ".mzML", full.names = T),
                                       list.files("MS2/tDDA/", pattern = ".mzML", full.names = T)),
                         targetFeatures = ft_dt_1a,
                         acquisitionMode = "DDA")

mSet <- PerformDDADeconvolution(mSet,
                                ppm1 = 10,
                                ppm2 = 25,
                                sn = 12,
                                filtering = 0,
                                window_size = 1,
                                intensity_thresh = 1.6e5,
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
                                 enableNL = FALSE,
                                 ncores = 8L)

mSet <- PerformResultsExport (mSet,
                              type = 3L,
                              topN = 10L,
                              ncores = 8L)
qs::qsave(mSet, file = "MS2/dda_res_optilcms_complete.qs")
