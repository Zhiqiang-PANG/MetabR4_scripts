setwd("/media/qiang/data2/MSV000089568/res/Negative")
library(OptiLCMS)

mSet1 <- ImportRawMSData(path = ".")
qs::qsave(mSet1, file = "mSet1.qs")
#load("mSet1.rda")
mSet1 <- PerformPeakProfiling(mSet1, Params = SetPeakParam(ppm = 25, 
                                                           bw = 3, 
                                                           mzdiff = 0.001, 
                                                           max_peakwidth = 35, 
                                                           min_peakwidth = 5, 
                                                           noise = 200, minFraction = 0.2), 
                              ncore = 4, 
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
qs::qsave(mSet1, file = "mSet1.qs")

library(OptiLCMS2ID)
mSet1 <- qs::qread("mSet1.qs")
mSet <- PerformMSnImport(mSet = mSet1,
                         filesPath = c(list.files(".", pattern = ".mzML", full.names = T, recursive = T)[c(1:6)]),
                         acquisitionMode = "DIA", 
                         SWATH_file = "DIA_SWATH_MS_experiment_file_neg.txt")
qs::qsave(mSet, file = "mSet2.qs")

mSet <- PerformDIADeconvolution(mSet,
                                min_width = 5,
                                span = 0.3,
                                ppm2 = 30,
                                sn = 12,
                                filtering = 0,
                                ncores = 1L)

mSet <- PerformSpectrumConsenus (mSet,
                                 ppm2 = 30,
                                 concensus_fraction = 0.25,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)
qs::qsave(mSet, file = "mSet2.qs")
mSet <- PerformDBSearchingBatch (mSet,
                                 ppm1 = 15,
                                 ppm2 = 30,
                                 rt_tol = 5,
                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_Bio.sqlite",
                                 use_rt = FALSE,
                                 enableNL = FALSE,
                                 ncores = 1L)

mSet <- PerformResultsExport (mSet,
                              type = 0L,
                              topN = 10L,
                              ncores = 1L)
qs::qsave(mSet, file = "dia_neg_optilcms_bio.qs")


####
mSet <- qs::qread("dia_neg_optilcms_bio.qs")
df <- OptiLCMS2ID::FormatMSnAnnotation(mSet, 5L, F)
mSet <- qs::qread("mSet1.qs")

dft0 <- cbind(mzmin = mSet@peakAnnotation[["camera_output"]][["mzmin"]],
             mzmax = mSet@peakAnnotation[["camera_output"]][["mzmax"]],
             rtmin = mSet@peakAnnotation[["camera_output"]][["rtmin"]],
             rtmax = mSet@peakAnnotation[["camera_output"]][["rtmax"]])

dft <- mSet@peakAnnotation[["camera_output"]][,-c(2,3,5,6, 35:37)]


pvalues <- vapply(1:nrow(dft), function(x){
  #cat(x, "\n")
  if(all(is.na(dft[x, c(3:12)])) & (!all(is.na(dft[x, c(13:30)])))){
    return(0.0)
  }
  if(!all(is.na(dft[x, c(3:12)])) & (all(is.na(dft[x, c(13:30)])))){
    return(0)
  }
  dft[x,is.na(dft[x,])] <- 0;
  return(t.test(as.numeric(dft[x, c(3:12)]), as.numeric(dft[x, c(13:30)]))$p.value)
}, double(1L))

tscores <- vapply(1:nrow(dft), function(x){
  if(all(is.na(dft[x, c(3:12)])) & (!all(is.na(dft[x, c(13:30)])))){
    return(10.0)
  }
  if(!all(is.na(dft[x, c(3:12)])) & (all(is.na(dft[x, c(13:30)])))){
    return(10)
  }
  
  dft[x,is.na(dft[x,])] <- 0;
  return(t.test(as.numeric(dft[x, c(3:12)]), as.numeric(dft[x, c(13:30)]))[["statistic"]][["t"]])
}, double(1L))

mzs <- vapply(1:nrow(dft), function(x){
  return(as.numeric(dft$mz[x]))
}, double(1L))

rts <- vapply(1:nrow(dft), function(x){
  return(as.numeric(dft$rt[x]))
}, double(1L))

df_res <- cbind(mz = round(mzs,5), rt = round(rts,digits = 2), t.score = round(tscores, 4), p.value = round(pvalues,4))

df_ms2 <- as.data.frame(matrix(nrow = length(rts), ncol=5))

dt_ms2 <- df
for(i in 1:nrow(df_ms2)){
  cat("Processing ==> ", i, "\n")
  for(j in 1:nrow(dt_ms2)){
    if((mzs[i] >= dt_ms2$mzmin[j]) & 
       (mzs[i] <= dt_ms2$mzmax[j]) & 
       (rts[i] >= dt_ms2$rtmin[j]) &
       (rts[i] <= dt_ms2$rtmax[j])){
      df_ms2[i,] <- c(dt_ms2$InchiKey_1[j], dt_ms2$InchiKey_2[j], dt_ms2$InchiKey_3[j], dt_ms2$InchiKey_4[j], dt_ms2$InchiKey_5[j])
    }
  }
}


write.table(df_ms2, file = "optilcms_mummichog_data/compound_ms2.txt", row.names = F, sep = ",")
write.table(df_res, file = "optilcms_mummichog_data/peaks_ms1.txt", row.names = F, sep = ",")










