### MS/MS
setwd("/media/qiang/data2/MTBLS2542/Metabolomics/POS/")
# load("mSet.rda")
# mSet1 <- mSet
# ft_dt_1a <- as.matrix(data.frame(mzmin = mSet1@peakAnnotation[["camera_output"]][["mzmin"]],
#                                  mzmax = mSet1@peakAnnotation[["camera_output"]][["mzmax"]],
#                                  rtmin = mSet1@peakAnnotation[["camera_output"]][["rtmin"]],
#                                  rtmax = mSet1@peakAnnotation[["camera_output"]][["rtmax"]]))
# qs::qsave(ft_dt_1a, file = "ft_dt_1a.qs")
ft_dt_1a <- qs::qread("ft_dt_1a.qs")
library(OptiLCMS2ID)

# mSet <- PerformMSnImport(filesPath = c(list.files("/ultraData/covid/MTBLS2542/mzML/POS/Metabolomics/", 
#                                                   pattern = ".mzML", full.names = T, recursive = T)),
#                          targetFeatures = ft_dt_1a,
#                          acquisitionMode = "DDA")
# 
# qs::qsave(mSet, file = "mSet_ms2.qs")
mSet <- qs::qread("mSet_ms2.qs")
system.time(mSet <- PerformDDADeconvolution(mSet,
                                ppm1 = 5,
                                ppm2 = 21,
                                sn = 12,
                                filtering = 0,
                                window_size = 1.5,
                                intensity_thresh = 1.6e5,
                                database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_Pathway.sqlite",
                                ncores = 1L))

mSet <- PerformSpectrumConsenus (mSet,
                                 ppm2 = 15,
                                 concensus_fraction = 0.2,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)

mSet <- PerformDBSearchingBatch (mSet,
                                 ppm1 = 10,
                                 ppm2 = 25,
                                 rt_tol = 5,
                                 database_path = "/data/COMPOUND_DBs/Curated_DB/MS2ID_Bio.sqlite",
                                 use_rt = FALSE,
                                 enableNL = FALSE,
                                 ncores = 1L)

mSet <- PerformResultsExport (mSet,
                              type = 0L,
                              topN = 10L,
                              ncores = 1L)
qs::qsave(mSet, file = "mSet_ms2.qs")

dtx <- FormatMSnAnnotation(mSet, 5L, F)
idx <- vapply(1:nrow(dtx), function(i){
  x <- dtx[i,c(1:4)]
  y <- which((as.numeric(mSet@MSnData[["peak_mtx"]][,1]) == as.numeric(x[1])) &
               (as.numeric(mSet@MSnData[["peak_mtx"]][,2]) == as.numeric(x[2])) &
               (as.numeric(mSet@MSnData[["peak_mtx"]][,3]) == as.numeric(x[3])) &
               (as.numeric(mSet@MSnData[["peak_mtx"]][,4]) == as.numeric(x[4])))
  return(y)
},integer(1L))
dt_ms1 <- read.csv(file = "peak_table_full.csv",header = T)
dt_ms1 <- dt_ms1[-1,]
labels <- dt_ms1[idx, 1]
dtx_ms2 <- cbind(labels, dtx)
write.csv(dtx_ms2, file = "ms2_compoundtable_pos.csv", row.names = F)

######
setwd("/media/qiang/data2/MTBLS2542/Metabolomics/POS/")
mSet <- qs::qread("mSet_ms2.qs")
df <- OptiLCMS2ID::FormatMSnAnnotation(mSet, 5L, F)
load("mSet.rda")

dft0 <- cbind(mzmin = mSet@peakAnnotation[["camera_output"]][["mzmin"]],
              mzmax = mSet@peakAnnotation[["camera_output"]][["mzmax"]],
              rtmin = mSet@peakAnnotation[["camera_output"]][["rtmin"]],
              rtmax = mSet@peakAnnotation[["camera_output"]][["rtmax"]])

dft <- mSet@peakAnnotation[["camera_output"]][,-c(2,3,5,6, 167:169)]


allmzfiles <- list.files("/ultraData/covid/MTBLS2542/mzML/POS/Metabolomics/", pattern = ".mzML", recursive = T)
allmzfiles_death <- sub(".mzML", "", basename(allmzfiles[grepl("Death", allmzfiles)]))
allmzfiles_mild <- sub(".mzML", "", basename(allmzfiles[grepl("Mild", allmzfiles)]))
idx_death <- vapply(allmzfiles_death, FUN = function(x){
  which(grepl(x, colnames(dft)))
}, integer(1L), USE.NAMES = F)
idx_mild <- vapply(allmzfiles_mild, FUN = function(x){
  which(grepl(x, colnames(dft)))
}, integer(1L), USE.NAMES = F)


pvalues <- vapply(1:nrow(dft), function(x){
  #cat(x, "\n")
  if(all(is.na(dft[x, idx_death])) & (!all(is.na(dft[x, idx_mild])))){
    return(0.0)
  }
  if(!all(is.na(dft[x, idx_death])) & (all(is.na(dft[x, idx_mild])))){
    return(0)
  }
  if(all(is.na(dft[x, idx_death])) & (all(is.na(dft[x, idx_mild])))){
    return(1)
  }
  dft[x,is.na(dft[x,])] <- 0;
  return(t.test(as.numeric(dft[x, idx_death]), as.numeric(dft[x, idx_mild]))$p.value)
}, double(1L))

tscores <- vapply(1:nrow(dft), function(x){
  if(all(is.na(dft[x, idx_death])) & (!all(is.na(dft[x, idx_mild])))){
    return(10.0)
  }
  if(!all(is.na(dft[x, idx_death])) & (all(is.na(dft[x, idx_mild])))){
    return(-10)
  }
  if(all(is.na(dft[x, idx_death])) & (all(is.na(dft[x, idx_mild])))){
    return(0)
  }
  dft[x,is.na(dft[x,])] <- 0;
  return(t.test(as.numeric(dft[x, idx_death]), as.numeric(dft[x, idx_mild]))[["statistic"]][["t"]])
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


