### WHOLE BLOOD - PLASMA -SERUM [R SCRIPTS]
setwd("/ultraData/new_wb/")
## DDA
### - C18 POS
{
  # 1. OptiLCMS
  rm(list = ls())
  mSet2 <- qs::qread("C18pos/mzML/MS2/dda_res_optilcms_complete.qs")
  library(OptiLCMS2ID);
  
  Inchikeys <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  mSet2 <- PerformResultsExport (mSet2,
                                 type = 1L,
                                 topN = 10L,
                                 ncores = 6L);
  
  cmpdsNMs <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  
  identifiedCMPD_dt <- data.frame(Inchikey=Inchikeys, CompoundName = cmpdsNMs)
  #identifiedCMPD_dt <- identifiedCMPD_dt[!is.na(identifiedCMPD_dt$Inchikey), ]
  
  mSet <- qs::qread("C18pos/OptiLCMS_ms2/ms1/mSet1.qs")
  dt <- mSet@dataSet
  idx_plasma <- which(dt[1,] == "Plasma")
  idx_serum <- which(dt[1,] == "Serum")
  idx_wb <- which(dt[1,] == "Whole_blood")
  idx_qc <- which(dt[1,] == "QC")
  feadt0 <- dt[-1,]
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& ((sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3) | is.na(sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3))
  })
  mzs <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][1])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  rts <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][2])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18pos/unique_ms1_features_idx_optilcms.qs")
  
  ft_dt_1a <- as.matrix(data.frame(mzmin = mSet@peakAnnotation[["camera_output"]][["mzmin"]],
                                   mzmax = mSet@peakAnnotation[["camera_output"]][["mzmax"]],
                                   rtmin = mSet@peakAnnotation[["camera_output"]][["rtmin"]],
                                   rtmax = mSet@peakAnnotation[["camera_output"]][["rtmax"]]))
  
  ft_dt_1b <- ft_dt_1a[mSet2@MSnResults[["Concensus_spec"]][[1]]+1,]
  Confirmed_idx <- vector(length = nrow(ft_dt_1b))
  for(i in 1:nrow(ft_dt_1b)){
    for(j in 1:length(mzs)){
      if((mzs[j] > ft_dt_1b[i,1]-0.002) &
         (mzs[j] < ft_dt_1b[i,2]+0.002) &
         (rts[j] > ft_dt_1b[i,3]-5) &
         (rts[j] < ft_dt_1b[i,4]+5)){
        Confirmed_idx[i] <- T;
        break;
      } 
    }
  }
  
  identifiedCMPD_dt0 <- identifiedCMPD_dt[Confirmed_idx, ]
  identifiedCMPD_dtx <- identifiedCMPD_dt0[!is.na(identifiedCMPD_dt0$CompoundName), ]
  save(identifiedCMPD_dtx, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18pos/Optilcms_resDT.rda")
  # 2. MSDIAL/MSFINDER
  rm(list = ls())
  resdt <- read.csv("C18pos/MSDIAL_ms2/dda_ms2/Structure result-2085.txt", sep = "\t")
  resdt <- resdt[resdt$Rank == 1, ]
  
  feadt <- read.csv("C18pos/MSDIAL_ms2/dda_ms1_Ms2/Height_0_2023111610.txt", sep = "\t")
  feadt <- feadt[-c(1:3),-c(1,4,5:32, ncol(feadt), ncol(feadt)-1)]
  idx <- grepl(pattern = "WB[0-9]+|PLASMA[0-9]+|SERUM[0-9]+|QC[0-9]+", x = feadt[1,])
  feadt1 <- feadt[,idx]
  feadt0 <- feadt1[-1,]
  colnms <- feadt1[1,]
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))
  
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18pos/unique_ms1_features_idx_msdial.qs")
  feadt_res <- feadt0[row_idx, ]
  feat_meta <- feadt[which(row_idx)+1,c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(resdt))
  for(k in 1:nrow(resdt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(resdt[k,5] - as.numeric(feat_meta$X.2[s])) < 0.005){
        this_rt <- as.numeric(strsplit(resdt$File.name[k], "_")[[1]][2])
        rt_val <-  as.numeric(feat_meta$X.1[s])
        if(abs(rt_val - this_rt) < 1/3){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  resdt_confirmed <- resdt[Confirmed_idx,]
  save(resdt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18pos/MSDIAL_resDT.rda")
  # 3. MZmine/SIRIUS
  rm(list = ls())
  dt <- read.csv("C18pos/SIRIUS_ms2/mzmine_dda_ms2_res/res/compound_identifications.tsv", sep = "\t")
  labels <- paste0(dt$ionMass, "__", dt$retentionTimeInSeconds)
  labelsu <- unique(labels)
  idxx <- vapply(labelsu, function(x){
    which(x == labels)[1]
  }, FUN.VALUE = integer(1L), USE.NAMES = F)
  dt <- dt[idxx, ]
  #dt <- dt[c(dt$formulaRank==1 & dt$X.predictedFPs==1),]
  dt <- dt[c(dt$formulaRank==1),]
  
  dt_ms1 <- read.csv("C18pos/SIRIUS_ms2/mzmine_dda_ms1_res/c18pos_ms1.csv")
  colnm_idx <- which(grepl("*.mzML.area", colnames(dt_ms1)))
  dt_ms1x <- dt_ms1[,c(7,9,colnm_idx)]
  
  colnms <- colnames(dt_ms1x)
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))
  
  row_idx <- apply(dt_ms1x, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18pos/unique_ms1_features_idx_mzmine.qs")
  feadt_res <- dt_ms1x[row_idx, c(1,2, idx_plasma, idx_serum, idx_wb, idx_qc)]
  feat_meta <- dt_ms1x[which(row_idx),c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(dt))
  for(k in 1:nrow(dt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(dt$ionMass[k] - as.numeric(feat_meta$mz[s])) < 0.005){
        this_rt <- as.numeric(dt$retentionTimeInSeconds[k])
        rt_val <-  as.numeric(feat_meta$rt[s])*60
        if(abs(rt_val - this_rt) < 20){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  dt_confirmed <- dt[Confirmed_idx,]
  save(dt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18pos/Mzmine_sirius_resDT.rda")
  
}

### - C18 NEG
{
  # 1. OptiLCMS
  rm(list = ls())
  mSet2 <- qs::qread("C18neg/mzML/MS2/dda_res_optilcms_complete.qs")
  library(OptiLCMS2ID);
  
  Inchikeys <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  mSet2 <- PerformResultsExport (mSet2,
                                 type = 1L,
                                 topN = 10L,
                                 ncores = 6L);
  
  cmpdsNMs <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  
  identifiedCMPD_dt <- data.frame(Inchikey=Inchikeys, CompoundName = cmpdsNMs)
  #identifiedCMPD_dt <- identifiedCMPD_dt[!is.na(identifiedCMPD_dt$Inchikey), ]
  
  mSet <- qs::qread("C18neg/OptiLCMS_ms2/ms1/mSet1.qs")
  dt <- mSet@dataSet
  idx_plasma <- which(dt[1,] == "Plasma")
  idx_serum <- which(dt[1,] == "Serum")
  idx_wb <- which(dt[1,] == "Whole_blood")
  idx_qc <- which(dt[1,] == "QC")
  feadt0 <- dt[-1,]
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& ((sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3) | is.na(sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3))
  })
  mzs <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][1])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  rts <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][2])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18neg/unique_ms1_features_idx_optilcms.qs")
  
  ft_dt_1a <- as.matrix(data.frame(mzmin = mSet@peakAnnotation[["camera_output"]][["mzmin"]],
                                   mzmax = mSet@peakAnnotation[["camera_output"]][["mzmax"]],
                                   rtmin = mSet@peakAnnotation[["camera_output"]][["rtmin"]],
                                   rtmax = mSet@peakAnnotation[["camera_output"]][["rtmax"]]))
  
  ft_dt_1b <- ft_dt_1a[mSet2@MSnResults[["Concensus_spec"]][[1]]+1,]
  Confirmed_idx <- vector(length = nrow(ft_dt_1b))
  for(i in 1:nrow(ft_dt_1b)){
    for(j in 1:length(mzs)){
      if((mzs[j] > ft_dt_1b[i,1]-0.002) &
         (mzs[j] < ft_dt_1b[i,2]+0.002) &
         (rts[j] > ft_dt_1b[i,3]-5) &
         (rts[j] < ft_dt_1b[i,4]+5)){
        Confirmed_idx[i] <- T;
        break;
      } 
    }
  }
  
  identifiedCMPD_dt0 <- identifiedCMPD_dt[Confirmed_idx, ]
  identifiedCMPD_dtx <- identifiedCMPD_dt0[!is.na(identifiedCMPD_dt0$CompoundName), ]
  save(identifiedCMPD_dtx, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18neg/Optilcms_resDT.rda")
  # 2. MSDIAL/MSFINDER
  rm(list = ls())
  resdt <- read.csv("C18neg/MSDIAL_ms2/dda_ms2_res/Structure result-2065.txt", sep = "\t")
  resdt <- resdt[resdt$Rank == 1, ]
  
  feadt <- read.csv("C18neg/MSDIAL_ms2/dda_ms1_ms2/Height_0_20231111553.txt", sep = "\t")
  feadt <- feadt[-c(1:3),-c(1,4,5:32, ncol(feadt), ncol(feadt)-1)]
  idx <- grepl(pattern = "WB[0-9]+|PLASMA[0-9]+|SERUM[0-9]+|QC[0-9]+", x = feadt[1,])
  feadt1 <- feadt[,idx]
  feadt0 <- feadt1[-1,]
  colnms <- feadt1[1,]
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))
  
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18neg/unique_ms1_features_idx_msdial.qs")
  feadt_res <- feadt0[row_idx, ]
  feat_meta <- feadt[which(row_idx)+1,c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(resdt))
  for(k in 1:nrow(resdt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(resdt[k,5] - as.numeric(feat_meta$X.2[s])) < 0.005){
        this_rt <- as.numeric(strsplit(resdt$File.name[k], "_")[[1]][2])
        rt_val <-  as.numeric(feat_meta$X.1[s])
        if(abs(rt_val - this_rt) < 1/3){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  resdt_confirmed <- resdt[Confirmed_idx,]
  save(resdt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18neg/MSDIAL_resDT.rda")
  # 3. MZmine/SIRIUS
  rm(list = ls())
  dt <- read.csv("C18neg/SIRIUS_ms2/mzmine_dda_ms2_res/res/compound_identifications.tsv", sep = "\t")
  labels <- paste0(dt$ionMass, "__", dt$retentionTimeInSeconds)
  labelsu <- unique(labels)
  idxx <- vapply(labelsu, function(x){
    which(x == labels)[1]
  }, FUN.VALUE = integer(1L), USE.NAMES = F)
  dt <- dt[idxx, ]
  #dt <- dt[c(dt$formulaRank==1 & dt$X.predictedFPs==1),]
  dt <- dt[c(dt$formulaRank==1),]
  
  dt_ms1 <- read.csv("C18neg/SIRIUS_ms2/mzmine_dda_ms1_res/c18neg_ms1.csv")
  colnm_idx <- which(grepl("*.mzML.area", colnames(dt_ms1)))
  colnm_idx2 <- which((colnames(dt_ms1) == "rt") | (colnames(dt_ms1) == "mz"))
  dt_ms1x <- dt_ms1[,c(colnm_idx2,colnm_idx)]
  
  colnms <- colnames(dt_ms1x)
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))
  
  row_idx <- apply(dt_ms1x, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18neg/unique_ms1_features_idx_mzmine.qs")
  feadt_res <- dt_ms1x[row_idx, c(1,2, idx_plasma, idx_serum, idx_wb, idx_qc)]
  feat_meta <- dt_ms1x[which(row_idx),c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(dt))
  for(k in 1:nrow(dt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(dt$ionMass[k] - as.numeric(feat_meta$mz[s])) < 0.005){
        this_rt <- as.numeric(dt$retentionTimeInSeconds[k])
        rt_val <-  as.numeric(feat_meta$rt[s])*60
        if(abs(rt_val - this_rt) < 20){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  dt_confirmed <- dt[Confirmed_idx,]
  save(dt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DDA/C18neg/Mzmine_sirius_resDT.rda")
}

### - HILIC POS
{
  # 1. OptiLCMS
  rm(list = ls())
  mSet2 <- qs::qread("HILICpos/mzML/MS2/dda_res_optilcms_complete.qs")
  library(OptiLCMS2ID);
  
  Inchikeys <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  mSet2 <- PerformResultsExport (mSet2,
                                 type = 1L,
                                 topN = 10L,
                                 ncores = 6L);
  
  cmpdsNMs <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  
  identifiedCMPD_dt <- data.frame(Inchikey=Inchikeys, CompoundName = cmpdsNMs)
  #identifiedCMPD_dt <- identifiedCMPD_dt[!is.na(identifiedCMPD_dt$Inchikey), ]
  
  mSet <- qs::qread("HILICpos/OptiLCMS_ms2/ms1/mSet1.qs")
  dt <- mSet@dataSet
  idx_plasma <- which(dt[1,] == "Plasma")
  idx_serum <- which(dt[1,] == "Serum")
  idx_wb <- which(dt[1,] == "Whole_blood")
  idx_qc <- which(dt[1,] == "QC")
  feadt0 <- dt[-1,]
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& ((sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3) | is.na(sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3))
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICpos/unique_ms1_features_idx_optilcms.qs")
  mzs <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][1])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  rts <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][2])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  
  ft_dt_1a <- as.matrix(data.frame(mzmin = mSet@peakAnnotation[["camera_output"]][["mzmin"]],
                                   mzmax = mSet@peakAnnotation[["camera_output"]][["mzmax"]],
                                   rtmin = mSet@peakAnnotation[["camera_output"]][["rtmin"]],
                                   rtmax = mSet@peakAnnotation[["camera_output"]][["rtmax"]]))
  
  ft_dt_1b <- ft_dt_1a[mSet2@MSnResults[["Concensus_spec"]][[1]]+1,]
  Confirmed_idx <- vector(length = nrow(ft_dt_1b))
  for(i in 1:nrow(ft_dt_1b)){
    for(j in 1:length(mzs)){
      if((mzs[j] > ft_dt_1b[i,1]-0.002) &
         (mzs[j] < ft_dt_1b[i,2]+0.002) &
         (rts[j] > ft_dt_1b[i,3]-5) &
         (rts[j] < ft_dt_1b[i,4]+5)){
        Confirmed_idx[i] <- T;
        break;
      } 
    }
  }
  
  identifiedCMPD_dt0 <- identifiedCMPD_dt[Confirmed_idx, ]
  identifiedCMPD_dtx <- identifiedCMPD_dt0[!is.na(identifiedCMPD_dt0$CompoundName), ]
  save(identifiedCMPD_dtx, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICpos/Optilcms_resDT.rda")
  # 2. MSDIAL/MSFINDER
  rm(list = ls())
  resdt <- read.csv("HILICpos/MSDIAL_MS2/dda_ms2/Structure result-2059.txt", sep = "\t")
  resdt <- resdt[resdt$Rank == 1, ]
  
  feadt <- read.csv("HILICpos/MSDIAL_MS2/dda_ms1_ms2/Height_0_2023111213.txt", sep = "\t")
  feadt <- feadt[-c(1:3),-c(1,4,5:32, ncol(feadt), ncol(feadt)-1)]
  idx <- grepl(pattern = "WB[0-9]+|PLASMA[0-9]+|SERUM[0-9]+|QC[0-9]+", x = feadt[1,])
  feadt1 <- feadt[,idx]
  feadt0 <- feadt1[-1,]
  colnms <- feadt1[1,]
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))
  
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICpos/unique_ms1_features_idx_msdial.qs")
  feadt_res <- feadt0[row_idx, ]
  feat_meta <- feadt[which(row_idx)+1,c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(resdt))
  for(k in 1:nrow(resdt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(resdt[k,5] - as.numeric(feat_meta$X.2[s])) < 0.005){
        this_rt <- as.numeric(strsplit(resdt$File.name[k], "_")[[1]][2])
        rt_val <-  as.numeric(feat_meta$X.1[s])
        if(abs(rt_val - this_rt) < 1/3){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  resdt_confirmed <- resdt[Confirmed_idx,]
  save(resdt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICpos/MSDIAL_resDT.rda")
  # 3. MZmine/SIRIUS
  rm(list = ls())
  dt <- read.csv("HILICpos/SIRIUS_ms2/mzmine_dda_ms2_res/res/compound_identifications.tsv", sep = "\t")
  labels <- paste0(dt$ionMass, "__", dt$retentionTimeInSeconds)
  labelsu <- unique(labels)
  idxx <- vapply(labelsu, function(x){
    which(x == labels)[1]
  }, FUN.VALUE = integer(1L), USE.NAMES = F)
  dt <- dt[idxx, ]
  #dt <- dt[c(dt$formulaRank==1 & dt$X.predictedFPs==1),]
  dt <- dt[c(dt$formulaRank==1),]
  
  dt_ms1 <- read.csv("HILICpos/SIRIUS_ms2/mzmine_dda_ms1_res/HILICpos_ms1.csv")
  colnm_idx <- which(grepl("*.mzML.area", colnames(dt_ms1)))
  colnm_idx2 <- which((colnames(dt_ms1) == "rt") | (colnames(dt_ms1) == "mz"))
  dt_ms1x <- dt_ms1[,c(colnm_idx2,colnm_idx)]
  
  colnms <- colnames(dt_ms1x)
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))
  
  row_idx <- apply(dt_ms1x, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICpos/unique_ms1_features_idx_mzmine.qs")
  feadt_res <- dt_ms1x[row_idx, c(1,2, idx_plasma, idx_serum, idx_wb, idx_qc)]
  feat_meta <- dt_ms1x[which(row_idx),c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(dt))
  for(k in 1:nrow(dt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(dt$ionMass[k] - as.numeric(feat_meta$mz[s])) < 0.005){
        this_rt <- as.numeric(dt$retentionTimeInSeconds[k])
        rt_val <-  as.numeric(feat_meta$rt[s])*60
        if(abs(rt_val - this_rt) < 20){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  dt_confirmed <- dt[Confirmed_idx,]
  save(dt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICpos/Mzmine_sirius_resDT.rda")
}

### - HILIC NEG
{
  # 1. OptiLCMS
  rm(list = ls())
  mSet2 <- qs::qread("HILICneg/mzML/MS2/dda_res_optilcms_complete.qs")
  library(OptiLCMS2ID);
  
  Inchikeys <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  mSet2 <- PerformResultsExport (mSet2,
                                 type = 1L,
                                 topN = 10L,
                                 ncores = 6L);
  
  cmpdsNMs <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  
  identifiedCMPD_dt <- data.frame(Inchikey=Inchikeys, CompoundName = cmpdsNMs)
  #identifiedCMPD_dt <- identifiedCMPD_dt[!is.na(identifiedCMPD_dt$Inchikey), ]
  
  mSet <- qs::qread("HILICneg/OptiLCMS_ms2/ms1/mSet1.qs")
  dt <- mSet@dataSet
  idx_plasma <- which(dt[1,] == "Plasma")
  idx_serum <- which(dt[1,] == "Serum")
  idx_wb <- which(dt[1,] == "Whole_blood")
  idx_qc <- which(dt[1,] == "QC")
  feadt0 <- dt[-1,]
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& ((sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3) | is.na(sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3))
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICneg/unique_ms1_features_idx_optilcms.qs")
  mzs <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][1])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  rts <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][2])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  
  ft_dt_1a <- as.matrix(data.frame(mzmin = mSet@peakAnnotation[["camera_output"]][["mzmin"]],
                                   mzmax = mSet@peakAnnotation[["camera_output"]][["mzmax"]],
                                   rtmin = mSet@peakAnnotation[["camera_output"]][["rtmin"]],
                                   rtmax = mSet@peakAnnotation[["camera_output"]][["rtmax"]]))
  
  ft_dt_1b <- ft_dt_1a[mSet2@MSnResults[["Concensus_spec"]][[1]]+1,]
  Confirmed_idx <- vector(length = nrow(ft_dt_1b))
  for(i in 1:nrow(ft_dt_1b)){
    for(j in 1:length(mzs)){
      if((mzs[j] > ft_dt_1b[i,1]-0.002) &
         (mzs[j] < ft_dt_1b[i,2]+0.002) &
         (rts[j] > ft_dt_1b[i,3]-5) &
         (rts[j] < ft_dt_1b[i,4]+5)){
        Confirmed_idx[i] <- T;
        break;
      } 
    }
  }
  
  identifiedCMPD_dt0 <- identifiedCMPD_dt[Confirmed_idx, ]
  identifiedCMPD_dtx <- identifiedCMPD_dt0[!is.na(identifiedCMPD_dt0$CompoundName), ]
  save(identifiedCMPD_dtx, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICneg/Optilcms_resDT.rda")
  # 2. MSDIAL/MSFINDER
  rm(list = ls())
  resdt <- read.csv("HILICneg/MSDIAL_MS2/dda_ms2/Structure result-2101.txt", sep = "\t")
  resdt <- resdt[resdt$Rank == 1, ]
  
  feadt <- read.csv("HILICneg/MSDIAL_MS2/dda_ms1_ms2/Height_0_20231111737.txt", sep = "\t")
  feadt <- feadt[-c(1:3),-c(1,4,5:32, ncol(feadt), ncol(feadt)-1)]
  idx <- grepl(pattern = "WB[0-9]+|PLASMA[0-9]+|SERUM[0-9]+|QC[0-9]+", x = feadt[1,])
  feadt1 <- feadt[,idx]
  feadt0 <- feadt1[-1,]
  colnms <- feadt1[1,]
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))
  
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICneg/unique_ms1_features_idx_msdial.qs")
  feadt_res <- feadt0[row_idx, ]
  feat_meta <- feadt[which(row_idx)+1,c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(resdt))
  for(k in 1:nrow(resdt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(resdt[k,5] - as.numeric(feat_meta$X.2[s])) < 0.005){
        this_rt <- as.numeric(strsplit(resdt$File.name[k], "_")[[1]][2])
        rt_val <-  as.numeric(feat_meta$X.1[s])
        if(abs(rt_val - this_rt) < 1/3){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  resdt_confirmed <- resdt[Confirmed_idx,]
  save(resdt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICneg/MSDIAL_resDT.rda")
  # 3. MZmine/SIRIUS
  rm(list = ls())
  dt <- read.csv("HILICneg/SIRIUS_ms2/mzmine_dda_ms2_res/res/compound_identifications.tsv", sep = "\t")
  labels <- paste0(dt$ionMass, "__", dt$retentionTimeInSeconds)
  labelsu <- unique(labels)
  idxx <- vapply(labelsu, function(x){
    which(x == labels)[1]
  }, FUN.VALUE = integer(1L), USE.NAMES = F)
  dt <- dt[idxx, ]
  #dt <- dt[c(dt$formulaRank==1 & dt$X.predictedFPs==1),]
  dt <- dt[c(dt$formulaRank==1),]
  
  dt_ms1 <- read.csv("HILICneg/SIRIUS_ms2/mzmine_dda_ms1_res/HILICneg_ms1.csv")
  colnm_idx <- which(grepl("*.mzML.area", colnames(dt_ms1)))
  colnm_idx2 <- which((colnames(dt_ms1) == "rt") | (colnames(dt_ms1) == "mz"))
  dt_ms1x <- dt_ms1[,c(colnm_idx2,colnm_idx)]
  
  colnms <- colnames(dt_ms1x)
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))
  
  row_idx <- apply(dt_ms1x, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4)))) 
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICneg/unique_ms1_features_idx_mzmine.qs")
  feadt_res <- dt_ms1x[row_idx, c(1,2, idx_plasma, idx_serum, idx_wb, idx_qc)]
  feat_meta <- dt_ms1x[which(row_idx),c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(dt))
  for(k in 1:nrow(dt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(dt$ionMass[k] - as.numeric(feat_meta$mz[s])) < 0.005){
        this_rt <- as.numeric(dt$retentionTimeInSeconds[k])
        rt_val <-  as.numeric(feat_meta$rt[s])*60
        if(abs(rt_val - this_rt) < 20){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  dt_confirmed <- dt[Confirmed_idx,]
  save(dt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DDA/HILICneg/Mzmine_sirius_resDT.rda")
}


## DIA
### - C18 POS
{
  # 1. OptiLCMS
  rm(list = ls())
  mSet2 <- qs::qread("C18pos/mzML/MS2/dia_res_optilcms_complete.qs")
  library(OptiLCMS2ID);

  Inchikeys <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  mSet2 <- PerformResultsExport (mSet2,
                                 type = 1L,
                                 topN = 10L,
                                 ncores = 6L);

  cmpdsNMs <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))

  identifiedCMPD_dt <- data.frame(Inchikey=Inchikeys, CompoundName = cmpdsNMs)
  #identifiedCMPD_dt <- identifiedCMPD_dt[!is.na(identifiedCMPD_dt$Inchikey), ]

  mSet <- qs::qread("C18pos/OptiLCMS_ms2/ms1/mSet1.qs")
  dt <- mSet@dataSet
  idx_plasma <- which(dt[1,] == "Plasma")
  idx_serum <- which(dt[1,] == "Serum")
  idx_wb <- which(dt[1,] == "Whole_blood")
  idx_qc <- which(dt[1,] == "QC")
  feadt0 <- dt[-1,]
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& ((sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3) | is.na(sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3))
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18pos/unique_ms1_features_idx_optilcms.qs")
  mzs <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][1])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  rts <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][2])}, FUN.VALUE = numeric(1L), USE.NAMES = F)

  ft_dt_1a <- as.matrix(data.frame(mzmin = mSet@peakAnnotation[["camera_output"]][["mzmin"]],
                                   mzmax = mSet@peakAnnotation[["camera_output"]][["mzmax"]],
                                   rtmin = mSet@peakAnnotation[["camera_output"]][["rtmin"]],
                                   rtmax = mSet@peakAnnotation[["camera_output"]][["rtmax"]]))

  ft_dt_1b <- ft_dt_1a[mSet2@MSnResults[["Concensus_spec"]][[1]]+1,]
  Confirmed_idx <- vector(length = nrow(ft_dt_1b))
  for(i in 1:nrow(ft_dt_1b)){
    for(j in 1:length(mzs)){
      if((mzs[j] > ft_dt_1b[i,1]-0.002) &
         (mzs[j] < ft_dt_1b[i,2]+0.002) &
         (rts[j] > ft_dt_1b[i,3]-5) &
         (rts[j] < ft_dt_1b[i,4]+5)){
        Confirmed_idx[i] <- T;
        break;
      }
    }
  }

  identifiedCMPD_dt0 <- identifiedCMPD_dt[Confirmed_idx, ]
  identifiedCMPD_dtx <- identifiedCMPD_dt0[!is.na(identifiedCMPD_dt0$CompoundName), ]
  save(identifiedCMPD_dtx, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18pos/Optilcms_resDT.rda")
  # 2. MSDIAL/MSFINDER
  rm(list = ls())
  resdt <- read.csv("C18pos/MSDIAL_ms2/dia_ms2/Structure result-2085.txt", sep = "\t")
  resdt <- resdt[resdt$Rank == 1, ]

  feadt <- read.csv("C18pos/MSDIAL_ms2/dda_ms1_Ms2/Height_0_2023111610.txt", sep = "\t")
  feadt <- feadt[-c(1:3),-c(1,4,5:32, ncol(feadt), ncol(feadt)-1)]
  idx <- grepl(pattern = "WB[0-9]+|PLASMA[0-9]+|SERUM[0-9]+|QC[0-9]+", x = feadt[1,])
  feadt1 <- feadt[,idx]
  feadt0 <- feadt1[-1,]
  colnms <- feadt1[1,]
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))

  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18pos/unique_ms1_features_idx_msdial.qs")
  feadt_res <- feadt0[row_idx, ]
  feat_meta <- feadt[which(row_idx)+1,c(1,2)]

  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(resdt))
  for(k in 1:nrow(resdt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(resdt[k,5] - as.numeric(feat_meta$X.2[s])) < 0.005){
        this_rt <- as.numeric(strsplit(resdt$File.name[k], "_")[[1]][2])
        rt_val <-  as.numeric(feat_meta$X.1[s])
        if(abs(rt_val - this_rt) < 1/3){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }

  resdt_confirmed <- resdt[Confirmed_idx,]
  save(resdt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18pos/MSDIAL_resDT.rda")
  # 3. XCMS/SIRIUS
  rm(list = ls())
  dt <- read.csv("C18pos/SIRIUS_ms2/xcms_dia_ms2_res/res/compound_identifications.tsv", sep = "\t")
  labels <- paste0(dt$ionMass, "__", dt$retentionTimeInSeconds)
  labelsu <- unique(labels)
  idxx <- vapply(labelsu, function(x){
    which(x == labels)[1]
  }, FUN.VALUE = integer(1L), USE.NAMES = F)
  dt <- dt[idxx, ]
  #dt <- dt[c(dt$formulaRank==1 & dt$X.predictedFPs==1),]
  dt <- dt[c(dt$formulaRank==1),]

  obj_ms1 <- qs::qread("C18pos/SIRIUS_ms2/xcms_dia_ms1_res/xcms_res.qs")
  dt_ms1 <- obj_ms1@assays@data@listData[["raw"]]
  dt_ms1x <- dt_ms1

  colnms <- colnames(dt_ms1x)
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))

  row_idx <- apply(dt_ms1x, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18pos/unique_ms1_features_idx_xcms.qs")
  feadt_res <- dt_ms1x[row_idx, c(idx_plasma, idx_serum, idx_wb, idx_qc)]
  feat_meta <- data.frame(mz = obj_ms1@elementMetadata@listData[["mzmed"]],
                          rt = obj_ms1@elementMetadata@listData[["rtmed"]])#dt_ms1x[which(row_idx),c(1,2)]

  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(dt))
  for(k in 1:nrow(dt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(dt$ionMass[k] - as.numeric(feat_meta$mz[s])) < 0.005){
        this_rt <- as.numeric(dt$retentionTimeInSeconds[k])
        rt_val <-  as.numeric(feat_meta$rt[s])
        if(abs(rt_val - this_rt) < 20){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }

  dt_confirmed <- dt[Confirmed_idx,]
  save(dt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18pos/XCMS_sirius_resDT.rda")

}

### - C18 NEG
{
  # 1. OptiLCMS
  rm(list = ls())
  mSet2 <- qs::qread("C18neg/mzML/MS2/dia_res_optilcms_complete.qs")
  library(OptiLCMS2ID);

  Inchikeys <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  mSet2 <- PerformResultsExport (mSet2,
                                 type = 1L,
                                 topN = 10L,
                                 ncores = 8L);

  cmpdsNMs <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))

  identifiedCMPD_dt <- data.frame(Inchikey=Inchikeys, CompoundName = cmpdsNMs)
  #identifiedCMPD_dt <- identifiedCMPD_dt[!is.na(identifiedCMPD_dt$Inchikey), ]

  mSet <- qs::qread("C18neg/OptiLCMS_ms2/ms1/mSet1.qs")
  dt <- mSet@dataSet
  idx_plasma <- which(dt[1,] == "Plasma")
  idx_serum <- which(dt[1,] == "Serum")
  idx_wb <- which(dt[1,] == "Whole_blood")
  idx_qc <- which(dt[1,] == "QC")
  feadt0 <- dt[-1,]
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& ((sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3) | is.na(sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3))
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18neg/unique_ms1_features_idx_optilcms.qs")
  mzs <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][1])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  rts <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][2])}, FUN.VALUE = numeric(1L), USE.NAMES = F)

  ft_dt_1a <- as.matrix(data.frame(mzmin = mSet@peakAnnotation[["camera_output"]][["mzmin"]],
                                   mzmax = mSet@peakAnnotation[["camera_output"]][["mzmax"]],
                                   rtmin = mSet@peakAnnotation[["camera_output"]][["rtmin"]],
                                   rtmax = mSet@peakAnnotation[["camera_output"]][["rtmax"]]))

  ft_dt_1b <- ft_dt_1a[mSet2@MSnResults[["Concensus_spec"]][[1]]+1,]
  Confirmed_idx <- vector(length = nrow(ft_dt_1b))
  for(i in 1:nrow(ft_dt_1b)){
    for(j in 1:length(mzs)){
      if((mzs[j] > ft_dt_1b[i,1]-0.002) &
         (mzs[j] < ft_dt_1b[i,2]+0.002) &
         (rts[j] > ft_dt_1b[i,3]-5) &
         (rts[j] < ft_dt_1b[i,4]+5)){
        Confirmed_idx[i] <- T;
        break;
      }
    }
  }

  identifiedCMPD_dt0 <- identifiedCMPD_dt[Confirmed_idx, ]
  identifiedCMPD_dtx <- identifiedCMPD_dt0[!is.na(identifiedCMPD_dt0$CompoundName), ]
  save(identifiedCMPD_dtx, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18neg/Optilcms_resDT.rda")
  # 2. MSDIAL/MSFINDER
  rm(list = ls())
  resdt <- read.csv("C18neg/MSDIAL_ms2/dia_ms2_res/Structure result-2066.txt", sep = "\t")
  resdt <- resdt[resdt$Rank == 1, ]

  feadt <- read.csv("C18neg/MSDIAL_ms2/dda_ms1_ms2/Height_0_20231111553.txt", sep = "\t")
  feadt <- feadt[-c(1:3),-c(1,4,5:32, ncol(feadt), ncol(feadt)-1)]
  idx <- grepl(pattern = "WB[0-9]+|PLASMA[0-9]+|SERUM[0-9]+|QC[0-9]+", x = feadt[1,])
  feadt1 <- feadt[,idx]
  feadt0 <- feadt1[-1,]
  colnms <- feadt1[1,]
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))

  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18neg/unique_ms1_features_idx_msdial.qs")
  feadt_res <- feadt0[row_idx, ]
  feat_meta <- feadt[which(row_idx)+1,c(1,2)]

  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(resdt))
  for(k in 1:nrow(resdt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(resdt[k,5] - as.numeric(feat_meta$X.2[s])) < 0.005){
        this_rt <- as.numeric(strsplit(resdt$File.name[k], "_")[[1]][2])
        rt_val <-  as.numeric(feat_meta$X.1[s])
        if(abs(rt_val - this_rt) < 1/3){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }

  resdt_confirmed <- resdt[Confirmed_idx,]
  save(resdt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18neg/MSDIAL_resDT.rda")
  # 3. XCMS/SIRIUS
  rm(list = ls())
  dt <- read.csv("C18neg/SIRIUS_ms2/xcms_dia_ms2_res/res/compound_identifications.tsv", sep = "\t")
  labels <- paste0(dt$ionMass, "__", dt$retentionTimeInSeconds)
  labelsu <- unique(labels)
  idxx <- vapply(labelsu, function(x){
    which(x == labels)[1]
  }, FUN.VALUE = integer(1L), USE.NAMES = F)
  dt <- dt[idxx, ]
  #dt <- dt[c(dt$formulaRank==1 & dt$X.predictedFPs==1),]
  dt <- dt[c(dt$formulaRank==1),]
  
  obj_ms1 <- qs::qread("C18neg/SIRIUS_ms2/xcms_dia_ms1_res/xcms_res.qs")
  dt_ms1 <- obj_ms1@assays@data@listData[["raw"]]
  dt_ms1x <- dt_ms1
  
  colnms <- colnames(dt_ms1x)
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))
  
  row_idx <- apply(dt_ms1x, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18neg/unique_ms1_features_idx_xcms.qs")
  feadt_res <- dt_ms1x[row_idx, c(idx_plasma, idx_serum, idx_wb, idx_qc)]
  feat_meta <- data.frame(mz = obj_ms1@elementMetadata@listData[["mzmed"]],
                          rt = obj_ms1@elementMetadata@listData[["rtmed"]])#dt_ms1x[which(row_idx),c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(dt))
  for(k in 1:nrow(dt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(dt$ionMass[k] - as.numeric(feat_meta$mz[s])) < 0.005){
        this_rt <- as.numeric(dt$retentionTimeInSeconds[k])
        rt_val <-  as.numeric(feat_meta$rt[s])
        if(abs(rt_val - this_rt) < 20){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  dt_confirmed <- dt[Confirmed_idx,]
  save(dt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DIA/C18neg/XCMS_sirius_resDT.rda")
  
}

### - HILIC POS
{
  # 1. OptiLCMS
  rm(list = ls())
  mSet2 <- qs::qread("HILICpos/mzML/MS2/dia_res_optilcms_complete.qs")
  library(OptiLCMS2ID);

  Inchikeys <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  mSet2 <- PerformResultsExport (mSet2,
                                 type = 1L,
                                 topN = 10L,
                                 ncores = 8L);

  cmpdsNMs <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))

  identifiedCMPD_dt <- data.frame(Inchikey=Inchikeys, CompoundName = cmpdsNMs)
  #identifiedCMPD_dt <- identifiedCMPD_dt[!is.na(identifiedCMPD_dt$Inchikey), ]

  mSet <- qs::qread("HILICpos/OptiLCMS_ms2/ms1/mSet1.qs")
  dt <- mSet@dataSet
  idx_plasma <- which(dt[1,] == "Plasma")
  idx_serum <- which(dt[1,] == "Serum")
  idx_wb <- which(dt[1,] == "Whole_blood")
  idx_qc <- which(dt[1,] == "QC")
  feadt0 <- dt[-1,]
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& ((sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3) | is.na(sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3))
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICpos/unique_ms1_features_idx_optilcms.qs")
  mzs <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][1])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  rts <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][2])}, FUN.VALUE = numeric(1L), USE.NAMES = F)

  ft_dt_1a <- as.matrix(data.frame(mzmin = mSet@peakAnnotation[["camera_output"]][["mzmin"]],
                                   mzmax = mSet@peakAnnotation[["camera_output"]][["mzmax"]],
                                   rtmin = mSet@peakAnnotation[["camera_output"]][["rtmin"]],
                                   rtmax = mSet@peakAnnotation[["camera_output"]][["rtmax"]]))

  ft_dt_1b <- ft_dt_1a[mSet2@MSnResults[["Concensus_spec"]][[1]]+1,]
  Confirmed_idx <- vector(length = nrow(ft_dt_1b))
  for(i in 1:nrow(ft_dt_1b)){
    for(j in 1:length(mzs)){
      if((mzs[j] > ft_dt_1b[i,1]-0.002) &
         (mzs[j] < ft_dt_1b[i,2]+0.002) &
         (rts[j] > ft_dt_1b[i,3]-5) &
         (rts[j] < ft_dt_1b[i,4]+5)){
        Confirmed_idx[i] <- T;
        break;
      }
    }
  }

  identifiedCMPD_dt0 <- identifiedCMPD_dt[Confirmed_idx, ]
  identifiedCMPD_dtx <- identifiedCMPD_dt0[!is.na(identifiedCMPD_dt0$CompoundName), ]
  save(identifiedCMPD_dtx, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICpos/Optilcms_resDT.rda")
  # 2. MSDIAL/MSFINDER
  rm(list = ls())
  resdt <- read.csv("HILICpos/MSDIAL_MS2/dia_ms2/Structure result-2089.txt", sep = "\t")
  resdt <- resdt[resdt$Rank == 1, ]

  feadt <- read.csv("HILICpos/MSDIAL_MS2/dda_ms1_ms2/Height_0_2023111213.txt", sep = "\t")
  feadt <- feadt[-c(1:3),-c(1,4,5:32, ncol(feadt), ncol(feadt)-1)]
  idx <- grepl(pattern = "WB[0-9]+|PLASMA[0-9]+|SERUM[0-9]+|QC[0-9]+", x = feadt[1,])
  feadt1 <- feadt[,idx]
  feadt0 <- feadt1[-1,]
  colnms <- feadt1[1,]
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))

  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICpos/unique_ms1_features_idx_msdial.qs")
  feadt_res <- feadt0[row_idx, ]
  feat_meta <- feadt[which(row_idx)+1,c(1,2)]

  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(resdt))
  for(k in 1:nrow(resdt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(resdt[k,5] - as.numeric(feat_meta$X.2[s])) < 0.005){
        this_rt <- as.numeric(strsplit(resdt$File.name[k], "_")[[1]][2])
        rt_val <-  as.numeric(feat_meta$X.1[s])
        if(abs(rt_val - this_rt) < 1/3){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }

  resdt_confirmed <- resdt[Confirmed_idx,]
  save(resdt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICpos/MSDIAL_resDT.rda")
  # 3. XCMS/SIRIUS
  rm(list = ls())
  dt <- read.csv("HILICpos/SIRIUS_ms2/xcms_dia_ms2_res/res/compound_identifications.tsv", sep = "\t")
  labels <- paste0(dt$ionMass, "__", dt$retentionTimeInSeconds)
  labelsu <- unique(labels)
  idxx <- vapply(labelsu, function(x){
    which(x == labels)[1]
  }, FUN.VALUE = integer(1L), USE.NAMES = F)
  dt <- dt[idxx, ]
  #dt <- dt[c(dt$formulaRank==1 & dt$X.predictedFPs==1),]
  dt <- dt[c(dt$formulaRank==1),]
  
  obj_ms1 <- qs::qread("HILICpos/SIRIUS_ms2/xcms_dia_ms1_res/xcms_res.qs")
  dt_ms1 <- obj_ms1@assays@data@listData[["raw"]]
  dt_ms1x <- dt_ms1
  
  colnms <- colnames(dt_ms1x)
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))
  
  row_idx <- apply(dt_ms1x, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICpos/unique_ms1_features_idx_xcms.qs")
  feadt_res <- dt_ms1x[row_idx, c(idx_plasma, idx_serum, idx_wb, idx_qc)]
  feat_meta <- data.frame(mz = obj_ms1@elementMetadata@listData[["mzmed"]],
                          rt = obj_ms1@elementMetadata@listData[["rtmed"]])#dt_ms1x[which(row_idx),c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(dt))
  for(k in 1:nrow(dt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(dt$ionMass[k] - as.numeric(feat_meta$mz[s])) < 0.005){
        this_rt <- as.numeric(dt$retentionTimeInSeconds[k])
        rt_val <-  as.numeric(feat_meta$rt[s])
        if(abs(rt_val - this_rt) < 20){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  dt_confirmed <- dt[Confirmed_idx,]
  save(dt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICpos/XCMS_sirius_resDT.rda")
}

### - HILIC NEG
{
  # 1. OptiLCMS
  rm(list = ls())
  mSet2 <- qs::qread("HILICneg/mzML/MS2/dia_res_optilcms_complete.qs")
  library(OptiLCMS2ID);

  Inchikeys <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))
  mSet2 <- PerformResultsExport (mSet2,
                                 type = 1L,
                                 topN = 10L,
                                 ncores = 6L);

  cmpdsNMs <- vapply(1:length(mSet2@MSnResults[["DBAnnoteRes"]]), function(x){
    mSet2@MSnResults[["DBAnnoteRes"]][[x]][[1]][[2]][1]
  }, FUN.VALUE = character(1L))

  identifiedCMPD_dt <- data.frame(Inchikey=Inchikeys, CompoundName = cmpdsNMs)
  #identifiedCMPD_dt <- identifiedCMPD_dt[!is.na(identifiedCMPD_dt$Inchikey), ]

  mSet <- qs::qread("HILICneg/OptiLCMS_ms2/ms1/mSet1.qs")
  dt <- mSet@dataSet
  idx_plasma <- which(dt[1,] == "Plasma")
  idx_serum <- which(dt[1,] == "Serum")
  idx_wb <- which(dt[1,] == "Whole_blood")
  idx_qc <- which(dt[1,] == "QC")
  feadt0 <- dt[-1,]
  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& ((sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3) | is.na(sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3))
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICneg/unique_ms1_features_idx_optilcms.qs")
  mzs <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][1])}, FUN.VALUE = numeric(1L), USE.NAMES = F)
  rts <- vapply(feadt0[row_idx, 1], function(x){as.numeric(strsplit(x, "__")[[1]][2])}, FUN.VALUE = numeric(1L), USE.NAMES = F)

  ft_dt_1a <- as.matrix(data.frame(mzmin = mSet@peakAnnotation[["camera_output"]][["mzmin"]],
                                   mzmax = mSet@peakAnnotation[["camera_output"]][["mzmax"]],
                                   rtmin = mSet@peakAnnotation[["camera_output"]][["rtmin"]],
                                   rtmax = mSet@peakAnnotation[["camera_output"]][["rtmax"]]))

  ft_dt_1b <- ft_dt_1a[mSet2@MSnResults[["Concensus_spec"]][[1]]+1,]
  Confirmed_idx <- vector(length = nrow(ft_dt_1b))
  for(i in 1:nrow(ft_dt_1b)){
    for(j in 1:length(mzs)){
      if((mzs[j] > ft_dt_1b[i,1]-0.002) &
         (mzs[j] < ft_dt_1b[i,2]+0.002) &
         (rts[j] > ft_dt_1b[i,3]-5) &
         (rts[j] < ft_dt_1b[i,4]+5)){
        Confirmed_idx[i] <- T;
        break;
      }
    }
  }

  identifiedCMPD_dt0 <- identifiedCMPD_dt[Confirmed_idx, ]
  identifiedCMPD_dtx <- identifiedCMPD_dt0[!is.na(identifiedCMPD_dt0$CompoundName), ]
  save(identifiedCMPD_dtx, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICneg/Optilcms_resDT.rda")
  # 2. MSDIAL/MSFINDER
  rm(list = ls())
  resdt <- read.csv("HILICneg/MSDIAL_MS2/dia_ms2/Structure result-2085.txt", sep = "\t")
  resdt <- resdt[resdt$Rank == 1, ]

  feadt <- read.csv("HILICneg/MSDIAL_MS2/dda_ms1_ms2/Height_0_20231111737.txt", sep = "\t")
  feadt <- feadt[-c(1:3),-c(1,4,5:32, ncol(feadt), ncol(feadt)-1)]
  idx <- grepl(pattern = "WB[0-9]+|PLASMA[0-9]+|SERUM[0-9]+|QC[0-9]+", x = feadt[1,])
  feadt1 <- feadt[,idx]
  feadt0 <- feadt1[-1,]
  colnms <- feadt1[1,]
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))

  row_idx <- apply(feadt0, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICneg/unique_ms1_features_idx_msdial.qs")
  feadt_res <- feadt0[row_idx, ]
  feat_meta <- feadt[which(row_idx)+1,c(1,2)]

  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(resdt))
  for(k in 1:nrow(resdt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(resdt[k,5] - as.numeric(feat_meta$X.2[s])) < 0.005){
        this_rt <- as.numeric(strsplit(resdt$File.name[k], "_")[[1]][2])
        rt_val <-  as.numeric(feat_meta$X.1[s])
        if(abs(rt_val - this_rt) < 1/3){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }

  resdt_confirmed <- resdt[Confirmed_idx,]
  save(resdt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICneg/MSDIAL_resDT.rda")
  # 3. XCMS/SIRIUS
  rm(list = ls())
  dt <- read.csv("HILICneg/SIRIUS_ms2/xcms_dia_ms2_res/res/compound_identifications.tsv", sep = "\t")
  labels <- paste0(dt$ionMass, "__", dt$retentionTimeInSeconds)
  labelsu <- unique(labels)
  idxx <- vapply(labelsu, function(x){
    which(x == labels)[1]
  }, FUN.VALUE = integer(1L), USE.NAMES = F)
  dt <- dt[idxx, ]
  #dt <- dt[c(dt$formulaRank==1 & dt$X.predictedFPs==1),]
  dt <- dt[c(dt$formulaRank==1),]

  #dt_ms1 <- read.csv("HILICneg/SIRIUS_ms2/mzmine_dda_ms1_res/HILICneg_ms1.csv")
  #colnm_idx <- which(grepl("*.mzML.area", colnames(dt_ms1)))
  #dt_ms1x <- dt_ms1[,c(7,9,colnm_idx)]
  obj_ms1 <- qs::qread("HILICneg/SIRIUS_ms2/xcms_dia_ms1_res/xcms_res.qs")
  dt_ms1 <- obj_ms1@assays@data@listData[["raw"]]
  dt_ms1x <- dt_ms1

  colnms <- colnames(dt_ms1x)
  idx_plasma <- which(grepl("PLASMA[0-9]+", colnms))
  idx_serum <- which(grepl("SERUM[0-9]+", colnms))
  idx_wb <- which(grepl("WB[0-9]+", colnms))
  idx_qc <- which(grepl("QC[0-9]+", colnms))

  row_idx <- apply(dt_ms1x, 1, function(x){
    # plasma unique, serum unique, wb unique
    (((length(which(x[idx_plasma] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_serum] > 0)) > 10) & ((length(which(x[idx_plasma] > 0)) < 4) | (length(which(x[idx_wb] > 0)) < 4))) |
       ((length(which(x[idx_wb] > 0)) > 10) & ((length(which(x[idx_serum] > 0)) < 4) | (length(which(x[idx_plasma] > 0)) < 4))))
    #& (sd(as.numeric(x[idx_qc]))/mean(as.numeric(x[idx_qc])) < 0.3)
  })
  qs::qsave(row_idx, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICneg/unique_ms1_features_idx_xcms.qs")
  feadt_res <- dt_ms1x[row_idx, c(idx_plasma, idx_serum, idx_wb, idx_qc)]
  feat_meta <- data.frame(mz = obj_ms1@elementMetadata@listData[["mzmed"]],
                          rt = obj_ms1@elementMetadata@listData[["rtmed"]])#dt_ms1x[which(row_idx),c(1,2)]
  
  Identified_idx <- vector(length = length(row_idx))
  Confirmed_idx <- vector(length = nrow(dt))
  for(k in 1:nrow(dt)){
    for(s in 1:nrow(feat_meta)){
      if(abs(dt$ionMass[k] - as.numeric(feat_meta$mz[s])) < 0.005){
        this_rt <- as.numeric(dt$retentionTimeInSeconds[k])
        rt_val <-  as.numeric(feat_meta$rt[s])
        if(abs(rt_val - this_rt) < 20){
          feat_meta <- feat_meta[-s,]
          Identified_idx[s] <- T
          Confirmed_idx[k] <- T
          break;
        }
      }
    }
  }
  
  dt_confirmed <- dt[Confirmed_idx,]
  save(dt_confirmed, file = "/data/ms2_benchmark/whole_blood/results/DIA/HILICneg/Mzmine_sirius_resDT.rda")
}
