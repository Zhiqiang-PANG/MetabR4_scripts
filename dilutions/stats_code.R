### match ms1 -> ms2; 
rm(list = ls())
setwd("/data/ms2_benchmark/whole_blood/dilutions/")

# MetaboAnalyst (OptiLCMS)
{
  setwd("/data/ms2_benchmark/whole_blood/dilutions/1_optilcms")
  rm(list = ls())
  cors_list <- list()
  Concentration <- c(1/(1024+1), 1/(256+1), 1/(64+1), 1/(16+1), 1/(4+1), 1/(1+1), 4/(4+1), 16/(16+1), 64/(64+1), 256/(256+1), 1024/(1024+1))
  Concentrations <- vapply(Concentration, function(x){
    rep(x,3)
  }, FUN.VALUE = double(length = 3L))
  dim(Concentrations) <- c(1,33)
  Concentrations <- as.double(Concentrations)
  
  files <- c("C18pos_ms1/mSet1_c18pos.qs", 
             "C18neg_ms1/mSet1_c18neg.qs",
             "HILICpos_ms1/mSet1_HILICpos.qs", 
             "HILICneg_ms1/mSet1_HILICneg.qs")
  dda_files <- c("C18pos_ms2/DDA/dda_res_optilcms_complete.qs", 
                 "C18neg_ms2/DDA/dda_res_optilcms_complete.qs",
                 "HILICpos_ms2/DDA/dda_res_optilcms_complete.qs", 
                 "HILICneg_ms2/DDA/dda_res_optilcms_complete.qs")
  dia_files1 <- c("C18pos_ms2/DIA1/dia_res_optilcms_complete.qs", 
                  "C18neg_ms2/DIA1/dia_res_optilcms_complete.qs",
                  "HILICpos_ms2/DIA1/dia_res_optilcms_complete.qs", 
                  "HILICneg_ms2/DIA1/dia_res_optilcms_complete.qs")
  dia_files2 <- c("C18pos_ms2/DIA2/dia_res_optilcms_complete.qs", 
                  "C18neg_ms2/DIA2/dia_res_optilcms_complete.qs",
                  "HILICpos_ms2/DIA2/dia_res_optilcms_complete.qs", 
                  "HILICneg_ms2/DIA2/dia_res_optilcms_complete.qs")

  for(i in 1:4){
    f <- files[i]
    cat("==f ->", f, "\n")
    mSet <- qs::qread(f)
    df <- mSet@dataSet
    labels <- df$Sample[-1]
    df <- df[-1,grepl("^S[A-Z]+", colnames(df))]

    dfx <- df[,grepl("^SL|^SM", colnames(df))]
    idxr <- apply(dfx, 1, function(x){
      ((length(which(x[1:3] == 0)) > 1) & (length(which(x[4:6] == 0)) <2)) |
        (length(which(x[4:6] == 0)) > 1 & (length(which(x[1:3] == 0)) <2))
    })

    df <- df[idxr,]
    lbls <- labels[idxr]

    cor_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["estimate"]][["cor"]]
    })

    cat("% -->", length(which(abs(cor_vecs)>0.9))/length(cor_vecs)*100, "\n")
    #cors_list <- c(cors_list, list(cor_vecs))
    
    # read in mset2 -- DDA
    fd <- dda_files[i]
    mSet2 <- qs::qread(fd)
    peak_idx0 <- mSet2@MSnResults[["Concensus_spec"]][[1]] + 1
    peak_idx0 <- peak_idx0[vapply(mSet2@MSnResults[["DBAnnoteRes"]], function(x){
      length(x[[1]][["IDs"]]) > 0
      }, FUN.VALUE = logical(1L))]
    idxs_dda <- vapply(which(idxr), function(x){
      x %in% peak_idx0
    }, logical(1L))
    identified_mtx <- data.frame(cors = cor_vecs, identified = idxs_dda)
    qs::qsave(identified_mtx, sub("dda_res_optilcms_complete.qs", "dda_res_optilcms_cors_df.qs", fd))
    
    # read in mset2 -- DIA1
    fdi1 <- dia_files1[i]
    mSet2i1 <- qs::qread(fdi1)
    peak_idx0 <- mSet2i1@MSnResults[["Concensus_spec"]][[1]] + 1
    peak_idx0 <- peak_idx0[vapply(mSet2@MSnResults[["DBAnnoteRes"]], function(x){
      if(length(x[[1]][["IDs"]]) == 0){
        return(FALSE)
      } else {
        return(length(x[[1]][["IDs"]]) > 0 & max(x[[1]][["Dot_Similarity"]], na.rm = T)>0.25)
      }
    }, FUN.VALUE = logical(1L))]
    idxs_dia1 <- vapply(which(idxr), function(x){
      x %in% peak_idx0
    }, logical(1L))
    # read in mset2 -- DIA2
    fdi2 <- dia_files2[i]
    mSet2i2 <- qs::qread(fdi2)
    peak_idx0 <- mSet2i2@MSnResults[["Concensus_spec"]][[1]] + 1
    peak_idx0 <- peak_idx0[vapply(mSet2@MSnResults[["DBAnnoteRes"]], function(x){
      if(length(x[[1]][["IDs"]]) == 0){
        return(FALSE)
      } else {
        return(length(x[[1]][["IDs"]]) > 0 & max(x[[1]][["Dot_Similarity"]], na.rm = T)>0.25)
      }
    }, FUN.VALUE = logical(1L))]
    idxs_dia2 <- vapply(which(idxr), function(x){
      x %in% peak_idx0
    }, logical(1L))
    identified_mtx1 <- data.frame(cors = cor_vecs, identified = (idxs_dia1 | idxs_dia2))
    
    qs::qsave(identified_mtx1, sub("dia_res_optilcms_complete.qs", "dia_res_optilcms_cors_df.qs", fdi1))
  }
  
  
}

# MSDIAL
{
  setwd("/data/ms2_benchmark/whole_blood/dilutions/2_MSDIAL/")
  rm(list = ls())
  cors_list <- list()
  Concentration <- c(1/(1024+1), 1/(256+1), 1/(64+1), 1/(16+1), 1/(4+1), 1/(1+1), 4/(4+1), 16/(16+1), 64/(64+1), 256/(256+1), 1024/(1024+1))
  Concentrations <- vapply(Concentration, function(x){
    rep(x,3)
  }, FUN.VALUE = double(length = 3L))
  dim(Concentrations) <- c(1,33)
  Concentrations <- as.double(Concentrations)
  
  files <- c("C18pos_ms1/Height_0_2023319415.txt", 
             "C18neg_ms1/Height_0_2023319516.txt",
             "HILICpos_ms1/Height_0_20233191319.txt", 
             "HILICneg_ms1/Height_0_20233191423.txt")
  dda_files <- c("C18pos_ms2/dda_res/Structure result-2094.txt", 
                 "C18neg_ms2/dda_res/Structure result-2068.txt",
                 "HILICpos_ms2/dda_res/Structure result-2086.txt", 
                 "HILICneg_ms2/dda_res/Structure result-2082.txt")
  dia_files1 <- c("C18pos_ms2/dia_res1/Structure result-2071.txt", 
                  "C18neg_ms2/dia_res1/Structure result-2073.txt",
                  "HILICpos_ms2/dia_res1/Structure result-2085.txt", 
                  "HILICneg_ms2/dia_res1/Structure result-2115.txt")
  dia_files2 <- c("C18pos_ms2/dia_res2/Structure result-2113.txt", 
                  "C18neg_ms2/dia_res2/Structure result-2071.txt",
                  "HILICpos_ms2/dia_res2/Structure result-2084.txt", 
                  "HILICneg_ms2/dia_res2/Structure result-2114.txt")
  
  for(i in 1:4){
    f <- files[i]
    cat("==f ->", f, "\n")
    dt <- read.csv(f, sep = "\t")
    dt <- dt[-c(1:3),]
    colnames(dt) <- dt[1,]
    dt <- dt[-1,]
    dt <- dt[,-c(1,4:32)]
    
    dt_mz_rt <- dt[,c(1:2)]
    
    df <- dt[-1,grepl("^S[A-Z]+", colnames(dt))]
    
    dfx <- df[,grepl("^SL|^SM", colnames(df))]
    idxr <- apply(dfx, 1, function(x){
      ((length(which(x[1:3] == 0)) > 1) & (length(which(x[4:6] == 0)) <2)) |
        ((length(which(x[4:6] == 0)) > 1) & (length(which(x[1:3] == 0)) <2))
    })
    
    df <- df[idxr,]
    dt_mz_rt <- dt_mz_rt[idxr,]
    
    cor_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["estimate"]][["cor"]]
    })
    
    cat("-->", length(which(abs(cor_vecs)>0.9))/length(cor_vecs)*100, "\n")

    #cors_list <- c(cors_list, list(cor_vecs))
    
    # read in mset2 -- DDA
    fd <- dda_files[i]
    ms2_msfinder <- read.csv(fd, sep = "\t")
    ms2_msfinder <- ms2_msfinder[ms2_msfinder$Rank == 1,]
    
    secs <- as.double(dt_mz_rt$`Average Rt(min)`)*60
    mzs <- as.double(dt_mz_rt$`Average Mz`)
    msfinderrts <- vapply(1:nrow(ms2_msfinder), function(x) {as.double(strsplit(ms2_msfinder$File.name[x], "_")[[1]][2])*60}, double(1L))
    msfindermzs <- vapply(1:nrow(ms2_msfinder), function(x) {as.double(ms2_msfinder$Precursor.mz[x])}, double(1L))
    idxs_dda <- vapply(1:length(secs), function(x){
      thismz <- mzs[x]
      thisrt <- secs[x];
      any((abs(thismz - msfindermzs) < 0.001) & (abs(thisrt - msfinderrts) < 1))
    }, FUN.VALUE = logical(1L))
    
    identified_mtx <- data.frame(cors = cor_vecs, identified = idxs_dda)
    qs::qsave(identified_mtx, sub("Structure result-[0-9]+.txt", "dda_res_msdial_cors_df.qs", fd))
    
    # read in mset2 -- DIA1
    fdi1 <- dia_files1[i]
    ms2_msfinder <- read.csv(fdi1, sep = "\t")
    ms2_msfinder <- ms2_msfinder[ms2_msfinder$Rank == 1,]
    
    secs <- as.double(dt_mz_rt$`Average Rt(min)`)*60
    mzs <- as.double(dt_mz_rt$`Average Mz`)
    msfinderrts <- vapply(1:nrow(ms2_msfinder), function(x) {as.double(strsplit(ms2_msfinder$File.name[x], "_")[[1]][2])*60}, double(1L))
    msfindermzs <- vapply(1:nrow(ms2_msfinder), function(x) {as.double(ms2_msfinder$Precursor.mz[x])}, double(1L))
    idxs_dia1 <- vapply(1:length(secs), function(x){
      thismz <- mzs[x]
      thisrt <- secs[x];
      any((abs(thismz - msfindermzs) < 0.01) & (abs(thisrt - msfinderrts) < 15))
    }, FUN.VALUE = logical(1L))
    
    # read in mset2 -- DIA2
    fdi2 <- dia_files2[i]
    ms2_msfinder <- read.csv(fdi2, sep = "\t")
    ms2_msfinder <- ms2_msfinder[ms2_msfinder$Rank == 1,]
    
    secs <- as.double(dt_mz_rt$`Average Rt(min)`)*60
    mzs <- as.double(dt_mz_rt$`Average Mz`)
    msfinderrts <- vapply(1:nrow(ms2_msfinder), function(x) {as.double(strsplit(ms2_msfinder$File.name[x], "_")[[1]][2])*60}, double(1L))
    msfindermzs <- vapply(1:nrow(ms2_msfinder), function(x) {as.double(ms2_msfinder$Precursor.mz[x])}, double(1L))
    idxs_dia2 <- vapply(1:length(secs), function(x){
      thismz <- mzs[x]
      thisrt <- secs[x];
      any((abs(thismz - msfindermzs) < 0.01) & (abs(thisrt - msfinderrts) < 15))
    }, FUN.VALUE = logical(1L))
    
    identified_mtx1 <- data.frame(cors = cor_vecs, identified = (idxs_dia1 | idxs_dia2))
    
    qs::qsave(identified_mtx1, sub("Structure result-[0-9]+.txt", "dia_res_msdial_cors_df.qs", fdi1))
    
  }
  
  
}

# MzMine - DDA
{
  
  setwd("/data/ms2_benchmark/whole_blood/dilutions/3_mzMine/")
  rm(list = ls())
  cors_list <- list()
  Concentration <- c(1/(1024+1), 1/(256+1), 1/(64+1), 1/(16+1), 1/(4+1), 1/(1+1), 4/(4+1), 16/(16+1), 64/(64+1), 256/(256+1), 1024/(1024+1))
  Concentrations <- vapply(Concentration, function(x){
    rep(x,3)
  }, FUN.VALUE = double(length = 3L))
  dim(Concentrations) <- c(1,33)
  Concentrations <- as.double(Concentrations)
  
  files <- c("C18pos_ms1/ms1_dt.csv", 
             "C18neg_ms1/ms1_dt.csv",
             "HILICpos_ms1/ms1_dt.csv", 
             "HILICneg_ms1/ms1_dt.csv")
  dda_files <- c("C18pos_ms2/res/compound_identifications.tsv", 
                 "C18neg_ms2/res/compound_identifications.tsv",
                 "HILICpos_ms2/res/compound_identifications.tsv", 
                 "HILICneg_ms2/res/compound_identifications.tsv")

  for(i in 1:4){
    f <- files[i]
    cat("==f ->", f, "\n")
    dt <- read.csv(f, sep = ",")
    dtx <- dt[ ,grepl(".area", colnames(dt))]
    rt_mz_col_idx <- which(colnames(dt) == "rt" | colnames(dt) == "mz")
    dt_mz_rt <- dt[, rt_mz_col_idx]
      
    df <- dtx[ ,grepl("^datafile.S[A-Z]+", colnames(dtx))]
    idxc <- vector()
    for(m in LETTERS[1:13]){
      idxc <- c(idxc, which(grepl(paste0("^datafile.S",m), colnames(df))))
    }
    df <- df[,idxc]
    df[is.na(df)] <- 0
    dfx <- df[,grepl("^datafile.SL|^datafile.SM", colnames(df))]
    idxr <- apply(dfx, 1, function(x){
      ((length(which(x[1:3] == 0)) > 1) & (length(which(x[4:6] == 0)) <2)) |
        (length(which(x[4:6] == 0)) > 1 & (length(which(x[1:3] == 0)) <2))
    })
    
    df <- df[idxr,]
    dt_mz_rt <- dt_mz_rt[idxr, ]
    
    cor_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["estimate"]][["cor"]]
    })

    cat("-->", length(which(abs(cor_vecs)>0.9))/length(cor_vecs)*100, "\n")

    #cors_list <- c(cors_list, list(cor_vecs))
    
    # read in mset2 -- DDA
    fd <- dda_files[i]
    ms2_sirius <- read.csv(fd, sep = "\t")
    
    secs <- as.double(dt_mz_rt$rt)*60
    mzs <- as.double(dt_mz_rt$mz)
    mssiriusrts <- vapply(1:nrow(ms2_sirius), function(x) {as.double(ms2_sirius$retentionTimeInSeconds[x])}, double(1L))
    mssiriusmzs <- vapply(1:nrow(ms2_sirius), function(x) {as.double(ms2_sirius$ionMass[x])}, double(1L))
    idxs_dda <- vapply(1:length(secs), function(x){
      thismz <- mzs[x]
      thisrt <- secs[x];
      any((abs(thismz - mssiriusmzs) < 0.01) & (abs(thisrt - mssiriusrts) < 25))
    }, FUN.VALUE = logical(1L))
    print(length(which(idxs_dda)))
    
    identified_mtx <- data.frame(cors = cor_vecs, identified = idxs_dda)
    qs::qsave(identified_mtx, sub("res/compound_identifications.tsv", "dda_res_mzmine_sirius_cors_df.qs", fd))
    

  }
  
}

# XCMS -DIA
{
  
  setwd("/data/ms2_benchmark/whole_blood/dilutions/4_xcms/")
  rm(list = ls())
  cors_list <- list()
  Concentration <- c(1/(1024+1), 1/(256+1), 1/(64+1), 1/(16+1), 1/(4+1), 1/(1+1), 4/(4+1), 16/(16+1), 64/(64+1), 256/(256+1), 1024/(1024+1))
  Concentrations <- vapply(Concentration, function(x){
    rep(x,3)
  }, FUN.VALUE = double(length = 3L))
  dim(Concentrations) <- c(1,33)
  Concentrations <- as.double(Concentrations)
  
  files <- c("C18pos_ms1/xcms_res.qs", 
             "C18neg_ms1/xcms_res.qs",
             "HILICpos_ms1/xcms_res.qs", 
             "HILICneg_ms1/xcms_res.qs")

  dia_files1 <- c("C18pos_ms2/DIA1/compound_identifications.tsv", 
                  "C18neg_ms2/DIA1/compound_identifications.tsv",
                  "HILICpos_ms2/DIA1/compound_identifications.tsv", 
                  "HILICneg_ms2/DIA1/compound_identifications.tsv")
  dia_files2 <- c("C18pos_ms2/DIA2/compound_identifications.tsv", 
                  "C18neg_ms2/DIA2/res/compound_identifications.tsv",
                  "HILICpos_ms2/DIA2/compound_identifications.tsv", 
                  "HILICneg_ms2/DIA2/compound_identifications.tsv")
  
  
  for(i in 1:4){
    f <- files[i]
    cat("==f ->", f, "\n")
    dt <- qs::qread(f)
    dtx <- dt@assays@data@listData[["raw"]]
    dt_mz_rt <- data.frame(mz = dt@elementMetadata@listData[["mzmed"]],
                           rt = dt@elementMetadata@listData[["rtmed"]])
    
    df <- dtx[ ,grepl("^S[A-Z]+", colnames(dtx))]
    df[is.na(df)] <- 0
    dfx <- df[,grepl("^SL|^SM", colnames(df))]
    idxr <- apply(dfx, 1, function(x){
      ((length(which(x[1:3] == 0)) > 1) & (length(which(x[4:6] == 0)) <2)) |
        (length(which(x[4:6] == 0)) > 1 & (length(which(x[1:3] == 0)) <2))
    })
    
    df <- df[idxr,]
    dt_mz_rt <- dt_mz_rt[idxr, ]
    
    cor_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["estimate"]][["cor"]]
    })

    cat("-->", length(which(abs(cor_vecs)>0.9))/length(cor_vecs)*100, "\n")

    # read in mset2 -- DIA1
    fdi1 <- dia_files1[i]
    ms2_sirius <- read.csv(fdi1, sep = "\t")
    
    secs <- as.double(dt_mz_rt$rt)
    mzs <- as.double(dt_mz_rt$mz)
    mssiriusrts <- vapply(1:nrow(ms2_sirius), function(x) {as.double(ms2_sirius$retentionTimeInSeconds[x])}, double(1L))
    mssiriusmzs <- vapply(1:nrow(ms2_sirius), function(x) {as.double(ms2_sirius$ionMass[x])}, double(1L))
    idxs_dia1 <- vapply(1:length(secs), function(x){
      thismz <- mzs[x]
      thisrt <- secs[x];
      any((abs(thismz - mssiriusmzs) < 0.01) & (abs(thisrt - mssiriusrts) < 25))
    }, FUN.VALUE = logical(1L))
    print(length(which(idxs_dia1)))
    
    # read in mset2 -- DIA2
    fdi2 <- dia_files2[i]
    ms2_sirius <- read.csv(fdi2, sep = "\t")
    
    secs <- as.double(dt_mz_rt$rt)
    mzs <- as.double(dt_mz_rt$mz)
    mssiriusrts <- vapply(1:nrow(ms2_sirius), function(x) {as.double(ms2_sirius$retentionTimeInSeconds[x])}, double(1L))
    mssiriusmzs <- vapply(1:nrow(ms2_sirius), function(x) {as.double(ms2_sirius$ionMass[x])}, double(1L))
    idxs_dia2 <- vapply(1:length(secs), function(x){
      thismz <- mzs[x]
      thisrt <- secs[x];
      any((abs(thismz - mssiriusmzs) < 0.01) & (abs(thisrt - mssiriusrts) < 25))
    }, FUN.VALUE = logical(1L))
    print(length(which(idxs_dia2)))
    
    identified_mtx1 <- data.frame(cors = cor_vecs, identified = (idxs_dia1 | idxs_dia2))
    print(length(which(identified_mtx1$identified)))
    
    qs::qsave(identified_mtx1, sub("compound_identifications.tsv", "dia_res_mzmine_sirius_cors_df.qs", fdi1))
  }
  
  
}


# Heatmap : [stats + plotting]
{
  setwd("/data/ms2_benchmark/whole_blood/dilutions/1_optilcms")
  
  files <- c("C18pos_ms1/mSet1_c18pos.qs", "C18neg_ms1/mSet1_c18neg.qs",
             "HILICpos_ms1/mSet1_HILICpos.qs", "HILICneg_ms1/mSet1_HILICneg.qs")
  for(f in files){
    cat("==f ->", f, "\n")
    mSet <- qs::qread(f)
    df <- mSet@dataSet
    df <- df[,grepl("^S[A-Z]+|Sample", colnames(df))]
    dm <- vapply(1:11, function(x){rep(x,3)}, integer(3L))
    df[1,] <- c("Label", as.character(dm), as.character(rep("0",3)), as.character(rep("12",3)))
    rownames(df) <- NULL
    write.csv(df, file = sub(".qs","_peaktable4heatmap.csv",f), row.names = F)
  }
  
  for(k in c("C18pos_ms1", "C18neg_ms1", "HILICpos_ms1", "HILICneg_ms1")){
    cat("Processing ==> ", k ,"\n")
    setwd(paste0("/data/ms2_benchmark/whole_blood/dilutions/1_optilcms/", k ,"/heatmap/"))
    rm(list = ls()[ls() != "k"])
    library(MetaboAnalystR)
    mSet<-InitDataObjects("pktable", "stat", FALSE)
    mSet<-Read.TextData(mSet, list.files("..", pattern = "*_peaktable4heatmap.csv", full.names = T), "colu", "disc");
    mSet[["dataSet"]][["cls"]] <- mSet[["dataSet"]][["orig.cls"]] <- 
      factor(mSet[["dataSet"]][["cls"]], levels = c("0", "1", "2", "3","4","5", "6", "7", "8", "9", "10", "11", "12"))
    mSet<-SanityCheckData(mSet)
    mSet<-ReplaceMin(mSet);
    mSet<-SanityCheckData(mSet)
    #mSet<-FilterVariable(mSet, "none", "F", 25)
    mSet<-PreparePrenormData(mSet)
    mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
    mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
    mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
    mSet<-PlotHeatMap(mSet, "heatmap_0_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", T, T, NULL, T, F, T, T, T)
    mSet<-PlotHeatMap(mSet, "heatmap_1_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", F, T, NULL, T, F, T, T, T)
    mSet<-PlotHeatMap(mSet, "heatmap_0_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", T, T, NULL, T, F, T, T, T)
    mSet<-PlotHeatMap(mSet, "heatmap_1_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", F, T, NULL, T, F, T, T, T)
    
    # Filtering
    mSet<-FilterVariable(mSet, "none", "F", 0)
    mSet<-PreparePrenormData(mSet)
    mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
    mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
    mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
    mSet<-PlotHeatMap(mSet, "heatmap_0_filtered_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", T, T, NULL, T, F, T, T, T)
    mSet<-PlotHeatMap(mSet, "heatmap_1_filtered_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", F, T, NULL, T, F, T, T, T)
    mSet<-PlotHeatMap(mSet, "heatmap_0_filtered_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", T, T, NULL, T, F, T, T, T)
    mSet<-PlotHeatMap(mSet, "heatmap_1_filtered_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", F, T, NULL, T, F, T, T, T)
    
  }
  
    
}



