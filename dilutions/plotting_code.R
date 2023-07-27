
library(ggplot2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(viridis)
library(hrbrthemes)

# Dilution design
{
  Concentration_serum <- c(0, 1/(1024+1), 1/(256+1), 1/(64+1), 1/(16+1), 1/(4+1), 1/(1+1), 4/(4+1), 16/(16+1), 64/(64+1), 256/(256+1), 1024/(1024+1), 1)
  Concentration_urine <- rev(c(0, 1/(1024+1), 1/(256+1), 1/(64+1), 1/(16+1), 1/(4+1), 1/(1+1), 4/(4+1), 16/(16+1), 64/(64+1), 256/(256+1), 1024/(1024+1), 1))
  
  Ratio <- c("0:100", "1:1024", "1:256", "1:64", "1:16", "1:4", "1:1", "4:1", "16:1", "64:1", "256:1", "1024:1", "100:0")
  Types <- c(rep("Serum", 13), rep("Urine", 13))
  
  df <- data.frame(Concentration = (c(Concentration_serum, Concentration_urine)),
                   Ratio = rep(Ratio, 2),
                   Samples = Types)
  df$Ratio <- factor(df$Ratio, levels = c("0:100","1:1024", "1:256", "1:64", "1:16", "1:4", "1:1", "4:1", "16:1", "64:1", "256:1", "1024:1", "100:0"))

  
  fill <- c("#5F9EA0", "#E1B378")
  fill <- c("#56B4E9", "#F0E442")
  fill <- c("#40b8d0", "#b2d183")
  fill <- c("#40b8d0", "#E1B378")
  
  p1 <- ggplot(data = df, aes(y = Concentration, x = Ratio, group = Samples)) + 
    geom_line(aes(colour = Samples)) + 
    geom_point(aes(colour = Samples)) +
    ggtitle("Dilution Series Design") + theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10)) #+ coord_flip()
  p1 <- p1 + scale_colour_manual(values=fill) +
    xlab("Ratio (Serum to Urine)") + ylab("Concentration (Relative)")
  p1
  
  
  setwd("/data/ms2_benchmark/whole_blood/")
  
  Cairo::Cairo(630, 500,file = paste0("figures/dilutions/dilution_dsesign_new.png"),dpi = 90,bg = "white")
  print(p1)
  dev.off()
  
  Cairo::CairoPDF(width = 6.3, height = 5,file = paste0("figures/dilutions/dilution_dsesign_new.pdf"))
  print(p1)
  dev.off()
  
}

## Dilution Stats - MS1 peak Number
{
  Concentration <- c(1/(1024+1), 1/(256+1), 1/(64+1), 1/(16+1), 1/(4+1), 1/(1+1), 4/(4+1), 16/(16+1), 64/(64+1), 256/(256+1), 1024/(1024+1))
  Concentrations <- vapply(Concentration, function(x){
    rep(x,3)
  }, FUN.VALUE = double(length = 3L))
  dim(Concentrations) <- c(1,33)
  Concentrations <- as.double(Concentrations)
  Ratio <- c("1:1024", "1:256", "1:64", "1:16", "1:4", "1:1", "4:1", "16:1", "64:1", "256:1", "1024:1")
  
  ## OptiLCMS 
  setwd("/data/ms2_benchmark/whole_blood/dilutions/1_optilcms")
  cors_list <- pval_list <- list()
  files <- c("C18pos_ms1/mSet1_c18pos.qs", "C18neg_ms1/mSet1_c18neg.qs",
             "HILICpos_ms1/mSet1_HILICpos.qs", "HILICneg_ms1/mSet1_HILICneg.qs")
  for(f in files){
    cat("==f ->", f, "\n")
    mSet <- qs::qread(f)
    df <- mSet@dataSet
    df <- df[-1,grepl("^S[A-Z]+", colnames(df))]
    
    dfx <- df[,grepl("^SL|^SM", colnames(df))]
    idxr <- apply(dfx, 1, function(x){
      ((length(which(x[1:3] == 0)) > 1) & (length(which(x[4:6] == 0)) <2)) |
        (length(which(x[4:6] == 0)) > 1 & (length(which(x[1:3] == 0)) <2))
    })
    
    df <- df[idxr,]
    
    cor_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["estimate"]][["cor"]]
    })
    pval_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["p.value"]]
    })
    cat("-->", length(which(abs(cor_vecs)>0.9))/length(cor_vecs)*100, "\n")
    cors_list <- c(cors_list, list(cor_vecs))
    pval_list <- c(pval_list, list(pval_vecs))
  }
  
  save(cors_list, pval_list, file = "results_stats/peaks_cors_pvals_optilcms.rda")
  
  ## MSDIAL 
  setwd("/data/ms2_benchmark/whole_blood/dilutions/2_MSDIAL/")
  cors_list <- pval_list <- list()
  files <- c("C18pos_ms1/Height_0_2023319415.txt", "C18neg_ms1/Height_0_2023319516.txt",
             "HILICpos_ms1/Height_0_20233191319.txt", "HILICneg_ms1/Height_0_20233191423.txt")
  for(f in files){
    cat("==f ->", f, "\n")
    dt <- read.csv(f, sep = "\t")
    dt <- dt[-c(1:3),]
    colnames(dt) <- dt[1,]
    dt <- dt[-1,]
    dt <- dt[,-c(1,4:32)]
    
    df <- dt[-1,grepl("^S[A-Z]+", colnames(dt))]
    
    dfx <- df[,grepl("^SL|^SM", colnames(df))]
    idxr <- apply(dfx, 1, function(x){
      ((length(which(x[1:3] == 0)) > 1) & (length(which(x[4:6] == 0)) <2)) |
        (length(which(x[4:6] == 0)) > 1 & (length(which(x[1:3] == 0)) <2))
    })
    
    df <- df[idxr,]
    
    cor_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["estimate"]][["cor"]]
    })
    pval_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["p.value"]]
    })
    cat("-->", length(which(abs(cor_vecs)>0.9))/length(cor_vecs)*100, "\n")
    cors_list <- c(cors_list, list(cor_vecs))
    pval_list <- c(pval_list, list(pval_vecs))
  }
  
  save(cors_list, pval_list, file = "results_stats/peaks_cors_pvals_msdial.rda")
  
  ## MZmine
  setwd("/data/ms2_benchmark/whole_blood/dilutions/3_mzMine/")
  cors_list <- pval_list <- list()
  files <- c("C18pos_ms1/ms1_dt.csv", "C18neg_ms1/ms1_dt.csv",
             "HILICpos_ms1/ms1_dt.csv", "HILICneg_ms1/ms1_dt.csv")
  for(f in files){
    cat("==f ->", f, "\n")
    dt <- read.csv(f, sep = ",")
    dtx <- dt[ ,grepl(".area", colnames(dt))]
    
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
    
    cor_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["estimate"]][["cor"]]
    })
    pval_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["p.value"]]
    })
    cat("-->", length(which(abs(cor_vecs)>0.9))/length(cor_vecs)*100, "\n")
    cors_list <- c(cors_list, list(cor_vecs))
    pval_list <- c(pval_list, list(pval_vecs))
  }
  
  save(cors_list, pval_list, file = "results_stats/peaks_cors_pvals_mzMine.rda")
  
  ## XCMS
  setwd("/data/ms2_benchmark/whole_blood/dilutions/4_xcms/")
  cors_list <- pval_list <- list()
  files <- c("C18pos_ms1/xcms_res.qs", "C18neg_ms1/xcms_res.qs",
             "HILICpos_ms1/xcms_res.qs", "HILICneg_ms1/xcms_res.qs")
  for(f in files){
    cat("==f ->", f, "\n")
    dt <- qs::qread(f)
    dtx <- dt@assays@data@listData[["raw"]]
    
    df <- dtx[ ,grepl("^S[A-Z]+", colnames(dtx))]
    df[is.na(df)] <- 0
    dfx <- df[,grepl("^SL|^SM", colnames(df))]
    idxr <- apply(dfx, 1, function(x){
      ((length(which(x[1:3] == 0)) > 1) & (length(which(x[4:6] == 0)) <2)) |
        (length(which(x[4:6] == 0)) > 1 & (length(which(x[1:3] == 0)) <2))
    })
    
    df <- df[idxr,]
    
    cor_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["estimate"]][["cor"]]
    })
    pval_vecs <- apply(df, 1, FUN = function(x){
      x_vec <- as.double(x[1:33])
      cor_res <- cor.test(x_vec, Concentrations)
      cor_res[["p.value"]]
    })
    cat("-->", length(which(abs(cor_vecs)>0.9))/length(cor_vecs)*100, "\n")
    cors_list <- c(cors_list, list(cor_vecs))
    pval_list <- c(pval_list, list(pval_vecs))
  }
  
  save(cors_list, pval_list, file = "results_stats/peaks_cors_pvals_xcms.rda")
  
}

## Plotting dot plots [cors] - MS1 
{
  rm(list = ls())
  setwd("/data/ms2_benchmark/whole_blood/dilutions/")
  load("1_optilcms/results_stats/peaks_cors_pvals_optilcms.rda")
  optilcms_res <- cors_list
  load("2_MSDIAL/results_stats/peaks_cors_pvals_msdial.rda")
  msdial_res <- cors_list
  load("3_mzMine/results_stats/peaks_cors_pvals_mzMine.rda")
  mzmine_res <- cors_list
  load("4_xcms/results_stats/peaks_cors_pvals_xcms.rda")
  xcms_res <- cors_list
  
  mode_vec <- c("C18pos", "C18neg", "HILICpos", "HILICneg")
  
  for(i in 1:4){
    df1 <- data.frame(cors = optilcms_res[[i]], Tools = "MetaboAnalyst")
    df2 <- data.frame(cors = msdial_res[[i]], Tools = "MSDIAL")
    df3 <- data.frame(cors = mzmine_res[[i]], Tools = "MzMine")
    df4 <- data.frame(cors = xcms_res[[i]], Tools = "XCMS")
    df <- rbind(df1, df2, df3, df4)
    df$cors <- abs(df$cors)
    df$Tools <- factor(df$Tools, levels = c("MetaboAnalyst", "XCMS", "MSDIAL", "MzMine"))
    p2 <- ggplot(df, aes(x=Tools, y = cors)) + 
      geom_violin(trim = T)+
      geom_dotplot(dotsize=0.2, binaxis = "y", stackdir = "center", binwidth = 0.001, stackratio=2)
    p3 <- p2 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                            geom="pointrange", color="blue") + ggtitle("Statistics of Dilution Series MS Features") + theme_ipsum() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10)) + ylim(c(0.3,1.1))
    
    setwd("/data/ms2_benchmark/whole_blood/")
    Cairo::Cairo(620, 500,file = paste0("figures/dilutions/dilution_ms1_stats_", mode_vec[i], ".png"),dpi = 90,bg = "white")
    print(p3)
    dev.off()
    
    Cairo::CairoPDF(width = 6.2, height = 5,file = paste0("figures/dilutions/dilution_ms1_stats_", mode_vec[i], ".pdf"))
    print(p3)
    dev.off()
  }
  
}

## Plotting dot plots [pvalues] - MS1 
{
  rm(list = ls())
  setwd("/data/ms2_benchmark/whole_blood/dilutions/")
  load("1_optilcms/results_stats/peaks_cors_pvals_optilcms.rda")
  optilcms_res <- pval_list
  load("2_MSDIAL/results_stats/peaks_cors_pvals_msdial.rda")
  msdial_res <- pval_list
  load("3_mzMine/results_stats/peaks_cors_pvals_mzMine.rda")
  mzmine_res <- pval_list
  load("4_xcms/results_stats/peaks_cors_pvals_xcms.rda")
  xcms_res <- pval_list
  
  mode_vec <- c("C18pos", "C18neg", "HILICpos", "HILICneg")
  
  for(i in 1:4){
    df1 <- data.frame(pvals = optilcms_res[[i]], Tools = "MetaboAnalyst")
    df2 <- data.frame(pvals = msdial_res[[i]], Tools = "MSDIAL")
    df3 <- data.frame(pvals = mzmine_res[[i]], Tools = "MzMine")
    df4 <- data.frame(pvals = xcms_res[[i]], Tools = "XCMS")
    df <- rbind(df1, df2, df3, df4)
    df$pvals <- abs(-log10(df$pvals))
    df$Tools <- factor(df$Tools, levels = c("MetaboAnalyst", "XCMS", "MSDIAL", "MzMine"))
    p2 <- ggplot(df, aes(x=Tools, y = pvals)) + 
      geom_violin(trim = T)+
      geom_dotplot(dotsize=0.2, binaxis = "y", stackdir = "center", binwidth = 0.075, stackratio=2.5)
    p3 <- p2 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                            geom="pointrange", color="blue") + ggtitle("Statistics of Dilution Series MS Features") + theme_ipsum() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10)) + ylab("-Log(p)")
    
    setwd("/data/ms2_benchmark/whole_blood/")
    Cairo::Cairo(620, 500,file = paste0("figures/dilutions/dilution_ms1_stats_pvals_", mode_vec[i], ".png"),dpi = 90,bg = "white")
    print(p3)
    dev.off()
    
    Cairo::CairoPDF(width = 6.2, height = 5,file = paste0("figures/dilutions/dilution_ms1_stats_pvals_", mode_vec[i], ".pdf"))
    print(p3)
    dev.off()
  }
  
}

## Plotting Stats of MS2 annotation results [DDA]
{
  rm(list=ls())
  setwd("/data/ms2_benchmark/whole_blood/dilutions/")
  
  num_cum <- function(f){
    dt <- qs::qread(f)
    dt$cors <- abs(dt$cors)
    cor_vec <- seq(1, 0, by=-0.02)
    res <- vapply(cor_vec, FUN = function(x){
      rs <- dt$identified[dt$cors >= x]
      length(which(rs))
    }, integer(1L))
    return(res)
  }
  
  # prepare compound number
  Optilcms_files_dda <- c("1_optilcms/C18pos_ms2/DDA/dda_res_optilcms_cors_df.qs",
                          "1_optilcms/C18neg_ms2/DDA/dda_res_optilcms_cors_df.qs",
                          "1_optilcms/HILICpos_ms2/DDA/dda_res_optilcms_cors_df.qs",
                          "1_optilcms/HILICneg_ms2/DDA/dda_res_optilcms_cors_df.qs")
  
  msdial_files_dda <- c("2_MSDIAL/C18pos_ms2/dda_res/dda_res_msdial_cors_df.qs",
                        "2_MSDIAL/C18neg_ms2/dda_res/dda_res_msdial_cors_df.qs",
                        "2_MSDIAL/HILICpos_ms2/dda_res/dda_res_msdial_cors_df.qs",
                        "2_MSDIAL/HILICneg_ms2/dda_res/dda_res_msdial_cors_df.qs")
  
  mzmine_files_dda <- c("3_mzMine/C18pos_ms2/dda_res_mzmine_sirius_cors_df.qs",
                        "3_mzMine/C18neg_ms2/dda_res_mzmine_sirius_cors_df.qs",
                        "3_mzMine/HILICpos_ms2/dda_res_mzmine_sirius_cors_df.qs",
                        "3_mzMine/HILICneg_ms2/dda_res_mzmine_sirius_cors_df.qs")
  
  
  optilcms_list <- lapply(Optilcms_files_dda, num_cum)
  msdial_list <- lapply(msdial_files_dda, num_cum)
  mzmine_list <- lapply(mzmine_files_dda, num_cum)
  
  optilcms_vec <- vapply(1:51, function(x){optilcms_list[[1]][x] + optilcms_list[[2]][x] + optilcms_list[[3]][x] + optilcms_list[[4]][x]}, integer(1L))
  msdial_vec <- vapply(1:51, function(x){msdial_list[[1]][x] + msdial_list[[2]][x] + msdial_list[[3]][x] + msdial_list[[4]][x]}, integer(1L))
  mzmine_vec <- vapply(1:51, function(x){mzmine_list[[1]][x] + mzmine_list[[2]][x] + mzmine_list[[3]][x] + mzmine_list[[4]][x]}, integer(1L))
  
  df <- data.frame(Compound_num = c(optilcms_vec, msdial_vec, mzmine_vec),
                   Tools = c(rep("MetaboAnalyst", 51), rep("MSDIAL/MSFinder", 51), rep("MzMine/SIRIUS", 51)),
                   Cors_r = rep(round(seq(1, 0, by=-0.02),2), 3))
  #df$Cors_r <- factor(df$Cors_r, levels = round(seq(1, 0, by=-0.02),2))
  df$Tools <- factor(df$Tools, levels = c("MetaboAnalyst", "MSDIAL/MSFinder", "MzMine/SIRIUS"))
  
  p5 <- ggplot(data = df, aes(y = Compound_num, x = Cors_r, group = Tools)) + 
    geom_line(aes(colour = Tools)) + 
    geom_point(aes(colour = Tools), size = 0.5) +
    ggtitle("Identified Compounds and Correlations (DDA)") + theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10)) #+ coord_flip()
  p5 <- p5 + xlab("Correlation Coefficient (r)") + ylab("Number of Identified Compounds") +
    scale_x_reverse(breaks = round(seq(1, 0, by = -0.05),1))
  
  setwd("/data/ms2_benchmark/whole_blood/")
  Cairo::Cairo(620, 500,file = paste0("figures/dilutions/dilution_ms2_dda_compare.png"),dpi = 90,bg = "white")
  print(p5)
  dev.off()
  
  Cairo::CairoPDF(width = 6.2, height = 5,file = paste0("figures/dilutions/dilution_ms2_dda_compare.pdf"))
  print(p5)
  dev.off()
  
}

## Plotting Stats of MS2 annotation results [DIA]
{
  rm(list=ls())
  setwd("/data/ms2_benchmark/whole_blood/dilutions/")
  
  num_cum <- function(f){
    dt <- qs::qread(f)
    dt$cors <- abs(dt$cors)
    cor_vec <- seq(1, 0, by=-0.02)
    res <- vapply(cor_vec, FUN = function(x){
      rs <- dt$identified[dt$cors >= x]
      length(which(rs))
    }, integer(1L))
    return(res)
  }
  
  # prepare compound number
  Optilcms_files_dia <- c("1_optilcms/C18pos_ms2/DIA1/dia_res_optilcms_cors_df.qs",
                          "1_optilcms/C18neg_ms2/DIA1/dia_res_optilcms_cors_df.qs",
                          "1_optilcms/HILICpos_ms2/DIA1/dia_res_optilcms_cors_df.qs",
                          "1_optilcms/HILICneg_ms2/DIA1/dia_res_optilcms_cors_df.qs")
  
  msdial_files_dia <- c("2_MSDIAL/C18pos_ms2/dia_res1/dia_res_msdial_cors_df.qs",
                        "2_MSDIAL/C18neg_ms2/dia_res1/dia_res_msdial_cors_df.qs",
                        "2_MSDIAL/HILICpos_ms2/dia_res1/dia_res_msdial_cors_df.qs",
                        "2_MSDIAL/HILICneg_ms2/dia_res1/dia_res_msdial_cors_df.qs")
  
  xcms_files_dia <- c("4_xcms/C18pos_ms2/DIA1/dia_res_mzmine_sirius_cors_df.qs",
                      "4_xcms/C18neg_ms2/DIA1/dia_res_mzmine_sirius_cors_df.qs",
                      "4_xcms/HILICpos_ms2/DIA1/dia_res_mzmine_sirius_cors_df.qs",
                      "4_xcms/HILICneg_ms2/DIA1/dia_res_mzmine_sirius_cors_df.qs")
  
  
  optilcms_list <- lapply(Optilcms_files_dia, num_cum)
  msdial_list <- lapply(msdial_files_dia, num_cum)
  xcms_list <- lapply(xcms_files_dia, num_cum)
  
  optilcms_vec <- vapply(1:51, function(x){optilcms_list[[1]][x] + optilcms_list[[2]][x] + optilcms_list[[3]][x] + optilcms_list[[4]][x]}, integer(1L))
  msdial_vec <- vapply(1:51, function(x){msdial_list[[1]][x] + msdial_list[[2]][x] + msdial_list[[3]][x] + msdial_list[[4]][x]}, integer(1L))
  xcms_vec <- vapply(1:51, function(x){xcms_list[[1]][x] + xcms_list[[2]][x] + xcms_list[[3]][x] + xcms_list[[4]][x]}, integer(1L))
  
  df <- data.frame(Compound_num = c(optilcms_vec, msdial_vec, xcms_vec),
                   Tools = c(rep("MetaboAnalyst", 51), rep("MSDIAL/MSFinder", 51), rep("XCMS/SIRIUS", 51)),
                   Cors_r = rep(round(seq(1, 0, by=-0.02),2), 3))
  #df$Cors_r <- factor(df$Cors_r, levels = round(seq(1, 0, by=-0.02),2))
  df$Tools <- factor(df$Tools, levels = c("MetaboAnalyst", "MSDIAL/MSFinder", "XCMS/SIRIUS"))
  
  p6 <- ggplot(data = df, aes(y = Compound_num, x = Cors_r, group = Tools)) + 
    geom_line(aes(colour = Tools)) + 
    geom_point(aes(colour = Tools), size = 0.5) +
    ggtitle("Identified Compounds and Correlations (SWATH-DIA)") + theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10)) #+ coord_flip()
  p6 <- p6 + xlab("Correlation Coefficient (r)") + ylab("Number of Identified Compounds") +
    scale_x_reverse(breaks = round(seq(1, 0, by = -0.05),1))
  
  setwd("/data/ms2_benchmark/whole_blood/")
  Cairo::Cairo(620, 500,file = paste0("figures/dilutions/dilution_ms2_DIA_compare.png"),dpi = 90,bg = "white")
  print(p6)
  dev.off()
  
  Cairo::CairoPDF(width = 6.2, height = 5,file = paste0("figures/dilutions/dilution_ms2_DIA_compare.pdf"))
  print(p6)
  dev.off()
  
}



