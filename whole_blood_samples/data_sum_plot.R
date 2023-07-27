library(ggplot2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(viridis)
library(hrbrthemes)

# SUM of all identical peaks
{
  setwd("/data/ms2_benchmark/whole_blood")
  rm(list = ls())
  {
    # optilcms - dda
    optilcms_ms1_c18pos_dda <- qs::qread("results/DDA/C18pos/unique_ms1_features_idx_optilcms.qs")
    optilcms_ms1_c18neg_dda <- qs::qread("results/DDA/C18neg/unique_ms1_features_idx_optilcms.qs")
    optilcms_ms1_HILICpos_dda <- qs::qread("results/DDA/HILICpos/unique_ms1_features_idx_optilcms.qs")
    optilcms_ms1_HILICneg_dda <- qs::qread("results/DDA/HILICneg/unique_ms1_features_idx_optilcms.qs")
    
    # msdial - dda
    msdial_ms1_c18pos_dda <- qs::qread("results/DDA/C18pos/unique_ms1_features_idx_msdial.qs")
    msdial_ms1_c18neg_dda <- qs::qread("results/DDA/C18neg/unique_ms1_features_idx_msdial.qs")
    msdial_ms1_HILICpos_dda <- qs::qread("results/DDA/HILICpos/unique_ms1_features_idx_msdial.qs")
    msdial_ms1_HILICneg_dda <- qs::qread("results/DDA/HILICneg/unique_ms1_features_idx_msdial.qs")
    
    # mzmine - dda
    mzmine_ms1_c18pos_dda <- qs::qread("results/DDA/C18pos/unique_ms1_features_idx_mzmine.qs")
    mzmine_ms1_c18neg_dda <- qs::qread("results/DDA/C18neg/unique_ms1_features_idx_mzmine.qs")
    mzmine_ms1_HILICpos_dda <- qs::qread("results/DDA/HILICpos/unique_ms1_features_idx_mzmine.qs")
    mzmine_ms1_HILICneg_dda <- qs::qread("results/DDA/HILICneg/unique_ms1_features_idx_mzmine.qs")
    
    Peak_num <- c(length(which(optilcms_ms1_c18pos_dda)), length(which(!optilcms_ms1_c18pos_dda)),
                  length(which(optilcms_ms1_c18neg_dda)), length(which(!optilcms_ms1_c18neg_dda)),
                  length(which(optilcms_ms1_HILICpos_dda)), length(which(!optilcms_ms1_HILICpos_dda)),
                  length(which(optilcms_ms1_HILICneg_dda)), length(which(!optilcms_ms1_HILICneg_dda)),
                  length(which(msdial_ms1_c18pos_dda)), length(which(!msdial_ms1_c18pos_dda)),
                  length(which(msdial_ms1_c18neg_dda)), length(which(!msdial_ms1_c18neg_dda)),
                  length(which(msdial_ms1_HILICpos_dda)), length(which(!msdial_ms1_HILICpos_dda)),
                  length(which(msdial_ms1_HILICneg_dda)), length(which(!msdial_ms1_HILICneg_dda)),
                  length(which(mzmine_ms1_c18pos_dda)), length(which(!mzmine_ms1_c18pos_dda)),
                  length(which(mzmine_ms1_c18neg_dda)), length(which(!mzmine_ms1_c18neg_dda)),
                  length(which(mzmine_ms1_HILICpos_dda)), length(which(!mzmine_ms1_HILICpos_dda)),
                  length(which(mzmine_ms1_HILICneg_dda)), length(which(!mzmine_ms1_HILICneg_dda)))
    ft_char <- rep(c('Unique Features','Generic Features'), 12)
    Tools <- c(rep("MetaboAnalyst", 8), rep("MSDIAL", 8),rep("MzMine", 8))
    Modes <- c(rep(c("C18pos", "C18pos", "C18neg", "C18neg", "HILICpos", "HILICpos", "HILICneg", "HILICneg"), 3))
    
    df <- data.frame(MS_Feature_Number = Peak_num,
                     Feature_characteristics = ft_char,
                     Tools = Tools,
                     Modes = Modes, stringsAsFactors = T)
    
    df$Tools <- factor(df$Tools,
                       levels= c("MetaboAnalyst", "MSDIAL", "MzMine"))
    
    df$Feature_characteristics <- factor(df$Feature_characteristics,
                                         levels= c("Generic Features", "Unique Features"))
    
    
    # fill <- c("#5F9EA0", "#E1B378")
    # fill <- c("#56B4E9", "#F0E442")
    fill <- c("#40b8d0", "#b2d183")
    
    p1 <- ggplot(data = df) + 
      geom_bar(aes(y = MS_Feature_Number, x = Modes, fill = Feature_characteristics), stat="identity") + 
      facet_grid(~Tools, scales = "free", space = "free") + 
      ggtitle("MetaboAnalyst vs. MSDIAL vs. MzMine") + #theme_ipsum() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8.5)) #+ coord_flip()
    p1 <- p1 + scale_fill_manual(values=fill) 
    p1
    
    setwd("/data/ms2_benchmark/whole_blood/")
    
    Cairo::Cairo(620, 500,file = paste0("figures/ms_feature_dda.png"),dpi = 90,bg = "white")
    print(p1)
    dev.off()
    
    Cairo::CairoPDF(width = 6.2, height = 5,file = paste0("figures/ms_feature_dda.pdf"))
    print(p1)
    dev.off()
    
  }

  #--------#
  rm(list = ls())
  {
    # optilcms - dia
    optilcms_ms1_c18pos_dia <- qs::qread("results/DIA/C18pos/unique_ms1_features_idx_optilcms.qs")
    optilcms_ms1_c18neg_dia <- qs::qread("results/DIA/C18neg/unique_ms1_features_idx_optilcms.qs")
    optilcms_ms1_HILICpos_dia <- qs::qread("results/DIA/HILICpos/unique_ms1_features_idx_optilcms.qs")
    optilcms_ms1_HILICneg_dia <- qs::qread("results/DIA/HILICneg/unique_ms1_features_idx_optilcms.qs")
    
    # msdial - dia
    msdial_ms1_c18pos_dia <- qs::qread("results/DIA/C18pos/unique_ms1_features_idx_msdial.qs")
    msdial_ms1_c18neg_dia <- qs::qread("results/DIA/C18neg/unique_ms1_features_idx_msdial.qs")
    msdial_ms1_HILICpos_dia <- qs::qread("results/DIA/HILICpos/unique_ms1_features_idx_msdial.qs")
    msdial_ms1_HILICneg_dia <- qs::qread("results/DIA/HILICneg/unique_ms1_features_idx_msdial.qs")
    
    # xcms - dia
    xcms_ms1_c18pos_dia <- qs::qread("results/DIA/C18pos/unique_ms1_features_idx_xcms.qs")
    xcms_ms1_c18neg_dia <- qs::qread("results/DIA/C18neg/unique_ms1_features_idx_xcms.qs")
    xcms_ms1_HILICpos_dia <- qs::qread("results/DIA/HILICpos/unique_ms1_features_idx_xcms.qs")
    xcms_ms1_HILICneg_dia <- qs::qread("results/DIA/HILICneg/unique_ms1_features_idx_xcms.qs")
    
    
    Peak_num <- c(length(which(optilcms_ms1_c18pos_dia)), length(which(!optilcms_ms1_c18pos_dia)),
                  length(which(optilcms_ms1_c18neg_dia)), length(which(!optilcms_ms1_c18neg_dia)),
                  length(which(optilcms_ms1_HILICpos_dia)), length(which(!optilcms_ms1_HILICpos_dia)),
                  length(which(optilcms_ms1_HILICneg_dia)), length(which(!optilcms_ms1_HILICneg_dia)),
                  length(which(msdial_ms1_c18pos_dia)), length(which(!msdial_ms1_c18pos_dia)),
                  length(which(msdial_ms1_c18neg_dia)), length(which(!msdial_ms1_c18neg_dia)),
                  length(which(msdial_ms1_HILICpos_dia)), length(which(!msdial_ms1_HILICpos_dia)),
                  length(which(msdial_ms1_HILICneg_dia)), length(which(!msdial_ms1_HILICneg_dia)),
                  length(which(xcms_ms1_c18pos_dia)), length(which(!xcms_ms1_c18pos_dia)),
                  length(which(xcms_ms1_c18neg_dia)), length(which(!xcms_ms1_c18neg_dia)),
                  length(which(xcms_ms1_HILICpos_dia)), length(which(!xcms_ms1_HILICpos_dia)),
                  length(which(xcms_ms1_HILICneg_dia)), length(which(!xcms_ms1_HILICneg_dia)))
    ft_char <- rep(c('Unique Features','Generic Features'), 12)
    Tools <- c(rep("MetaboAnalyst", 8), rep("MSDIAL", 8),rep("XCMS", 8))
    Modes <- c(rep(c("C18pos", "C18pos", "C18neg", "C18neg", "HILICpos", "HILICpos", "HILICneg", "HILICneg"), 3))
    
    df <- data.frame(MS_Feature_Number = Peak_num,
                     Feature_characteristics = ft_char,
                     Tools = Tools,
                     Modes = Modes, stringsAsFactors = T)
    
    df$Tools <- factor(df$Tools,
                       levels= c("MetaboAnalyst", "MSDIAL", "XCMS"))
    
    df$Feature_characteristics <- factor(df$Feature_characteristics,
                                         levels= c("Generic Features", "Unique Features"))
    
    
    fill <- c("#5F9EA0", "#E1B378")
    # fill <- c("#56B4E9", "#F0E442")
    #fill <- c("#40b8d0", "#b2d183")
    
    p2 <- ggplot(data = df) + 
      geom_bar(aes(y = MS_Feature_Number, x = Modes, fill = Feature_characteristics), stat="identity") + 
      facet_grid(~Tools, scales = "free", space = "free") + 
      ggtitle("MetaboAnalyst vs. MSDIAL vs. MzMine") + #theme_ipsum() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8.5)) #+ coord_flip()
    p2 <- p2 + scale_fill_manual(values=fill) 
    p2
    
    setwd("/data/ms2_benchmark/whole_blood/")
    
    Cairo::Cairo(620, 500,file = paste0("figures/ms_feature_dia.png"),dpi = 90,bg = "white")
    print(p2)
    dev.off()
    
    Cairo::CairoPDF(width = 6.2, height = 5,file = paste0("figures/ms_feature_dia.pdf"))
    print(p2)
    dev.off()
  }

  
  ## merge DDA and DIA
  rm(list = ls())
  {
    
    # optilcms - dda
    optilcms_ms1_c18pos_dda <- qs::qread("results/DDA/C18pos/unique_ms1_features_idx_optilcms.qs")
    optilcms_ms1_c18neg_dda <- qs::qread("results/DDA/C18neg/unique_ms1_features_idx_optilcms.qs")
    optilcms_ms1_HILICpos_dda <- qs::qread("results/DDA/HILICpos/unique_ms1_features_idx_optilcms.qs")
    optilcms_ms1_HILICneg_dda <- qs::qread("results/DDA/HILICneg/unique_ms1_features_idx_optilcms.qs")
    
    # msdial - dda
    msdial_ms1_c18pos_dda <- qs::qread("results/DDA/C18pos/unique_ms1_features_idx_msdial.qs")
    msdial_ms1_c18neg_dda <- qs::qread("results/DDA/C18neg/unique_ms1_features_idx_msdial.qs")
    msdial_ms1_HILICpos_dda <- qs::qread("results/DDA/HILICpos/unique_ms1_features_idx_msdial.qs")
    msdial_ms1_HILICneg_dda <- qs::qread("results/DDA/HILICneg/unique_ms1_features_idx_msdial.qs")
    
    # mzmine - dda
    mzmine_ms1_c18pos_dda <- qs::qread("results/DDA/C18pos/unique_ms1_features_idx_mzmine.qs")
    mzmine_ms1_c18neg_dda <- qs::qread("results/DDA/C18neg/unique_ms1_features_idx_mzmine.qs")
    mzmine_ms1_HILICpos_dda <- qs::qread("results/DDA/HILICpos/unique_ms1_features_idx_mzmine.qs")
    mzmine_ms1_HILICneg_dda <- qs::qread("results/DDA/HILICneg/unique_ms1_features_idx_mzmine.qs")
    
    # xcms - dia
    xcms_ms1_c18pos_dia <- qs::qread("results/DIA/C18pos/unique_ms1_features_idx_xcms.qs")
    xcms_ms1_c18neg_dia <- qs::qread("results/DIA/C18neg/unique_ms1_features_idx_xcms.qs")
    xcms_ms1_HILICpos_dia <- qs::qread("results/DIA/HILICpos/unique_ms1_features_idx_xcms.qs")
    xcms_ms1_HILICneg_dia <- qs::qread("results/DIA/HILICneg/unique_ms1_features_idx_xcms.qs")
    
    
    Peak_num <- c(length(which(optilcms_ms1_c18pos_dda)), length(which(!optilcms_ms1_c18pos_dda)),
                  length(which(optilcms_ms1_c18neg_dda)), length(which(!optilcms_ms1_c18neg_dda)),
                  length(which(optilcms_ms1_HILICpos_dda)), length(which(!optilcms_ms1_HILICpos_dda)),
                  length(which(optilcms_ms1_HILICneg_dda)), length(which(!optilcms_ms1_HILICneg_dda)),
                  length(which(msdial_ms1_c18pos_dda)), length(which(!msdial_ms1_c18pos_dda)),
                  length(which(msdial_ms1_c18neg_dda)), length(which(!msdial_ms1_c18neg_dda)),
                  length(which(msdial_ms1_HILICpos_dda)), length(which(!msdial_ms1_HILICpos_dda)),
                  length(which(msdial_ms1_HILICneg_dda)), length(which(!msdial_ms1_HILICneg_dda)),
                  length(which(mzmine_ms1_c18pos_dda)), length(which(!mzmine_ms1_c18pos_dda)),
                  length(which(mzmine_ms1_c18neg_dda)), length(which(!mzmine_ms1_c18neg_dda)),
                  length(which(mzmine_ms1_HILICpos_dda)), length(which(!mzmine_ms1_HILICpos_dda)),
                  length(which(mzmine_ms1_HILICneg_dda)), length(which(!mzmine_ms1_HILICneg_dda)),
                  length(which(xcms_ms1_c18pos_dia)), length(which(!xcms_ms1_c18pos_dia)),
                  length(which(xcms_ms1_c18neg_dia)), length(which(!xcms_ms1_c18neg_dia)),
                  length(which(xcms_ms1_HILICpos_dia)), length(which(!xcms_ms1_HILICpos_dia)),
                  length(which(xcms_ms1_HILICneg_dia)), length(which(!xcms_ms1_HILICneg_dia)))
    
    ft_char <- rep(c('Unique Features','Generic Features'), 16)
    Tools <- c(rep("MetaboAnalyst", 8), rep("MSDIAL", 8),rep("MZmine", 8), rep("XCMS", 8))
    Modes <- c(rep(c("C18 ESI+", "C18 ESI+", "C18 ESI-", "C18 ESI-", "HILIC ESI+", "HILIC ESI+", "HILIC ESI-", "HILIC ESI-"), 4))
    
    df <- data.frame(MS_Feature_Number = Peak_num,
                     Feature_characteristics = ft_char,
                     Tools = Tools,
                     Modes = Modes, stringsAsFactors = T)
    
    df$Tools <- factor(df$Tools,
                       levels= c("MetaboAnalyst", "MSDIAL", "MZmine", "XCMS"))
    
    df$Feature_characteristics <- factor(df$Feature_characteristics,
                                         levels= c("Generic Features", "Unique Features"))
    
    
    # fill <- c("#5F9EA0", "#E1B378")
    # fill <- c("#56B4E9", "#F0E442")
    fill <- c("#40b8d0", "#b2d183")
    fill <- c("#5F9EA0", "#E1B378")
    
    p1 <- ggplot(data = df) + 
      geom_bar(aes(y = MS_Feature_Number, x = Modes, fill = Feature_characteristics), stat="identity") + 
      facet_grid(~Tools, scales = "free", space = "free") + 
      ggtitle("MS Features stats in blood samples of different tools") + #theme_ipsum() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8.5)) #+ coord_flip()
    p1 <- p1 + scale_fill_manual(values=fill) 
    p1
    
    setwd("/data/ms2_benchmark/whole_blood/")
    
    Cairo::Cairo(620, 500,file = paste0("figures/ms_feature_stats.png"),dpi = 90,bg = "white")
    print(p1)
    dev.off()
    
    Cairo::CairoPDF(width = 6.2, height = 5,file = paste0("figures/ms_feature_stats.pdf"))
    print(p1)
    dev.off()
    
    
    
    
  }

  
}

# SUM of all compounds identified
{
  setwd("/data/ms2_benchmark/whole_blood")
  rm(list = ls())
  res_optilcms_dda <- c(
    "results/DDA/C18neg/Optilcms_resDT.rda",
    "results/DDA/C18pos/Optilcms_resDT.rda",
    "results/DDA/HILICneg/Optilcms_resDT.rda",
    "results/DDA/HILICpos/Optilcms_resDT.rda"
  )
  res_msdial_dda <- c(
    "results/DDA/C18neg/MSDIAL_resDT.rda",
    "results/DDA/C18pos/MSDIAL_resDT.rda",
    "results/DDA/HILICneg/MSDIAL_resDT.rda",
    "results/DDA/HILICpos/MSDIAL_resDT.rda"
  )
  res_sirius_dda <- c(
    "results/DDA/C18neg/Mzmine_sirius_resDT.rda",
    "results/DDA/C18pos/Mzmine_sirius_resDT.rda",
    "results/DDA/HILICneg/Mzmine_sirius_resDT.rda",
    "results/DDA/HILICpos/Mzmine_sirius_resDT.rda"
  )
  
  res_optilcms_dia <- c(
    "results/DIA/C18neg/Optilcms_resDT.rda",
    "results/DIA/C18pos/Optilcms_resDT.rda",
    "results/DIA/HILICneg/Optilcms_resDT.rda",
    "results/DIA/HILICpos/Optilcms_resDT.rda"
  )
  res_msdial_dia <- c(
    "results/DIA/C18neg/MSDIAL_resDT.rda",
    "results/DIA/C18pos/MSDIAL_resDT.rda",
    "results/DIA/HILICneg/MSDIAL_resDT.rda",
    "results/DIA/HILICpos/MSDIAL_resDT.rda"
  )
  res_sirius_dia <- c(
    "results/DIA/C18neg/XCMS_sirius_resDT.rda",
    "results/DIA/C18pos/XCMS_sirius_resDT.rda",
    "results/DIA/HILICneg/XCMS_sirius_resDT.rda",
    "results/DIA/HILICpos/XCMS_sirius_resDT.rda"
  )
  
  identified_cmpd_num <- integer()
  for(f in res_optilcms_dda){load(f); identified_cmpd_num <- c(identified_cmpd_num, nrow(identifiedCMPD_dtx))}
  for(f in res_msdial_dda){load(f); identified_cmpd_num <- c(identified_cmpd_num, nrow(resdt_confirmed))}
  for(f in res_sirius_dda){load(f); identified_cmpd_num <- c(identified_cmpd_num, nrow(dt_confirmed))}
  for(f in res_optilcms_dia){load(f); identified_cmpd_num <- c(identified_cmpd_num, nrow(identifiedCMPD_dtx))}
  for(f in res_msdial_dia){load(f); identified_cmpd_num <- c(identified_cmpd_num, nrow(resdt_confirmed))}
  for(f in res_sirius_dia){load(f); identified_cmpd_num <- c(identified_cmpd_num, nrow(dt_confirmed))}
  
  Tools <- c(rep("MetaboAnalyst", 4), rep("MSDIAL/MSFINDER", 4), rep("MzMine/SIRIUS", 4),
             rep("MetaboAnalyst", 4), rep("MSDIAL/MSFINDER", 4), rep("XCMS/SIRIUS", 4))
  
  Modes <- c(rep(c("C18_ESI-", "C18_ESI+", "HILIC_ESI-", "HILIC_ESI+"), 6))
  
  Acquisitions <- c(rep("DDA", 12), rep("SWATH-DIA", 12))
  
  df <- data.frame(Compound_number = identified_cmpd_num,
                   Tools = Tools,
                   Modes = Modes,
                   Acquisitions = Acquisitions)
  
  df$Tools <- factor(df$Tools,
                     levels= c("MetaboAnalyst", "MSDIAL/MSFINDER", "MzMine/SIRIUS", "XCMS/SIRIUS"))
  df$Acquisitions <- factor(df$Acquisitions,
                            levels= c("DDA", "SWATH-DIA"))
  
  #fill <- c("#5F9EA0", "#E1B378", "#40b8d0", "#b2d183")
  # fill <- c("#56B4E9", "#F0E442")
  # fill <- c("#40b8d0", "#b2d183")
  fill <- c("#bbbbbb", "#888888", "#666666", "#444444")
  p3 <- ggplot(data = df) + 
    geom_bar(aes(y = Compound_number, x = Tools, fill = Modes), stat="identity", width = 0.6) + 
    facet_grid(~Acquisitions, scales = "free", space = "free") + 
    ggtitle("MetaboAnalyst vs. MSDIAL vs. MzMine") + #theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8.5)) #+ coord_flip()
  p3 <- p3 + scale_fill_manual(values=fill) 
  p3
  
  setwd("/data/ms2_benchmark/whole_blood/")
  
  Cairo::Cairo(520, 500,file = paste0("figures/cmpd_stats.png"),dpi = 90,bg = "white")
  print(p3)
  dev.off()
  
  Cairo::CairoPDF(width = 5.2, height = 5,file = paste0("figures/cmpd_stats.pdf"))
  print(p3)
  dev.off()
  
}

# SUM of all inchikeys - venn diagram
{
  setwd("/data/ms2_benchmark/whole_blood")
  rm(list = ls())
  res_optilcms_dda <- c(
    "results/DDA/C18neg/Optilcms_resDT.rda",
    "results/DDA/C18pos/Optilcms_resDT.rda",
    "results/DDA/HILICneg/Optilcms_resDT.rda",
    "results/DDA/HILICpos/Optilcms_resDT.rda"
  )
  res_msdial_dda <- c(
    "results/DDA/C18neg/MSDIAL_resDT.rda",
    "results/DDA/C18pos/MSDIAL_resDT.rda",
    "results/DDA/HILICneg/MSDIAL_resDT.rda",
    "results/DDA/HILICpos/MSDIAL_resDT.rda"
  )
  res_sirius_dda <- c(
    "results/DDA/C18neg/Mzmine_sirius_resDT.rda",
    "results/DDA/C18pos/Mzmine_sirius_resDT.rda",
    "results/DDA/HILICneg/Mzmine_sirius_resDT.rda",
    "results/DDA/HILICpos/Mzmine_sirius_resDT.rda"
  )
  
  res_optilcms_dia <- c(
    "results/DIA/C18neg/Optilcms_resDT.rda",
    "results/DIA/C18pos/Optilcms_resDT.rda",
    "results/DIA/HILICneg/Optilcms_resDT.rda",
    "results/DIA/HILICpos/Optilcms_resDT.rda"
  )
  res_msdial_dia <- c(
    "results/DIA/C18neg/MSDIAL_resDT.rda",
    "results/DIA/C18pos/MSDIAL_resDT.rda",
    "results/DIA/HILICneg/MSDIAL_resDT.rda",
    "results/DIA/HILICpos/MSDIAL_resDT.rda"
  )
  res_sirius_dia <- c(
    "results/DIA/C18neg/XCMS_sirius_resDT.rda",
    "results/DIA/C18pos/XCMS_sirius_resDT.rda",
    "results/DIA/HILICneg/XCMS_sirius_resDT.rda",
    "results/DIA/HILICpos/XCMS_sirius_resDT.rda"
  )
  
  identified_cmpd_num <- list()
  for(f in res_optilcms_dda){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(identifiedCMPD_dtx$Inchikey))}
  for(f in res_msdial_dda){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(resdt_confirmed$InChIKey))}
  for(f in res_sirius_dda){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(dt_confirmed$InChIkey2D))}
  
  for(f in res_optilcms_dia){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(identifiedCMPD_dtx$Inchikey))}
  for(f in res_msdial_dia){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(resdt_confirmed$InChIKey))}
  for(f in res_sirius_dia){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(dt_confirmed$InChIkey2D))}
  
  identified_cmpd_all <- lapply(identified_cmpd_num, function(x){
    vapply(x, function(i){strsplit(i,"-")[[1]][1]}, FUN.VALUE = character(length = 1L), USE.NAMES = F)
  })
  
  Optilcms_dda_inchikeys <- unique(unlist(identified_cmpd_all[1:4]))
  msdial_dda_inchikeys <- unique(unlist(identified_cmpd_all[5:8]))
  sirius_dda_inchikeys <- unique(unlist(identified_cmpd_all[9:12]))
  
  Optilcms_dia_inchikeys <- unique(unlist(identified_cmpd_all[13:16]))
  msdial_dia_inchikeys <- unique(unlist(identified_cmpd_all[17:20]))
  sirius_dia_inchikeys <- unique(unlist(identified_cmpd_all[21:24]))
  
  ### Venn plots of different tools
  require("VennDiagram");
  vennPlot3 <- function(dataList, image){
    venn.plot <- venn.diagram(
      dataList,
      filename = NULL,
      euler.d = FALSE,
      scaled = FALSE,
      col = "transparent",
      fill = c("red", "blue", "green"),
      alpha = 0.35,
      label.col = c("darkred", "white", "darkblue", "white",
                             "white", "white", "darkgreen"),
                             cex = 1.5,
      fontfamily = "serif",
      fontface = "bold",
      cat.col = c("darkred", "darkblue", "darkgreen"),
      cat.cex = 1.5,
      cat.pos = c(-27, 27, 175),
      cat.dist = c(0.055, 0.055, 0.055),
      cat.fontfamily = "serif",
      rotation.degree = 0,
      margin = 0.2
    );
    #grid.draw(venn.plot)
    #dev.off()
    Cairo::Cairo(640, 660,file = paste0(image,".png"),dpi = 90,bg = "white")
    grid.draw(venn.plot)
    dev.off()
    
    Cairo::CairoPDF(width = 6.2, height = 6.3,file = paste0(image,".pdf"))
    grid.draw(venn.plot)
    dev.off()
  }
  
  # DDA DATASET
  dataList = list(
    MSDIAL_MSFINDER = msdial_dda_inchikeys[!is.na(msdial_dda_inchikeys)],
    MzMine_SIRIUS = sirius_dda_inchikeys[!is.na(sirius_dda_inchikeys)],
    MetaboAnalyst = Optilcms_dda_inchikeys[!is.na(Optilcms_dda_inchikeys)]
  )
  vennPlot3(dataList, image = "figures/DDA_inchikeys_3methods")
  
  # DIA DATASET
  dataList = list(
    MSDIAL_MSFINDER = msdial_dia_inchikeys[!is.na(msdial_dia_inchikeys)],
    MzMine_SIRIUS = sirius_dia_inchikeys[!is.na(sirius_dia_inchikeys)],
    MetaboAnalyst = Optilcms_dia_inchikeys[!is.na(Optilcms_dia_inchikeys)]
  )
  vennPlot3(dataList, image = "figures/DIA_inchikeys_3methods")
  
}

# HDMB database matching results - ratio
{
  setwd("/data/ms2_benchmark/whole_blood")
  rm(list = ls())
  hdmbdb <- qs::qread("/ultraData/HMDB/Dec_12/hmdb_metabolites_named.qs")
  
  blood_inchikeys <- vapply(1:217920, function(x){ #920
    #cat("x => ", x, "\n")
    if(!is.null(hdmbdb[[x]][["biological_properties"]][["biospecimen_locations"]])){
      vec <- unlist(hdmbdb[[x]][["biological_properties"]][["biospecimen_locations"]])
      if("Blood" %in% vec){
        if(is.null(hdmbdb[[x]][["inchikey"]])){
          ""
        } else {
          hdmbdb[[x]][["inchikey"]]
        }
      } else {
        ""
      }
    } else {
      ""
    }
  }, FUN.VALUE = character(1L))
  blood_inchikeys <- blood_inchikeys[blood_inchikeys !=""];
  
  all_inchikeys <- vapply(1:217920, function(x){ 
    if(is.null(hdmbdb[[x]][["inchikey"]])){
      ""
      } else {
        hdmbdb[[x]][["inchikey"]]
      }
  }, FUN.VALUE = character(1L))
  all_inchikeys <- all_inchikeys[all_inchikeys != ""]

  ##
  res_optilcms_dda <- c(
    "results/DDA/C18neg/Optilcms_resDT.rda",
    "results/DDA/C18pos/Optilcms_resDT.rda",
    "results/DDA/HILICneg/Optilcms_resDT.rda",
    "results/DDA/HILICpos/Optilcms_resDT.rda"
  )
  res_msdial_dda <- c(
    "results/DDA/C18neg/MSDIAL_resDT.rda",
    "results/DDA/C18pos/MSDIAL_resDT.rda",
    "results/DDA/HILICneg/MSDIAL_resDT.rda",
    "results/DDA/HILICpos/MSDIAL_resDT.rda"
  )
  res_sirius_dda <- c(
    "results/DDA/C18neg/Mzmine_sirius_resDT.rda",
    "results/DDA/C18pos/Mzmine_sirius_resDT.rda",
    "results/DDA/HILICneg/Mzmine_sirius_resDT.rda",
    "results/DDA/HILICpos/Mzmine_sirius_resDT.rda"
  )
  
  res_optilcms_dia <- c(
    "results/DIA/C18neg/Optilcms_resDT.rda",
    "results/DIA/C18pos/Optilcms_resDT.rda",
    "results/DIA/HILICneg/Optilcms_resDT.rda",
    "results/DIA/HILICpos/Optilcms_resDT.rda"
  )
  res_msdial_dia <- c(
    "results/DIA/C18neg/MSDIAL_resDT.rda",
    "results/DIA/C18pos/MSDIAL_resDT.rda",
    "results/DIA/HILICneg/MSDIAL_resDT.rda",
    "results/DIA/HILICpos/MSDIAL_resDT.rda"
  )
  res_sirius_dia <- c(
    "results/DIA/C18neg/XCMS_sirius_resDT.rda",
    "results/DIA/C18pos/XCMS_sirius_resDT.rda",
    "results/DIA/HILICneg/XCMS_sirius_resDT.rda",
    "results/DIA/HILICpos/XCMS_sirius_resDT.rda"
  )
  
  identified_cmpd_num <- list()
  for(f in res_optilcms_dda){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(identifiedCMPD_dtx$Inchikey))}
  for(f in res_msdial_dda){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(resdt_confirmed$InChIKey))}
  for(f in res_sirius_dda){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(dt_confirmed$InChIkey2D))}
  
  for(f in res_optilcms_dia){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(identifiedCMPD_dtx$Inchikey))}
  for(f in res_msdial_dia){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(resdt_confirmed$InChIKey))}
  for(f in res_sirius_dia){load(f); identified_cmpd_num <- c(identified_cmpd_num, list(dt_confirmed$InChIkey2D))}
  
  identified_cmpd_all <- lapply(identified_cmpd_num, function(x){
    vapply(x, function(i){strsplit(i,"-")[[1]][1]}, FUN.VALUE = character(length = 1L), USE.NAMES = F)
  })
  
  Optilcms_dda_inchikeys <- unique(unlist(identified_cmpd_all[1:4]))
  msdial_dda_inchikeys <- unique(unlist(identified_cmpd_all[5:8]))
  sirius_dda_inchikeys <- unique(unlist(identified_cmpd_all[9:12]))
  
  Optilcms_dia_inchikeys <- unique(unlist(identified_cmpd_all[13:16]))
  msdial_dia_inchikeys <- unique(unlist(identified_cmpd_all[17:20]))
  sirius_dia_inchikeys <- unique(unlist(identified_cmpd_all[21:24]))
  
  blood_inchikeys_x <- vapply(blood_inchikeys, function(i){strsplit(i,"-")[[1]][1]}, FUN.VALUE = character(length = 1L), USE.NAMES = F)
  all_inchikeys_x <- vapply(all_inchikeys, function(i){strsplit(i,"-")[[1]][1]}, FUN.VALUE = character(length = 1L), USE.NAMES = F)
  
  ratio_inHMDB <- ratio_inHMDB_blood <- double(6L)
  ## DDA
  # - OptiLCMS
  Optilcms_dda_inDB <- vapply(Optilcms_dda_inchikeys, function(i){
    i %in% blood_inchikeys_x
  }, logical(1L))
  ratio_inHMDB_blood[1] <- length(which(Optilcms_dda_inDB))/length(Optilcms_dda_inDB)
  Optilcms_dda_inDB_all <- vapply(Optilcms_dda_inchikeys, function(i){
    i %in% all_inchikeys_x
  }, logical(1L))
  ratio_inHMDB[1] <- length(which(Optilcms_dda_inDB_all))/length(Optilcms_dda_inDB_all)
  # - MSDIAL
  msdial_dda_inDB <- vapply(msdial_dda_inchikeys, function(i){
    i %in% blood_inchikeys_x
  }, logical(1L))
  ratio_inHMDB_blood[2] <- length(which(msdial_dda_inDB))/length(msdial_dda_inDB)
  msdial_dda_inDB_all <- vapply(msdial_dda_inchikeys, function(i){
    i %in% all_inchikeys_x
  }, logical(1L))
  ratio_inHMDB[2] <- length(which(msdial_dda_inDB_all))/length(msdial_dda_inDB_all)
  # - SIRIUS
  sirius_dda_inDB <- vapply(sirius_dda_inchikeys, function(i){
    i %in% blood_inchikeys_x
  }, logical(1L))
  ratio_inHMDB_blood[3] <- length(which(sirius_dda_inDB))/length(sirius_dda_inDB)
  sirius_dda_inDB_all <- vapply(sirius_dda_inchikeys, function(i){
    i %in% all_inchikeys_x
  }, logical(1L))
  ratio_inHMDB[3] <- length(which(sirius_dda_inDB_all))/length(sirius_dda_inDB_all)
  
  # DIA
  # - OptiLCMS
  Optilcms_dia_inDB <- vapply(Optilcms_dia_inchikeys, function(i){
    i %in% blood_inchikeys_x
  }, logical(1L))
  ratio_inHMDB_blood[4] <- length(which(Optilcms_dia_inDB))/length(Optilcms_dia_inchikeys)
  Optilcms_dia_inDB_all <- vapply(Optilcms_dia_inchikeys, function(i){
    i %in% all_inchikeys_x
  }, logical(1L))
  ratio_inHMDB[4] <- length(which(Optilcms_dia_inDB_all))/length(Optilcms_dia_inchikeys)
  
  # - MSDIAL
  msdial_dia_inDB <- vapply(msdial_dia_inchikeys, function(i){
    i %in% blood_inchikeys_x
  }, logical(1L))
  ratio_inHMDB_blood[5] <- length(which(msdial_dia_inDB))/length(msdial_dia_inchikeys)
  msdial_dia_inDB_all <- vapply(msdial_dia_inchikeys, function(i){
    i %in% all_inchikeys_x
  }, logical(1L))
  ratio_inHMDB[5] <- length(which(msdial_dia_inDB_all))/length(msdial_dia_inchikeys)
  
  # - SIRIUS
  sirius_dia_inDB <- vapply(sirius_dia_inchikeys, function(i){
    i %in% blood_inchikeys_x
  }, logical(1L))
  ratio_inHMDB_blood[6] <-  length(which(sirius_dia_inDB))/length(sirius_dia_inchikeys)
  sirius_dia_inDB_all <- vapply(sirius_dia_inchikeys, function(i){
    i %in% all_inchikeys_x
  }, logical(1L))
  ratio_inHMDB[6] <- length(which(sirius_dia_inDB_all))/length(sirius_dia_inchikeys)
  
  
  #### ratio values
  ratio_values <- c(ratio_inHMDB_blood, ratio_inHMDB)
  ratio_values <- ratio_values*100
  databases <- c(rep("HMDB_blood", 6), rep("HMDB_all", 6))
  Tools <- rep(c("MetaboAnalyst", "MSFINDER", "SIRIUS"), 4)
  Modes <- c(rep("DDA", 3), rep("SWATH-DIA", 3), rep("DDA", 3), rep("SWATH-DIA", 3))
  df <- data.frame(Percentage = ratio_values,
                   Databases = databases,
                   Tools = Tools,
                   Modes = Modes)
  
  df$Tools <- factor(df$Tools,
                     levels= c("MetaboAnalyst", "MSFINDER", "SIRIUS"))
  df$Modes <- factor(df$Modes,
                     levels= c("DDA", "SWATH-DIA"))
  
  save(df, file = "hmdb_results/hmdb_matching_ratio.rda")
  #fill <- c("#5F9EA0", "#E1B378", "#40b8d0", "#b2d183")
  # fill <- c("#56B4E9", "#F0E442")
  # fill <- c("#40b8d0", "#b2d183")
  fill <- c("#aaaaaa", "#666666", "#444444")
  p4 <- ggplot(data = df) + 
    geom_bar(aes(y = Percentage, x = Databases, fill = Tools), stat="identity", width = 0.6, position = "dodge") + 
    ggtitle("MetaboAnalyst vs. MSFINDER vs. SIRIUS") + #theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8.5)) #+ coord_flip()
  p4 <- p4 + scale_fill_manual(values=fill) 
  p4
  
  setwd("/data/ms2_benchmark/whole_blood/")
  
  Cairo::Cairo(520, 500,file = paste0("figures/hmdb_matching_ratio.png"),dpi = 90,bg = "white")
  print(p4)
  dev.off()
  
  Cairo::CairoPDF(width = 5.2, height = 5,file = paste0("figures/hmdb_matching_ratio.pdf"))
  print(p4)
  dev.off()
}