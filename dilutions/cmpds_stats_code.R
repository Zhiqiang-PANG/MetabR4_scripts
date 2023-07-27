setwd("/data/ms2_benchmark/whole_blood/dilutions")
library(OptiLCMS2ID)
{
  res_dda_c18neg <- qs::qread("1_optilcms/C18neg_ms2/DDA/dda_res_optilcms_complete.qs")
  dt_c18neg1 <- OptiLCMS2ID::FormatMSnAnnotation(res_dda_c18neg, topN = 5L)
  res_dia_c18neg <- qs::qread("1_optilcms/C18neg_ms2/DIA1/dia_res_optilcms_complete.qs")
  res_dia_c18neg <- OptiLCMS2ID::PerformResultsExport(mSet = res_dia_c18neg, type = 0L, topN = 5L, ncores = 8L)
  dt_c18neg2 <- OptiLCMS2ID::FormatMSnAnnotation(res_dia_c18neg, topN = 5L)
  
  res_dia_c18neg2 <- qs::qread("1_optilcms/C18neg_ms2/DIA2/dia_res_optilcms_complete.qs")
  dt_c18neg3 <- OptiLCMS2ID::FormatMSnAnnotation(res_dia_c18neg2, topN = 5L)
  
  dt_c18neg <- rbind(dt_c18neg1, dt_c18neg2, dt_c18neg3)
  qs::qsave(dt_c18neg, file = "1_optilcms/C18neg_ms2/all_cmpd_dt_c18neg.qs")
  
}

{
  rm(list = ls())
  res_dda_c18pos <- qs::qread("1_optilcms/C18pos_ms2/DDA/dda_res_optilcms_complete.qs")
  dt_c18pos1 <- OptiLCMS2ID::FormatMSnAnnotation(res_dda_c18pos, topN = 5L)
  res_dia_c18pos <- qs::qread("1_optilcms/C18pos_ms2/DIA1/dia_res_optilcms_complete.qs")
  res_dia_c18pos <- OptiLCMS2ID::PerformResultsExport(mSet = res_dia_c18pos, type = 0L, topN = 5L, ncores = 8L)
  dt_c18pos2 <- OptiLCMS2ID::FormatMSnAnnotation(res_dia_c18pos, topN = 5L)
  
  res_dia_c18pos2 <- qs::qread("1_optilcms/C18pos_ms2/DIA2/dia_res_optilcms_complete.qs")
  res_dia_c18pos2 <- OptiLCMS2ID::PerformResultsExport(mSet = res_dia_c18pos2, type = 0L, topN = 5L, ncores = 8L)
  dt_c18pos3 <- OptiLCMS2ID::FormatMSnAnnotation(res_dia_c18pos2, topN = 5L)
  
  dt_c18pos <- rbind(dt_c18pos1, dt_c18pos2, dt_c18pos3)
  qs::qsave(dt_c18pos, file = "1_optilcms/C18pos_ms2/all_cmpd_dt_c18pos.qs")
  
}

{
  res_dda_HILICneg <- qs::qread("1_optilcms/HILICneg_ms2/DDA/dda_res_optilcms_complete.qs")
  dt_HILICneg1 <- OptiLCMS2ID::FormatMSnAnnotation(res_dda_HILICneg, topN = 5L)
  res_dia_HILICneg <- qs::qread("1_optilcms/HILICneg_ms2/DIA1/dia_res_optilcms_complete.qs")
  res_dia_HILICneg <- OptiLCMS2ID::PerformResultsExport(mSet = res_dia_HILICneg, type = 0L, topN = 5L, ncores = 8L)
  dt_HILICneg2 <- OptiLCMS2ID::FormatMSnAnnotation(res_dia_HILICneg, topN = 5L)
  
  res_dia_HILICneg2 <- qs::qread("1_optilcms/HILICneg_ms2/DIA2/dia_res_optilcms_complete.qs")
  dt_HILICneg3 <- OptiLCMS2ID::FormatMSnAnnotation(res_dia_HILICneg2, topN = 5L)
  
  dt_HILICneg <- rbind(dt_HILICneg1, dt_HILICneg2, dt_HILICneg3)
  qs::qsave(dt_HILICneg, file = "1_optilcms/HILICneg_ms2/all_cmpd_dt_HILICneg.qs")
  
}

{
  res_dda_HILICpos <- qs::qread("1_optilcms/HILICpos_ms2/DDA/dda_res_optilcms_complete.qs")
  dt_HILICpos1 <- OptiLCMS2ID::FormatMSnAnnotation(res_dda_HILICpos, topN = 5L)
  res_dia_HILICpos <- qs::qread("1_optilcms/HILICpos_ms2/DIA1/dia_res_optilcms_complete.qs")
  res_dia_HILICpos <- OptiLCMS2ID::PerformResultsExport(mSet = res_dia_HILICpos, type = 0L, topN = 5L, ncores = 8L)
  dt_HILICpos2 <- OptiLCMS2ID::FormatMSnAnnotation(res_dia_HILICpos, topN = 5L)
  
  res_dia_HILICpos2 <- qs::qread("1_optilcms/HILICpos_ms2/DIA2/dia_res_optilcms_complete.qs")
  dt_HILICpos3 <- OptiLCMS2ID::FormatMSnAnnotation(res_dia_HILICpos2, topN = 5L)
  
  dt_HILICpos <- rbind(dt_HILICpos1, dt_HILICpos2, dt_HILICpos3)
  qs::qsave(dt_HILICpos, file = "1_optilcms/HILICpos_ms2/all_cmpd_dt_HILICpos.qs")
  
}

dt_c18neg <- qs::qread("1_optilcms/C18neg_ms2/all_cmpd_dt_c18neg.qs")
dt_c18pos <- qs::qread("1_optilcms/C18pos_ms2/all_cmpd_dt_c18pos.qs")
dt_hilicneg <- qs::qread("1_optilcms/HILICneg_ms2/all_cmpd_dt_HILICneg.qs")
dt_hilicpos <- qs::qread("1_optilcms/HILICpos_ms2/all_cmpd_dt_HILICpos.qs")

dt_all <- rbind(dt_c18neg, dt_c18pos, dt_hilicneg, dt_hilicpos)
dt_all <- dt_all[(dt_all$Score_1 > 80) & (!grepl("Compound_", dt_all$Compound_1)), ]

