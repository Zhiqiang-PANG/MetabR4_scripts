# load("/data/ms2_benchmark/whole_blood/results/DDA/C18pos/Optilcms_resDT.rda")
# dt1 <- identifiedCMPD_dtx
# load("/data/ms2_benchmark/whole_blood/results/DDA/C18neg/Optilcms_resDT.rda")
# dt2 <- identifiedCMPD_dtx
# load("/data/ms2_benchmark/whole_blood/results/DDA/HILICneg/Optilcms_resDT.rda")
# dt3 <- identifiedCMPD_dtx
# load("/data/ms2_benchmark/whole_blood/results/DDA/HILICpos/Optilcms_resDT.rda")
# dt4 <- identifiedCMPD_dtx
# 
# dt_DDA <- rbind(dt1, dt2, dt3, dt4)
# 
# 
# load("/data/ms2_benchmark/whole_blood/results/DIA/C18pos/Optilcms_resDT.rda")
# dt1 <- identifiedCMPD_dtx
# load("/data/ms2_benchmark/whole_blood/results/DIA/C18neg/Optilcms_resDT.rda")
# dt2 <- identifiedCMPD_dtx
# load("/data/ms2_benchmark/whole_blood/results/DIA/HILICneg/Optilcms_resDT.rda")
# dt3 <- identifiedCMPD_dtx
# load("/data/ms2_benchmark/whole_blood/results/DIA/HILICpos/Optilcms_resDT.rda")
# dt4 <- identifiedCMPD_dtx
# 
# dt_DIA <- rbind(dt1, dt2, dt3, dt4)
# 
# dt <- rbind(dt_DIA, dt_DDA)
# write.table(unique(dt$Inchikey), file = "/data/ms2_benchmark/whole_blood/results/all_cmpd_inchikeys_optilcms.txt", quote = F, row.names = F, col.names = F)


res_dda <- qs::qread("/media/qiang/My Passport/Ultradata_backup/new_wb/C18pos/mzML/MS2/dda_res_optilcms_complete.qs")


