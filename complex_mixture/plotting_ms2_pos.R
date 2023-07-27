load("/data/ms2_benchmark/MTBLS2207/IROA/optilcms_pos_dda/all_score_raw_new.rda")
t.test(all_score_deco, all_score_noDeco, paired = T, alternative = "greater")
# p-value = 1.106e-10

# Load ggplot2
library(ggplot2)
library(dplyr)
library(scales)
library(hrbrthemes)
library(viridis)
# Data
data <- data.frame(Algorithm = c(rep("Deco", 99),rep("noDeco", 99)), Score = c(all_score_deco, all_score_noDeco))
# Calculates mean, sd, se and IC
my_sum <- data %>%
  group_by(Algorithm) %>%
  summarise( 
    n=n(),
    score=mean(Score),
    sd=sd(Score)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))


# Standard Error
p <- ggplot(my_sum) +
  geom_bar( aes(x=Algorithm, y=score), stat="identity", fill="forestgreen",width = 0.5, alpha=0.5) +
  scale_y_continuous(limits=c(65,75),oob = rescale_none) +
  geom_errorbar(aes(x=Algorithm, ymin=score-se, ymax=score+se), width=0.2, colour="orange", alpha=0.9, size=1) +
  ggtitle("Matching Scoring Comparison (DDA, ESI+)")+ 
  theme_ipsum()

Cairo::Cairo(file = "/data/ms2_benchmark/MTBLS2207/IROA/optilcms_pos_dda/all_score_raw_new.png",  unit="in", dpi=72, width=6, height=6, type="png", bg="white");
print(p)
dev.off()


Cairo::Cairo(600, 700,file = paste0("/data/ms2_benchmark/MTBLS2207/IROA/Figures/bar_deco_nodeco/score_bar_deco_nodeco_pos.png"),
             dpi = 90,bg = "white")
print(p)
dev.off()

Cairo::CairoPDF(width = 6, height = 7,
                file = paste0("/data/ms2_benchmark/MTBLS2207/IROA/Figures/bar_deco_nodeco/score_bar_deco_nodeco_pos.pdf"))
print(p)
dev.off()



###### Now let's plot venn diagram
require("VennDiagram");
rm(list = ls())
setwd("/data/ms2_benchmark/MTBLS2207/IROA/")
vennPlot <- function(dataList, image){
  venn.plot <- venn.diagram(
    dataList,
    filename = NULL,
    euler.d = FALSE,
    scaled = FALSE,
    col = "transparent",
    #fill = brewer.pal(12, "Set3")[c(5,7,12,10)],
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    alpha = 0.50,
    label.col = c("orange", "white", "darkorchid4", "white", 
                  "white", "white", "white", "white", "darkblue", "white", 
                  "white", "white", "white", "darkgreen", "white"),
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.07,
    cat.fontfamily = "serif",
    rotation.degree = 0,
    margin = 0.2
  );
  
  #dev.off()
  Cairo::Cairo(640, 665,file = paste0(image,".png"),dpi = 90,bg = "white")
  grid.draw(venn.plot)
  dev.off()
  
  Cairo::CairoPDF(width = 6, height = 6.5,file = paste0(image,".pdf"))
  grid.draw(venn.plot)
  dev.off()
}
load("msdial_pos_dda/matched_inchis_msdial_dda_pos.rda")
load("mzmine_sirius_pos_dda/matched_inchis_mzmine_pos.rda")
load("optilcms_pos_dda/matched_inchis_pos_dda.rda")
load("optilcms_pos_dda/matched_inchis_pos_noDeco_dda.rda")

dataList = list(
  MSDIAL = matched_inchis_msdial_dda_pos,
  MzMine = matched_inchis_mzmine_pos,
  MetaboAnalyst = matched_inchis_pos_dda,
  Metabo_noDeco = matched_inchis_noDeco_dda
)
vennPlot(dataList, image = "Figures/venn_IROA_1/DDA_pos_inchikeys_4methods")

require("VennDiagram");
rm(list = ls())
setwd("/data/ms2_benchmark/MTBLS2207/IROA/")
vennPlot2 <- function(dataList, image){
  
  venn.plot <- venn.diagram(
    dataList,
    filename = NULL,
    euler.d = FALSE,
    scaled = FALSE,
    col = "transparent",
    #fill = brewer.pal(12, "Set3")[c(5,7,12,10)],
    fill = metan::ggplot_color(length(dataList)),
    alpha = 0.50,
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.07,
    cat.fontfamily = "serif",
    rotation.degree = 0,
    margin = 0.2
  );
  
  #dev.off()
  Cairo::Cairo(640, 660,file = paste0(image,".png"),dpi = 90,bg = "white")
  grid.draw(venn.plot)
  dev.off()
  
  Cairo::CairoPDF(width = 6.4, height = 6.5,file = paste0(image,".pdf"))
  grid.draw(venn.plot)
  dev.off()
}
load("optilcms_pos_dda/matched_inchis_pos_dda.rda")
load("optilcms_pos_dda/matched_inchis_bioDB_dda.rda")
dataList = list(
  All_DB = matched_inchis_pos_dda,
  Bio_DB = matched_inchis_bioDB_dda
)
vennPlot2(dataList, image = "Figures/venn_IROA_1/DDA_pos_inchikeys_2DBs")


###### venn diagram -- DIA
require("VennDiagram");
rm(list = ls())
setwd("/data/ms2_benchmark/MTBLS2207/IROA/")
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
load("msdial_pos_dia/matched_inchis_msdial_dia_pos.rda")
load("xcms_sirius_dia/matched_inchis_xcms_pos.rda")
load("optilcms_pos_dia/matched_inchis_pos_dia.rda")

dataList = list(
  MSDIAL = matched_inchis_msdial_dia_pos,
  XCMS = matched_inchis_xcms_pos,
  MetaboAnalyst = matched_inchis_pos_dia
)
vennPlot3(dataList, image = "Figures/venn_IROA_1/DIA_pos_inchikeys_3methods")

require("VennDiagram");
rm(list = ls())
setwd("/data/ms2_benchmark/MTBLS2207/IROA/")
vennPlot2 <- function(dataList, image){
  
  venn.plot <- venn.diagram(
    dataList,
    filename = NULL,
    euler.d = FALSE,
    scaled = FALSE,
    col = "transparent",
    #fill = brewer.pal(12, "Set3")[c(5,7,12,10)],
    fill = metan::ggplot_color(length(dataList)),
    alpha = 0.50,
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.07,
    cat.fontfamily = "serif",
    rotation.degree = 0,
    margin = 0.2
  );
  
  #dev.off()
  Cairo::Cairo(640, 660,file = paste0(image,".png"),dpi = 90,bg = "white")
  grid.draw(venn.plot)
  dev.off()
  
  Cairo::CairoPDF(width = 6.4, height = 6.5,file = paste0(image,".pdf"))
  grid.draw(venn.plot)
  dev.off()
}
load("optilcms_pos_dia/matched_inchis_pos_dia.rda")
load("optilcms_pos_dia/matched_inchis_biodb_pos_dia.rda")
dataList = list(
  All_DB = matched_inchis_pos_dia,
  Bio_DB = matched_inchis_BioDB_pos_dia
)
vennPlot2(dataList, image = "Figures/venn_IROA_1/DIA_pos_inchikeys_2DBs")






