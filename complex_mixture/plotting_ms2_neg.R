load("/data/ms2_benchmark/MTBLS2207/IROA/optilcms_neg_dda/all_scores.rda")
t.test(all_score_deco, all_score_noDeco, paired = T, alternative = "greater")
# p-value = 1.106e-10

# Load ggplot2
library(ggplot2)
library(dplyr)
library(scales)
library(hrbrthemes)
library(viridis)
# Data
data <- data.frame(Algorithm = c(rep("Deco", 90),rep("noDeco", 90)), Score = c(all_score_deco, all_score_noDeco))
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
  scale_y_continuous(limits=c(60,70),oob = rescale_none) +
  geom_errorbar(aes(x=Algorithm, ymin=score-se, ymax=score+se), width=0.2, colour="orange", alpha=0.9, size=1) +
  ggtitle("Matching Scoring Comparison (DDA, ESI-)")+ 
  theme_ipsum()

Cairo::Cairo(file = "/data/ms2_benchmark/MTBLS2207/IROA/optilcms_neg_dda/score_bar_deco_nodeco.png",  unit="in", 
             dpi=72, width=6, height=6, type="png", bg="white");
print(p)
dev.off()

Cairo::Cairo(600, 700,file = paste0("/data/ms2_benchmark/MTBLS2207/IROA/Figures/bar_deco_nodeco/score_bar_deco_nodeco_neg.png"),
             dpi = 90,bg = "white")
print(p)
dev.off()

Cairo::CairoPDF(width = 6, height = 7,
                file = paste0("/data/ms2_benchmark/MTBLS2207/IROA/Figures/bar_deco_nodeco/score_bar_deco_nodeco_neg.pdf"))
print(p)
dev.off()


###### Now let's plot venn diagram -- DDA
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
load("msdial_neg_dda/matched_inchis_msdial_dda_neg.rda")
load("mzmine_sirius_neg_dda/matched_inchis_mzmine_neg.rda")
load("optilcms_neg_dda/matched_inchis_neg_dda.rda")
load("optilcms_neg_dda/matched_inchis_neg_noDeco_dda.rda")

dataList = list(
  MSDIAL = matched_inchis_msdial_dda_neg,
  MzMine = matched_inchis_mzmine_neg,
  MetaboAnalyst = matched_inchis_neg_dda,
  Metabo_noDeco = matched_inchis_neg_noDeco_dda
)
vennPlot(dataList, image = "Figures/venn_IROA_1/DDA_neg_inchikeys_4methods")

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
load("optilcms_neg_dda/matched_inchis_neg_dda.rda")
load("optilcms_neg_dda/matched_inchis_neg_BioDB_dda.rda")
dataList = list(
  All_DB = matched_inchis_neg_dda,
  Bio_DB = matched_inchis_neg_BioDB_dda
)
vennPlot2(dataList, image = "Figures/venn_IROA_1/DDA_neg_inchikeys_2DBs")

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
load("msdial_neg_dia/matched_inchis_msdial_dia_neg.rda")
load("xcms_sirius_dia/matched_inchis_xcms_neg.rda")
load("optilcms_neg_dia/matched_inchis_neg_dia.rda")

dataList = list(
  MSDIAL = matched_inchis_msdial_dia_neg,
  XCMS = matched_inchis_xcms_neg,
  MetaboAnalyst = matched_inchis_neg_dia
)
vennPlot3(dataList, image = "Figures/venn_IROA_1/DIA_neg_inchikeys_3methods")

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
load("optilcms_neg_dia/matched_inchis_neg_dia.rda")
load("optilcms_neg_dia/matched_inchis_neg_Biodb_dia.rda")
dataList = list(
  All_DB = matched_inchis_neg_dia,
  Bio_DB = matched_inchis_neg_Biodb_dia
)
vennPlot2(dataList, image = "Figures/venn_IROA_1/DIA_neg_inchikeys_2DBs")



