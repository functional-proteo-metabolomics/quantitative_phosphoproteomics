######################################
# The script plots Venn diagram of the analysed data.
# Input data is available either in /Dataset directory in RData format or in ProteomeXchange PXD030450 dataset.
# December 20, 2021, Yang Zhang

# clean all variables
rm(list=ls())

# the scripts start here
if (!require('eulerr')) install.packages('eulerr')
library(eulerr)
if (!require('wesanderson')) install.packages('wesanderson');
library(wesanderson)

load(file = "Dataset\\SiteTableV_cell.rdata")

SiteTableV_cell <- x 
a1 <-  plot(euler(list(LFQ = SiteTableV_cell$LFQ, siSILAC = SiteTableV_cell$siSILAC, TMT = SiteTableV_cell$TMT), shape = "ellipse"), quantities = TRUE, 
            fills = list(fill = wes_palette("GrandBudapest2", n = 3, type = "continuous"), alpha = 0.8),
            labels = list(col = "blue", cex=1), cat.fontfamily = "Arial", lty = 3, reverse = TRUE)
ggsave(file="VennD of Phosphosites_Cellline.svg", plot=a1, width=3.2, height=3.2)

load(file = "Dataset\\SiteTableV_tissue.rdata")
SiteTableV_tissue <- x 
a2 <-  plot(euler(list(LFQ = SiteTableV_tissue$LFQ, siSILAC = SiteTableV_tissue$siSILAC, TMT =SiteTableV_tissue$TMT), shape = "ellipse"), quantities = TRUE, 
            fills = list(fill = wes_palette("GrandBudapest2", n = 3, type = "continuous"), alpha = 0.8),
            labels = list(col = "blue", cex=1), cat.fontfamily = "Arial", lty = 3, reverse = TRUE)
ggsave(file="VennD of Phosphosites_Tissue.svg", plot=a2, width=3.2, height=3.2)
