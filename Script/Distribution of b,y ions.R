######################################
# The script calculates the distribution of b and y ion for the identifed phosphopeptides and peptides.
# Input data is available either in /Dataset directory in RData format or in ProteomeXchange PXD030450 dataset.
# December 20, 2021, Yang Zhang

# clean all variables
rm(list=ls())

# library
if (!require('stringr')) install.packages('stringr')
library(stringr)
if (!require('reshape2')) install.packages('reshape2')
library(reshape2)
if (!require('ggplot2')) install.packages('ggplot2')
library(ggplot2)

# bunion
countb <- function(x) {
  str_count(as.character(unname(unlist(x))), "b")
}

county <- function(x) {
  str_count(as.character(unname(unlist(x))), "y")
}

countAll <- function(x) {
  str_count(as.character(unname(unlist(x))), ";") + 1
}

# load the Rdata
load(file = "Dataset\\MSMS.Rdata")

# phosphopeptides
# Cell line
# SILAC
PhoMSMS_C_siSILAC <- PhoMSMS_C_siSILAC[which(PhoMSMS_C_siSILAC$Raw.file == "190920_ID_59_YZ_cellline_SILAC_C1_R1_enriched"),]
PhoMSMS_C_siSILAC <- PhoMSMS_C_siSILAC[, c("Sequence", "Matches")]

# LFQ
PhoMSMS_C_LFQ <- PhoMSMS_C_LFQ[which(PhoMSMS_C_LFQ$Raw.file == "190918_ID_23_YZ_cellline_LFQ_C1_R1_enriched"),]
PhoMSMS_C_LFQ <- PhoMSMS_C_LFQ[, c("Sequence", "Matches")]

# TMT
PhoMSMS_C_TMT <- PhoMSMS_C_TMT[which(PhoMSMS_C_TMT$Raw.file == "190417_ID_3_YZ_cellline_TMT_pair1_enriched"),]
PhoMSMS_C_TMT <- PhoMSMS_C_TMT[, c("Sequence", "Matches")]

# calculate the percentage of b and y ions
# TMT 
PhoMSMS_C_TMT$bn <- countb(PhoMSMS_C_TMT$Matches)
PhoMSMS_C_TMT$yn <- county(PhoMSMS_C_TMT$Matches)
PhoMSMS_C_TMT$All <- countAll(PhoMSMS_C_TMT$Matches)

PhoMSMS_C_TMT$b <- PhoMSMS_C_TMT$bn/PhoMSMS_C_TMT$All *100
PhoMSMS_C_TMT$y <- PhoMSMS_C_TMT$yn/PhoMSMS_C_TMT$All *100

# siSILAC
PhoMSMS_C_siSILAC$bn <- countb(PhoMSMS_C_siSILAC$Matches)
PhoMSMS_C_siSILAC$yn <- county(PhoMSMS_C_siSILAC$Matches)
PhoMSMS_C_siSILAC$All <- countAll(PhoMSMS_C_siSILAC$Matches)

PhoMSMS_C_siSILAC$b <- PhoMSMS_C_siSILAC$bn/PhoMSMS_C_siSILAC$All *100
PhoMSMS_C_siSILAC$y <- PhoMSMS_C_siSILAC$yn/PhoMSMS_C_siSILAC$All *100

# LFQ
PhoMSMS_C_LFQ$bn <- countb(PhoMSMS_C_LFQ$Matches)
PhoMSMS_C_LFQ$yn <- county(PhoMSMS_C_LFQ$Matches)
PhoMSMS_C_LFQ$All <- countAll(PhoMSMS_C_LFQ$Matches)

PhoMSMS_C_LFQ$b <- PhoMSMS_C_LFQ$bn/PhoMSMS_C_LFQ$All *100
PhoMSMS_C_LFQ$y <- PhoMSMS_C_LFQ$yn/PhoMSMS_C_LFQ$All *100

# prepare the dataset for plotting
PhoMSMS_C_LFQ <- melt(PhoMSMS_C_LFQ[, c(1, 6, 7)], id.vars = "Sequence")
PhoMSMS_C_siSILAC <- melt(PhoMSMS_C_siSILAC[, c(1, 6, 7)], id.vars = "Sequence")
PhoMSMS_C_TMT <- melt(PhoMSMS_C_TMT[, c(1, 6, 7)], id.vars = "Sequence")

PhoMSMS_C_LFQ <- data.frame(Method = "LFQ", PhoMSMS_C_LFQ[,])
PhoMSMS_C_siSILAC <- data.frame(Method = "siSILAC", PhoMSMS_C_siSILAC[,])
PhoMSMS_C_TMT <- data.frame(Method = "TMT", PhoMSMS_C_TMT[,])

PhoMSMS_C <- rbind(PhoMSMS_C_LFQ, PhoMSMS_C_siSILAC, PhoMSMS_C_TMT)
colnames(PhoMSMS_C)[3] <- "fragment ion type"

# plot 
D2 <- ggplot(PhoMSMS_C, aes(x = value, color = `fragment ion type`))+
  geom_density(size = 0.5, show_guide=FALSE)+
  stat_density(aes(x = value, color = `fragment ion type`),
               geom="line",position="identity")+
  facet_wrap(~Method)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        strip.background = element_rect(fill="lightblue", size = 3))+
  theme(plot.title = element_text(size = 8))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size=8),
        legend.key = element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("frequency")+
  xlab("percentage")+
  labs(color = "fragment\nion type")+ 
  ggtitle("tumor tissue")+
  theme(plot.title = element_text(lineheight=.8, face="bold"))
ggsave(file="The distribution of b and y ion_pho_cell line.svg", plot=D2, width=3.2, height=2.8)

# tumor tissue
# SILAC
PhoMSMS_T_siSILAC <- PhoMSMS_T_siSILAC[which(PhoMSMS_T_siSILAC$Raw.file == "190925_ID_135_YZ_tissue_SILAC_C1_R1_enriched"),]
PhoMSMS_T_siSILAC <- PhoMSMS_T_siSILAC[, c("Sequence", "Matches")]

# LFQ
PhoMSMS_T_LFQ <- PhoMSMS_T_LFQ[which(PhoMSMS_T_LFQ$Raw.file == "190923_ID_99_YZ_tissue_LFQ_C1_R1_enriched"),]
PhoMSMS_T_LFQ <- PhoMSMS_T_LFQ[, c("Sequence", "Matches")]

# TMT
PhoMSMS_T_TMT <- PhoMSMS_T_TMT[which(PhoMSMS_T_TMT$Raw.file == "190417_ID_79_YZ_tissue_TMT_pair1_enriched"),]
PhoMSMS_T_TMT <- PhoMSMS_T_TMT[, c("Sequence", "Matches")]

# calculate the percentage of b and y ions
# TMT 
PhoMSMS_T_TMT$bn <- countb(PhoMSMS_T_TMT$Matches)
PhoMSMS_T_TMT$yn <- county(PhoMSMS_T_TMT$Matches)
PhoMSMS_T_TMT$All <- countAll(PhoMSMS_T_TMT$Matches)

PhoMSMS_T_TMT$b <- PhoMSMS_T_TMT$bn/PhoMSMS_T_TMT$All *100
PhoMSMS_T_TMT$y <- PhoMSMS_T_TMT$yn/PhoMSMS_T_TMT$All *100

# siSILAC
PhoMSMS_T_siSILAC$bn <- countb(PhoMSMS_T_siSILAC$Matches)
PhoMSMS_T_siSILAC$yn <- county(PhoMSMS_T_siSILAC$Matches)
PhoMSMS_T_siSILAC$All <- countAll(PhoMSMS_T_siSILAC$Matches)

PhoMSMS_T_siSILAC$b <- PhoMSMS_T_siSILAC$bn/PhoMSMS_T_siSILAC$All *100
PhoMSMS_T_siSILAC$y <- PhoMSMS_T_siSILAC$yn/PhoMSMS_T_siSILAC$All *100

# LFQ
PhoMSMS_T_LFQ$bn <- countb(PhoMSMS_T_LFQ$Matches)
PhoMSMS_T_LFQ$yn <- county(PhoMSMS_T_LFQ$Matches)
PhoMSMS_T_LFQ$All <- countAll(PhoMSMS_T_LFQ$Matches)

PhoMSMS_T_LFQ$b <- PhoMSMS_T_LFQ$bn/PhoMSMS_T_LFQ$All *100
PhoMSMS_T_LFQ$y <- PhoMSMS_T_LFQ$yn/PhoMSMS_T_LFQ$All *100

# prepare the dataset for plotting
PhoMSMS_T_LFQ <- melt(PhoMSMS_T_LFQ[, c(1, 6, 7)], id.vars = "Sequence")
PhoMSMS_T_siSILAC <- melt(PhoMSMS_T_siSILAC[, c(1, 6, 7)], id.vars = "Sequence")
PhoMSMS_T_TMT <- melt(PhoMSMS_T_TMT[, c(1, 6, 7)], id.vars = "Sequence")

PhoMSMS_T_LFQ <- data.frame(Method = "LFQ", PhoMSMS_T_LFQ[,])
PhoMSMS_T_siSILAC <- data.frame(Method = "siSILAC", PhoMSMS_T_siSILAC[,])
PhoMSMS_T_TMT <- data.frame(Method = "TMT", PhoMSMS_T_TMT[,])

PhoMSMS_T <- rbind(PhoMSMS_T_LFQ, PhoMSMS_T_siSILAC, PhoMSMS_T_TMT)
colnames(PhoMSMS_T)[3] <- "fragment ion type"

# plot 
D2 <- ggplot(PhoMSMS_T, aes(x = value, color = `fragment ion type`))+
  geom_density(size = 0.5, show_guide=FALSE)+
  stat_density(aes(x = value, color = `fragment ion type`),
               geom="line",position="identity")+
  facet_wrap(~Method)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        strip.background = element_rect(fill="lightblue", size = 3))+
  theme(plot.title = element_text(size = 8))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size=8),
        legend.key = element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("frequency")+
  xlab("percentage")+
  labs(color = "fragment\nion type")+ 
  ggtitle("tumor tissue")+
  theme(plot.title = element_text(lineheight=.8, face="bold"))
ggsave(file="The distribution of b and y ion_pho_tumor tissue.svg", plot=D2, width=3.2, height=2.8)

# peptides
# Cell line
# SILAC
MSMS_C_siSILAC <- MSMS_C_siSILAC[which(MSMS_C_siSILAC$Raw.file == "190508_ID_41_YZ_cellline_SILAC_C1_R1_total"),]
MSMS_C_siSILAC <- MSMS_C_siSILAC[, c("Sequence", "Matches")]

# LFQ
MSMS_C_LFQ <- MSMS_C_LFQ[which(MSMS_C_LFQ$Raw.file == "190503_ID_5_YZ_cellline_LFQ_C1_R1_total"),]
MSMS_C_LFQ <- MSMS_C_LFQ[, c("Sequence", "Matches")]

# TMT
MSMS_C_TMT <- MSMS_C_TMT[which(MSMS_C_TMT$Raw.file == "190417_ID_1_YZ_cellline_TMT_pair1_total"),]
MSMS_C_TMT <- MSMS_C_TMT[, c("Sequence", "Matches")]

# calculate the percentage of b and y ions
# TMT
MSMS_C_TMT$bn <- countb(MSMS_C_TMT$Matches)
MSMS_C_TMT$yn <- county(MSMS_C_TMT$Matches)
MSMS_C_TMT$All <- countAll(MSMS_C_TMT$Matches)

MSMS_C_TMT$b <- MSMS_C_TMT$bn/MSMS_C_TMT$All *100
MSMS_C_TMT$y <- MSMS_C_TMT$yn/MSMS_C_TMT$All *100

# siSILAC
MSMS_C_siSILAC$bn <- countb(MSMS_C_siSILAC$Matches)
MSMS_C_siSILAC$yn <- county(MSMS_C_siSILAC$Matches)
MSMS_C_siSILAC$All <- countAll(MSMS_C_siSILAC$Matches)

MSMS_C_siSILAC$b <- MSMS_C_siSILAC$bn/MSMS_C_siSILAC$All *100
MSMS_C_siSILAC$y <- MSMS_C_siSILAC$yn/MSMS_C_siSILAC$All *100

# LFQ 
MSMS_C_LFQ$bn <- countb(MSMS_C_LFQ$Matches)
MSMS_C_LFQ$yn <- county(MSMS_C_LFQ$Matches)
MSMS_C_LFQ$All <- countAll(MSMS_C_LFQ$Matches)

MSMS_C_LFQ$b <- MSMS_C_LFQ$bn/MSMS_C_LFQ$All *100
MSMS_C_LFQ$y <- MSMS_C_LFQ$yn/MSMS_C_LFQ$All *100

# prepare the dataset for plotting
MSMS_C_LFQ <- melt(MSMS_C_LFQ[, c(1, 6, 7)], id.vars = "Sequence")
MSMS_C_siSILAC <- melt(MSMS_C_siSILAC[, c(1, 6, 7)], id.vars = "Sequence")
MSMS_C_TMT <- melt(MSMS_C_TMT[, c(1, 6, 7)], id.vars = "Sequence")

MSMS_C_LFQ <- data.frame(Method = "LFQ", MSMS_C_LFQ[,])
MSMS_C_siSILAC <- data.frame(Method = "siSILAC", MSMS_C_siSILAC[,])
MSMS_C_TMT <- data.frame(Method = "TMT", MSMS_C_TMT[,])

MSMS_C <- rbind(MSMS_C_LFQ, MSMS_C_siSILAC, MSMS_C_TMT)
colnames(MSMS_C)[3] <- "fragment ion type"

# plot 
D2 <- ggplot(MSMS_C, aes(x = value, color = `fragment ion type`))+
  geom_density(size = 0.5, show_guide=FALSE)+
  stat_density(aes(x = value, color = `fragment ion type`),
               geom="line",position="identity")+
  facet_wrap(~Method)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        strip.background = element_rect(fill="lightblue", size = 3))+
  theme(plot.title = element_text(size = 8))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size=8),
        legend.key = element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("frequency")+
  xlab("percentage")+
  labs(color = "fragment\nion type")+ 
  ggtitle("tumor tissue")+
  theme(plot.title = element_text(lineheight=.8, face="bold"))
ggsave(file="The distribution of b and y ion_cell line.svg", plot=D2, width=3.2, height=2.8)

# tumor tissue
# SILAC
MSMS_T_siSILAC <- MSMS_T_siSILAC[which(MSMS_T_siSILAC$Raw.file == "190930_117_YZ_tissue_SILAC_C1_R1_total"),]
MSMS_T_siSILAC <- MSMS_T_siSILAC[, c("Sequence", "Matches")]

# LFQ
MSMS_T_LFQ <- MSMS_T_LFQ[which(MSMS_T_LFQ$Raw.file == "190927_ID_81_YZ_tissue_LFQ_C1_R1_total"),]
MSMS_T_LFQ <- MSMS_T_LFQ[, c("Sequence", "Matches")]

# TMT
MSMS_T_TMT <- MSMS_T_TMT[which(MSMS_T_TMT$Raw.file == "190417_ID_77_YZ_tissue_TMT_pair1_total"),]
MSMS_T_TMT <- MSMS_T_TMT[, c("Sequence", "Matches")]

# calculate the percentage of b and y ions
MSMS_T_TMT$bn <- countb(MSMS_T_TMT$Matches)
MSMS_T_TMT$yn <- county(MSMS_T_TMT$Matches)
MSMS_T_TMT$All <- countAll(MSMS_T_TMT$Matches)

MSMS_T_TMT$b <- MSMS_T_TMT$bn/MSMS_T_TMT$All *100
MSMS_T_TMT$y <- MSMS_T_TMT$yn/MSMS_T_TMT$All *100

MSMS_T_siSILAC$bn <- countb(MSMS_T_siSILAC$Matches)
MSMS_T_siSILAC$yn <- county(MSMS_T_siSILAC$Matches)
MSMS_T_siSILAC$All <- countAll(MSMS_T_siSILAC$Matches)

MSMS_T_siSILAC$b <- MSMS_T_siSILAC$bn/MSMS_T_siSILAC$All *100
MSMS_T_siSILAC$y <- MSMS_T_siSILAC$yn/MSMS_T_siSILAC$All *100

MSMS_T_LFQ$bn <- countb(MSMS_T_LFQ$Matches)
MSMS_T_LFQ$yn <- county(MSMS_T_LFQ$Matches)
MSMS_T_LFQ$All <- countAll(MSMS_T_LFQ$Matches)

MSMS_T_LFQ$b <- MSMS_T_LFQ$bn/MSMS_T_LFQ$All *100
MSMS_T_LFQ$y <- MSMS_T_LFQ$yn/MSMS_T_LFQ$All *100

# prepare the dataset for plotting
MSMS_T_LFQ <- melt(MSMS_T_LFQ[, c(1, 6, 7)], id.vars = "Sequence")
MSMS_T_siSILAC <- melt(MSMS_T_siSILAC[, c(1, 6, 7)], id.vars = "Sequence")
MSMS_T_TMT <- melt(MSMS_T_TMT[, c(1, 6, 7)], id.vars = "Sequence")

MSMS_T_LFQ <- data.frame(Method = "LFQ", MSMS_T_LFQ[,])
MSMS_T_siSILAC <- data.frame(Method = "siSILAC", MSMS_T_siSILAC[,])
MSMS_T_TMT <- data.frame(Method = "TMT", MSMS_T_TMT[,])

MSMS_T <- rbind(MSMS_T_LFQ, MSMS_T_siSILAC, MSMS_T_TMT)
colnames(MSMS_T)[3] <- "fragment ion type"

# plot 
D2 <- ggplot(MSMS_T, aes(x = value, color = `fragment ion type`))+
  geom_density(size = 0.5, show_guide=FALSE)+
  stat_density(aes(x = value, color = `fragment ion type`),
               geom="line",position="identity")+
  facet_wrap(~Method)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        strip.background = element_rect(fill="lightblue", size = 3))+
  theme(plot.title = element_text(size = 8))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size=8),
        legend.key = element_rect(fill = "white", colour = "white"))+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  ylab("frequency")+
  xlab("percentage")+
  labs(color = "fragment\nion type")+ 
  ggtitle("tumor tissue")+
  theme(plot.title = element_text(lineheight=.8, face="bold"))
ggsave(file="The distribution of b and y ion_tumor tissue.svg", plot=D2, width=3.2, height=2.8)

