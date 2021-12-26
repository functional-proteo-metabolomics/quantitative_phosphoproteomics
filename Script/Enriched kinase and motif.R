######################################
# The script plots the kinase enrichment analysis and sequence motif analysis.
# Input data is available either in /Dataset directory in RData format or in ProteomeXchange PXD030450 dataset.
# December 20, 2021, Yang Zhang

# clean all variables
rm(list=ls())

# this script is used to plotted the kinase enrichment analysis and sequence motif analysis

# library
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require('clusterProfiler')) BiocManager::install('clusterProfiler')
library(clusterProfiler)
if (!require('tidyverse')) install.packages('tidyverse')
library(tidyverse)

# kinase enrichment analysis
# load the phosphosites which is quantified as significantly differential 
load(file = "Dataset\\ORA_data.rdata")
Phosphosites_Q <- x

# read the kinase-substrate database 
KS_database <- read.table("Kinase_Substrate_Dataset_PhosphositesPlus.txt", header=T, sep = "\t", quote = "", fill = TRUE)
KS_database <- KS_database %>% 
  filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human")

KS_database$Phosphosites <- unlist(lapply(paste(KS_database$SUB_ACC_ID, KS_database$SUB_MOD_RSD, sep = '_'),"[",1))

KS_database <- KS_database %>% 
  dplyr::select(KINASE, Phosphosites)  %>% 
  dplyr::rename(term = KINASE, gene = Phosphosites)

Kinases_enrich <- vector(mode = "list", length = 10)
Kinases_enrich
class(Kinases_enrich)

for (i in 1: length(Phosphosites_Q)) {
  Kinases_enrich[[i]] <- enricher(Phosphosites_Q[[i]] , pvalueCutoff = 0.05, pAdjustMethod = "BH", 
         minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE = KS_database)@result 
  names(Kinases_enrich)[i]<- names(Phosphosites_Q)[i] 
}
 

for (i in 1: length(Kinases_enrich)) {
  Kinases_enrich[[i]] <- Kinases_enrich[[i]][which(Kinases_enrich[[i]]$ID == "CDK2" | Kinases_enrich[[i]]$ID == "ERK1" | Kinases_enrich[[i]]$ID == "ERK2" | Kinases_enrich[[i]]$ID == "Akt1" | 
                       Kinases_enrich[[i]]$ID == "AMPKA1" | Kinases_enrich[[i]]$ID == "PAK1" | Kinases_enrich[[i]]$ID == "CDK5" 
                       |Kinases_enrich[[i]]$ID == "mTOR" | Kinases_enrich[[i]]$ID == "Chk1" | Kinases_enrich[[i]]$ID ==  "AMPKA1"), c(1, 6, 9)]
}


Kinases_enrich[[1]] <- data.frame(method = "LFQ", Kinases_enrich[[1]][,])
Kinases_enrich[[2]] <- data.frame(method = "LFQ+MBR", Kinases_enrich[[2]][,])
Kinases_enrich[[3]] <- data.frame(method = "LFQ", Kinases_enrich[[3]][,])
Kinases_enrich[[4]] <- data.frame(method = "LFQ+MBR", Kinases_enrich[[4]][,])

Kinases_enrich[[5]] <- data.frame(method = "siSILAC+MBRreq", Kinases_enrich[[5]][,])
Kinases_enrich[[6]] <- data.frame(method = "siSILAC", Kinases_enrich[[6]][,])
Kinases_enrich[[7]] <- data.frame(method = "siSILAC+MBRreq", Kinases_enrich[[7]][,])
Kinases_enrich[[8]] <- data.frame(method = "siSILAC", Kinases_enrich[[8]][,])

Kinases_enrich[[9]] <- data.frame(method = "TMT", Kinases_enrich[[9]][,])
Kinases_enrich[[10]] <- data.frame(method = "TMT", Kinases_enrich[[10]][,])

cellline <- rbind(Kinases_enrich[[1]], Kinases_enrich[[2]], Kinases_enrich[[5]], Kinases_enrich[[6]], Kinases_enrich[[9]])
tissue <- rbind(Kinases_enrich[[3]], Kinases_enrich[[4]], Kinases_enrich[[7]], Kinases_enrich[[8]], Kinases_enrich[[10]])


tissue$method <- factor(tissue$method, levels=c("LFQ+MBR", "LFQ", "siSILAC+MBRreq", "siSILAC", "TMT"))
cellline$method <- factor(cellline$method, levels=c("LFQ+MBR", "LFQ", "siSILAC+MBRreq", "siSILAC", "TMT"))
tissue$ID <- factor(tissue$ID, levels=c("Chk1", "mTOR", "CDK5", "PAK1", "AMPKA1", "Akt1", "ERK2", "ERK1", "CDK2"))
cellline$ID <- factor(cellline$ID, levels=c("Chk1", "mTOR", "CDK5", "PAK1", "AMPKA1", "Akt1", "ERK2", "ERK1", "CDK2"))

# plot()
N1 <- ggplot(tissue) + 
  aes(x=method, y=ID, size=Count, col=-log10(p.adjust)) + 
  geom_point()+
  scale_colour_gradient(low = "blue", high = "red") +
  theme_light()+
  theme(axis.text=element_text(size=9))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size=8)) +
  xlab("")+
  ylab("")
ggsave(file="ORA_tumor tissue.svg", plot=N1, width=6, height=6)

N2 <- ggplot(cellline) + 
  aes(x=method, y=ID, size=Count, col=-log10(p.adjust)) + 
  geom_point()+
  scale_colour_gradient(low = "blue", high = "red") +
  theme_light()+
  theme(axis.text=element_text(size=9))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size=8)) +
  xlab("")+
  ylab("")
ggsave(file="ORA_cell line.svg", plot=N2, width=6, height=6)


# motif analysis, senquence log figure of the phosphosites 
# get the uniprot ID and Sequence.
ProteinDatabase <- read.fasta(file = "20201111_uniprot-human-filtered-reviewed_yes+AND+organism__Homo+sapiens+(Human)--.fasta", seqtype = c("AA"), as.string = TRUE)
ProteinSequence = as.data.frame(unlist(ProteinDatabase))
ProteinSequence$ProteinID <- row.names(ProteinSequence)
ProteinSequence$ProteinID <- unlist(lapply(strsplit(as.vector(ProteinSequence$ProteinID),"\\|"), "[", 2))
row.names(ProteinSequence) <- NULL
colnames(ProteinSequence) <-c("Sequence", "uniprotID")
ProteinSequence <- ProteinSequence[c(2, 1)]
class(ProteinSequence)
str(ProteinSequence)
ProteinSequence$Sequence <- as.character(ProteinSequence$Sequence)

for (i in 1: length(Phosphosites_Q)) {
  Phosphosites_Q[[i]] <- as.data.frame(Phosphosites_Q[[i]])
  Phosphosites_Q[[i]]$uniprotID <- unlist(lapply(strsplit(as.vector(Phosphosites_Q[[i]]$`Phosphosites_Q[[i]]`),"\\_"),"[",1))
  Phosphosites_Q[[i]]$localization <- unlist(lapply(strsplit(as.vector(Phosphosites_Q[[i]]$`Phosphosites_Q[[i]]`),"\\_"),"[",2))
  Phosphosites_Q[[i]]$localization <- substring(Phosphosites_Q[[i]]$localization, 2)
  Phosphosites_Q[[i]] <- merge(Phosphosites_Q[[i]], ProteinSequence, by="uniprotID")
  Phosphosites_Q[[i]]$Sequence1 <- str_sub(Phosphosites_Q[[i]]$Sequence, start = as.numeric(Phosphosites_Q[[i]]$localization)-5, end = as.numeric(Phosphosites_Q[[i]]$localization)+5)
  Phosphosites_Q[[i]]$Number <- nchar(Phosphosites_Q[[i]]$Sequence1)
  Phosphosites_Q[[i]] <- subset(Phosphosites_Q[[i]], Phosphosites_Q[[i]]$Number > 10)
}
  
N1 <- 
ggseqlogo(Phosphosites_Q[[1]]$Sequence1, method = 'prob', seq_type='aa' )
ggsave(file="seqlogo_LFQ_cellline.svg", plot=N1, width=7, height=3)


N2 <- 
ggseqlogo(Phosphosites_Q[[2]]$Sequence1, method = 'prob', seq_type='aa' )
ggsave(file="seqlogo_LFQ with MBR_cellline.svg", plot=N2, width=7, height=3)

N3 <- 
ggseqlogo(Phosphosites_Q[[3]]$Sequence1, method = 'prob', seq_type='aa' )
ggsave(file="seqlogo_LFQ_tissue.svg", plot=N3, width=7, height=3)
 
N4 <- 
ggseqlogo(Phosphosites_Q[[4]]$Sequence1, method = 'prob', seq_type='aa' )
ggsave(file="seqlogo_LFQ with MBR_tissue.svg", plot=N4, width=7, height=3)
 
N5 <- 
ggseqlogo(Phosphosites_Q[[5]]$Sequence1, method = 'prob', seq_type='aa' )
ggsave(file="seqlogo_siSILAC_cellline.svg", plot=N5, width=7, height=3)


N6 <- 
ggseqlogo(Phosphosites_Q[[6]]$Sequence1, method = 'prob', seq_type='aa' )
ggsave(file="seqlogo_siSILAC with MBRreq_cellline.svg", plot=N6, width=7, height=3)

N7 <- 
ggseqlogo(Phosphosites_Q[[7]]$Sequence1, method = 'prob', seq_type='aa' )
ggsave(file="seqlogo_siSILAC_tissue.svg", plot=N7, width=7, height=3)
 
N8 <- 
ggseqlogo(Phosphosites_Q[[8]]$Sequence1, method = 'prob', seq_type='aa' )
ggsave(file="seqlogo_siSILAC with MBRreq_tissue.svg", plot=N8, width=7, height=3)
 
N9 <- 
ggseqlogo(Phosphosites_Q[[9]]$Sequence1, method = 'prob', seq_type='aa' )
ggsave(file="seqlogo_TMT_cellline.svg", plot=N9, width=7, height=3)


N10 <- 
ggseqlogo(Phosphosites_Q[[10]]$Sequence1, method = 'prob', seq_type='aa' )
ggsave(file="seqlogo_TMT_tissue.svg", plot=N10, width=7, height=3)
