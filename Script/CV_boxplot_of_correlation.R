######################################
# The script plots the density plot of our CV values and the boxplot of the Pearson correlation parameters.
# Input data is available either in /Dataset directory in RData format or in ProteomeXchange PXD030450 dataset.
# December 20, 2021, Yang Zhang

# clean all variables
rm(list=ls())

# the density plot of the CV 
# library
if (!require('reshape')) install.packages('reshape')
library(reshape)
if (!require('ggplot2')) install.packages('ggplot2')
library(ggplot2)

# functions 
NAcount <- function(x){
  sum(is.na(x))
} 

CV <- function(x){
  sd <- sd(x, na.rm = TRUE)
  mean <- mean(x, na.rm=TRUE)
  y <- (sd/mean)*100
  return(y)
}

gg_color_hue <- function(n) {
  hcl(h = n, l = 65, c = 100)[1:n]
}

# load the dataset 
load(file = "Dataset\\Quantificationlist_Psites.rdata")
Quantificationlist_Psites <- x

CV_phosphosites <- vector(mode = "list", length = 10)
CV_phosphosites
class(CV_phosphosites)

for (i in 1:length(Quantificationlist_Psites)) {
  Quantificationlist_Psites[[i]]$NA_Sample6X <- apply(Quantificationlist_Psites[[i]][ , c(2:7)], 1, NAcount)
  Quantificationlist_Psites[[i]]$NA_Sample3X <- apply(Quantificationlist_Psites[[i]][ , c(8:13)], 1, NAcount)
  Quantificationlist_Psites[[i]]$NA_Sample1X <- apply(Quantificationlist_Psites[[i]][ , c(14:19)], 1, NAcount)
  
  Quantificationlist_Psites[[i]] <- Quantificationlist_Psites[[i]][which(Quantificationlist_Psites[[i]]$NA_Sample6X < 4 
                                                                         & Quantificationlist_Psites[[i]]$NA_Sample3X < 4 
                                                                         & Quantificationlist_Psites[[i]]$NA_Sample1X < 4), ]
  
  Quantificationlist_Psites[[i]]$"sample-6X" <- apply(Quantificationlist_Psites[[i]][ , c(2:7)], 1, CV)
  Quantificationlist_Psites[[i]]$"sample-3X" <- apply(Quantificationlist_Psites[[i]][ , c(8:13)], 1, CV)
  Quantificationlist_Psites[[i]]$"sample-1X" <- apply(Quantificationlist_Psites[[i]][ , c(14:19)], 1, CV)
  
  CV_phosphosites[[i]] <- Quantificationlist_Psites[[i]][, -c(2:22)]
  CV_phosphosites[[i]] <- melt(CV_phosphosites[[i]], id.vars = "Phosphosites", mearsure.vars = c("sample-6X", "sample-3X", "sample-1X"))
}

names(CV_phosphosites) <- names(Quantificationlist_Psites)

CV_phosphosites[[1]] <- data.frame(method = "LFQ_MBR", CV_phosphosites[[1]][,])
CV_phosphosites[[2]] <- data.frame(method = "LFQ", CV_phosphosites[[2]][,])
CV_phosphosites[[3]] <- data.frame(method = "siSILAC_MBRreq", CV_phosphosites[[3]][,])
CV_phosphosites[[4]] <- data.frame(method = "siSILAC", CV_phosphosites[[4]][,])
CV_phosphosites[[5]] <- data.frame(method = "TMT", CV_phosphosites[[5]][,])

CV_phosphosites[[6]] <- data.frame(method = "LFQ_MBR", CV_phosphosites[[6]][,])
CV_phosphosites[[7]] <- data.frame(method = "LFQ", CV_phosphosites[[7]][,])
CV_phosphosites[[8]] <- data.frame(method = "siSILAC_MBRreq", CV_phosphosites[[8]][,])
CV_phosphosites[[9]] <- data.frame(method = "siSILAC", CV_phosphosites[[9]][,])
CV_phosphosites[[10]] <- data.frame(method = "TMT", CV_phosphosites[[10]][,])


# plot the CV density plot for LFQ with MBR, siSILAC with MBRreq and TMT 
# tumor tissue
CV_tissue1 <- rbind(CV_phosphosites[[6]], CV_phosphosites[[8]], CV_phosphosites[[10]])
colnames(CV_tissue1)[3] <- "sample group"

CV_tissue1$method <- factor(CV_tissue1$method,
                             levels = c('LFQ_MBR','siSILAC_MBRreq','TMT'), ordered = TRUE)

a1 <- ggplot(CV_tissue1, aes(x = value, color = method, linetype = `sample group`)) + 
  geom_density(size = 0.5, show_guide = FALSE)+
  stat_density(aes(x = value, color = method, linetype = `sample group`),
               geom = "line",position = "identity")+
  scale_linetype_manual(values = c("solid", "dotted", "twodash"))+ 
  xlab("CV")+
  ylab("frequency")+
  xlim(0,120)+
  scale_color_manual(values = gg_color_hue(3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        strip.background = element_rect(fill ="lightblue", size = 2))+
  theme(plot.title = element_text(size = 10))+
  theme(legend.text = element_text(size = 6), legend.title = element_text(size = 8),
        legend.key = element_rect(fill = "white", colour = "white"))+
  theme(legend.position = c(0.75, 0.55))+
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  ggtitle("tumor tissue")+
  theme(plot.title = element_text(lineheight = .8, face = "bold"))

ggsave(file="CV_Phosphosites_tumor.svg", plot = a1, width = 2.8, height = 2.8)

# cell line
CV_celline1 <- rbind(CV_phosphosites[[1]], CV_phosphosites[[3]], CV_phosphosites[[5]])
colnames(CV_celline1)[3] <- "sample group"

CV_celline1$method <- factor(CV_celline1$method,
                             levels = c('LFQ_MBR','siSILAC_MBRreq','TMT'), ordered = TRUE)

a1 <- ggplot(CV_celline1, aes(x = value, color = method, linetype = `sample group`)) + 
  geom_density(size = 0.5, show_guide = FALSE) +
  stat_density(aes(x = value, color = method, linetype = `sample group`),
               geom = "line",position = "identity") +
  scale_linetype_manual(values = c("solid", "dotted", "twodash")) + 
  xlab("CV") +
  ylab("frequency") +
  xlim(0,120) +
  scale_color_manual(values = gg_color_hue(3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        strip.background = element_rect(fill ="lightblue", size = 2)) +
  theme(plot.title = element_text(size = 10)) +
  theme(legend.text = element_text(size = 6), legend.title = element_text(size = 8),
        legend.key = element_rect(fill = "white", colour = "white")) +
  theme(legend.position = c(0.75, 0.55)) +
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) + 
  ggtitle("tumor celline") +
  theme(plot.title = element_text(lineheight = .8, face = "bold"))

ggsave(file="CV_Phosphosites_cellline.svg", plot = a1, width = 2.8, height = 2.8)

# the CV density plot for the  

#### #### 




