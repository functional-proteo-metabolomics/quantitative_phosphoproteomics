These scripts were used for analyzing the phosphoproteome data from Ovarian cancer tissue and SKOV3 cell line.


# Requirements:
These scripts depend on the following packages, which will be installed once first execution of the script:
reshape
ggplot2
eulerr
wesanderson
clusterProfiler-
tidyverse
stringr
ggseqlogo



# Before running these scripts some things have to be specified:
All the raw and MaxQuant processed  files has been submitted to ProteomeXchange with PXD030450 identifier. 
The phosphosites identification and quantification information are extrated from the "Phospho (STY)Sites.txt" file.
The phospho-peptides identification information are extrated from the "modificationSpecificPeptides.txt" file.
The b and y ion information are extracted from the "msms.txt" file.

# Scripts are in Scripts directory
	
* The [*"Distribution of b,y ions.R"*](Script/Distribution of b,y ions.R) script calculates the distribution of b and y ion for the identifed phosphopeptides and peptides.
	
* The [*"Enriched kinase and motif.R"*](Script/Enriched kinase and motif.R) script plots the kinase enrichment analysis and sequence motif analysis.
	
* The [*"CV plot and boxplot of the correlation parameters.R"*](Script/CV plot and boxplot of the correlation parameters.R) script plots the density plot of our CV values and the boxplot of the Pearson correlation parameters.
	
* The [*"VennDiagram.R"*](Script/VennDiagram.R) script plots Venn diagram of the analysed data.
