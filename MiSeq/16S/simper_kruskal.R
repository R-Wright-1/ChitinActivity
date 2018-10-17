## Gabriel Erni Cassola
## 28.02.2018
## MiSeq samples (PhD)
#______________________#

# following SOP for mothur processed MiSeq data available from: http://rpubs.com/dillmcfarlan/R_microbiotaSOP

#install.packages("ape")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("gplots")
#install.packages("lme4")
#install.packages("phangorn")
#install.packages("plotly")
#install.packages("tidyr")
#install.packages("vegan")
#install.packages("VennDiagram")
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

#packages
library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(phyloseq)

OTU2 <- read.table(file = "/Users/u1560915/Documents/OneDrive/PhD_PlasticOceans/Experiments/MiSeq/r_and_python_chitin/12_LS_percent_r.csv", header = TRUE, sep = ",")
#0_percent_r.csv, 2_daily_r_percent_r.csv, 3_daily_percent_r.csv, 4_LG_percent_r.csv, 5_RG_percent_r.csv, 6_beg_end_percent_r.csv, 7_beginning_percent_r.csv, 8_beg_end_percent_r.csv, 9_long_GR_percent_r.csv, 10_short_GR_percent_r.csv, 11_beg_end_S_percent_r.csv, 12_LS_percent_r.csv

tax2 <- read.table(file = "/Users/u1560915/Documents/OneDrive/PhD_PlasticOceans/Experiments/MiSeq/r_and_python_chitin/stability.final.pick.taxonomy", header = TRUE, sep = "\t")

meta2 <- read.table(file = "/Users/u1560915/Documents/OneDrive/PhD_PlasticOceans/Experiments/MiSeq/r_and_python_chitin/12_LS_meta.csv", header = TRUE, row.names = 1, sep = ",")
#0_meta.csv, 2_daily_r_meta.csv, 3_daily_meta.csv, 4_LG_meta.csv, 5_RG_meta.csv, 6_beg_end_meta.csv, 7_beginning_meta.csv, 8_beg_end_meta.csv, 9_long_GR_meta.csv, 10_short_GR_meta.csv, 11_beg_end_S_meta.csv, 12_LS_meta.csv

row.names(OTU2) <- OTU2$Group
OTU2.clean <- OTU2[, -which(names(OTU2) %in% c("label", "numOtus", "Group"))]
OTU2.clean <- OTU2.clean[order(row.names(OTU2.clean)),]

row.names(tax2) <- tax2$OTU
tax.2 <- tax2[row.names(tax2) %in% colnames(OTU2.clean),]
tax.clean.a <- separate(tax.2, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = "\\;")
tax.clean <- tax.clean.a[, -which(names(tax.clean.a) %in% c("Size", "Strain", "OTU"))]

source("/Users/u1560915/Documents/OneDrive/PhD_PlasticOceans/Experiments/MiSeq/r_and_python_chitin/simper_pretty.r")
source("/Users/u1560915/Documents/OneDrive/PhD_PlasticOceans/Experiments/MiSeq/r_and_python_chitin/R_krusk.r")

simper.pretty(OTU2.clean, meta2, c("treatment"), perc_cutoff = 1, low_cutoff = "y", low_val = 0.01, "trt")

simper.results <- data.frame(read.csv("/Users/u1560915/Documents/OneDrive/PhD_PlasticOceans/Experiments/MiSeq/r_and_python_chitin/trt_clean_simper.csv"))
kruskal.pretty(OTU2.clean, meta2, simper.results, c("treatment"), "trt", tax2)

KW.results <- data.frame(read.csv("/Users/u1560915/Documents/OneDrive/PhD_PlasticOceans/Experiments/MiSeq/r_and_python_chitin/trt_krusk_simper.csv"))
KW.results.signif <- KW.results[KW.results$fdr_krusk_p.val < 0.05,] 			#extract significant results
KW.results.signif <- KW.results.signif[with(KW.results.signif, order(OTU)),]	#order OTUs by number
head(KW.results.signif)









