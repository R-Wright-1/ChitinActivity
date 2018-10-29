## Robyn Wright
## MiSeq samples
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
#library(ape)
#library(dplyr)
#library(ggplot2)
#library(gplots)
#library(lme4)
#library(phangorn)
#library(plotly)
#library(tidyr)
#library(VennDiagram)
#library(phyloseq)

library(vegan)


OTU2 <- read.table(file = "/Users/u1560915/Documents/GitHub/CommunityAnalysis/InitialAnalysis_18S/treatment1_daily_samples_for_simper_r.csv", header = TRUE, sep = ",")
# treatment2_all_rg_ls_samples_for_simper_r.csv

meta2 <- read.table(file = "/Users/u1560915/Documents/GitHub/CommunityAnalysis/InitialAnalysis_18S/treatment1_daily_meta.csv", header = TRUE, row.names = 1, sep = ",")
# treatment2_all_rg_ls_meta.csv

row.names(OTU2) <- OTU2$Group
OTU2.clean <- OTU2[, -which(names(OTU2) %in% c("label", "numOtus", "Group"))]
OTU2.clean <- OTU2.clean[order(row.names(OTU2.clean)),]

source("/Users/u1560915/Documents/GitHub/CommunityAnalysis/InitialAnalysis_18S/simper_pretty.r")

simper.pretty(OTU2.clean, meta2, c("treatment"), perc_cutoff = 1, low_cutoff = "y", low_val = 0.01, "treatment1_daily")









