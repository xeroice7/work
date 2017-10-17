### Packages 
install.packages("scales")
install.packages("devtools")
install.packages("cowplot")
install.packages("data.table")
install.packages("tidyverse")

library(scales)
library(reshape2)
library(devtools)
library(cowplot)
library(data.table)
library(tidyverse)
library(VennDiagram)

#iPS to MCF7 Comparisons (No Somatic Filters - can include somatic hits)
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsHUViPS/")
mcf7huvips_overlap <- data.frame(read.csv(paste(file.path, "HUViPSMCF7OverlapGENE.csv", sep=""), header = FALSE, sep = ",", stringsAsFactors = FALSE))

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/MCF7vsFiPS/")
mcf7fips_overlap <- data.frame(read.csv(paste(file.path, "FiPSMCF7OverlapGENE.csv", sep=""), header = FALSE, sep = ",", stringsAsFactors = FALSE))

mcf7ips_overlap <- intersect(mcf7fips_overlap$V1, mcf7huvips_overlap$V1)

write.csv(mcf7ips_overlap, "MCF7SomaticiPSOverlap.csv")
