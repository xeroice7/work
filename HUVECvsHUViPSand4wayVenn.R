### Packages 
install.packages("dplyr")
install.packages("ggplot2")
install.packages("stringr")
install.packages("scales")
install.packages("devtools")
install.packages("cowplot")
install.packages("data.table")
install.packages("tidyverse")
install.packages("VennDiagram")
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(reshape2)
library(devtools)
library(cowplot)
library(data.table)
library(tidyverse)
library(VennDiagram)

#Functions needed
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

# HUVEC/HUViPS comparison
#Tidying Data
huv <- data.frame(read.csv(paste(file.path, "HUVECHUViPSSummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(huv$X, " ", 5)
huv <- cbind.data.frame(splitcolumns, huv)
names(huv)[names(huv) == '1'] <- 'Line'  #Change the first column to "Line"
names(huv)[names(huv) == '3'] <- 'PassNum'  #Change the second column to "PassNum"
names(huv)[names(huv) == '4'] <- 'Biotin'  #Change the first column to "Line"
huv <- subset(huv, select=-c(`2`, `5`, X))
huv <- filter(huv, Species == "HUMAN") #Filtering only human hits
splitcolumns <- str_split_fixed(huv$Gene.Names," ", 2)
huv <- cbind.data.frame(splitcolumns, huv)
names(huv)[names(huv) == '1'] <- 'Gene'  #Change the first column to "Line"
names(huv)[names(huv) == '2'] <- 'Synonyms'  #Change the second column to "PassNum"
huv <- subset(huv, select=-Gene.Names)
huvslim <- data.frame(huv$Line, huv$Gene, huv$Synonyms, huv$Biotin, huv$Unused, huv$Total, huv$Gene.Ontology)
colnames(huvslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

#Analysis of HUVEC
huvecslim <- filter(huvslim, Line == "HUVEC")  #Keep only HUVEC
huvecbiotin <- filter(huvecslim, Biotin == "+") #Keep only +Biotin for HUVEC
colnames(huvecbiotin)[colnames(huvecbiotin) == 'Unused'] <- 'UnusedHUVECBiotin'
colnames(huvecbiotin)[colnames(huvecbiotin) == 'Total'] <- 'TotalHUVECBiotin'

huvecnobiotin <- filter(huvecslim, Biotin == "No") #Keep only No Biotin for HUVEC
colnames(huvecnobiotin)[colnames(huvecnobiotin) == 'Unused'] <- 'UnusedHUVECNoBiotin'
colnames(huvecnobiotin)[colnames(huvecnobiotin) == 'Total'] <- 'TotalHUVECNoBiotin'

huvec <- merge.data.frame(huvecbiotin, huvecnobiotin, by = "Gene", all = TRUE)
huvec$UnusedHUVECBiotin <- as.numeric(as.character(huvec$UnusedHUVECBiotin))
huvec$UnusedHUVECNoBiotin <- as.numeric(as.character(huvec$UnusedHUVECNoBiotin))
huvec$TotalHUVECBiotin <- as.numeric(as.character(huvec$TotalHUVECBiotin))
huvec$TotalHUVECNoBiotin <- as.numeric(as.character(huvec$TotalHUVECNoBiotin))

huvec$UnusedHUVECBiotin <- na.zero(huvec$UnusedHUVECBiotin)
huvec$UnusedHUVECNoBiotin <- na.zero(huvec$UnusedHUVECNoBiotin)
huvec$TotalHUVECBiotin <- na.zero(huvec$TotalHUVECBiotin)
huvec$TotalHUVECNoBiotin <- na.zero(huvec$TotalHUVECNoBiotin)

huvecdf1 <- data.frame(huvecnobiotin$Gene, huvecnobiotin$Line)
colnames(huvecdf1) <- c("Gene", "Line")
huvecno1 <- nrow(huvecdf1)
huvecdf2 <- data.frame(huvecbiotin$Gene, huvecbiotin$Line)
colnames(huvecdf2) <- c("Gene", "Line")
huvecno2 <- nrow(huvecdf2)

huvecdf3 <- setdiff(huvecdf1, huvecdf2)
huvecdf4 <- setdiff(huvecdf2, huvecdf1)
huvecdf5 <- semi_join(huvecdf1, huvecdf2)
huvecno5 <- nrow(huvecdf5)

grid.newpage()
draw.pairwise.venn(huvecno1, huvecno2, huvecno5, category = c("HUVEC No Biotin", "HUVEC + Biotin"), lty = rep("solid", 2), 
                  lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("sienna1", "slateblue1"), alpha = rep(0.5, 1), 
                  cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

huvecbiotinunique <- anti_join(huvecdf2, huvecdf1)
write.csv(huvecbiotinunique, "UniqueHUVECBiotin.csv")

#Analysis of HUViPS
huvipsslim <- filter(huvslim, Line == "HUViPS4F1")  #Keep only HUViPS4F1
huvipsbiotin <- filter(huvipsslim, Biotin == "+") #Keep only +Biotin for HUViPS4F1
colnames(huvipsbiotin)[colnames(huvipsbiotin) == 'Unused'] <- 'UnusedHUViPSBiotin'
colnames(huvipsbiotin)[colnames(huvipsbiotin) == 'Total'] <- 'TotalHUViPSBiotin'

huvipsnobiotin <- filter(huvipsslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(huvipsnobiotin)[colnames(huvipsnobiotin) == 'Unused'] <- 'UnusedHUViPSNoBiotin'
colnames(huvipsnobiotin)[colnames(huvipsnobiotin) == 'Total'] <- 'TotalHUViPSNoBiotin'

huvips <- merge.data.frame(huvipsbiotin, huvipsnobiotin, by = "Gene", all = TRUE)
huvips$UnusedHUViPSBiotin <- as.numeric(as.character(huvips$UnusedHUViPSBiotin))
huvips$UnusedHUViPSNoBiotin <- as.numeric(as.character(huvips$UnusedHUViPSNoBiotin))
huvips$TotalHUViPSBiotin <- as.numeric(as.character(huvips$TotalHUViPSBiotin))
huvips$TotalHUViPSNoBiotin <- as.numeric(as.character(huvips$TotalHUViPSNoBiotin))

huvips$UnusedHUViPSBiotin <- na.zero(huvips$UnusedHUViPSBiotin)
huvips$UnusedHUViPSNoBiotin <- na.zero(huvips$UnusedHUViPSNoBiotin)
huvips$TotalHUViPSBiotin <- na.zero(huvips$TotalHUViPSBiotin)
huvips$TotalHUViPSNoBiotin <- na.zero(huvips$TotalHUViPSNoBiotin)

huvdf1 <- data.frame(huvipsnobiotin$Gene, huvipsnobiotin$Line)
colnames(huvdf1) <- c("Gene", "Line")
huvno1 <- nrow(huvdf1)
huvdf2 <- data.frame(huvipsbiotin$Gene, huvipsbiotin$Line)
colnames(huvdf2) <- c("Gene", "Line")
huvno2 <- nrow(huvdf2)

huvdf3 <- setdiff(huvdf1, huvdf2)
huvdf4 <- setdiff(huvdf2, huvdf1)

huvdf5 <- semi_join(huvdf1, huvdf2)
huvno5 <- nrow(huvdf5)

grid.newpage()
draw.pairwise.venn(huvno1, huvno2, huvno5, category = c("HUViPS4F5 No Biotin", "HUViPS4F5 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("yellow", "dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

huvipsbiotinunique <- anti_join(huvdf2, huvdf1)
write.csv(huvipsbiotinunique, "UniqueHUViPSBiotin.csv")

# To compare together
huvectotal <- filter(huvec, TotalHUVECNoBiotin < TotalHUVECBiotin)
huvipstotal <- filter(huvips, TotalHUViPSNoBiotin < TotalHUViPSBiotin)

huvecdf6 <- data.frame(huvectotal$Gene)
colnames(huvecdf6) <- c("Gene")
huvecno6 <- nrow(huvecdf6)
huvdf6 <- data.frame(huvipstotal$Gene)
colnames(huvdf6) <- c("Gene")
huvno6 <- nrow(huvdf6)

huvhuvecdf1 <- setdiff(huvecdf6, huvdf6)
huvhuvecdf2 <- setdiff(huvdf6, huvecdf6)
huvhuvecdf3 <- semi_join(huvecdf6, huvdf6)
huvhuvecno3 <- nrow(huvhuvecdf3)

grid.newpage()
draw.pairwise.venn(huvecno6, huvno6, huvhuvecno3, category = c("HUVEC", "HUViPS4F1"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("slateblue1","dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

huvecunique <- anti_join(huvecdf6, huvdf6)
huvipsunique <- anti_join(huvdf6, huvecdf6)
huvechuvips <- intersect(huvecdf6$Gene, huvdf6$Gene)
write.csv(huvecunique, "UniqueHUVECtoHUViPS.csv")
write.csv(huvipsunique, "UniqueHUViPStoHUVEC.csv")
write.csv(huvechuvips, "HUViPSHUVECOverlap.csv")
write.csv(huvectotal$Gene, "HUVECReal.csv")
write.csv(huvipstotal$Gene, "HUViPSReal.csv")

#To compare which is bigger TOTAL
huvecvshuvips <- merge.data.frame(huvipstotal, huvectotal, by = "Gene", all = TRUE)

huvecvshuvips$UnusedHUViPSBiotin <- na.zero(huvecvshuvips$UnusedHUViPSBiotin)
huvecvshuvips$UnusedHUViPSNoBiotin <- na.zero(huvecvshuvips$UnusedHUViPSNoBiotin)
huvecvshuvips$TotalHUViPSBiotin <- na.zero(huvecvshuvips$TotalHUViPSBiotin)
huvecvshuvips$TotalHUViPSNoBiotin <- na.zero(huvecvshuvips$TotalHUViPSNoBiotin)
huvecvshuvips$UnusedHUVECBiotin <- na.zero(huvecvshuvips$UnusedHUVECBiotin)
huvecvshuvips$UnusedHUVECNoBiotin <- na.zero(huvecvshuvips$UnusedHUVECNoBiotin)
huvecvshuvips$TotalHUVECBiotin <- na.zero(huvecvshuvips$TotalHUVECBiotin)
huvecvshuvips$TotalHUVECNoBiotin <- na.zero(huvecvshuvips$TotalHUVECNoBiotin)

huvipsup <- filter(huvecvshuvips, TotalHUViPSBiotin > TotalHUVECBiotin)
huvecup <- filter(huvecvshuvips, TotalHUViPSBiotin < TotalHUVECBiotin)
write.csv(huvipsup$Gene, "HUViPSUp.csv")
write.csv(huvecup$Gene, "HUVECUp.csv")

#Getting rid of duplicates
huvecdf1 <- data.frame(unique(huvecdf1$Gene))
huvecdf2 <- data.frame(unique(huvecdf2$Gene))
huvdf1 <- data.frame(unique(huvdf1$Gene))
huvdf2 <- data.frame(unique(huvdf2$Gene))
colnames(huvecdf1) <- c("Gene")
colnames(huvecdf2) <- c("Gene")
colnames(huvdf1) <- c("Gene")
colnames(huvdf2) <- c("Gene")

#Defining the areas
area1 <- nrow(huvecdf1)
area2 <- nrow(huvecdf2)
area3 <- nrow(huvdf1)
area4 <- nrow(huvdf2)

n12_df <- semi_join(huvecdf1, huvecdf2, by = "Gene")
n12 <- nrow(n12_df)
write.csv(n12_df, "n12.csv")

n13_df <- semi_join(huvecdf1, huvdf1, by = "Gene")
n13 <- nrow(n13_df)
write.csv(n13_df, "n13.csv")

n14_df <- semi_join(huvecdf1, huvdf2, by = "Gene")
n14 <- nrow(n14_df)
write.csv(n14_df, "n14.csv")

n23_df <- semi_join(huvecdf2, huvdf1, by = "Gene")
n23 <- nrow(n23_df)
write.csv(n23_df, "n23.csv")

n24_df <- semi_join(huvecdf2, huvdf2, by = "Gene")
n24 <- nrow(n24_df)
write.csv(n24_df, "n24.csv")

n34_df <- semi_join(huvdf1, huvdf2, by = "Gene")
n34 <- nrow(n34_df)
write.csv(n34_df, "n34.csv")

n123_df <- semi_join(n12_df, huvdf1, by = "Gene")
n123 <- nrow(n123_df)
write.csv(n123_df, "n123.csv")

n124_df <- semi_join(n12_df, huvdf2, by = "Gene")
n124 <- nrow(n124_df)
write.csv(n124_df, "n124.csv")

n134_df <- semi_join(huvecdf1, n34_df, by = "Gene")
n134 <- nrow(n134_df)
write.csv(n134_df, "n134.csv")

n234_df <- semi_join(huvecdf2, n34_df, by = "Gene")
n234 <- nrow(n234_df)
write.csv(n234_df, "n234.csv")

n1234_df <- semi_join(n12_df, n34_df, by = "Gene")
n1234 <- nrow(n1234_df)
write.csv(n1234_df, "n1234.csv")

#Writing the data frames for the exact gene number
true_n1 <- anti_join(huvecdf1, n13_df, by = "Gene")
true_n1 <- anti_join(true_n1, n14_df, by = "Gene")
true_n1 <- anti_join(true_n1, n12_df, by = "Gene")
write.csv(true_n1, "TRUE_N1.csv")

true_n2 <- anti_join(huvecdf2, n23_df, by = "Gene")
true_n2 <- anti_join(true_n2, n24_df, by = "Gene")
true_n2 <- anti_join(true_n2, n12_df, by = "Gene")
write.csv(true_n2, "TRUE_N2.csv")

true_n3 <- anti_join(huvdf1, n13_df, by = "Gene")
true_n3 <- anti_join(true_n3, n23_df, by = "Gene")
true_n3 <- anti_join(true_n3, n34_df, by = "Gene")
write.csv(true_n3, "TRUE_N3.csv")

true_n4 <- anti_join(huvdf2, n14_df, by = "Gene")
true_n4 <- anti_join(true_n4, n24_df, by = "Gene")
true_n4 <- anti_join(true_n4, n34_df, by = "Gene")
write.csv(true_n4, "TRUE_N4.csv")

true_n12 <- anti_join(n12_df, n1234_df, by = "Gene")
true_n12 <- anti_join(true_n12, n124_df, by = "Gene")
true_n12 <- anti_join(true_n12, n123_df, by = "Gene")
write.csv(true_n12, "TRUE_N12.csv")

true_n13 <- anti_join(n13_df, n1234_df, by = "Gene")
true_n13 <- anti_join(true_n13, n134_df, by = "Gene")
true_n13 <- anti_join(true_n13, n123_df, by = "Gene")
write.csv(true_n13, "TRUE_N13.csv")

true_n34 <- anti_join(n34_df, n1234_df, by = "Gene")
true_n34 <- anti_join(true_n34, n134_df, by = "Gene")
true_n34 <- anti_join(true_n34, n234_df, by = "Gene")
write.csv(true_n34, "TRUE_N34.csv")

true_n24 <- anti_join(n24_df, n1234_df, by = "Gene")
true_n24 <- anti_join(true_n24, n234_df, by = "Gene")
true_n24 <- anti_join(true_n24, n124_df, by = "Gene")
write.csv(true_n24, "TRUE_N24.csv")

true_n14 <- anti_join(n14_df, n1234_df, by = "Gene")
true_n14 <- anti_join(true_n14, n134_df, by = "Gene")
true_n14 <- anti_join(true_n14, n124_df, by = "Gene")
write.csv(true_n14, "TRUE_N14.csv")

true_n23 <- anti_join(n23_df, n1234_df, by = "Gene")
true_n23 <- anti_join(true_n23, n123_df, by = "Gene")
true_n23 <- anti_join(true_n23, n234_df, by = "Gene")
write.csv(true_n23, "TRUE_N23.csv")

true_n124 <- anti_join(n124_df, n1234_df, by = "Gene")
write.csv(true_n124, "TRUE_N124.csv")

true_n123 <- anti_join(n123_df, n1234_df, by = "Gene")
write.csv(true_n123, "TRUE_N123.csv")

true_n134 <- anti_join(n134_df, n1234_df, by = "Gene")
write.csv(true_n134, "TRUE_N134.csv")

true_n234 <- anti_join(n234_df, n1234_df, by = "Gene")
write.csv(true_n234, "TRUE_N234.csv")

# 4 way Venn
grid.newpage()
draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234, 
               category = c("HUVEC No Biotin", "HUVEC + Biotin", "HUViPS4F1 No Biotin", "HUViPS4F1 + Biotin"), 
               lwd = rep(2, 4), 
               lty = rep("solid", 4), 
               col = rep("black", 4), 
               fill = c("red", "darkred", "yellow", "goldenrod"))


