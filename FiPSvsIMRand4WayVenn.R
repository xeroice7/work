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

# Functions needed
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

####### FiPS/IMR comparison

# Tidying Data
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
fips <- data.frame(read.csv(paste(file.path, "IMRFiPSSummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(fips$X, " ", 4)
fips <- cbind.data.frame(splitcolumns, fips)
names(fips)[names(fips) == '1'] <- 'Line'  #Change the first column to "Line"
names(fips)[names(fips) == '2'] <- 'PassNum'  #Change the second column to "PassNum"
names(fips)[names(fips) == '3'] <- 'Biotin'  #Change the first column to "Line"
fips <- subset(fips, select=-c(`4`, X))
fips <- filter(fips, Species == "HUMAN") #Filtering only human hits
splitcolumns <- str_split_fixed(fips$Gene.Names," ", 2)
fips <- cbind.data.frame(splitcolumns, fips)
names(fips)[names(fips) == '1'] <- 'Gene'  #Change the first column to "Line"
names(fips)[names(fips) == '2'] <- 'Synonyms'  #Change the second column to "PassNum"
fips <- subset(fips, select=-Gene.Names)
fipsslim <- data.frame(fips$Line, fips$Gene, fips$Synonyms, fips$Biotin, fips$Unused, fips$Total, fips$Gene.Ontology)
colnames(fipsslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

# Analysis of IMR
imrslim <- filter(fipsslim, Line == "IMR90")  #Keep only IMR
imrbiotin <- filter(imrslim, Biotin == "+") #Keep only +Biotin for IMR
colnames(imrbiotin)[colnames(imrbiotin) == 'Unused'] <- 'UnusedIMRBiotin'
colnames(imrbiotin)[colnames(imrbiotin) == 'Total'] <- 'TotalIMRBiotin'

imrnobiotin <- filter(imrslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(imrnobiotin)[colnames(imrnobiotin) == 'Unused'] <- 'UnusedIMRNoBiotin'
colnames(imrnobiotin)[colnames(imrnobiotin) == 'Total'] <- 'TotalIMRNoBiotin'

imr <- merge.data.frame(imrbiotin, imrnobiotin, by = "Gene", all = TRUE)
imr$UnusedIMRBiotin <- as.numeric(as.character(imr$UnusedIMRBiotin))
imr$UnusedIMRNoBiotin <- as.numeric(as.character(imr$UnusedIMRNoBiotin))
imr$TotalIMRBiotin <- as.numeric(as.character(imr$TotalIMRBiotin))
imr$TotalIMRNoBiotin <- as.numeric(as.character(imr$TotalIMRNoBiotin))

imr$UnusedIMRBiotin <- na.zero(imr$UnusedIMRBiotin)
imr$TotalIMRBiotin <- na.zero(imr$TotalIMRBiotin)
imr$UnusedIMRNoBiotin <- na.zero(imr$UnusedIMRNoBiotin)
imr$TotalIMRNoBiotin <- na.zero(imr$TotalIMRNoBiotin)

imrdf1 <- data.frame(imrnobiotin$Gene, imrnobiotin$Line)
colnames(imrdf1) <- c("Gene", "Line")
imrno1 <- nrow(imrdf1)
imrdf2 <- data.frame(imrbiotin$Gene, imrbiotin$Line)
colnames(imrdf2) <- c("Gene", "Line")
imrno2 <- nrow(imrdf2)

imrdf3 <- setdiff(imrdf1, imrdf2)
imrdf4 <- setdiff(imrdf2, imrdf1)
imrdf5 <- semi_join(imrdf1, imrdf2)
imrno5 <- nrow(imrdf5)

grid.newpage()
draw.pairwise.venn(imrno1, imrno2, imrno5, category = c("IMR90 No Biotin", "IMR90 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("red", "blue"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

imrbiotinunique <- anti_join(imrdf2, imrdf1)
write.csv(imrbiotinunique, "UniqueIMRBiotin.csv")

# Analysis of FiPS
fipsslim <- filter(fipsslim, Line == "FiPS4F5")  #Keep only FiPS4F5
fipsbiotin <- filter(fipsslim, Biotin == "+") #Keep only +Biotin for FiPS4F5
colnames(fipsbiotin)[colnames(fipsbiotin) == 'Unused'] <- 'UnusedFiPSBiotin'
colnames(fipsbiotin)[colnames(fipsbiotin) == 'Total'] <- 'TotalFiPSBiotin'

fipsnobiotin <- filter(fipsslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(fipsnobiotin)[colnames(fipsnobiotin) == 'Unused'] <- 'UnusedFiPSNoBiotin'
colnames(fipsnobiotin)[colnames(fipsnobiotin) == 'Total'] <- 'TotalFiPSNoBiotin'

fips <- merge.data.frame(fipsbiotin, fipsnobiotin, by = "Gene", all = TRUE)
fips$UnusedFiPSBiotin <- as.numeric(as.character(fips$UnusedFiPSBiotin))
fips$UnusedFiPSNoBiotin <- as.numeric(as.character(fips$UnusedFiPSNoBiotin))
fips$TotalFiPSBiotin <- as.numeric(as.character(fips$TotalFiPSBiotin))
fips$TotalFiPSNoBiotin <- as.numeric(as.character(fips$TotalFiPSNoBiotin))

fips$UnusedFiPSBiotin <- na.zero(fips$UnusedFiPSBiotin)
fips$UnusedFiPSNoBiotin <- na.zero(fips$UnusedFiPSNoBiotin)
fips$TotalFiPSBiotin <- na.zero(fips$TotalFiPSBiotin)
fips$TotalFiPSNoBiotin <- na.zero(fips$TotalFiPSNoBiotin)

fipsdf1 <- data.frame(fipsnobiotin$Gene, fipsnobiotin$Line)
colnames(fipsdf1) <- c("Gene", "Line")
fipsno1 <- nrow(fipsdf1)
fipsdf2 <- data.frame(fipsbiotin$Gene, fipsbiotin$Line)
colnames(fipsdf2) <- c("Gene", "Line")
fipsno2 <- nrow(fipsdf2)

fipsdf3 <- setdiff(fipsdf1, fipsdf2)
fipsdf4 <- setdiff(fipsdf2, fipsdf1)
fipsdf5 <- semi_join(fipsdf1, fipsdf2)
fipsno5 <- nrow(fipsdf5)

grid.newpage()
draw.pairwise.venn(fipsno1, fipsno2, fipsno5, category = c("FiPS4F5 No Biotin", "FiPS4F5 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("purple", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

fipsbiotinunique <- anti_join(fipsdf2, fipsdf1)
write.csv(fipsbiotinunique, "UniqueFiPSBiotin.csv")
write.csv(fipsdf5, "FiPSBiotinNoBiotinOverlap.csv")


# Combined Analysis 
imrtotal <- filter(imr, TotalIMRNoBiotin < TotalIMRBiotin)
fipstotal <- filter(fips, TotalFiPSNoBiotin < TotalFiPSBiotin)

imrdf6 <- data.frame(imrtotal$Gene)
colnames(imrdf6) <- c("Gene")
imrno6 <- nrow(imrdf6)
fipsdf6 <- data.frame(fipstotal$Gene)
colnames(fipsdf6) <- c("Gene")
fipsno6 <- nrow(fipsdf6)

fipsimrdf1 <- setdiff(imrdf6, fipsdf6)
fipsimrdf2 <- setdiff(fipsdf6, imrdf6)
fipsimrdf3 <- semi_join(imrdf6, fipsdf6)
fipsimrno3 <- nrow(fipsimrdf3)

grid.newpage()
draw.pairwise.venn(imrno6, fipsno6, fipsimrno3, category = c("IMR90", "FiPS4F5"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("blue", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

imrunique <- anti_join(imrdf6, fipsdf6)
fipsunique <- anti_join(fipsdf6, imrdf6)
imrfips <- intersect(imrdf6$Gene, fipsdf6$Gene)
write.csv(imrunique, "UniqueIMRtoFiPS.csv")
write.csv(fipsunique, "UniqueFiPStoIMR.csv")
write.csv(imrfips, "FiPSIMROverlap.csv")
write.csv(imrtotal$Gene, "IMRReal.csv")
write.csv(fipstotal$Gene, "FiPSReal.csv")

#To compare which is bigger unused
imrvsfips <- merge.data.frame(fipstotal, imrtotal, by = "Gene", all = TRUE)

imrvsfips$UnusedFiPSBiotin <- na.zero(imrvsfips$UnusedFiPSBiotin)
imrvsfips$UnusedFiPSNoBiotin <- na.zero(imrvsfips$UnusedFiPSNoBiotin)
imrvsfips$TotalFiPSBiotin <- na.zero(imrvsfips$TotalFiPSBiotin)
imrvsfips$TotalFiPSNoBiotin <- na.zero(imrvsfips$TotalFiPSNoBiotin)
imrvsfips$UnusedIMRBiotin <- na.zero(imrvsfips$UnusedIMRBiotin)
imrvsfips$UnusedIMRNoBiotin <- na.zero(imrvsfips$UnusedIMRNoBiotin)
imrvsfips$TotalIMRBiotin <- na.zero(imrvsfips$TotalIMRBiotin)
imrvsfips$TotalIMRNoBiotin <- na.zero(imrvsfips$TotalIMRNoBiotin)

fipsup <- filter(imrvsfips, TotalFiPSBiotin > TotalIMRBiotin)
imrup <- filter(imrvsfips, TotalFiPSBiotin < TotalIMRBiotin)
write.csv(fipsup$Gene, "FiPSUp.csv")
write.csv(imrup$Gene, "IMRUp.csv")


#Getting rid of duplicates
imrdf1 <- data.frame(unique(imrdf1$Gene))
imrdf2 <- data.frame(unique(imrdf2$Gene))
fipsdf1 <- data.frame(unique(fipsdf1$Gene))
fipsdf2 <- data.frame(unique(fipsdf2$Gene))
colnames(imrdf1) <- c("Gene")
colnames(imrdf2) <- c("Gene")
colnames(fipsdf1) <- c("Gene")
colnames(fipsdf2) <- c("Gene")

#Defining the areas
area1 <- nrow(imrdf1)
area2 <- nrow(imrdf2)
area3 <- nrow(fipsdf1)
area4 <- nrow(fipsdf2)

n12_df <- semi_join(imrdf1, imrdf2, by = "Gene")
n12 <- nrow(n12_df)
write.csv(n12_df, "n12.csv")

n13_df <- semi_join(imrdf1, fipsdf1, by = "Gene")
n13 <- nrow(n13_df)
write.csv(n13_df, "n13.csv")

n14_df <- semi_join(imrdf1, fipsdf2, by = "Gene")
n14 <- nrow(n14_df)
write.csv(n14_df, "n14.csv")

n23_df <- semi_join(imrdf2, fipsdf1, by = "Gene")
n23 <- nrow(n23_df)
write.csv(n23_df, "n23.csv")

n24_df <- semi_join(imrdf2, fipsdf2, by = "Gene")
n24 <- nrow(n24_df)
write.csv(n24_df, "n24.csv")

n34_df <- semi_join(fipsdf1, fipsdf2, by = "Gene")
n34 <- nrow(n34_df)
write.csv(n34_df, "n34.csv")

n123_df <- semi_join(n12_df, fipsdf1, by = "Gene")
n123 <- nrow(n123_df)
write.csv(n123_df, "n123.csv")

n124_df <- semi_join(n12_df, fipsdf2, by = "Gene")
n124 <- nrow(n124_df)
write.csv(n124_df, "n124.csv")

n134_df <- semi_join(imrdf1, n34_df, by = "Gene")
n134 <- nrow(n134_df)
write.csv(n134_df, "n134.csv")

n234_df <- semi_join(imrdf2, n34_df, by = "Gene")
n234 <- nrow(n234_df)
write.csv(n234_df, "n234.csv")

n1234_df <- semi_join(n12_df, n34_df, by = "Gene")
n1234 <- nrow(n1234_df)
write.csv(n1234_df, "n1234.csv")

#Writing the data frames for the exact gene number
true_n1 <- anti_join(imrdf1, n13_df, by = "Gene")
true_n1 <- anti_join(true_n1, n14_df, by = "Gene")
true_n1 <- anti_join(true_n1, n12_df, by = "Gene")
write.csv(true_n1, "TRUE_N1.csv")

true_n2 <- anti_join(imrdf2, n23_df, by = "Gene")
true_n2 <- anti_join(true_n2, n24_df, by = "Gene")
true_n2 <- anti_join(true_n2, n12_df, by = "Gene")
write.csv(true_n2, "TRUE_N2.csv")

true_n3 <- anti_join(fipsdf1, n13_df, by = "Gene")
true_n3 <- anti_join(true_n3, n23_df, by = "Gene")
true_n3 <- anti_join(true_n3, n34_df, by = "Gene")
write.csv(true_n3, "TRUE_N3.csv")

true_n4 <- anti_join(fipsdf2, n14_df, by = "Gene")
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
               category = c("IMR No Biotin", "IMR + Biotin", "FiPS4F5 No Biotin", "FiPS4F5 + Biotin"), 
               lwd = rep(2, 4), 
               lty = rep("solid", 4), 
               col = rep("black", 4), 
               fill = c("lightgreen", "green", "lightblue", "blue"))
