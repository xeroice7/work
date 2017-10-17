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
install.packages("Vennerable", repos="http://R-Forge.R-project.org") 
install.packages("reshape")
install.packages("graph")
install.packages("RBGL")
install.packages("gtools")
install.packages("xtable")
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
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("graph")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("RBGL")

# Define Function for converting NAs to 0s when comparing value scores
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

#HUViPS
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
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

#MCF7
mcf7 <- data.frame(read.csv(paste(file.path, "MCF7SummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(mcf7$X, "/", 2)
mcf7 <- cbind.data.frame(splitcolumns, mcf7)
names(mcf7)[names(mcf7) == '1'] <- 'Line'  #Change the first column to "Line"
names(mcf7)[names(mcf7) == '2'] <- 'YN'  #Change the second column to "YN"
splitcolumns <- str_split_fixed(mcf7$Line, " ", 3)
mcf7 <- cbind.data.frame(splitcolumns, mcf7)
names(mcf7)[names(mcf7) == '1'] <- 'Cell'  #Change the first column to "Cell"
names(mcf7)[names(mcf7) == '3'] <- 'Dox'  #Change the first column to "Line"
splitcolumns <- str_split_fixed(mcf7$YN, " ", 2)
splitcolumns <- splitcolumns[,1]
mcf7 <- cbind.data.frame(splitcolumns, mcf7)
names(mcf7)[names(mcf7) == "splitcolumns"] <- "Biotin"
mcf7 <- subset(mcf7, select=-c(`2`, Line, YN, X))
mcf7 <- filter(mcf7, Species == "HUMAN") #Filtering only human hits
mcf7 <- filter(mcf7, Dox == "No Dox") #Filtering only No Dox hits
splitcolumns <- str_split_fixed(mcf7$Gene.Names," ", 2)
mcf7 <- cbind.data.frame(splitcolumns, mcf7)
names(mcf7)[names(mcf7) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(mcf7)[names(mcf7) == '2'] <- 'Synonyms'  #Change the first column to "Line"
mcf7 <- subset(mcf7, select=-Gene.Names)
mcf7slim <- data.frame(mcf7$Cell, mcf7$Gene, mcf7$Synonyms, mcf7$Biotin, mcf7$Unused, mcf7$Total, mcf7$Gene.Ontology)
colnames(mcf7slim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

mcf7biotin <- filter(mcf7slim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(mcf7biotin)[colnames(mcf7biotin) == 'Unused'] <- 'UnusedMCF7Biotin'
colnames(mcf7biotin)[colnames(mcf7biotin) == 'Total'] <- 'TotalMCF7Biotin'

mcf7nobiotin <- filter(mcf7slim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(mcf7nobiotin)[colnames(mcf7nobiotin) == 'Unused'] <- 'UnusedMCF7NoBiotin'
colnames(mcf7nobiotin)[colnames(mcf7nobiotin) == 'Total'] <- 'TotalMCF7NoBiotin'

mcf7 <- merge.data.frame(mcf7biotin, mcf7nobiotin, by = "Gene", all = TRUE)
mcf7$UnusedMCF7Biotin <- as.numeric(as.character(mcf7$UnusedMCF7Biotin))
mcf7$UnusedMCF7NoBiotin <- as.numeric(as.character(mcf7$UnusedMCF7NoBiotin))
mcf7$TotalMCF7Biotin <- as.numeric(as.character(mcf7$TotalMCF7Biotin))
mcf7$TotalMCF7NoBiotin <- as.numeric(as.character(mcf7$TotalMCF7NoBiotin))

mcf7$UnusedMCF7Biotin <- na.zero(mcf7$UnusedMCF7Biotin)
mcf7$UnusedMCF7NoBiotin <- na.zero(mcf7$UnusedMCF7NoBiotin)
mcf7$TotalMCF7Biotin <- na.zero(mcf7$TotalMCF7Biotin)
mcf7$TotalMCF7NoBiotin <- na.zero(mcf7$TotalMCF7NoBiotin)

mcf7df1 <- data.frame(mcf7nobiotin$Gene, mcf7nobiotin$Line)
colnames(mcf7df1) <- c("Gene", "Line")
mcf7no1 <- nrow(mcf7df1)
mcf7df2 <- data.frame(mcf7biotin$Gene, mcf7biotin$Line)
colnames(mcf7df2) <- c("Gene", "Line")
mcf7no2 <- nrow(mcf7df2)

mcf7df3 <- setdiff(mcf7df1, mcf7df2)
mcf7df4 <- setdiff(mcf7df2, mcf7df2)
mcf7df5 <- semi_join(mcf7df1, mcf7df2)
mcf7no5 <- nrow(mcf7df5)

grid.newpage()
draw.pairwise.venn(mcf7no1, mcf7no2, mcf7no5, category = c("MCF7 No Biotin", "MCF7 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("deepskyblue", "firebrick1"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mcf7biotinunique <- anti_join(mcf7df2, mcf7df1)
write.csv(mcf7biotinunique, "UniqueMCF7Biotin.csv")


# To compare together
mcf7total <- filter(mcf7, TotalMCF7NoBiotin < TotalMCF7Biotin)
huvipstotal <- filter(huvips, TotalHUViPSNoBiotin < TotalHUViPSBiotin)

huvipsvsmcf7 <- merge.data.frame(mcf7total, huvipstotal, by = "Gene", all = TRUE)

huvipsvsmcf7$UnusedMCF7Biotin <- na.zero(huvipsvsmcf7$UnusedMCF7Biotin)
huvipsvsmcf7$UnusedMCF7NoBiotin <- na.zero(huvipsvsmcf7$UnusedMCF7NoBiotin)
huvipsvsmcf7$TotalMCF7Biotin <- na.zero(huvipsvsmcf7$TotalMCF7Biotin)
huvipsvsmcf7$TotalMCF7NoBiotin <- na.zero(huvipsvsmcf7$TotalMCF7NoBiotin)
huvipsvsmcf7$UnusedHUViPSBiotin <- na.zero(huvipsvsmcf7$UnusedHUViPSBiotin)
huvipsvsmcf7$UnusedHUViPSNoBiotin <- na.zero(huvipsvsmcf7$UnusedHUViPSNoBiotin)
huvipsvsmcf7$TotalHUViPSBiotin <- na.zero(huvipsvsmcf7$TotalHUViPSBiotin)
huvipsvsmcf7$TotalHUViPSNoBiotin <- na.zero(huvipsvsmcf7$TotalHUViPSNoBiotin)

mcf7df6 <- data.frame(unique(mcf7total$Gene))
colnames(mcf7df6) <- c("Gene")
mcf7no6 <- nrow(mcf7df6)
huvdf6 <- data.frame(unique(huvipstotal$Gene))
colnames(huvdf6) <- c("Gene")
huvno6 <- nrow(huvdf6)

huvmcf7df1 <- setdiff(mcf7df6, huvdf6)
huvmcf7df2 <- setdiff(huvdf6, mcf7df6)
huvmcf7df3 <- semi_join(mcf7df6, huvdf6)
huvmcf7no3 <- nrow(unique(huvmcf7df3))

grid.newpage()
draw.pairwise.venn(mcf7no6, huvno6, huvmcf7no3, category = c("MCF7", "HUViPS4F1"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("firebrick","dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mcf7uniquetohuvips <- anti_join(mcf7df6, huvdf6)
huvipsuniquetomcf7 <- anti_join(huvdf6, mcf7df6)
mcf7huvips <- intersect(mcf7df6$Gene, huvdf6$Gene)
write.csv(mcf7uniquetohuvips, "UniqueMCF7toHUVIPS.csv")
write.csv(huvipsuniquetomcf7, "UniqueHUViPStoMCF7.csv")
write.csv(mcf7huvips, "HUViPSMCF7Overlap.csv")
write.csv(mcf7total$Gene, "MCF7Real.csv")
write.csv(huvipstotal$Gene, "HUViPSReal.csv")

#######
# 4 way Venn Diagram
#######

#Getting rid of duplicates
mcf7df1 <- data.frame(unique(mcf7df1$Gene))
mcf7df2 <- data.frame(unique(mcf7df2$Gene))
huvdf1 <- data.frame(unique(huvdf1$Gene))
huvdf2 <- data.frame(unique(huvdf2$Gene))
colnames(mcf7df1) <- c("Gene")
colnames(mcf7df2) <- c("Gene")
colnames(huvdf1) <- c("Gene")
colnames(huvdf2) <- c("Gene")

#Defining the areas
area1 <- nrow(mcf7df1)
area2 <- nrow(mcf7df2)
area3 <- nrow(huvdf1)
area4 <- nrow(huvdf2)

n12_df <- semi_join(mcf7df1, mcf7df2, by = "Gene")
n12 <- nrow(n12_df)
write.csv(n12_df, "n12.csv")

n13_df <- semi_join(mcf7df1, huvdf1, by = "Gene")
n13 <- nrow(n13_df)
write.csv(n13_df, "n13.csv")

n14_df <- semi_join(mcf7df1, huvdf2, by = "Gene")
n14 <- nrow(n14_df)
write.csv(n14_df, "n14.csv")

n23_df <- semi_join(mcf7df2, huvdf1, by = "Gene")
n23 <- nrow(n23_df)
write.csv(n23_df, "n23.csv")

n24_df <- semi_join(mcf7df2, huvdf2, by = "Gene")
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

n134_df <- semi_join(mcf7df1, n34_df, by = "Gene")
n134 <- nrow(n134_df)
write.csv(n134_df, "n134.csv")

n234_df <- semi_join(mcf7df2, n34_df, by = "Gene")
n234 <- nrow(n234_df)
write.csv(n234_df, "n234.csv")

n1234_df <- semi_join(n12_df, n34_df, by = "Gene")
n1234 <- nrow(n1234_df)
write.csv(n1234_df, "n1234.csv")

#Writing the data frames for the exact gene number
true_n1 <- anti_join(mcf7df1, n13_df, by = "Gene")
true_n1 <- anti_join(true_n1, n14_df, by = "Gene")
true_n1 <- anti_join(true_n1, n12_df, by = "Gene")
write.csv(true_n1, "TRUE_N1.csv")

true_n2 <- anti_join(mcf7df2, n23_df, by = "Gene")
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
               category = c("MCF7 No Biotin", "MCF7 + Biotin", "HUViPS4F1 No Biotin", "HUViPS4F1 + Biotin"), 
               lwd = rep(2, 4), 
               lty = rep("solid", 4), 
               col = rep("black", 4), 
               fill = c("lightyellow", "darkgoldenrod1", "red", "firebrick"))
