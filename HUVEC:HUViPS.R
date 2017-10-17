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


### Bring in the data
# HUVEC/HUViPS comparison
huv <- data.frame(read.csv(paste(file.path, "HUVECHUViPSSummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(huv$X, " ", 5)
huv <- cbind.data.frame(splitcolumns, huv)
names(huv)[names(huv) == '1'] <- 'Line'  #Change the first column to "Line"
names(huv)[names(huv) == '3'] <- 'PassNum'  #Change the second column to "PassNum"
names(huv)[names(huv) == '4'] <- 'Biotin'  #Change the first column to "Line"
huv <- subset(huv, select=-c(`2`, `5`, X))
huv <- filter(huv, Species == "HUMAN") #Filtering only human hits
splitcolumns <- colsplit(huv$Gene.Names," ",c("Gene","Synonyms"))
huv <- cbind.data.frame(splitcolumns, huv)
huv <- subset(huv, select=-Gene.Names)
huvslim <- data.frame(huv$Line, huv$Gene, huv$Synonyms, huv$Biotin, huv$Unused, huv$Total, huv$Gene.Ontology)
colnames(huvslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")

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

na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}
huvec$UnusedHUVECBiotin <- na.zero(huvec$UnusedHUVECBiotin)
huvec$UnusedHUVECNoBiotin <- na.zero(huvec$UnusedHUVECNoBiotin)
huvec$TotalHUVECBiotin <- na.zero(huvec$TotalHUVECBiotin)
huvec$TotalHUVECNoBiotin <- na.zero(huvec$TotalHUVECNoBiotin)

hdf1 <- data.frame(huvecnobiotin$Gene, huvecnobiotin$Line)
colnames(hdf1) <- c("Gene", "Line")
hdf2 <- data.frame(huvecbiotin$Gene, huvecbiotin$Line)
colnames(hdf2) <- c("Gene", "Line")

hdf3 <- setdiff(hdf1, hdf2)
no1 <- nrow(hdf1)
hdf4 <- setdiff(hdf2, hdf1)
no2 <- nrow(hdf2)
hdf5 <- semi_join(hdf1, hdf2)
no5 <- nrow(hdf5)

grid.newpage()
draw.pairwise.venn(no1, no2, no5, category = c("HUVEC No Biotin", "HUVEC + Biotin"), lty = rep("solid", 2), 
                  lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("yellow", "dark green"), alpha = rep(0.5, 1), 
                  cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

huvecbiotinunique <- anti_join(hdf2, hdf1)
write.csv(huvecbiotinunique, "UniqueHUVECBiotin.csv")

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

hdf6 <- data.frame(huvipsnobiotin$Gene, huvipsnobiotin$Line)
colnames(hdf6) <- c("Gene", "Line")
hdf7 <- data.frame(huvipsbiotin$Gene, huvipsbiotin$Line)
colnames(hdf7) <- c("Gene", "Line")

hdf8 <- setdiff(hdf6, hdf7)
no6 <- nrow(hdf6)
hdf9 <- setdiff(hdf7, hdf6)
no7 <- nrow(hdf7)
hdf10 <- semi_join(hdf6, hdf7)
no10 <- nrow(hdf10)

grid.newpage()
draw.pairwise.venn(no6, no7, no10, category = c("HUViPS4F5 No Biotin", "HUViPS4F5 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("yellow", "dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)



huvipsbiotinunique <- anti_join(hdf7, hdf6)
write.csv(huvipsbiotinunique, "UniqueHUViPSBiotin.csv")


huvectotal <- filter(huvec, TotalHUVECNoBiotin < TotalHUVECBiotin)
huvipstotal <- filter(huvips, TotalHUViPSNoBiotin < TotalHUViPSBiotin)

hdf11 <- data.frame(huvipstotal$Gene)
colnames(hdf11) <- c("Gene")
hdf12 <- data.frame(huvectotal$Gene)
colnames(hdf12) <- c("Gene")

hdf13 <- setdiff(hdf11, hdf12)
no11 <- nrow(hdf11)
hdf14 <- setdiff(hdf12, hdf11)
no12 <- nrow(hdf12)
hdf15 <- semi_join(hdf11, hdf12)
no15 <- nrow(hdf15)

grid.newpage()
draw.pairwise.venn(no11, no12, no15, category = c("HUViPS4F5", "HUVEC"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("yellow", "dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

huvecunique <- anti_join(hdf12, hdf11)
huvipsunique <- anti_join(hdf11, hdf12)
huvechuvips <- intersect(hdf11, hdf12)
write.csv(huvecunique, "UniqueHUVEC.csv")
write.csv(huvipsunique, "UniqueHUViPS.csv")
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



