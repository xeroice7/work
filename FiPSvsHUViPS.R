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
huvslim <- data.frame(huv$Line, huv$Gene, huv$Synonyms, huv$Biotin, huv$Unused, huv$Total, huv$Gene.Ontology, huv$Entry)
colnames(huvslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology", "ProteinID")

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

huvdf1 <- data.frame(huvipsnobiotin$Gene, huvipsnobiotin$Line, huvipsnobiotin$ProteinID)
colnames(huvdf1) <- c("Gene", "Line", "ProteinID")
huvno1 <- nrow(huvdf1)
huvdf2 <- data.frame(huvipsbiotin$Gene, huvipsbiotin$Line, huvipsbiotin$ProteinID)
colnames(huvdf2) <- c("Gene", "Line", "ProteinID")
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

#FiPS
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
fipsslim <- data.frame(fips$Line, fips$Gene, fips$Synonyms, fips$Biotin, fips$Unused, fips$Total, fips$Gene.Ontology, fips$Entry)
colnames(fipsslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology", "ProteinID")

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

fipsdf1 <- data.frame(fipsnobiotin$Gene, fipsnobiotin$Line, fipsnobiotin$ProteinID)
colnames(fipsdf1) <- c("Gene", "Line", "ProteinID")
fipsno1 <- nrow(fipsdf1)
fipsdf2 <- data.frame(fipsbiotin$Gene, fipsbiotin$Line, fipsbiotin$ProteinID)
colnames(fipsdf2) <- c("Gene", "Line", "ProteinID")
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


# To compare together
fipstotal <- filter(fips, TotalFiPSNoBiotin < TotalFiPSBiotin)
huvipstotal <- filter(huvips, TotalHUViPSNoBiotin < TotalHUViPSBiotin)

huvipsvsfips <- merge(fipstotal, huvipstotal, by = "Gene", all = TRUE)

huvipsvsfips$UnusedFiPSBiotin <- na.zero(huvipsvsfips$UnusedFiPSBiotin)
huvipsvsfips$UnusedFiPSNoBiotin <- na.zero(huvipsvsfips$UnusedFiPSNoBiotin)
huvipsvsfips$TotalFiPSBiotin <- na.zero(huvipsvsfips$TotalFiPSBiotin)
huvipsvsfips$TotalFiPSNoBiotin <- na.zero(huvipsvsfips$TotalFiPSNoBiotin)
huvipsvsfips$UnusedHUViPSBiotin <- na.zero(huvipsvsfips$UnusedHUViPSBiotin)
huvipsvsfips$UnusedHUViPSNoBiotin <- na.zero(huvipsvsfips$UnusedHUViPSNoBiotin)
huvipsvsfips$TotalHUViPSBiotin <- na.zero(huvipsvsfips$TotalHUViPSBiotin)
huvipsvsfips$TotalHUViPSNoBiotin <- na.zero(huvipsvsfips$TotalHUViPSNoBiotin)

huvdf6 <- data.frame(huvipstotal$Gene)
colnames(huvdf6) <- c("Gene")
huvno6 <- nrow(huvdf6)
fipsdf6 <- data.frame(fipstotal$Gene)
colnames(fipsdf6) <- c("Gene")
fipsno6 <- nrow(fipsdf6)

fipshuvdf1 <- setdiff(huvdf6, fipsdf6)
fipshuvdf2 <- setdiff(fipsdf6, huvdf6)
fipshuvdf3 <- semi_join(huvdf6, fipsdf6)
fipshuvno3 <- nrow(fipshuvdf3)

grid.newpage()
draw.pairwise.venn(huvno6, fipsno6, fipshuvno3, category = c("HUViPS4F1", "FiPS4F5"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark green", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

huvunique <- anti_join(huvdf6, fipsdf6)
fipsunique <- anti_join(fipsdf6, huvdf6)
huvfips <- intersect(huvdf6$Gene, fipsdf6$Gene)
write.csv(fipsunique, "UniqueFiPStoHUVIPS.csv")
write.csv(huvunique, "UniqueHUViPStoFiPS.csv")
write.csv(huvfips, "HUViPSFiPSOverlap.csv")
write.csv(fipstotal$Gene, "FiPSReal.csv")
write.csv(huvipstotal$Gene, "HUViPSReal.csv")


#Compare what is unique to just both iPS from their somatic sources
fips1 <- data.frame(read.csv(paste(file.path, "UniqueFiPStoIMR.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
huvips1 <- data.frame(read.csv(paste(file.path, "UniqueHUViPStoHUVEC.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
iPS1 <- data.frame(fips1$Gene)
colnames(iPS1) <- c("Gene")
iPS2 <- data.frame(huvips1$Gene)
colnames(iPS2) <- c("Gene")
iPS1$Gene <- as.character(iPS1$Gene)
iPS2$Gene <- as.character(iPS2$Gene)

df16 <- anti_join(iPS1, iPS2)
no16 <- nrow(iPS1)
df17 <- anti_join(iPS2, iPS1)
no17 <- nrow(iPS2)
df18 <- semi_join(iPS1, iPS2)
no18 <- nrow(df18)


grid.newpage()
draw.pairwise.venn(no16, no17, no18, category = c("FiPS4F5", "HUViPS4F1"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("pink", "dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

write.csv(df16, "TotalUniqueFiPS.csv")
write.csv(df17, "TotalUniqueHUViPS.csv")
write.csv(df18, "TotalUniqueiPSOverlapGENE.csv")

########
### Comparing both unique iPS hits to MDA

uniqueipsmdadf1 <- setdiff(mdadf6, df18)
uniqueipsmdadf2 <- setdiff(df18, mdadf6)
uniqueipsmdadf3 <- semi_join(mdadf6, df18)
ipsmdano3 <- nrow(uniqueipsmdadf3)

grid.newpage()
draw.pairwise.venn(mdano6, no18, ipsmdano3, category = c("MDA-MB-231", "Unique iPS"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("light green", "dark gray"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mdatoipsunique <- anti_join(mdadf6, df18)
ipstomdaunique <- anti_join(df18, mdadf6)
mdauniqueips <- intersect(mdadf6$Gene, df18$Gene)
write.csv(mdatoipsunique, "UniqueMDAtoiPS.csv")
write.csv(ipstomdaunique, "UniqueiPStoMDA.csv")
write.csv(mdauniqueips, "iPSMDAOverlap.csv")


### Comparing both unique iPS hits to MCF7

uniqueipsmcf7df1 <- setdiff(mcf7df6, df18)
uniqueipsmcf7df2 <- setdiff(df18, mcf7df6)
uniqueipsmcf7df3 <- semi_join(mcf7df6, df18)
ipsmcf7no3 <- nrow(uniqueipsmcf7df3)

grid.newpage()
draw.pairwise.venn(mcf7no6, no18, ipsmcf7no3, category = c("MCF7", "Unique iPS"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("firebrick1", "dark gray"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mcf7toipsunique <- anti_join(mcf7df6, df18)
ipstomcf7unique <- anti_join(df18, mcf7df6)
mcf7uniqueips <- intersect(mcf7df6$Gene, df18$Gene)
write.csv(mcf7toipsunique, "UniqueMCF7toiPS.csv")
write.csv(ipstomcf7unique, "UniqueiPStoMCF7.csv")
write.csv(mcf7uniqueips, "iPSMCF7Overlap.csv")

### Comparing both unique iPS hits to A549

uniqueipsa549df1 <- setdiff(a549df6, df18)
uniqueipsa549df2 <- setdiff(df18, a549df6)
uniqueipsa549df3 <- semi_join(a549df6, df18)
ipsa549no3 <- nrow(uniqueipsa549df3)

grid.newpage()
draw.pairwise.venn(a549no6, no18, ipsa549no3, category = c("A549", "Unique iPS"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark blue", "dark gray"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

a549toipsunique <- anti_join(a549df6, df18)
ipstoa549unique <- anti_join(df18, a549df6)
a549uniqueips <- intersect(a549df6$Gene, df18$Gene)
write.csv(a549toipsunique, "UniqueA549toiPS.csv")
write.csv(ipstoa549unique, "UniqueiPStoA549.csv")
write.csv(a549uniqueips, "iPSA549Overlap.csv")

### Comparing both unique iPS hits to NCI-H1299

uniqueipsncidf1 <- setdiff(ncidf6, df18)
uniqueipsncidf2 <- setdiff(df18, ncidf6)
uniqueipsncidf3 <- semi_join(ncidf6, df18)
ipsncino3 <- nrow(uniqueipsncidf3)

grid.newpage()
draw.pairwise.venn(ncino6, no18, ipsncino3, category = c("NCI-H1299", "Unique iPS"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("orange", "dark gray"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

ncitoipsunique <- anti_join(ncidf6, df18)
ipstonciunique <- anti_join(df18, ncidf6)
nciuniqueips <- intersect(ncidf6$Gene, df18$Gene)
write.csv(ncitoipsunique, "UniqueNCItoiPS.csv")
write.csv(ipstonciunique, "UniqueiPStoNCI.csv")
write.csv(nciuniqueips, "iPSNCIOverlap.csv")

### Comparisons of all cell lines:

z <- semi_join(ncidf6, df18)
z <- semi_join(z, mcf7df6)
z <- semi_join(z, a549df6)
z <- semi_join(z, mdadf6)

##################
##################

fipsup <- filter(imrvsfips, TotalFiPSBiotin > TotalIMRBiotin)
imrup <- filter(imrvsfips, TotalFiPSBiotin < TotalIMRBiotin)
write.csv(fipsup$Gene, "FiPSUp.csv")
write.csv(imrup$Gene, "IMRUp.csv")

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

fipsunique <- anti_join(hdf12, hdf11)
huvipsunique <- anti_join(hdf11, hdf12)
huvechuvips <- intersect(hdf11, hdf12)

write.csv(huvipsunique, "UniqueHUViPS.csv")
write.csv(huvechuvips, "HUViPSHUVECOverlap.csv")

write.csv(huvipstotal$Gene, "HUViPSReal.csv")

#To compare which is bigger TOTAL
huvecvshuvips <- merge.data.frame(huvipstotal, huvectotal, by = "Gene", all = TRUE)

huvecvshuvips$UnusedHUViPSBiotin <- na.zero(huvecvshuvips$UnusedHUViPSBiotin)
huvecvshuvips$UnusedHUViPSNoBiotin <- na.zero(huvecvshuvips$UnusedHUViPSNoBiotin)
huvecvshuvips$TotalHUViPSBiotin <- na.zero(huvecvshuvips$TotalHUViPSBiotin)
huvecvshuvips$TotalHUViPSNoBiotin <- na.zero(huvecvshuvips$TotalHUViPSNoBiotin)

huvipsup <- filter(huvecvshuvips, TotalHUViPSBiotin > TotalHUVECBiotin)

write.csv(huvipsup$Gene, "HUViPSUp.csv")



############
############
############

imrtotal <- filter(imr, TotalIMRNoBiotin < TotalIMRBiotin)
fipstotal <- filter(fips, TotalFiPSNoBiotin < TotalFiPSBiotin)

df11 <- data.frame(fipstotal$Gene)
colnames(df11) <- c("Gene")
df12 <- data.frame(imrtotal$Gene)
colnames(df12) <- c("Gene")

df13 <- setdiff(df11, df12)
no11 <- nrow(df11)
df14 <- setdiff(df12, df11)
no12 <- nrow(df12)
df15 <- semi_join(df11, df12)
no15 <- nrow(df15)

grid.newpage()
draw.pairwise.venn(no11, no12, no15, category = c("FiPS4F5", "IMR90"), lty = rep("blank", 2), 
                   fill = c("blue", "red"), alpha = rep(0.5, 2), cat.pos = c(0, 0),
                   cat.dist = rep(0.025, 2), scaled = TRUE)

imrunique <- anti_join(df12, df11)
fipsunique <- anti_join(df11, df12)
imrfips <- intersect(df11, df12)
write.csv(imrunique, "UniqueIMR.csv")
write.csv(fipsunique, "UniqueFiPS.csv")
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

########
# 4 way Venn
########

#Getting rid of duplicates
fipsdf1 <- data.frame(unique(fipsdf1$Gene))
fipsdf2 <- data.frame(unique(fipsdf2$Gene))
huvdf1 <- data.frame(unique(huvdf1$Gene))
huvdf2 <- data.frame(unique(huvdf2$Gene))
colnames(fipsdf1) <- c("Gene")
colnames(fipsdf2) <- c("Gene")
colnames(huvdf1) <- c("Gene")
colnames(huvdf2) <- c("Gene")

#Defining the areas
area1 <- nrow(fipsdf1)
area2 <- nrow(fipsdf2)
area3 <- nrow(huvdf1)
area4 <- nrow(huvdf2)

n12_df <- semi_join(fipsdf1, fipsdf2, by = "Gene")
n12 <- nrow(n12_df)
write.csv(n12_df, "n12.csv")

n13_df <- semi_join(fipsdf1, huvdf1, by = "Gene")
n13 <- nrow(n13_df)
write.csv(n13_df, "n13.csv")

n14_df <- semi_join(fipsdf1, huvdf2, by = "Gene")
n14 <- nrow(n14_df)
write.csv(n14_df, "n14.csv")

n23_df <- semi_join(fipsdf2, huvdf1, by = "Gene")
n23 <- nrow(n23_df)
write.csv(n23_df, "n23.csv")

n24_df <- semi_join(fipsdf2, huvdf2, by = "Gene")
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

n134_df <- semi_join(fipsdf1, n34_df, by = "Gene")
n134 <- nrow(n134_df)
write.csv(n134_df, "n134.csv")

n234_df <- semi_join(fipsdf2, n34_df, by = "Gene")
n234 <- nrow(n234_df)
write.csv(n234_df, "n234.csv")

n1234_df <- semi_join(n12_df, n34_df, by = "Gene")
n1234 <- nrow(n1234_df)
write.csv(n1234_df, "n1234.csv")

#Writing the data frames for the exact gene number
true_n1 <- anti_join(fipsdf1, n13_df, by = "Gene")
true_n1 <- anti_join(true_n1, n14_df, by = "Gene")
true_n1 <- anti_join(true_n1, n12_df, by = "Gene")
write.csv(true_n1, "TRUE_N1.csv")

true_n2 <- anti_join(fipsdf2, n23_df, by = "Gene")
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
               category = c("FiPS4F5 No Biotin", "FiPS4F5 + Biotin", "HUViPS4F1 No Biotin", "HUViPS4F1 + Biotin"), 
               lwd = rep(2, 4), 
               lty = rep("solid", 4), 
               col = rep("black", 4), 
               fill = c("lightgreen", "green", "red", "firebrick"))
