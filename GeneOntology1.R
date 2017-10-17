install.packages("stringi")
install.packages("splitstackshape")
install.packages("tidyverse")
install.packages("VennDiagram")
library(stringi)
library(splitstackshape)
library(tidyverse)
library(VennDiagram)

# Define Function for converting NAs to 0s when comparing value scores
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

######
# HUViPS
######
# Bring in data
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
names(huv)[names(huv) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(huv)[names(huv) == '2'] <- 'Synonyms'  #Change the first column to "Line"
huv <- subset(huv, select=-Gene.Names)
huvslim <- data.frame(huv$Line, huv$Gene, huv$Synonyms, huv$Biotin, huv$Unused, huv$Total, huv$Gene.Ontology, huv$Entry)
colnames(huvslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology", "Entry")

# Separate Biotin from No Biotin
huvipsslim <- filter(huvslim, Line == "HUViPS4F1")  #Keep only HUViPS4F1
huvipsbiotin <- filter(huvipsslim, Biotin == "+") #Keep only +Biotin for HUViPS4F1
colnames(huvipsbiotin)[colnames(huvipsbiotin) == 'Unused'] <- 'UnusedHUViPSBiotin'
colnames(huvipsbiotin)[colnames(huvipsbiotin) == 'Total'] <- 'TotalHUViPSBiotin'

huvipsnobiotin <- filter(huvipsslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(huvipsnobiotin)[colnames(huvipsnobiotin) == 'Unused'] <- 'UnusedHUViPSNoBiotin'
colnames(huvipsnobiotin)[colnames(huvipsnobiotin) == 'Total'] <- 'TotalHUViPSNoBiotin'

# Compare Biotin to No Biotin
huvdf1 <- data.frame(huvipsnobiotin$Gene, huvipsnobiotin$`Gene Ontology`, huvipsnobiotin$Entry)
colnames(huvdf1) <- c("Gene", "GO", "ProteinID")
huvno1 <- nrow(huvdf1)
huvdf2 <- data.frame(huvipsbiotin$Gene, huvipsbiotin$`Gene Ontology`, huvipsbiotin$Entry)
colnames(huvdf2) <- c("Gene", "GO", "ProteinID")
huvno2 <- nrow(huvdf2)

huvdf3 <- setdiff(huvdf1, huvdf2)
huvdf4 <- setdiff(huvdf2, huvdf1)

huvdf5 <- semi_join(huvdf1, huvdf2)
huvno5 <- nrow(huvdf5)

# Plot Venn Diagram
grid.newpage()
draw.pairwise.venn(huvno1, huvno2, huvno5, category = c("HUViPS4F5 No Biotin", "HUViPS4F5 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("yellow", "dark green"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

# Gene Ontology for HUViPS
huvipsbiotinunique <- anti_join(huvdf2, huvdf1)
splitGO <- stri_split_fixed(huvipsbiotinunique$GO, ";", simplify = TRUE) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA
huvBioGO <- as.data.frame(table(unlist(splitGO)))
huvBioGO <- filter(huvBioGO, huvBioGO$Freq >= 5)
write.csv(huvipsbiotinunique, "UniqueHUViPSBiotin.csv")
write.csv(huvBioGO, "UniqueHUViPSBiotinGO.csv")
huvBioPID <- as.data.frame(huvipsbiotinunique$ProteinID)
write.csv(huvBioPID, "UniqueHUViPSBiotinPIDs.csv")

huvipsnobiotinunique <- anti_join(huvdf1, huvdf2)
splitGO <- stri_split_fixed(huvipsnobiotinunique$GO, ";", simplify = TRUE) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA
huvNoBioGO <- as.data.frame(table(unlist(splitGO)))
write.csv(huvipsnobiotinunique, "UniqueHUViPSNOBiotin.csv")
write.csv(huvNoBioGO, "UniqueHUViPSNOBiotinGO.csv")

huvipscommon <- semi_join(huvdf1, huvdf2)
splitGO <- stri_split_fixed(huvipsbiotinunique$GO, ";", simplify = TRUE) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA
huvComGO <- as.data.frame(table(unlist(splitGO)))
write.csv(huvipscommon, "CommonHUViPS.csv")
write.csv(huvBioGO, "CommonHUViPSGO.csv")

ggplot(huvBioGO) + 
  geom_bar(aes(x=Var1, y=Freq), stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))

# Solidifying HUViPS 
huvips <- merge.data.frame(huvipsbiotin, huvipsnobiotin, by = "Gene", all = TRUE)
huvips$UnusedHUViPSBiotin <- as.numeric(as.character(huvips$UnusedHUViPSBiotin))
huvips$UnusedHUViPSNoBiotin <- as.numeric(as.character(huvips$UnusedHUViPSNoBiotin))
huvips$TotalHUViPSBiotin <- as.numeric(as.character(huvips$TotalHUViPSBiotin))
huvips$TotalHUViPSNoBiotin <- as.numeric(as.character(huvips$TotalHUViPSNoBiotin))

huvips$UnusedHUViPSBiotin <- na.zero(huvips$UnusedHUViPSBiotin)
huvips$UnusedHUViPSNoBiotin <- na.zero(huvips$UnusedHUViPSNoBiotin)
huvips$TotalHUViPSBiotin <- na.zero(huvips$TotalHUViPSBiotin)
huvips$TotalHUViPSNoBiotin <- na.zero(huvips$TotalHUViPSNoBiotin)

huvipstotal <- filter(huvips, TotalHUViPSNoBiotin < TotalHUViPSBiotin)
huvipstotal <- unite(huvipstotal, Entry, c(Entry.x, Entry.y), sep='')
huvipstotal$Entry <- substr(huvipstotal$Entry, 0, 6)
write.csv(c(huvipstotal$Gene, huvipstotal$Entry), "HUViPSRealSomatics.csv")

######
# FiPS
######
# Bring in data
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
colnames(fipsslim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology", "Entry")

# Separate Biotin from No Biotin
fipsslim <- filter(fipsslim, Line == "FiPS4F5")  #Keep only FiPS4F5
fipsbiotin <- filter(fipsslim, Biotin == "+") #Keep only +Biotin for FiPS4F5
colnames(fipsbiotin)[colnames(fipsbiotin) == 'Unused'] <- 'UnusedFiPSBiotin'
colnames(fipsbiotin)[colnames(fipsbiotin) == 'Total'] <- 'TotalFiPSBiotin'

fipsnobiotin <- filter(fipsslim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(fipsnobiotin)[colnames(fipsnobiotin) == 'Unused'] <- 'UnusedFiPSNoBiotin'
colnames(fipsnobiotin)[colnames(fipsnobiotin) == 'Total'] <- 'TotalFiPSNoBiotin'

# Compare Biotin to No Biotin
fipsdf1 <- data.frame(fipsnobiotin$Gene, fipsnobiotin$`Gene Ontology`, fipsnobiotin$Entry)
colnames(fipsdf1) <- c("Gene", "GO", "ProteinID")
fipsno1 <- nrow(fipsdf1)
fipsdf2 <- data.frame(fipsbiotin$Gene, fipsbiotin$`Gene Ontology`, fipsbiotin$Entry)
colnames(fipsdf2) <- c("Gene", "GO", "ProteinID")
fipsno2 <- nrow(fipsdf2)

fipsdf3 <- setdiff(fipsdf1, fipsdf2)
fipsdf4 <- setdiff(fipsdf2, fipsdf1)
fipsdf5 <- semi_join(fipsdf1, fipsdf2)
fipsno5 <- nrow(fipsdf5)

# Plot Venn Diagram
grid.newpage()
draw.pairwise.venn(fipsno1, fipsno2, fipsno5, category = c("FiPS4F5 No Biotin", "FiPS4F5 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("purple", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

# Gene Ontology for FiPS
fipsbiotinunique <- anti_join(fipsdf2, fipsdf1)
splitGO <- stri_split_fixed(fipsbiotinunique$GO, ";", simplify = TRUE) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA
fipsBioGO <- as.data.frame(table(unlist(splitGO)))
fipsBioGO <- filter(fipsBioGO, fipsBioGO$Freq >= 5)
write.csv(fipsbiotinunique, "UniqueFiPSBiotin.csv")
write.csv(fipsBioGO, "UniqueFiPSBiotinGO.csv")
fipsBioPID <- as.data.frame(fipsbiotinunique$ProteinID)
write.csv(fipsBioPID, "UniqueFiPSBiotinPIDs.csv")

fipsnobiotinunique <- anti_join(fipsdf1, fipsdf2)
splitGO <- stri_split_fixed(fipsnobiotinunique$GO, ";", simplify = TRUE) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA
fipsNoBioGO <- as.data.frame(table(unlist(splitGO)))
write.csv(fipsnobiotinunique, "UniqueFiPSNOBiotin.csv")
write.csv(fipsNoBioGO, "UniqueFiPSNOBiotinGO.csv")

fipscommon <- semi_join(fipsdf1, fipsdf2)
splitGO <- stri_split_fixed(fipsbiotinunique$GO, ";", simplify = TRUE) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA
fipsComGO <- as.data.frame(table(unlist(splitGO)))
write.csv(fipscommon, "CommonFiPS.csv")
write.csv(fipsBioGO, "CommonFiPSGO.csv")

ggplot(fipsBioGO) + 
  geom_bar(aes(x=Var1, y=Freq), stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))

# Solidifying FiPS
fips <- merge.data.frame(fipsbiotin, fipsnobiotin, by = "Gene", all = TRUE)
fips$UnusedFiPSBiotin <- as.numeric(as.character(fips$UnusedFiPSBiotin))
fips$UnusedFiPSNoBiotin <- as.numeric(as.character(fips$UnusedFiPSNoBiotin))
fips$TotalFiPSBiotin <- as.numeric(as.character(fips$TotalFiPSBiotin))
fips$TotalFiPSNoBiotin <- as.numeric(as.character(fips$TotalFiPSNoBiotin))

fips$UnusedFiPSBiotin <- na.zero(fips$UnusedFiPSBiotin)
fips$UnusedFiPSNoBiotin <- na.zero(fips$UnusedFiPSNoBiotin)
fips$TotalFiPSBiotin <- na.zero(fips$TotalFiPSBiotin)
fips$TotalFiPSNoBiotin <- na.zero(fips$TotalFiPSNoBiotin)

fipstotal <- filter(fips, TotalFiPSNoBiotin < TotalFiPSBiotin)
fipstotal <- unite(fipstotal, Entry, c(Entry.x, Entry.y), sep='')
fipstotal$Entry <- substr(fipstotal$Entry, 0, 6)
write.csv(c(fipstotal$Gene, fipstotal$Entry), "FiPSRealSomatics.csv")

######
# MCF7
######
# Bring in data
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
mcf7slim <- data.frame(mcf7$Cell, mcf7$Gene, mcf7$Synonyms, mcf7$Biotin, mcf7$Unused, mcf7$Total, mcf7$Gene.Ontology, mcf7$Entry)
colnames(mcf7slim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology", "Entry")

# Separate Biotin from No Biotin
mcf7biotin <- filter(mcf7slim, Biotin == "+") #Keep only +Biotin for MCF7
colnames(mcf7biotin)[colnames(mcf7biotin) == 'Unused'] <- 'UnusedMCF7Biotin'
colnames(mcf7biotin)[colnames(mcf7biotin) == 'Total'] <- 'TotalMCF7Biotin'

mcf7nobiotin <- filter(mcf7slim, Biotin == "No") #Keep only No Biotin for FiPS
colnames(mcf7nobiotin)[colnames(mcf7nobiotin) == 'Unused'] <- 'UnusedMCF7NoBiotin'
colnames(mcf7nobiotin)[colnames(mcf7nobiotin) == 'Total'] <- 'TotalMCF7NoBiotin'

# Compare Biotin to No Biotin
mcf7df1 <- data.frame(mcf7nobiotin$Gene, mcf7nobiotin$`Gene Ontology`, mcf7nobiotin$Entry)
colnames(mcf7df1) <- c("Gene", "GO", "ProteinID")
mcf7no1 <- nrow(mcf7df1)
mcf7df2 <- data.frame(mcf7biotin$Gene, mcf7biotin$`Gene Ontology`, mcf7biotin$Entry)
colnames(mcf7df2) <- c("Gene", "GO", "ProteinID")
mcf7no2 <- nrow(mcf7df2)

mcf7df3 <- setdiff(mcf7df1, mcf7df2)
mcf7df4 <- setdiff(mcf7df2, mcf7df1)
mcf7df5 <- semi_join(mcf7df1, mcf7df2)
mcf7no5 <- nrow(mcf7df5)

# Plot Venn Diagram
grid.newpage()
draw.pairwise.venn(mcf7no1, mcf7no2, mcf7no5, category = c("MCF7 No Biotin", "MCF7 + Biotin"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("deepskyblue", "firebrick1"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

# Gene Ontology for MCF7
mcf7biotinunique <- anti_join(mcf7df2, mcf7df1)
splitGO <- stri_split_fixed(mcf7biotinunique$GO, ";", simplify = TRUE) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA
mcf7BioGO <- as.data.frame(table(unlist(splitGO)))
mcf7BioGO <- filter(mcf7BioGO, mcf7BioGO$Freq >= 5)
write.csv(mcf7biotinunique, "UniqueMCF7Biotin.csv")
write.csv(mcf7BioGO, "UniqueMCF7BiotinGO.csv")
mcf7BioPID <- as.data.frame(mcf7biotinunique$ProteinID)
write.csv(mcf7BioPID, "UniqueMCF7BiotinPIDs.csv")

mcf7nobiotinunique <- anti_join(mcf7df1, mcf7df2)
splitGO <- stri_split_fixed(mcf7nobiotinunique$GO, ";", simplify = TRUE) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA
mcf7NoBioGO <- as.data.frame(table(unlist(splitGO)))
write.csv(mcf7nobiotinunique, "UniqueMCF7NOBiotin.csv")
write.csv(mcf7NoBioGO, "UniqueMCF7NOBiotinGO.csv")

mcf7common <- semi_join(mcf7df1, mcf7df2)
splitGO <- stri_split_fixed(mcf7biotinunique$GO, ";", simplify = TRUE) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA
mcf7ComGO <- as.data.frame(table(unlist(splitGO)))
write.csv(mcf7common, "CommonMCF7.csv")
write.csv(mcf7BioGO, "CommonMCF7GO.csv")

ggplot(mcf7BioGO) + 
  geom_bar(aes(x=Var1, y=Freq), stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))

# Solidifying MCF7
mcf7df <- merge.data.frame(mcf7biotin, mcf7nobiotin, by = "Gene", all = TRUE)
mcf7df$UnusedMCF7Biotin <- as.numeric(as.character(mcf7df$UnusedMCF7Biotin))
mcf7df$UnusedMCF7NoBiotin <- as.numeric(as.character(mcf7df$UnusedMCF7NoBiotin))
mcf7df$TotalMCF7Biotin <- as.numeric(as.character(mcf7df$TotalMCF7Biotin))
mcf7df$TotalMCF7NoBiotin <- as.numeric(as.character(mcf7df$TotalMCF7NoBiotin))

mcf7df$UnusedMCF7Biotin <- na.zero(mcf7df$UnusedMCF7Biotin)
mcf7df$UnusedMCF7NoBiotin <- na.zero(mcf7df$UnusedMCF7NoBiotin)
mcf7df$TotalMCF7Biotin <- na.zero(mcf7df$TotalMCF7Biotin)
mcf7df$TotalMCF7NoBiotin <- na.zero(mcf7df$TotalMCF7NoBiotin)

mcf7totaldf <- filter(mcf7df, TotalMCF7NoBiotin < TotalMCF7Biotin)
mcf7totaldf <- unite(mcf7totaldf, Entry, c(Entry.x, Entry.y), sep='')
mcf7totaldf$Entry <- substr(mcf7totaldf$Entry, 0, 6)
write.csv(c(mcf7totaldf$Gene, mcf7totaldf$Entry), "MCF7Real.csv")

######
# Combining all 3 cell lines together
######

huvdf6 <- data.frame(huvipstotal$Gene, huvipstotal$Entry)
colnames(huvdf6) <- c("Gene", "ProteinID")
huvno6 <- nrow(huvdf6)
fipsdf6 <- data.frame(fipstotal$Gene, fipstotal$Entry)
colnames(fipsdf6) <- c("Gene", "ProteinID")
fipsno6 <- nrow(fipsdf6)
mcf7df6 <- data.frame(mcf7totaldf$Gene, mcf7totaldf$Entry)
colnames(mcf7df6) <- c("Gene", "ProteinID")
mcf7no6 <- nrow(mcf7df6)

fipshuvdf1 <- anti_join(huvdf6, fipsdf6)
fipshuvdf2 <- anti_join(fipsdf6, huvdf6)
fipshuvdf3 <- semi_join(huvdf6, fipsdf6)
fipshuvno3 <- nrow(fipshuvdf3)

fipsmcf7df1 <- anti_join(mcf7df6, fipsdf6)
fipsmcf7df2 <- anti_join(fipsdf6, mcf7df6)
fipsmcf7df3 <- semi_join(mcf7df6, fipsdf6)
fipsmcf7no3 <- nrow(unique(fipsmcf7df3))

mcf7huvdf1 <- anti_join(huvdf6, mcf7df6)
mcf7huvdf2 <- anti_join(mcf7df6, huvdf6)
mcf7huvdf3 <- semi_join(huvdf6, mcf7df6)
mcf7huvno3 <- nrow(mcf7huvdf3)

# Combined Venn Diagrams
grid.newpage()
draw.pairwise.venn(huvno6, fipsno6, fipshuvno3, category = c("HUViPS4F1", "FiPS4F5"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark green", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

grid.newpage()
draw.pairwise.venn(mcf7no6, fipsno6, fipsmcf7no3, category = c("MCF7", "FiPS4F5"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark green", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

grid.newpage()
draw.pairwise.venn(huvno6, mcf7no6, mcf7huvno3, category = c("HUViPS4F1", "MCF7"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark green", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

# HUViPS4F1 / FiPS comparisons (can contain somatic proteins)
huvuniquetofipssomatic <- anti_join(huvdf6, fipsdf6)
write.csv(huvuniquetofipssomatic, "UniqueHUViPStoFiPSSomatics.csv")
huvipsPIDtofipssomatic <- as.data.frame(huvuniquetofipssomatic$ProteinID)
write.csv(huvipsPIDtofipssomatic, "UniqueHUViPStoFiPSPIDSomatics.csv")

fipsuniquetohuvipssomatic <- anti_join(fipsdf6, huvdf6)
write.csv(fipsuniquetohuvipssomatic, "UniqueFiPStoHUVIPSSomatics.csv")
fipsPIDtohuvipssomatic <- as.data.frame(fipsuniquetohuvipssomatic$ProteinID)
write.csv(fipsPIDtohuvipssomatic, "UniqueFiPStoHUViPSPIDSomatics.csv")

huvfipsGENEsomatic <- intersect(huvdf6$Gene, fipsdf6$Gene)
huvfipsPIDsomatic <- intersect(huvdf6$ProteinID, fipsdf6$ProteinID)
write.csv(huvfipsGENEsomatic, "HUViPSFiPSOverlapGENESomatics.csv")
write.csv(huvfipsPIDsomatic, "HUViPSFiPSOverlapPIDSomatics.csv")

# HUViPS4F1 / MCF7 comparisons (can contain somatic proteins)
huvuniquetomcf7somatic <- anti_join(huvdf6, mcf7df6)
write.csv(huvuniquetomcf7somatic, "UniqueHUViPStoMCF7Somatics.csv")
huvipsPIDtomcf7somatic <- as.data.frame(huvuniquetomcf7somatic$ProteinID)
write.csv(huvipsPIDtomcf7somatic, "UniqueHUViPStoMCF7PIDSomatics.csv")

mcf7uniquetohuvipssomatic <- anti_join(mcf7df6, huvdf6)
write.csv(mcf7uniquetohuvipssomatic, "UniqueMCF7toHUViPSSomatics.csv")
mcf7PIDhuvipssomatic <- as.data.frame(mcf7uniquetohuvipssomatic$ProteinID)
write.csv(mcf7PIDhuvipssomatic, "UniqueMCF7toHUViPSPIDSomatics.csv")

huvmcf7GENEsomatic <- intersect(huvdf6$Gene, mcf7df6$Gene)
huvmcf7PIDsomatic <- intersect(huvdf6$ProteinID, mcf7df6$ProteinID)
write.csv(huvmcf7GENEsomatic, "HUViPSMCF7OverlapGENESomatics.csv")
write.csv(huvmcf7PIDsomatic, "HUViPSMCF7OverlapPIDSomatics.csv")

# FiPS4F5 / MCF7 comparisons (can contain somatic proteins)
fipsuniquetomcf7somatic <- anti_join(fipsdf6, mcf7df6)
write.csv(fipsuniquetomcf7somatic, "UniqueFiPStoMCF7Somatics.csv")
fipsPIDtomcf7somatic <- as.data.frame(fipsuniquetomcf7somatic$ProteinID)
write.csv(fipsPIDtomcf7somatic, "UniqueFiPStoMCF7PIDSomatics.csv")

mcf7uniquetofipssomatic <- anti_join(mcf7df6, fipsdf6)
write.csv(mcf7uniquetofipssomatic, "UniqueMCF7toFiPSSomatics.csv")
mcf7PIDfipssomatic <- as.data.frame(mcf7uniquetofipssomatic$ProteinID)
write.csv(mcf7PIDfipssomatic, "UniqueMCF7toFiPSPIDSomatics.csv")

fipsmcf7GENEsomatic <- intersect(fipsdf6$Gene, mcf7df6$Gene)
fipsmcf7PIDsomatic <- intersect(fipsdf6$ProteinID, mcf7df6$ProteinID)
write.csv(fipsmcf7GENEsomatic, "FiPSMCF7OverlapGENESomatics.csv")
write.csv(fipsmcf7PIDsomatic, "FiPSMCF7OverlapPIDSomatics.csv")

# Bring in data that is iPS unique
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
write.csv(df18, "TotalUniqueiPSOverlap.csv")

fips2 <- data.frame(read.csv(paste(file.path, "IMRFiPSSummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(fips2$X, " ", 4)
fips2 <- cbind.data.frame(splitcolumns, fips2)
names(fips2)[names(fips2) == '1'] <- 'Line'  #Change the first column to "Line"
names(fips2)[names(fips2) == '2'] <- 'PassNum'  #Change the second column to "PassNum"
names(fips2)[names(fips2) == '3'] <- 'Biotin'  #Change the first column to "Line"
splitcolumns <- str_split_fixed(fips2$Gene.Names," ", 2)
fips2 <- cbind.data.frame(splitcolumns, fips2)
names(fips2)[names(fips2) == '1'] <- 'Gene'  #Change the first column to "Line"
names(fips2)[names(fips2) == '2'] <- 'Synonyms'  #Change the second column to "PassNum"
fips2 <- filter(fips2, Species == "HUMAN")
fips2 <- filter(fips2, Line == "FiPS4F5")
fips2 <- filter(fips2, Biotin == "+")
fips3 <- filter(fips2, Gene %in% df16$Gene)
write.csv(fips3$Entry, "TotalUniqueFiPSPID.csv")

huv2 <- data.frame(read.csv(paste(file.path, "HUVECHUViPSSummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(huv2$X, " ", 5)
huv2 <- cbind.data.frame(splitcolumns, huv2)
names(huv2)[names(huv2) == '1'] <- 'Line'  #Change the first column to "Line"
names(huv2)[names(huv2) == '2'] <- 'P'
names(huv2)[names(huv2) == '3'] <- 'PassNum'  #Change the second column to "PassNum"
names(huv2)[names(huv2) == '4'] <- 'Biotin'  #Change the first column to "Line"
splitcolumns <- str_split_fixed(huv2$Gene.Names," ", 2)
huv2 <- cbind.data.frame(splitcolumns, huv2)
names(huv2)[names(huv2) == '1'] <- 'Gene'  #Change the first column to "Cell"
huv2 <- filter(huv2, Species == "HUMAN")
huv2 <- filter(huv2, Line == "HUViPS4F1")
huv2 <- filter(huv2, Biotin == "+")
huv3 <- filter(huv2, Gene %in% df17$Gene)
write.csv(huv3$Entry, "TotalUniqueHUViPSPID.csv")

huv4 <- filter(huv2, Gene %in% df18$Gene)
write.csv(huv4$Entry, "TotalUniqueiPSPID.csv")

########
# To compare together
mcf7total <- filter(mcf7, TotalMCF7NoBiotin < TotalMCF7Biotin)
fipstotal <- filter(fips, TotalFiPSNoBiotin < TotalFiPSBiotin)

mcf7df6 <- data.frame(unique(mcf7total$Gene))
colnames(mcf7df6) <- c("Gene")
mcf7no6 <- nrow(mcf7df6)
fipsdf6 <- data.frame(unique(fipstotal$Gene))
colnames(fipsdf6) <- c("Gene")
fipsno6 <- nrow(fipsdf6)

grid.newpage()
draw.pairwise.venn(mcf7no6, fipsno6, fipsmcf7no3, category = c("MCF7", "FiPS4F5"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("firebrick1", "pink"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)

mcf7uniquetofips <- anti_join(mcf7df6, fipsdf6)
fipsuniquetomcf7 <- anti_join(fipsdf6, mcf7df6)
mcf7fips <- intersect(mcf7df6$Gene, fipsdf6$Gene)
write.csv(mcf7uniquetofips, "UniqueMCF7toFIPS.csv")
write.csv(fipsuniquetomcf7, "UniqueFiPStoMCF7.csv")
write.csv(mcf7fips, "FiPSMCF7Overlap.csv")
write.csv(mcf7total$Gene, "MCF7Real.csv")
write.csv(fipstotal$Gene, "FiPSReal.csv")

### Comparing both unique iPS hits to MCF7

uniqueipsmcf7df1 <- setdiff(mcf7df6, df18)
uniqueipsmcf7df2 <- setdiff(df18, mcf7df6)
uniqueipsmcf7df3 <- semi_join(mcf7df6, df18)
ipsmcf7no3 <- nrow(uniqueipsmcf7df3)

grid.newpage()
draw.pairwise.venn(mcf7no6, no18, ipsmcf7no3, category = c("MCF7", "Unique iPS"), lty = rep("solid", 2), 
                   lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("firebrick1", "dark gray"), alpha = rep(0.5, 1), 
                   cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
colnames(mcf7df6) <- c("Gene")
mcf7toipsunique <- anti_join(mcf7df6, df18)
ipstomcf7unique <- anti_join(df18, mcf7df6)
mcf7uniqueips <- intersect(mcf7df6$Gene, df18$Gene)
write.csv(mcf7toipsunique, "UniqueMCF7toiPS.csv")
write.csv(ipstomcf7unique, "UniqueiPStoMCF7.csv")
write.csv(mcf7uniqueips, "iPSMCF7Overlap.csv")

fips4 <- filter(fips2, Gene %in% mcf7uniqueips)
write.csv(fips4$Entry, "iPSMCF7OverlapPID.csv")

mcf7toipsunique <- as.character(mcf7toipsunique$Gene)
mcf7.1 <- filter(mcf7.2, Gene %in% mcf7toipsunique)
write.csv(mcf7.1$Entry, "UniqueMCF7toiPSPID.csv")

ipstomcf7unique <- as.character(ipstomcf7unique$Gene)
fips6 <- filter(fips2, Gene %in% ipstomcf7unique)
write.csv(fips6$Entry, "UniqueiPStoMCF7PID.csv")

# MCF7 to HUViPS4F1
mcf7uniquetohuvips <- anti_join(mcf7df6, huvdf6)
huvipsuniquetomcf7 <- anti_join(huvdf6, mcf7df6)
mcf7huvips <- intersect(mcf7df6$Gene, huvdf6$Gene)

mcf7.2 <- data.frame(read.csv(paste(file.path, "MCF7SummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
splitcolumns <- str_split_fixed(mcf7.2$X, "/", 2)
mcf7.2 <- cbind.data.frame(splitcolumns, mcf7.2)
names(mcf7.2)[names(mcf7.2) == '1'] <- 'Line'  #Change the first column to "Line"
names(mcf7.2)[names(mcf7.2) == '2'] <- 'YN'  #Change the second column to "YN"
splitcolumns <- str_split_fixed(mcf7.2$Line, " ", 3)
mcf7.2 <- cbind.data.frame(splitcolumns, mcf7.2)
names(mcf7.2)[names(mcf7.2) == '1'] <- 'Cell'  #Change the first column to "Cell"
names(mcf7.2)[names(mcf7.2) == '3'] <- 'Dox'  #Change the first column to "Line"
splitcolumns <- str_split_fixed(mcf7.2$YN, " ", 2)
splitcolumns <- splitcolumns[,1]
mcf7.2 <- cbind.data.frame(splitcolumns, mcf7.2)
names(mcf7.2)[names(mcf7.2) == "splitcolumns"] <- "Biotin"
mcf7.2 <- subset(mcf7.2, select=-c(`2`, Line, YN, X))
splitcolumns <- str_split_fixed(mcf7.2$Gene.Names," ", 2)
mcf7.2 <- cbind.data.frame(splitcolumns, mcf7.2)
names(mcf7.2)[names(mcf7.2) == '1'] <- 'Gene'  #Change the first column to "Cell"
names(mcf7.2)[names(mcf7.2) == '2'] <- 'Synonyms'  #Change the first column to "Line"
mcf7.2 <- filter(mcf7.2, Species == "HUMAN") #Filtering only human hits
mcf7.2 <- filter(mcf7.2, Dox == "No Dox") #Filtering only No Dox hits
mcf7.2 <- filter(mcf7.2, Biotin == "+")

#Overlap MCF7/HUViPS
mcf7.3 <- filter(mcf7.2, Gene %in% mcf7huvips)
mcf7.3 <- subset(mcf7.3, !duplicated(mcf7.3$Entry))
write.csv(mcf7.3$Entry, "MCF7HUViPSOverlapPID.csv")

#HUViPS Specific
mcf7.4 <- filter(mcf7.2, Gene %in% huvipsuniquetomcf7$Gene)
mcf7.4 <- subset(mcf7.4, !duplicated(mcf7.4$Entry))
write.csv(mcf7.4$Entry, "UniqueHUViPStoMCF7PID.csv")

#MCF7 Specific
mcf7.5 <- filter(mcf7.2, Gene %in% mcf7uniquetohuvips$Gene)
mcf7.5 <- subset(mcf7.5, !duplicated(mcf7.5$Entry))
write.csv(mcf7.5$Entry, "UniqueMCF7toHUViPSPID.csv")

#Overlap MCF7/FiPS
mcf7.6 <- filter(mcf7.2, Gene %in% mcf7fips)
mcf7.6 <- subset(mcf7.6, !duplicated(mcf7.6$Entry))
write.csv(mcf7.6$Entry, "MCF7FiPSOverlapPID.csv")

#FiPS Specific
mcf7.7 <- filter(mcf7.2, Gene %in% fipsuniquetomcf7$Gene)
mcf7.7 <- subset(mcf7.7, !duplicated(mcf7.7$Entry))
write.csv(mcf7.7$Entry, "UniqueFiPStoMCF7PID.csv")

#MCF7 Specific
mcf7.8 <- filter(mcf7.2, Gene %in% mcf7uniquetofips$Gene)
mcf7.8 <- subset(mcf7.8, !duplicated(mcf7.8$Entry))
write.csv(mcf7.8$Entry, "UniqueMCF7toFiPSPID.csv")

