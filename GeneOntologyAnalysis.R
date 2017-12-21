ipsoverlap <- read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/HUViPSvsFiPS/HUViPSFiPSOverlap.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
ipsoverlap <- ipsoverlap$x

a549huvips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/A549vsHUViPS/HUViPSA549OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
a549fips_overlap <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/MS Analysis/iPS Comparisons/A549vsFiPS/FiPSA549OverlapGENE.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE))
a549ips_overlap <- intersect(a549fips_overlap$V1, a549huvips_overlap$V1)

huv <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/02-19-16/HUVECHUViPSSummaryFDRContaminants.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
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
huv <- subset(huv, Gene %in% a549ips_overlap)
huv <- filter(huv, Line == "HUViPS4F1")  #Keep only HUViPS4F1
huv <- filter(huv, Biotin == "+") #Keep only +Biotin for HUViPS4F1
huv <- subset(huv, select=c(Gene, Gene.Ontology))


splitGO <- stri_split_fixed(huv$Gene.Ontology, ";", simplify = TRUE) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA
huvGO <- as.data.frame(table(unlist(splitGO)))
huvGO <- arrange(huvGO, desc(Freq))
huvGO$Var1 <- factor(huvGO$Var1, levels = huvGO$Var1[order(huvGO$Freq)])
huvGO <- huvGO[1:50,]
#huvGO <- filter(huvGO, huvGO$Freq >= 25)

ggplot(huvGO) +
  geom_bar(aes(x=Var1, y=Freq), stat = 'identity') + 
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60)) + 
  theme_minimal() + 
  ylab("Frequency of GO Term") + 
  xlab("GO Term") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25, size = 10), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 10)) + coord_flip()

