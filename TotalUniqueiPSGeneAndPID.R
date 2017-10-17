library(stringi)
library(splitstackshape)

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
hits <- data.frame(read.csv(paste(file.path, "HUVECHUViPSSummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
hits <- filter(hits, Species == "HUMAN") #Filtering only human hits

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
pid<- data.frame(read.csv(paste(file.path, "TotalUniqueiPSPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

colnames(pid) <- c("No", "Entry")

joined <- left_join(pid, hits)

split <- str_split_fixed(joined$Gene.Names," ", 2)
joined <- cbind.data.frame(split, joined)
names(joined)[names(joined) == '1'] <- 'Gene'  #Change the first column to "Line"
names(joined)[names(joined) == '2'] <- 'Synonyms'  #Change the second column to "PassNum"
joined <- subset(joined, select=-Gene.Names)

joined1 <- data.frame(joined$Entry, joined$Gene)
joined1 <- unique(joined1)

colnames(joined1) <- c("PID", "Gene")
write.csv(joined1, "TotalUniqueiPSGENEandPID.csv")

joined2 <- data.frame(joined$Gene, joined$Entry, joined$Gene.Ontology)
colnames(joined2) <- c("Gene", "Entry", "GO")
splitGO <- data.frame(stri_split_fixed(joined2$GO, ";", simplify = TRUE)) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA

joinedGO <- as.data.frame(table(unlist(splitGO)))
joinedGO <- filter(joinedGO, joinedGO$Freq >= 5)

ggplot(joinedGO) + 
  geom_bar(aes(x=Var1, y=Freq), stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12))

#Membrane specific hits
membraneGO <- filter(joined2, grepl('cell surface|extracellular|plasma membrane', GO))
membraneGO <- unique(membraneGO)

splitGO <- data.frame(stri_split_fixed(membraneGO$GO, ";", simplify = TRUE)) #splits column into a variable amount of columns
splitGO[splitGO == ""] = NA # changes all "" to NA

surface_targets <- cbind.data.frame(membraneGO, splitGO)
surface_targets <- subset(surface_targets, select = -GO)
counts <- as.data.frame(table(unlist(surface_targets)))
counts <- filter(counts, counts$Freq >= 5)

ggplot(counts) + 
  geom_bar(aes(x=Var1, y=Freq), stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12))

write.csv(membraneGO, "MembraneGO.csv")
