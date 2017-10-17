install.packages('readxl')
install.packages('reshape2')
install.packages('tidyverse')
install.packages('stringr')

library(readxl)    
library(reshape2)
library(tidyverse)
library(stringr)

#Function to make read in all the sheets of a spreadsheet
read_excel_allsheets <- function(filename) {
  sheets <- readxl::excel_sheets(filename)
  x <-    lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  names(x) <- sheets
  x
}
#mysheets <- read_excel_allsheets("S1_File.xlsx")

# Define Function for converting NAs to 0s when comparing value scores
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
ips <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/")
atlas <- data.frame(read.csv(paste(file.path, "S1_File_TableB.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

names(atlas)[names(atlas) == 'ID_link.'] <- 'PID'

pos <- atlas[atlas$PID %in% ips$PID,]

pos <- melt(pos)

pos <- subset(pos, select=-c(Organisme., CD, X))

pos <- pos[pos$variable != "ENTREZ.geneID",]
pos <- pos[pos$variable != "Protein.count.",]

pos$value <- na.zero(pos$value)

split <- str_split_fixed(pos$CSPA.category," - ", 2)

pos <- cbind.data.frame(pos, split)

pos <- subset(pos, select=-c(CSPA.category, `1`))

colnames(pos) <- c("PID", "Line", "Presence", "Confidence")

pos$Presence[pos$Presence == 0] <- "Absent"
pos$Presence[pos$Presence == 1] <- "Present"

pos$Confidence = factor(pos$Confidence, levels = c("high confidence", "putative", "unspecific"))

colors <- c("blue","red","yellow")

names(colors) = c("high confidence", "putative", "unspecific")

pos <- merge(pos, ips[, c("PID", "Gene")], by="PID")

ggplot(pos) + 
  facet_wrap(~ Presence, scales = "free") +
  geom_tile(aes(x = Gene, y = Line, fill = Confidence)) +  
  scale_fill_manual(values=colors) +
  theme(axis.text.x = element_text(size = 8, angle = 90), 
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 8))


###### MOUSE
file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
m_pid <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSMouseGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/")
m_atlas <- data.frame(read.csv(paste(file.path, "S1_File_TableC.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

names(m_atlas)[names(m_atlas) == 'ID_link.'] <- 'PID'

m_pos <- m_atlas[m_atlas$PID %in% m_pid$PID,]

m_pos <- melt(m_pos)

m_pos <- subset(m_pos, select=-organisme.)

m_pos <- m_pos[m_pos$variable != "Protein.count.",]

m_pos$value <- na.zero(m_pos$value)

split <- str_split_fixed(m_pos$CSPA.category," - ", 2)

m_pos <- cbind.data.frame(m_pos, split)

m_pos <- subset(m_pos, select=-c(CSPA.category, `1`))

colnames(m_pos) <- c("PID", "Line", "Presence", "Confidence")

m_pos$Presence[m_pos$Presence == 0] <- "Absent"
m_pos$Presence[m_pos$Presence == 1] <- "Present"

m_pos$Confidence = factor(m_pos$Confidence, levels = c("high confidence", "putative", "unspecific"))

colors <- c("blue","red","yellow")

names(colors) = c("high confidence", "putative", "unspecific")

m_pos <- merge(m_pos, m_pid[, c("PID", "Gene")], by="PID")

ggplot(m_pos) + 
  facet_wrap(~ Presence, scales = "free") +
  geom_tile(aes(x = Gene, y = Line, fill = Confidence)) +  
  scale_fill_manual(values=colors) +
  theme(axis.text.x = element_text(size = 8, angle = 90), 
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 8))
