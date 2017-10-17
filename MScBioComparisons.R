file.path <- ("~/Desktop/")
cna <- data.frame(read.csv(paste(file.path, "CNA_Genes.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
gene <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)


hits <- cna[cna$Gene %in% gene,]


file.path <- ("~/Desktop/")
cnaHG <- data.frame(read.csv(paste(file.path, "CNAGenes_LungHG.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
mutatedHG <- data.frame(read.csv(paste(file.path, "MutatedGenes_LungHG.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
cnaLG <- data.frame(read.csv(paste(file.path, "CNAGenes_LungLG.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
mutatedLG <- data.frame(read.csv(paste(file.path, "MutatedGenes_LungLG.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

file.path <- ("~/Desktop/")
totalgene <- data.frame(read.csv(paste(file.path, "iPSGene.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
totalgene <- as.factor(totalgene$x)

totalhitsCNA <- cnaHG[cnaHG$Gene %in% totalgene,]
totalhitsMUT <- mutatedHG[mutatedHG$Gene %in% totalgene,]

diffCNA <- cnaHG[!(cnaHG$Gene %in% cnaLG$Gene),]
diffMUT <- mutatedHG[!(mutatedHG$Gene %in% mutatedLG$Gene),]
