library(tidyverse)
library(RColorBrewer)
library(reshape2)

nan.na <- function (x) {
  x[is.nan(x)] <- NA
  return(x)
}

file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/brca/tcga/")
proteinDF <- data.frame(read.csv(paste(file.path, "data_protein_quantification_Zscores.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
rnaDF <- read.csv(paste(file.path, "data_mRNA_median_Zscores.csv", sep=""), sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
gene <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

split.names <- str_split_fixed(proteinDF$Composite.Element.REF,"[|]", 2) #special character, so needed [ ]
proteinDF <- data.frame(proteinDF, split.names)
proteinDF <- subset(proteinDF, select = -c(Composite.Element.REF, `X2`))
names(proteinDF)[names(proteinDF) == 'X1'] <- 'Gene'

ips.targets <- proteinDF[proteinDF$Gene %in% gene,]
ips.targets <- melt(ips.targets)
ips.targets <- subset(ips.targets, select = -variable)
ips.targets$value <- nan.na(ips.targets$value)
ips.targets <- ips.targets %>%
  group_by(Gene) %>%
  mutate(mean = mean(value, na.rm = TRUE))

ips.hits <- rnaDF[rnaDF$Hugo_Symbol %in% gene,]
ips.hits <- subset(ips.hits, select = -Entrez_Gene_Id)
ips.hits <- melt(ips.hits)
ips.hits <- subset(ips.hits, select = -variable)
ips.hits$value <- nan.na(ips.hits$value)

ggplot(ips.hits) + 
  geom_point(aes(x=Hugo_Symbol, y=value)) + 
  ggtitle("Median Z-Scores (mRNA) vs Gene by Tumor Grade") + 
  xlab("Median Z-score") + 
  ylab("Gene") + 
  theme(axis.text.x = element_text(angle = 90))

ips.hits <- ips.hits %>%
  group_by(Hugo_Symbol) %>%
  mutate(mean = mean(value, na.rm = TRUE))



######
# cBio - TCGA Nature, 2012 Dataset
#####

file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/brca_tcga_pub2012/brca_tcga_pub/")
rnaDF <- read.csv(paste(file.path, "data_mRNA_median_Zscores.csv", sep=""), sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
patientDF <- data.frame(read.csv(paste(file.path, "data_clinical.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))

colnames(patientDF) <- lapply(patientDF[5,], as.character)
patientDF <- patientDF[-(1:5),]

file.path <- ("~/Desktop/Clay/Mass Spec Results/02-19-16/")
gene <- data.frame(read.csv(paste(file.path, "TotalUniqueiPSGENEandPID.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
gene <- as.factor(gene$Gene)

ips.hits <- rnaDF[rnaDF$Hugo_Symbol %in% gene,]
ips.hits <- subset(ips.hits, select = -Entrez_Gene_Id)
ips.hits <- melt(ips.hits)
ips.hits$tumor_grade <- patientDF$Tumor[match(ips.hits$variable, patientDF$PATIENT_ID)]
ips.hits$stage <- patientDF$`AJCC Stage`[match(ips.hits$variable, patientDF$PATIENT_ID)]
ips.hits$value <- nan.na(ips.hits$value)

#Plot of Z-scores vs Tumor Grade by Gene
ggplot(ips.hits) + 
  geom_point(aes(x=Hugo_Symbol, y=value, color = tumor_grade), position = position_jitter(width = 0.4)) + 
  ggtitle("Median Z-Scores (mRNA) vs Gene by Tumor Grade (TCGA, Nature 2012)") + 
  xlab("Gene") + 
  ylab("Median Z-score") + 
  theme_dark() + 
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))

#Plot of Z-scores vs Tumor Stage by Gene
ggplot(ips.hits) + 
  geom_point(aes(x=Hugo_Symbol, y=value, color = stage), position = position_jitter(width = 0.4)) + 
  ggtitle("Median Z-Scores (mRNA) vs Gene by Stage (TCGA, Nature 2012)") + 
  xlab("Gene") + 
  ylab("Median Z-score") + 
  theme_dark() + 
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))

#Plotting Z-Scores vs Gene by Tumor Grade ----- Faceted by TUMOR GRADE, colored by STAGE
ggplot(ips.hits) + 
  facet_wrap(~tumor_grade, scales = 'free') + 
  geom_point(aes(x=Hugo_Symbol, y=value, color = stage), position = position_jitter(width = 0.4)) + 
  ggtitle("Median Z-scores (mRNA) vs Gene by Tumor Grade (TCGA, Nature 2012)") + 
  xlab("Gene") + 
  ylab("Median Z-score") + 
  theme_dark() + 
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))

#Plotting Expression Level vs Gene by Tumor Stage  ----- Faceted by STAGE, colored by TUMOR GRADE
ggplot(ips.hits) + 
  facet_wrap(~stage, scales = 'free') +
  geom_point(aes(x=Hugo_Symbol, y=value, color = tumor_grade), position = position_jitter(width = 0.4)) + 
  ggtitle("Median Z-scores (mRNA) vs Gene by Tumor Stage (TCGA, Nature, 2012)") + 
  xlab("Gene") + 
  ylab("Median Z-score") + 
  theme_dark() + 
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))

hits <- subset(ips.hits, select = -variable)

ips.means <- hits %>%
  group_by(stage, tumor_grade, Hugo_Symbol) %>%
  mutate(avg = mean(value, na.rm = TRUE))

ips.means <- ips.means %>%
  group_by(stage, tumor_grade, Hugo_Symbol) %>%
  mutate(sd = sd(value, na.rm = TRUE))

ggplot(ips.means) + 
  facet_wrap(~stage, scales = 'free') + 
  geom_point(aes(x=Hugo_Symbol, y=avg, color=tumor_grade)) + 
  ggtitle("Mean Z-Scores vs Gene by Tumor Stage (TCGA, Nature, 2012)") + 
  xlab("Gene") + 
  ylab("Mean Z-Score") + 
  theme_dark() + 
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))

write.csv(ips.means, "means.csv")


######
# cBio - METABRIC Nature 2012 & Nat Com 2016 Dataset
#####

file.path <- ("~/Desktop/Clay/Mass Spec Results/WebData/brca_metabric/")
tumordata <- read.csv(paste(file.path, "data_clinical_supp_sample.csv", sep=""), sep = ",", stringsAsFactors = FALSE)
patientdata <- read.csv(paste(file.path, "data_clinical_supp_patient.csv", sep=""), sep = ",", stringsAsFactors = FALSE)
expressiondata <- read.csv(paste(file.path, "data_expression.csv", sep=""), sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

#Merging Data frames
patientDF <- merge(x = tumordata, y = patientdata, by = "PATIENT_ID", all = TRUE)

expression <- subset(expressiondata, select = -Entrez_Gene_Id)
expression <- melt(expression)

ips.hits <- expression[expression$Hugo_Symbol %in% gene,]
ips.hits$tumor_grade <- patientDF$GRADE[match(ips.hits$variable, patientDF$PATIENT_ID)]
ips.hits$stage <- patientDF$TUMOR_STAGE[match(ips.hits$variable, patientDF$PATIENT_ID)]
ips.hits$value <- nan.na(ips.hits$value)

#Plotting Expression Level vs Gene by Tumor Grade
ggplot(ips.hits) + 
  geom_point(aes(x=Hugo_Symbol, y=value, color = tumor_grade), position = position_jitter(width = 0.4)) + 
  ggtitle("Expression Level vs Gene by Tumor Grade (METABRIC, Nat 2012 & Nat Com 2016)") + 
  xlab("Gene") + 
  ylab("Expression Level") + 
  theme_dark() + 
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))

#Plotting Expression Level vs Gene by Tumor Stage
ggplot(ips.hits) + 
  geom_point(aes(x=Hugo_Symbol, y=value, color = stage), position = position_jitter(width = 0.4)) + 
  ggtitle("Expression Level vs Gene by Tumor Stage (METABRIC, Nat 2012 & Nat Com 2016)") + 
  xlab("Gene") + 
  ylab("Expression Level") + 
  theme_dark() + 
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))

#Plotting Expression Level vs Gene by Tumor Grade ----- Faceted by TUMOR GRADE, colored by STAGE
ggplot(ips.hits) + 
  facet_wrap(~tumor_grade, scales = 'free') + 
  geom_point(aes(x=Hugo_Symbol, y=value, color = stage), position = position_jitter(width = 0.4)) + 
  ggtitle("Expression Level vs Gene by Tumor Grade (METABRIC, Nat 2012 & Nat Com 2016)") + 
  xlab("Gene") + 
  ylab("Expression Level") + 
  theme_dark() + 
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))

#Plotting Expression Level vs Gene by Tumor Stage  ----- Faceted by STAGE, colored by TUMOR GRADE
ggplot(ips.hits) + 
  facet_wrap(~stage, scales = 'free') +
  geom_point(aes(x=Hugo_Symbol, y=value, color = tumor_grade), position = position_jitter(width = 0.4)) + 
  ggtitle("Expression Level vs Gene by Tumor Stage (METABRIC, Nat 2012 & Nat Com 2016)") + 
  xlab("Gene") + 
  ylab("Expression Level") + 
  theme_dark() + 
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))

hits <- subset(ips.hits, select = -variable)

ips.means <- hits %>%
  group_by(stage, tumor_grade, Hugo_Symbol) %>%
  mutate(avg = mean(value, na.rm = TRUE))

ips.means <- ips.means %>%
  group_by(stage, tumor_grade, Hugo_Symbol) %>%
  mutate(sd = sd(value, na.rm = TRUE))

ggplot(ips.means) + 
  facet_wrap(~stage, scales = 'free') + 
  geom_point(aes(x=Hugo_Symbol, y=avg, color=tumor_grade)) + 
  ggtitle("Mean Expression Level vs Gene by Tumor Stage (METABRIC, Nat 2012 & Nat Com 2016)") + 
  xlab("Gene") + 
  ylab("Mean Expression Level") + 
  theme_dark() + 
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))


######Hierarchial Clustering Analysis
install.packages("ISLR")
library(ISLR)
zzz <- expressiondata
zzz <- zzz[,-(1:2)]
colnames(zzz) <- NULL
#ncidat1 = t(NCI60$data)
ncidat <- t(zzz)
#colnames(ncidat1) = NCI60$labs
#unique(colnames(ncidat1))
#dim(ncidat)

#colnames(ncidat) <- lapply(ncidat[1,], as.character)
colnames(ncidat) <- expressiondata$Hugo_Symbol
#ncidat <- ncidat[-1,]  
#ncidat <- ncidat[,-1]  
X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
sv = t(X)
#is.na(X) <- do.call(cbind,lapply(X, is.infinite))
#length(which(!is.finite(as.matrix(X))))
w <- which(is.na(as.matrix(sv)))
#sv <- sv[which(!is.finite(as.matrix(sv)))] <- 0
sv[w] <- 0
sv = svd(sv)
U = sv$u
V = sv$v
D = sv$d


aa<- grep("grey",colors())
bb<- grep("green",colors())
cc<-  grep("red",colors())
gcol2<- colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

## use the genes that drive the first PC1. This is the first major patter in the data
k=1
ord1<- order(abs(V[,k]),decreasing=TRUE)
x1<- as.matrix(X[ord1[1:nrow(X)],])
heatmap(x1,col=gcol2)


expression$tumor_grade <- patientDF$GRADE[match(expression$variable, patientDF$PATIENT_ID)]
expression$stage <- patientDF$TUMOR_STAGE[match(expression$variable, patientDF$PATIENT_ID)]
#expression$value <- nan.na(ips.hits$value)

x <- subset(expression, select = -variable)

y <- dcast(x, Hugo_Symbol ~ stage, mean)
z <- dcast(x, Hugo_Symbol ~ stage, median)
#write.csv(y, "test.csv")

dy <- dist(as.matrix(y))   # find distance matrix 
dz <- dist(as.matrix(z))

hcy <- hclust(dy)                # apply hirarchical clustering 
hcz <- hclust(dz)

plot(hcy)

