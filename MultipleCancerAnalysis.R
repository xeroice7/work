##### To find most likely/intersting/promising candidates - BY STAGE:
#Total iPS Overlap - STAGE
metabric_st <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/Breast_Metabric_Stats_Stage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
metabric_gr <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_metabric/Analysis/Breast_Metabric_Stats_Grade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
brca_tcga_st <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_tcga_pub/brca_tcga_pub/Analysis/Breast_TCGA_Stats_Stage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
brca_tcga_gr <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Breast/brca_tcga_pub/brca_tcga_pub/Analysis/Breast_TCGA_Stats_Grade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

luad_tcga_st <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/Analysis/Lung_LUAD_TCGA_Stats_Stage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
luad_tcga_gr <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/luad_tcga_prov/tcga/Analysis/Lung_LUAD_TCGA_Stats_Grade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
lusc_tcga_st <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/lusc_tcga_prov/tcga/Analysis/Lung_LUSC_TCGA_Stats_Stage.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))
lusc_tcga_gr <- data.frame(read.csv("~/Desktop/Clay/Mass Spec Results/WebData/Lung/lusc_tcga_prov/tcga/Analysis/Lung_LUSC_TCGA_Stats_Grade.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE))

stage.br <- merge(metabric_st, brca_tcga_st, by = "Genes")
stage_sig.br <- filter(stage.br, pvalues.x < 0.05 & pvalues.y < 0.05)
stage_near.br <- filter(stage.br, pvalues.x < 0.1 & pvalues.y < 0.1)
stage_near.br <- filter(stage_near.br, !stage_near.br$Genes %in% stage_sig.br$Genes)

grade.br <- merge(metabric_gr, brca_tcga_gr, by = "Genes")
grade_sig.br <- filter(grade.br, pvalues.x < 0.05 & pvalues.y < 0.05)
grade_near.br <- filter(grade.br, pvalues.x < 0.1 & pvalues.y < 0.1)
grade_near.br <- filter(grade_near.br, !grade_near.br$Genes %in% grade_sig.br$Genes)

total_sig.br <- merge(stage_sig.br, grade_sig.br, by = "Genes")
total_near.br <- merge(stage_near.br, grade_near.br, by = "Genes")

write.csv(grade_sig.br, "Significant_Grade_BR.csv")
write.csv(stage_sig.br, "Significant_Stage_BR.csv")
write.csv(grade_near.br, "Interesting_Grade_BR.csv")
write.csv(stage_near.br, "Interesting_Stage_BR.csv")

stage.lu <- merge(luad_tcga_st, lusc_tcga_st, by = "Genes")
stage_sig.lu <- filter(stage.lu, pvalues.x < 0.05 & pvalues.y < 0.05)
stage_near.lu <- filter(stage.lu, pvalues.x < 0.1 & pvalues.y < 0.1)
stage_near.lu <- filter(stage_near.lu, !stage_near.lu$Genes %in% stage_sig.lu$Genes)

grade.lu <- merge(luad_tcga_gr, lusc_tcga_gr, by = "Genes")
grade_sig.lu <- filter(grade.lu, pvalues.x < 0.05 & pvalues.y < 0.05)
grade_near.lu <- filter(grade.lu, pvalues.x < 0.1 & pvalues.y < 0.1)
grade_near.lu <- filter(grade_near.lu, !grade_near.lu$Genes %in% grade_sig.lu$Genes)

total_sig.lu <- merge(stage_sig.lu, grade_sig.lu, by = "Genes")
total_near.lu <- merge(stage_near.lu, grade_near.lu, by = "Genes")

write.csv(grade_sig.lu, "Significant_Grade_LU.csv")
write.csv(stage_sig.lu, "Significant_Stage_LU.csv")
write.csv(grade_near.lu, "Interesting_Grade_LU.csv")
write.csv(stage_near.lu, "Interesting_Stage_LU.csv")