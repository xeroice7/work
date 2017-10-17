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
               fill = c("lightgreen", "green", "lightblue", "blue"))