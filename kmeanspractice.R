data1 <- data.frame(read.csv("~/Desktop/data_1024.csv", sep = ""), header = TRUE, sep = ",", stringsAsFactors = FALSE)
ggplot(data1) + 
  geom_point(aes(x = Distance_Feature, y = Speeding_Feature))

clusters <- kmeans(data1[, c(2:3)], 4, nstart = 20)

table(clusters$cluster, data1$Speeding_Feature)
