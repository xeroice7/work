
# Reveal number of observations and variables in two different ways
str(iris)
dim(iris)

# Show first and last observations in the iris data set
head(iris)
tail(iris)

# Summarize the iris data set
summary(iris)

# The Wage dataset is available

# Build Linear Model: lm_wage (coded already)
lm_wage <- lm(wage ~ age, data = Wage) #lm(response ~ predictor)

# Define data.frame: unseen (coded already)
unseen <- data.frame(age = 60)

# Predict the wage for a 60-year old worker
predict(lm_wage, unseen)

# The emails dataset is already loaded into your workspace

# Show the dimensions of emails
dim(emails)

# Inspect definition of spam_classifier()
spam_classifier <- function(x){
  prediction <- rep(NA, length(x)) # initialize prediction vector
  prediction[x > 4] <- 1
  prediction[x >= 3 & x <= 4] <- 0
  prediction[x >= 2.2 & x < 3] <- 1
  prediction[x >= 1.4 & x < 2.2] <- 0
  prediction[x > 1.25 & x < 1.4] <- 1
  prediction[x <= 1.25] <- 0
  return(prediction) # prediction is either 0 or 1
}

# Apply the classifier to the avg_capital_seq column: spam_pred
spam_pred <- spam_classifier(emails$avg_capital_seq)

# Compare spam_pred to emails$spam. Use ==
spam_pred == emails$spam

# linkedin is already available in your workspace

# Create the days vector
days <- c(1:21)

# Fit a linear model called on the linkedin views per day: linkedin_lm
linkedin_lm <- lm(linkedin ~ days)

# Predict the number of views for the next three days: linkedin_pred
future_days <- data.frame(days = 22:24)
linkedin_pred <- predict(linkedin_lm, future_days)

# Plot historical data and predictions
plot(linkedin ~ days, xlim = c(1, 24))
points(22:24, linkedin_pred, col = "green")

# Set random seed. Don't remove this line.
set.seed(1)

# Chop up iris in my_iris and species
my_iris <- iris[-5]
species <- iris$Species

# Perform k-means clustering on my_iris: kmeans_iris
kmeans_iris <- kmeans(my_iris, 3)

# Compare the actual Species to the clustering using table()
table(species, kmeans_iris$cluster)

# Plot Petal.Width against Petal.Length, coloring by cluster
plot(Petal.Length ~ Petal.Width, data = my_iris, col = kmeans_iris$cluster)


# Set random seed. Don't remove this line.
set.seed(1)

# Take a look at the iris dataset
str(iris)
summary(iris)

# A decision tree model has been built for you
tree <- rpart(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
              data = iris, method = "class")

# A dataframe containing unseen observations
unseen <- data.frame(Sepal.Length = c(5.3, 7.2),
                     Sepal.Width = c(2.9, 3.9),
                     Petal.Length = c(1.7, 5.4),
                     Petal.Width = c(0.8, 2.3))

# Predict the label of the unseen observations. Print out the result.
predict(tree, unseen, type = "class")

# The cars data frame is pre-loaded

# Set random seed. Don't remove this line.
set.seed(1)

# Explore the cars dataset
str(cars)
summary(cars)

# Group the dataset into two clusters: km_cars
km_cars <- kmeans(cars, 2)

# Print out the contents of each cluster
km_cars$cluster

# The cars data frame is pre-loaded

# Set random seed. Don't remove this line
set.seed(1)

# Group the dataset into two clusters: km_cars
km_cars <- kmeans(cars, 2)

# Add code: color the points in the plot based on the clusters
plot(cars, col = km_cars$cluster)

# Print out the cluster centroids
km_cars$centers

# Replace the ___ part: add the centroids to the plot
points(km_cars$centers, pch = 22, bg = c(1, 2), cex = 2)