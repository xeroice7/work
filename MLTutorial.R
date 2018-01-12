# Part 1: Load the data ===================

# Read in `iris` data
iris <- read.csv(url("http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data"), 
                 header = FALSE) 

# Print first lines
head(iris)

# Add column names
names(iris) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width", "Species")

# Check the result
iris

# Part 2: Know the data =====================

# Load in `ggvis`
library(ggvis)

# Iris scatter plot
iris %>% ggvis(~Sepal.Length, ~Sepal.Width, fill = ~Species) %>% layer_points()

iris %>% ggvis(~Petal.Length, ~Petal.Width, fill = ~Species) %>% layer_points()


# Overall correlation `Petal.Length` and `Petal.Width`
cor(iris$Petal.Length, iris$Petal.Width)

# Return values of `iris` levels 
x=levels(iris$Species)

# Print Setosa correlation matrix
print(x[1])
cor(iris[iris$Species==x[1],1:4])

# Print Versicolor correlation matrix
print(x[2])
cor(iris[iris$Species==x[2],1:4])

# Print Virginica correlation matrix
print(x[3])
cor(iris[iris$Species==x[3],1:4])

# Return all `iris` data
iris

# Return first 5 lines of `iris`
head(iris)

# Return structure of `iris`
str(iris)

# Division of `Species`
table(iris$Species) # Displays breakdown of # of species for each set

# Percentual division of `Species`
round(prop.table(table(iris$Species)) * 100, digits = 1)  # Gives breakdown in percentages

# Summary overview of `iris`
summary(iris) 

# Refined summary overview
summary(iris[c("Petal.Width", "Sepal.Width")]) # Only summary of these two columns

# Part 3: Which attributes to choose? ==================

#Need to decide on the use cases that would be relevant for your data set
#Think abolut what your data set may teach you or information you can get from it
#Then, come up with algorithms that you need to apply to get the desired results

# Part 4: Prepare the workspace =============================

#Most algorithms are not built into R, so grab the required packages

library(class)

#To see if any package is installed
#any(grepl("<name of your package>", installed.packages()))

# Part 5: Prepare the data =================================

#We want to normalize the data when the data has very different ranges for 
#variables, ie if your dataset has just two attributes, X and Y, and X has values 
#that range from 1 to 1000, while Y has values that only go from 1 to 100, then Y’s 
#influence on the distance function will usually be overpowered by X’s influence.
#When you normalize, you actually adjust the range of all features, so that distances 
#between variables with larger ranges will not be over-emphasised.

# Build your own `normalize()` function
normalize <- function(x) {
  num <- x - min(x)
  denom <- max(x) - min(x)
  return (num/denom)
}

# Normalize the `iris` data
iris_norm <- as.data.frame(lapply(iris[1:4], normalize))
 
# Summarize `iris_norm`
summary(iris_norm)

set.seed(1234) #This is a # in R's random # generator, which is reproducible

#Take a sample with a size that is set as a # of rows of the dataset, and choose from a 
# vector of 2 elements and assign either 1 or 2 to the 150 rows with probability
# weights of 0.67 and 0.33, and replace = T means you assign a 1 or 2 to a certain row
# and then reset the vector of 2 to its original state, meaning for the next rows
# in the dataset, you can assign either a 1 or 2 each time again
# This essentially divides the training set to a random 2/3 of the data, and the test
# set to 1/3 of the data
ind <- sample(2, nrow(iris), replace=TRUE, prob=c(0.67, 0.33)) 

#We use the 1 and 2s to define our training and test sets

# Compose training set
iris.training <- iris[ind==1, 1:4] 

# Inspect training set
head(iris.training)

# Compose test set
iris.test <- iris[ind==2, 1:4]

# Inspect test set
head(iris.test)

# We want to predict the attribute "Species", but we do want to include it into 
#the KNN algorithm or there will not be a prediction for it, so we need to store 
#class variables in factor vectors and divde them over training and test sets

# Compose `iris` training labels
iris.trainLabels <- iris[ind==1,5]

# Inspect result
print(iris.trainLabels)

# Compose `iris` test labels
iris.testLabels <- iris[ind==2, 5]

# Inspect result
print(iris.testLabels)

#Part 6: The KNN Model ==============================

#Now we want to find the k nearest neighbors of the training set, so we take the 
# knn function and use some simple arguments, and use an odd k number ti avoid ties in 
# voting scores. The output of the knn function is a factor vector with the predicted classes
# for each row of the test data, but we do not want to insert test labels yet b/c we
# want to verify if our model is actually good

# Build the model
iris_pred <- knn(train = iris.training, test = iris.test, cl = iris.trainLabels, k=3)

# Inspect `iris_pred`
iris_pred

# Part 7: Evaluation of the model ============================

#Need to determine if our model is actually good at predicting results

# Put `iris.testLabels` in a data frame
irisTestLabels <- data.frame(iris.testLabels)

# Merge `iris_pred` and `iris.testLabels` 
merge <- data.frame(iris_pred, iris.testLabels) #This takes the predicted and test labels to determine whether our model is good 

# Specify column names for `merge`
names(merge) <- c("Predicted Species", "Observed Species")

# Inspect `merge` 
merge # Take a look at the results of our predictions

library(gmodels)

CrossTable(x = iris.testLabels, y = iris_pred, prop.chisq=FALSE)

#Now try making a model using the caret package in a similar way

library(caret)

# Create index to split based on labels  
index <- createDataPartition(iris$Species, p=0.75, list=FALSE)

# Subset training set with index
iris.training <- iris[index,]

# Subset test set with index
iris.test <- iris[-index,]

# Overview of algos supported by caret
names(getModelInfo())

# Train a model
model_knn <- train(iris.training[, 1:4], iris.training[, 5], method='knn')

#Possible to make other models simply by changing the method argument
model_cart <- train(iris.training[, 1:4], iris.training[, 5], method='rpart2')

# Predict the labels of the test set
predictions<-predict(object=model_knn,iris.test[,1:4])

# Evaluate the predictions
table(predictions)

# Confusion matrix 
confusionMatrix(predictions,iris.test[,5]) # Allows us to see the stats and results behind the predictions

# Now lets try preprocessing data with scaling and centering 
# Train the model with preprocessing
model_knn <- train(iris.training[, 1:4], iris.training[, 5], method='knn', preProcess=c("center", "scale"))

# Predict values
predictions<-predict.train(object=model_knn,iris.test[,1:4], type="raw")

# Confusion matrix
confusionMatrix(predictions,iris.test[,5])


