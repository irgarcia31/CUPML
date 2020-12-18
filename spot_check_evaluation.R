#' Spot-check algorithms and prediction
#' 
#' This function is thougth to be used inside a cross-validation loop, therefore,
#' no caret::trControl is used to specify cross-validation settings.
#' 
#' @param dataTrain data.frame with samples as rows and features as columns. 
#' There must be an extra column with class labels.
#' @param dataTest data.frame with samples as rows and features as columns. 
#' There must be an extra column with class labels.
#' @param label character specifying the column name with class labels.
#' @param metric A string that specifies what summary metric will be used to select the 
#' optimal model. By default, possible values are "RMSE" and "Rsquared" for regression 
#' and "Accuracy" and "Kappa" for classification.
#' Default: Accuracy
#' @param weights A numeric vector of case weights. This argument will only affect models 
#' that allow case weights.
#' Default: NULL
#' @param algorithms character vector with the machine learning algortihms to spot-check.
#' It is desirable to try a mixture of algorithms. The function uses "lda", "nb", 
#' "svmRadial", "rpart", "knn", "rf", "C5.0", "treebag" and "gbm" by default
#' @param return a list with $Accuracy and $Kappa information of the models specify in 
#' algorithms

spot.check <- function(dataTrain, 
                       dataTest, 
                       label,
                       metrics = "Accuracy", 
                       algorithms = c("lda", "nb", "svmRadial", "rpart", "knn", "rf", "C5.0", "treebag", "gbm"), 
                       weights = NULL) {
  
  # First, check if input data is in the correct format
  if (!is.data.frame(dataTrain) | !is.data.frame(dataTest)) {
    warning("Input must be a data.frame.")
    dataTrain <- as.data.frame(dataTrain)
    dataTest <- as.data.frame(dataTest)
  }
  if (!is.factor(dataTrain[, length(dataTrain)]) | !is.factor(dataTest[, length(dataTest)])) {
    warning("Labels columns must be a factor. Converting to factor...")
    dataTrain[, length(dataTrain)] <- as.factor(dataTrain[, length(dataTrain)])
    dataTest[, length(dataTest)] <- as.factor(dataTest[, length(dataTest)])
  }

  # Spot-check algorithms
  accuracy <- list()
  kappa <- list()
  for (i in algorithms){
    set.seed(seed)
    # Train the model with the training data
    fit <- caret::train(dataTrain[, -length(dataTrain)], 
                 dataTrain[, length(dataTrain)],
                 method = i, 
                 metric = metric,
                 weights = weights)
    # Make predictions with the testing data
    library(caret)
    predictions <- predict(fit, dataTest)
    matrix <- caret::confusionMatrix(predictions, dataTest[, label])
    # Save the metrics of each model
    accuracy[[i]] <- matrix$overall["Accuracy"]
    kappa[[i]] <- matrix$overall["Kappa"]
  }
  
  my_list <- list("Accuracy" = accuracy, "Kappa" = kappa)
  return(my_list)

}

# LDA
# Warning messages:
#   1: In lda.default(x, grouping, ...) : variables are collinear
# Multicollinearity means that your predictors are correlated. Why is this bad?
# Because LDA, like regression techniques involves computing a matrix inversion, 
# which is inaccurate if the determinant is close to 0 (i.e. two or more variables
# are almost a linear combination of each other).
# More importantly, it makes the estimated coefficients impossible to interpret. 
# If an increase in 𝑋1, say, is associated with an decrease in 𝑋2 and they both increase variable 𝑌, every change in 𝑋1 will be compensated by a change in 𝑋2 and you will underestimate the effect of 𝑋1 on 𝑌. In LDA, you would underestimate the effect of 𝑋1
# on the classification.
# If all you care for is the classification per se, and that after training your 
# model on half of the data and testing it on the other half you get 85-95% accuracy I'd say it is fine.

# NB
# Warning messages:
#   1: In FUN(X[[i]], ...) :
#   Numerical 0 probability for all classes with observation 1
