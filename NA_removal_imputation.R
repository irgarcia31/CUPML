#' Handle NA for Machine Learning Classification 
#' 
#' Many Machine Learning algorithms do not support missing values, therefore it is 
#' necessary to handle them. For this purpose, the following function remove probes
#' with too many NAs and impute the remaining missing values, avoiding data leakage.
#' 
#' @param dataTrain a data frame with training samples in rows and CpG in columns.
#' For "knn" and "median" methods, non-numeric data will not be pre-processed and 
#' their values will be in the data frame produced by the predict function. 
#' Therefore, it is possible to have additional columns with labels for ML 
#' classification.
#' For "mean" method it is necessary to exclude this column from the computation.
#' Therefore, it is necessary to specify the name column in @param class_label 
#' @param dataTest a data frame with testing samples in rows and CpG in columns. 
#' For "knn" and "median" methods, non-numeric data will not be pre-processed and 
#' their values will be in the data frame produced by the predict function. 
#' Therefore, it is possible to have additional columns with labels for ML 
#' classification.
#' For "mean" method it is necessary to exclude this column from the computation.
#' Therefore, it is necessary to specify the name column in @param class_label 
#' @param method a character value with imputation method
#' Default: "knn"
#' @param ProbeCutoff a numeric value of frecuency of missing values to filter out 
#' probes
#' Default: 0.1
#' @param k the number of nearest neighbours from the training set to use for 
#' imputation
#' Default: 5
#' @param class_label a character value corresponding to the column name that containing 
#' the labels use to classify
#' @param verbose logical value. Prints a log as the computations proceed
#' Default: FALSE
#' @details Three imputation criteria can be used: "knn", median" and "mean".
#' k-nearest neighbour imputation is carried out by finding the k closest samples 
#' (Euclidian distance) in the training set. This method is simple, accurate and 
#' accepts missing values, but it has much higher computational cost
#' Imputation via medians takes the median of each predictor in the training set, 
#' and uses them to fill missing values. This method is simple, fast, and accepts 
#' missing values, but treats each predictor independently, and may be inaccurate.
#' Imputation via means takes the mean of each predictor in the training set, 
#' and uses them to fill missing values. 
#' @return A list with the same training and testing data frame but without probes 
#' with too many NAs and missing value imputation

remove.impute.NA <- function(dataTrain, 
                            dataTest,
                            method = "median",
                            ProbeCutoff = 0.1,
                            k = 5,
                            class_label = "",
                            verbose = FALSE)
{
  if (verbose == TRUE){
    print("<<<STARTING REMOVAL AND IMPUTATION OF NA...>>>")
    print(paste("Number of NAs in the training set before removal:", sum(is.na(dataTrain))))
    print(paste("Number of NAs in the testing set before removal:", sum(is.na(dataTest))))
  }
  # Remove probes with too many NA values (> 10% by default)
  print(paste("(1) Removing probes with more than ", ProbeCutoff*100, "% of NAs..."))
  keep <- mclapply(dataTrain, function(x) {
    percentage_NA <-sum(is.na(x)) / dim(dataTrain)[1]
    return(ifelse(percentage_NA > ProbeCutoff, FALSE, TRUE))
  })
  if (verbose == TRUE){
    print(paste0("Number of probes to remove based on % of NAs: ", length(unlist(keep)) - sum(unlist(keep))))
  }
  dataTrain <- dataTrain[, unlist(keep)]
  dataTest <- dataTest[, unlist(keep)]
  if (method == "median") {
    print("(2) Imputing remaining NAs...")
    print("Method: Median")
    # Impute the remaining NAs with the mean 
    # transform the data frame to a numerical matrix to speed up the imputation
    class_label_train <- as.factor(dataTrain[, class_label])
    class_label_test <- as.factor(dataTest[, class_label])
    train_matrix <- dataTrain[, !names(dataTrain) %in% class_label]
    test_matrix <- dataTest[, !names(dataTest) %in% class_label]
    train_matrix <- t(train_matrix)
    test_matrix <- t(test_matrix)
    # extract the median value of each CpG to impute NA
    median_NA <- matrixStats::rowMedians(train_matrix, na.rm = TRUE)
    # create a logical matrix with NA information
    logical_matrix_NA_train <- is.na(train_matrix)
    logical_matrix_NA_test <- is.na(test_matrix)
    if (verbose == TRUE){
      print(paste("Number of NAs in the training set before imputation:", sum(is.na(train_matrix))))
      print(paste("Number of NAs in the testing set before imputation:", sum(is.na(test_matrix))))
    }
    # impute NA with median values from the training data
    train_matrix[logical_matrix_NA_train] <- median_NA[row(train_matrix)][logical_matrix_NA_train]
    test_matrix[logical_matrix_NA_test] <- median_NA[row(test_matrix)][logical_matrix_NA_test]
    if (verbose == TRUE){
      print(paste("Number of NAs in the training set after imputation:", sum(is.na(train_matrix))))
      print(paste("Number of NAs in the testing set after imputation:", sum(is.na(test_matrix))))
    }
    # transform the matrix to a transpose data frame with tumor tissue site information
    # with samples in rows and CpG in columns
    train_matrix_df <- as.data.frame(t(train_matrix))
    test_matrix_df <- as.data.frame(t(test_matrix))
    train_df_tumor <- cbind(train_matrix_df, as.factor(class_label_train))
    test_df_tumor <- cbind(test_matrix_df, as.factor(class_label_test))
    colnames(train_df_tumor)[dim(train_df_tumor)[2]] <- class_label
    colnames(test_df_tumor)[dim(test_df_tumor)[2]] <- class_label
    my_list <- list("dataTrain" = train_df_tumor, "dataTest" = test_df_tumor)
    return(my_list)
  }
  else if (method == "mean") {
    print("(2) Imputing remaining NAs...")
    print("Method: Mean")
    # Impute the remaining NAs with the mean 
    # transform the data frame to a numerical matrix to speed up the imputation
    class_label_train <- as.factor(dataTrain[, class_label])
    class_label_test <- as.factor(dataTest[, class_label])
    train_matrix <- dataTrain[, !names(dataTrain) %in% class_label]
    test_matrix <- dataTest[, !names(dataTest) %in% class_label]
    train_matrix <- t(train_matrix)
    test_matrix <- t(test_matrix)
    # extract the median value of each CpG to impute NA
    mean_NA <- rowMeans(train_matrix, na.rm = TRUE)
    # create a logical matrix with NA information
    logical_matrix_NA_train <- is.na(train_matrix)
    logical_matrix_NA_test <- is.na(test_matrix)
    if (verbose == TRUE){
      print(paste("Number of NAs in the training set before imputation:", sum(is.na(train_matrix))))
      print(paste("Number of NAs in the testing set before imputation:", sum(is.na(test_matrix))))
    }
    # impute NA with median values from the training data
    train_matrix[logical_matrix_NA_train] <- mean_NA[row(train_matrix)][logical_matrix_NA_train] # very slow
    test_matrix[logical_matrix_NA_test] <- mean_NA[row(test_matrix)][logical_matrix_NA_test]
    if (verbose == TRUE){
      print(paste("Number of NAs in the training set after imputation:", sum(is.na(train_matrix))))
      print(paste("Number of NAs in the testing set after imputation:", sum(is.na(test_matrix))))
    }
    # transform the matrix to a transpose data frame with tumor tissue site information
    # with samples in rows and CpG in columns
    train_matrix_df <- as.data.frame(t(train_matrix))
    test_matrix_df <- as.data.frame(t(test_matrix))
    train_df_tumor <- cbind(train_matrix_df, as.factor(class_label_train))
    test_df_tumor <- cbind(test_matrix_df, as.factor(class_label_test))
    colnames(train_df_tumor)[dim(train_df_tumor)[2]] <- class_label
    colnames(test_df_tumor)[dim(test_df_tumor)[2]] <- class_label
    my_list <- list("dataTrain" = train_df_tumor, "dataTest" = test_df_tumor)
    return(my_list)
  }
  print("<<< END OF REMOVAL AND IMPUTATION OF NA >>>")
}
