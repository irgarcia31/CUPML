# PACKAGES
library(data.table) # load dataset
library(RColorBrewer) # color palette
library(Rtsne) # data visualization
library(umap) # data visualization
library(factoextra) # data visualization
library(caret) # machine learning
library(Hmisc) # mean and median imputation
library(limma) # Differential Methylation Analysis
library(impute) # impute knn
library(doParallel) # paralellize

# GLOBAL VARIABLES
# Path where the data of the project is stored
DATA_DIR <- "/mnt/ElRaid/irodriguez/Nextcloud/share/20201005_IRodriguez_all_TCGA/all_betas/"
# Path where the results of the project are stored
RESULTS_DIR <- "/mnt/ElRaid/irodriguez/Nextcloud/share/20201129_IRodriguez_ML_TCGA/results/"
# Path where the code of the project is stored
CODE_DIR <- "/mnt/ElRaid/irodriguez/Nextcloud/share/20201129_IRodriguez_ML_TCGA/code/"
# CSV file with samples information, it MUST have a column with tissue site 
# information for supervised classification
SAMPLESHEET <- "samplesheet_ines.csv"
# Name of the column with tissue site information
SAMPLE_GROUP <- "Sample_Group"
# CSV file with b-values
BVALUES <- "all_samples.csv"
seed <- 1234
# If FALSE, no log will be printed
verbose <- FALSE
# Pearson's correlation threshold above which probes are considered to be higly 
# correlated
cutoff_correlated <- 0.8
# ML model evaluation metrics
metric <- "Accuracy"
# accuracy and kappa: default metrics on binary and multi-class datasets (kappa 
# is more suitable for unbalanced data).
# RMSE and R^2: default metrics on regression datasets.
# AUROC: only for binary datasets
# LogLoss: useful in multi-class datasets
# For Cross-Validation
# Select number of repeats for cross-validation
r <- 2
# Select number of folds for cross-validation
k <- 3
# Select number of features for limma selection
n <- 1000

# FUNCTIONS
source(file = file.path(CODE_DIR, "NA_removal_imputation.R"))
source(file = file.path(CODE_DIR, "limma_fs.R"))
source(file = file.path(CODE_DIR, "spot_check_evaluation.R"))

################################################################################### 
##STEP 1: LOAD DATA
###################################################################################
print("1. DATA LOADING")
# # for testing
# tbVals <- NULL
# for (l in list.files(DATA_DIR)) {
#   tbVals_project <- readRDS(paste0(DATA_DIR, l))
#   tbVals_30 <- tbVals_project[sample(nrow(tbVals_project),size=30), ]
#   tbVals <- rbind(tbVals, tbVals_30)                
# }

samplesheet <- read.csv(file.path(DATA_DIR, SAMPLESHEET), skip = 7) # there has to be a better way that "skip"
system.time(dataset <- fread(file.path(DATA_DIR, BVALUES))) # this function belongs to data.frame, is a faster way to load large datasets
Sample_Group <- c(list(SAMPLE_GROUP), as.list(samplesheet[, SAMPLE_GROUP]))
system.time(dataset <- rbindlist(list(dataset, Sample_Group)))
# BVALUES has samples as columns and CpG as rows, it is necessary to generate the 
# transpose
system.time(dataset <- transpose(dataset))
# impossible to create an additional row and do the transpose due a lack of 
# RAM memory 

################################################################################### 
##STEP 2: DATA EXPLORATION
###################################################################################
print("2. DATA EXPLORATION")

# Colour palette and shape for visualiazing 32 tumoral types
pal <- rep(brewer.pal(8, "Dark2"), 4)
pch <- c(rep(15,8), rep(16,8), rep(17, 8), rep(18,8))

## t-SNE = t-Distributed Stochastic Neighbor Embedding (of 50000 random CpGs)
## ---------------------------------------------------
# is a non-linear technique for dimensionality reduction that is particularly well 
# suited for the visualization of high-dimensional datasets.
print("Performing t-Distributed Stochastic Neighbor Embedding (t-SNE)...")
set.seed(1234)
tbVals_na <- tbVals[ , apply(tbVals, 2, function(x) !any(is.na(x)))]
tbVals_tsne_subset <- sample(tbVals_na[-length(tbVals_na)], 50000)
tsne <- Rtsne(as.matrix(tbVals_tsne_subset), 
              dims = 2, 
              perplexity=6, 
              verbose=TRUE, 
              max_iter = 5000,
              check_duplicates = FALSE)
png(file = file.path(RESULTS_DIR, paste0("tsne_final.png")),
    res = 300,
    height = 2700,
    width = 2700)
plot(tsne$Y, 
     pch = pch[factor(tbVals_na$tumor_tissue_site)], 
     cex = 1 , 
     main = paste0("t-sne"), 
     col = pal[factor(tbVals_na$tumor_tissue_site)],
     xlab = "First dimension", 
     ylab = "Second dimension")
legend("topright", legend=levels(factor(tbVals_na$tumor_tissue_site)), text.col=pal, col = pal,
       bg="white", pch = pch, cex=0.5)
dev.off()

dataset <- tbVals

################################################################################
## STEP 3: PREPARE DATA
################################################################################
print("3. DATA PREPARATION FOR ML")
# Split-out validation dataset
print("Splitting-out validation dataset...")
set.seed(1234)
trainIndex <- createDataPartition(dataset[, class_label], p=0.70, list=FALSE)
dataTrain <- dataset[ trainIndex,]
dataTest <- dataset[-trainIndex,]

################################################################################
## STEP 4: EVALUATE ALGORITHMS
################################################################################
dataTrain$ID <- c(1:(nrow(dataTrain)))

# repeated 3 times 10-fold cross-validation

dataTrain <- iris
dataTrain$ID <- c(1:(nrow(dataTrain)))

cl <- makePSOCKcluster(8)
registerDoParallel(cl)

# loop for doing 3 repetitions of cross-validation
for (j in 1:r) {
	
	# create lists for saving the final metrics
	final_accuracy <- list()
	final_kappa <- list()
	
	# loop for generating 10 folds from the training dataset and perform cross-validation
	for (i in 1:k) {
		
		print(paste(j, " x ", i))
    		
		# create lists for saving the metrics
    accuracy <- list()
    kappa <- list()
    
		# create 10 random folds from the training dataset
    print("Creating folds...")
   	set.seed(seed)
    folds <- createFolds(dataTrain[, class_label], k = k)
    
		# splitting the data into training and testing
		train <- dataTrain[dataTrain$ID %in% unlist(folds[-i]), ]
		test <- dataTrain[dataTrain$ID %in% folds[[i]], ]
		train$ID <- NULL
		test$ID <- NULL
    
		# A) FEATURE SELECTION
		print("Performing feature selection...")
    
		#   a.1) Handle NA
		print(" 1) NAs...")
		#   Some ML algorithm do not support missing values
		#   To adress this problem, a custom function has been created
		#   remove.inpute.NA function remove probes with more than "x" % of NAs and impute the rest of them
		set.seed(seed)
		data_list <- remove.impute.NA(dataTrain = train, dataTest = test, method = "median", class_label = class_label)
		train <- data_list$dataTrain
		test <- data_list$dataTest

		#   a.2) DMPs
		print(" 2) Limma...")
		set.seed(seed)
		keep <- limma_fs(train[, -length(train)], factor(train[, length(train)]), n)
		#   create a new dataset with DMPs
		keep <- colnames(train) %in% keep
		keep[length(keep)] <- TRUE
		print(paste("Features:", sum(keep)))
		train <- train[, keep]
		test <- test[, keep]
    
		#   a.3) Remove correlated attributes
		print(" 3) Correlated attributes...")
		#   find attributes that are highly correlated and remove them
		set.seed(seed)
		correlations <- cor(train[,-length(train)])
		highlyCorrelated <- findCorrelation(correlations, cutoff=cutoff_correlated)
		print(paste("Features:", dim(train)[2] - length(highlyCorrelated)))
		#   create a new dataset without highly correlated features
		train <- train[,-highlyCorrelated]
		test <- test[,-highlyCorrelated]
    
		#   a.4) Wrapped methods: Recursive Feature Selection
		print(" 4) RFE...")
		set.seed(seed)
		#   define the control using a random forest selection function
		control <- rfeControl(functions=rfFuncs, method="cv", number=10, repeats = 5)
		#   define a vector with the number of features that should be trained with
		subset <- seq(1, ncol(train), by = 1)
		#   run the RFE algorithm
		results <- rfe(train[,-length(train)],
		   as.factor(train[,length(train)]),
		   size = subset,
		   metric = "Kappa",
		   rfeControl = control)
		##   summarize the results
		#results
  	  	##   plot and save the results
		#png(filename = file.path(RESULTS_DIR, paste0("predictors_", j,"x", i, ".png")),
		#    res = 200,
    		#    width = 1000,
		#    height = 800)
		#plot(results, type=c("g", "o"))
		#dev.off()
		#write.csv(results$fit$importance, file = file.path(RESULTS_DIR, paste0("predictors_", j,"x", i, ".csv")))
		#write.csv(results$fit$confusion, file = file.path(RESULTS_DIR, paste0("predictors_confusion_matrix", j,"x", i, ".csv")))
		#   list the chosen features
		predictors <- predictors(results)
		keep <- colnames(train) %in% predictors
		keep[length(keep)] <- TRUE
		print(paste("Features:", sum(keep)))
		train <- train[, keep]
		test <- test[, keep]
    
		# Handle class imbalance
		model_weights <- NULL
		for (m in train[, class_label]){
			model_weights <- c(model_weights, (1/table(train[, class_label])[m]) * 0.5)
		}
    
		# B) SPOT-CHECK ALGORITHMS
		print("Performing spot-check...")
		set.seed(seed)
		models <- spot.check(train, test, label = class_label, weights = model_weights)
    
		# C) METRICS
		accuracy[[i]] <- models$Accuracy
		print(accuracy)
		# accuracy <- rbind(accuracy, models$Accuracy)
		kappa[[i]] <- models$Kappa
		print(kappa)
		# kappa <- rbind(kappa, models$Kappa)
	}
	
	final_accuracy[[j]] <- accuracy
	print(final_accuracy)
	# final_accuracy <- rbind(final_accuracy, accuracy)
	final_kappa[[j]] <- kappa
	print(final_kappa)
	# final_kappa <- rbind(final_kappa, kappa)
}

stopCluster(cl)

# summarize Best Model
print(fit.glmnet)

# ################################################################################
# ## STEP 5: IMPROVE ACCURACY
# ################################################################################
# # a) Algorithm Tuning
# # Grid Search
# trainControl <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
# set.seed(seed)
# tunegrid <- expand.grid(.alpha=seq(0, 1, by=0.05), .lambda=seq(0, 0.1, by=0.001))
# glmnetGrid <- train(tumor_tissue_site~., data=dataTrain, method="glmnet", metric=metric, tuneGrid=tunegrid,
#                 trControl=trainControl, na.action = na.omit)
# print(glmnetGrid)
# plot(glmnetGrid)
# # b) Ensembles
# # 10-fold cross validation with 3 repeats
# trainControl <- trainControl(method="repeatedcv", number=10, repeats=3)
# metric <- "Accuracy"
# # Bagged CART
# set.seed(seed)
# fit.treebag <- train(tumor_tissue_site~., data=dataTrain, method="treebag", metric=metric, preProc=c("center", "scale"),
#                      trControl=trainControl, na.action = na.omit)
# # Random Forest
# set.seed(seed)
# fit.rf <- train(tumor_tissue_site~., data=dataTrain, method="rf", metric=metric, preProc=c("center", "scale"),
#                 trControl=trainControl, na.action = na.omit)
# # Stochastic Gradient Boosting
# set.seed(seed)
# fit.gbm <- train(tumor_tissue_site~., data=dataTrain, method="gbm", metric=metric, preProc=c("center", "scale"),
#                  trControl=trainControl, verbose=FALSE, na.action = na.omit)
# # C5.0
# set.seed(seed)
# fit.c50 <- train(tumor_tissue_site~., data=dataTrain, method="C5.0", metric=metric, preProc=c("center", "scale"),
#                  trControl=trainControl, na.action = na.omit)
# # Warning messages:
# #   1: 'trials' should be <= 1 for this object. Predictions generated using 1 trials
# # Compare results
# ensembleResults <- resamples(list(BAG=fit.treebag, RF=fit.rf, GBM=fit.gbm, C50=fit.c50))
# summary(ensembleResults)
# png("results/ML/bVals_na_ensemble.png",
#     height = 500,
#     width = 400)
# dotplot(ensembleResults)
# dev.off()
# # summarize Best Model
# print(fit.treebag)
# ################################################################################
# ## STEP 6: FINALIZE MODEL
# ################################################################################
# # a)Predictions on validation dataset
# predictions <- predict(fit.glmnet, dataTest)
# confusionMatrix(predictions, dataTest$tumor_tissue_site)
# # b)Create standalone model on entire training dataset
# set.seed(seed)
# finalModel <- randomForest(Class~., training, mtry=2, ntree=2000)
# # c)Save model for later use
# saveRDS(finalModel, "./finalModel.rds")
# 
# 
# #===============================================================================

print("THE END!! :)")
