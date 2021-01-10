#!/usr/bin/env Rscript
# this allows you to execute the script directly with ./[name of the script]
#==============================================================================
# Machine Learning workflow with tidymodels package

# Details:

# Note:

# Author: Inés Rodríguez García
# Version: 0.0.1
# Date: 2020/12/18
#==============================================================================
## PACKAGES
library(tidymodels)
library(data.table)

#==============================================================================
### GLOBAL VARIABLES
# data folder with a b-values CSV file previously generated.
DATA_DIR <- "/mnt/ElRaid/irodriguez/Nextcloud/share/20201005_IRodriguez_all_TCGA/all_betas"
DATA_DIR_test <- "/mnt/ElRaid/irodriguez/Nextcloud/share/20201129_IRodriguez_ML_TCGA/data/tb_values/"
# folder where the results will be written, it should exists.
RESULTS_DIR <- "/mnt/ElRaid/irodriguez/Nextcloud/share/20201129_IRodriguez_ML_TCGA/results"
# b-values CSV file name
BVALUES <- "all_samples.csv"
# sample sheet file name
SAMPLESHEET <- "samplesheet_ines.csv"
# column name with class label for ML classification
class_label <- "tumor_tissue_site"
# Pearson's correlation threshold above which probes are considered to be higly 
# correlated
threshold <- 0.8

# you should take into account that tidymodels ML is performed with a data frame
# with samples as rows and features as columns, so it may be necessary to modify
# BVALUES object to fit this requirements
#==============================================================================
### FUNCTIONS


#==============================================================================
### MAIN
###############################################################################
## Step 1: data loading
###############################################################################
# To perform supervised ML it is necessary to have a class label column
targets <- file.path(DATA_DIR, SAMPLESHEET)

# Load b-values
# In this particular case, the data frame has samples as rows and features as 
# columns, therefore it is necessay to transpose the data frame
# Besides, it is necessary to add a class label column 

tcga <- fread(file.path(DATA_DIR, BVALUES))
tcga <- tcga %>% t %>% as.data.frame

# Besides, it is necessary to add a class label column 
tcga_data <- tcga %>% 
	add_column(Class = as.factor(targets[, class_label]))

# load sample for testing
tcga_data <- NULL
for (i in list.files(DATA_DIR_test)) {
  tbVals_project <- readRDS(paste0(DATA_DIR_test, i))
  tbVals_30 <- tbVals_project[sample(nrow(tbVals_project),size=30), ]
  tcga_data<- rbind(tcga_data, tbVals_30)                
}


###############################################################################
## Step 2: data splitting
###############################################################################
# Fix the random numbers by setting the seed 
# This enables the analysis to be reproducible when random numbers are used 
set.seed(1234)

# Put 0.7 of the data into the training set
data_split <- initial_split(tcga_data, prop = 0.7, strata = Class)
# strata argument is needed to ensure that the distributions are the same in 
# both training and testing datasets when creating the resamples

# if the data is too small, this function cannot stratify.
# for this purpose caret works better

# Create data frames for the two sets:
train_data <- training(data_split)
test_data  <- testing(data_split)

library(caret)
set.seed(1234)
trainIndex <- createDataPartition(tcga_data[, class_label], p=0.70, list=FALSE)
train_data <- tcga_data[ trainIndex,]
test_data <- tcga_data[-trainIndex,]

# Check same distribution of output
bind_rows(
  as.data.frame(round(prop.table(table(train_data$tumor_tissue_site)),4)) %>% mutate(Data = "Train"),
  as.data.frame(round(prop.table(table(test_data$tumor_tissue_site)),4)) %>% mutate(Data = "Test")
) %>%
  spread(Var1, Freq) %>%
  kable(style = "pandoc",
        align = c('c','c','c'))


###############################################################################
## Step 3: preprocess and feature selection
###############################################################################
# Before training the model, a recipe to preprocess and select features must be 
# created.
tcga_rec <- 
	recipe(Class ~ ., data = train_data) %in%

# NA imputation
ratio_recipe <- tcga_rec %>%
  step_knnimpute(all_predictors(), neighbors = 5) %>%
  step_corr(all_predictors(), threshold = .5)
# apperently step_knnimpute() can impute datasets where every sample has at least
# one missing value, unlike caret::knnImpute
if (method == "knn") {
  system.time(rec <- recipes::recipe(class_label ~ ., data = dataTrain ))
  %>%
    recipes::step_knnimpute(all_predictors(), neighbors = neighbors) %>%
    recipes::prep(training = dataTrain)
recipe <- recipes::prep(recipe, training = dataTrain)
dataTrain <- recipes::bake(rec, dataTrain)
dataTest <- recipes::bake(rec, dataTest)
my_list <- list("dataTrain" = dataTrain, "dataTest" = dataTest)
return(my_list)
}

# PROBLEM: couldn't alocate vector of 3... GB size

# maybe add check_missing() which stop pred if NAs are detected

# filter CpGs (correlation, statistical)
# some wrapped feature selection

# prep() For a recipe with at least one preprocessing operation, estimate the 
# required parameters from a training set that can be later applied to other 
# data sets.
# juice() As steps are estimated by prep, these operations are applied to the 
# training set. Rather than running bake() to duplicate this processing, this 
# function will return variables from the processed training set.
# bake() For a recipe with at least one preprocessing operation that has been 
# trained by prep.recipe(), apply the computations to new data.
###############################################################################
## Step 4: evaluate the model
###############################################################################


###############################################################################
## Step 5: tune model parameters
###############################################################################


###############################################################################
## Step 6: save the model
###############################################################################
