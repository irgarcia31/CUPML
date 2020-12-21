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

#==============================================================================
### GLOBAL VARIABLES
# data folder with a b-values CSV file previously generated.
DATA_DIR <- "/mnt/ElRaid/irodriguez/Nextcloud/share/20201005_IRodriguez_all_TCGA/all_betas"
# folder where the results will be written, it should exists.
RESULTS_DIR <- "/mnt/ElRaid/irodriguez/Nextcloud/share/20201129_IRodriguez_ML_TCGA/results"
# b-values CSV file name
BVALUES <- "all_samples.csv"
# sample sheet file name
SAMPLESHEET <- "samplesheet_ines.csv"
# column name with class label for ML classification
class_label <- "tumor_tissue_site"

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

tcga <- read_csv(file.path(DATA_DIR, BVALUES)) %>% t %>% as.data.frame

# Besides, it is necessary to add a class label column 
tcga_data <- tcga %>% 
	add_column(Class = as.factor(targets[, class_label]))


###############################################################################
## Step 2: data splitting
###############################################################################
# Fix the random numbers by setting the seed 
# This enables the analysis to be reproducible when random numbers are used 
set.seed(1234)

# Put 0.7 of the data into the training set
data_split <- initial_split(tcga_data, prop = 0.7)

# Create data frames for the two sets:
train_data <- training(data_split)
test_data  <- testing(data_split)


###############################################################################
## Step 3: preprocess and feature selection
###############################################################################
# Before training the model, a recipe to preprocess and select features must be 
# created.
tcga_rec <- 
	recipe(Class ~ ., data = train_data) %in%



###############################################################################
## Step 4: evaluate the model
###############################################################################


###############################################################################
## Step 5: tune model parameters
###############################################################################


###############################################################################
## Step 6: save the model
###############################################################################
