# Install packages
pacman::p_load("tidyverse", "caret", "rsample", "ranger", "e1071", "tidyverse")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the big dataset
dat <- read_csv("data/full50_allstats.csv")
dim(dat)
dat$topt
# Read in the protein properties
prop <- read_csv("data/50_protein_props.csv")

# Select desired columns
rawdat <- dat %>%
  dplyr::filter(!is.na(temperature_range)) %>%
  dplyr::mutate(temp_status = case_when(temperature_range == "Thermophilic" ~ "NP",
                                        temperature_range == "Mesophilic" ~ "NP",
                                        temperature_range == "Psychrophilic" ~ "P")) %>%
  dplyr::select(contains(colnames(prop)), 131:163, temp_status,
                -temperature_range, -newnams, -sqs, -acc) 

table(rawdat$temp_status)

# Only keep variables with nonzero variance
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] # didn't remove any!
which_rem

# Check for duplicates(safety check)
dat <- rawdat[!duplicated(rawdat),]

# Set random seed 
set.seed(1234)

# Split into test and training data
dat_split <- rsample::initial_split(rawdat, strata = "temp_status", prop = 0.8)
dat_train <- rsample::training(dat_split)
dat_test  <- rsample::testing(dat_split)
nrow(dat_test)
nrow(dat_train)/nrow(dat) # 82.9% of the data is in dat_train

# Define our response

# Independent variables
x_train <- dat_train[,!colnames(dat_train) %in% c("nams", "temp_status")]
x_test <- dat_test[,!colnames(dat_test) %in% c("nams", "temp_status")]

# Dependent variable
y_train <- dat_train$temp_status
y_test <- dat_test$temp_status

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nams)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nams)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$nams)

# Train a random forest classifier using the 'ranger' package
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")

# Confusion matrix
getTrainPerf(rf) # Training set accuracy is 60% for random forest

# Testing set accuracy
rf_pred <- predict(rf, newdata = form_test)
rf_pred
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$temp_status))
cm_rf # 71% 

## Challenge 1. Unfortunately right now we aren't getting very good performance...
# See if you can improve this first of all by increasing the dataset size using the 73 OleAs.
# If possible, try to get at least 30 in each category (NP and P)

## Challenge 2. Try varying the proportions of training/test splits to be
# 70%, 80%, and 90% (training set). Does it appear to dramatically change the outcome of the 
# training and testing scores? What do you think is happening?

## Challenge 3. Regardless of if accuracy improves, calculate variable importance. Which
# variables are most important? What do they mean?

## Challenge 4. Calculate and plot the area under the receiver operating characteristic curve

## Optional challenge 5. Try building a custom tuning grid like in the sample_machine_learning.r script


## Challenge 6. Try predicting for topt instead -> the topt predictions 
# using the tool TOMER (https://www.biorxiv.org/content/10.1101/2020.05.06.081737v1)
# Have been uploaded to your data folder as data/73_OleA_tomer_results.csv

## Challenge 7. Try feature selection using random forest to reduce the number
# of input variables in the model. Does it improve performance?

## Super challenge 8. Try putting all the code inside a for loop and doing
# 100 independent training-test splits. Store the training and testing accuracy scores
# for each iteration of the loop in a vector. Plot the distribution of these
# scores as a histogram.

## Super challenge 9. What other features could you try. For example...
# could you add columns for the channel residues? positive patch residues?





