# Install packages
pacman::p_load("tidyverse", "caret", "rsample", "ranger", "e1071")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/TARA-thiolases/")

# Read in the big dataset
dat <- read_csv("data/full50_allstats.csv")
colnames(dat)

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
# 50 rows and 89 columns

# Only keep variables with nonzero variance
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] # didn't remove any!

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
getTrainPerf(rf) # Training set accuracy is 64% for random forest

# Testing set
rf_pred <- predict(rf, newdata = form_test)
rf_pred
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$temp_status))
cm_rf

# See how bad it does on the entire dataset
x_dat <- dat %>%
  select(-nams, -temp_status)
rf <- train(
  x = x_dat,
  y = dat$temp_status,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")
getTrainPerf(rf)

