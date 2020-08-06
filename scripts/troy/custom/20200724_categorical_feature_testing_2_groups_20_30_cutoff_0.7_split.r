# Install packages
pacman::p_load("tidyverse", "caret", "rsample", "ranger", "e1071", "tidyverse")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in dataset
# For 84 dataset with physiochemical properties
dat <- read_csv("data/residue_extraction/channelAB_84_OleA_physical_aa_features_extracted.csv") %>%
  full_join(read_csv("data/84_OleA_temps_noNAs.csv"))
rawdat <- dat %>%
  dplyr::filter(!is.na(temperature_range)) %>%
  dplyr::mutate(temp_status = case_when(temperature >= 30 ~ "highT",
                                        temperature <= 20 ~ "lowT",
                                        TRUE ~ "midT"
                                        )) %>%
  dplyr::select(1:121, temp_status) %>% 
  dplyr::filter(temp_status != "midT")



# Check counts in each condition
table(rawdat$temp_status)

# Create list to populate in loop
rf_list <- list()

# Only keep variables with nonzero variance
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 
which_rem
rawdat1 <- rawdat[,!colnames(rawdat) %in% which_rem]

# Check for duplicates(safety check)
dat <- rawdat1[!duplicated(rawdat1),]

# Set random seed 
set.seed(1234)

# Loop for random forest regression
for(i in 1:10){

  # Split into test and training data
  dat_split <- rsample::initial_split(rawdat, strata = "temp_status", prop = 0.7)
  dat_train <- rsample::training(dat_split)
  dat_test  <- rsample::testing(dat_split)

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
    trControl = trainControl(method = "cv", number = 10,
                             verboseIter = T, classProbs = T,
                             savePredictions = "final"),
    verbose = TRUE,
    importance = "permutation")

  # Save information in list

  rf_pred <- predict(rf, newdata = form_test)
  cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$temp_status))
  tr_acc <- getTrainPerf(rf)[1] 
  rf_list[[i]] <- c("training_Accuracy" = tr_acc$TrainAccuracy, "oob_error" = rf$finalModel$prediction.error, "testing" = cm_rf$overall[1])
}

# Turn list into data frame
final_df <- do.call(rbind.data.frame, rf_list) 
colnames(final_df) <- c("training_accuracy", "oob_error", "testing_accuracy")

# Write to file
write_csv(final_df, "output/residue_extraction/channelAB_categorical_10_loops_2_groups_20_30_cutoff_0.7_split.csv")

#Plot of variable importance
rf_imp <- varImp(rf, scale = FALSE, 
                 surrogates = FALSE, 
                 competes = FALSE)
rf_imp
ggplot(rf_imp, top = 20) + 
  xlab("") +
  theme_classic()


