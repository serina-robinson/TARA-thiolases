# Install packages
pacman::p_load("tidyverse", "caret", "Metrics", "ggpmisc", "ggrepel")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the dataset
# for 84 dataset with physiochemical properties for channels
origdat <- read_csv("data/residue_extraction/channelAB_84_OleA_physical_aa_features_extracted.csv") %>%
  full_join(read_csv("data/84_OleA_temps_noNAs.csv"))
rawdat <- origdat %>%
  dplyr::filter(!is.na(temperature_range)) %>%
  dplyr::select(1:121, temperature)

# Create list to populate in loop
rf_list <- list()

# Only keep variables with nonzero variance
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 
which_rem

# Check for duplicates(safety check)
dat <- rawdat[!duplicated(rawdat),]

# Set random seed 
set.seed(1234)

# Loop for random forest regression
for(i in 1:10){
  
  # Split into test and training data  
  dat_split <- rsample::initial_split(rawdat, strata = "temperature", prop = 0.8)
  dat_train <- rsample::training(dat_split)
  dat_test  <- rsample::testing(dat_split)

  # Independent variables
  x_train <- dat_train[,!colnames(dat_train) %in% c("nams", "temperature")]
  x_test <- dat_test[,!colnames(dat_test) %in% c("nams", "temperature")]

  # Dependent variable
  y_train <- dat_train$temperature
  y_test <- dat_test$temperature

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
                             verboseIter = T, 
                             savePredictions = "final"),
    verbose = TRUE,
    importance = "permutation")

  # Save information in list
  # Training data
  rf_list[[i]] <- getTrainPerf(rf) # Training RMSE, r-squared, and MAE
  rf_list[[i]]$OOBerror <- sqrt(rf$finalModel$prediction.error) # Out-of-bag error
  
  # Testing data
  rf_pred <- predict(rf, form_test) # Use model to predict test data
  rf_list[[i]]$TestRMSE <- Metrics::rmse(rf_pred, y_test) # Testing RMSE
  rf_df <- data.frame(cbind(rf_pred, y_test))
  my.formula <- y ~ x
  summ <- summary(lm(rf_df$y_test ~ rf_pred))
  rf_list[[i]]$TestRsquared <- summ$r.squared # Testing r-squared
}

# Turn list into data frame and select useful columns
final_df <- do.call(rbind.data.frame, rf_list) %>% 
  select(c(1:2, 5:7))
# Write to file
# write_csv(final_df, "output/residue_extraction/channelAB_regression_10_loops_0.8_split.csv")

# Plot observed vs. predicted

ggplot(rf_df, aes(x = y_test, y = rf_pred)) +
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", method = "lm",
              lty = 2, lwd = 1, alpha = .5) +
  theme_classic() +
  xlab("Observed temperature") +
  ylab("Predicted temperature") +
  ggpmisc::stat_poly_eq(formula = my.formula,
                        aes(label = paste(..rr.label.., sep = "~~~")),
                        parse = TRUE)
