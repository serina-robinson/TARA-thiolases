# Install packages
pacman::p_load("tidyverse", "caret", "Metrics", "ggpmisc", "ggrepel")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the dataset
###### for 84 dataset with physiochemical properties for channels
origdat <- read_csv("data/residue_extraction/channelAB_84_OleA_physical_aa_features_extracted.csv") %>%
  full_join(read_csv("data/84_OleA_temps_noNAs.csv"))
rawdat <- origdat %>%
  dplyr::filter(!is.na(temperature_range)) %>%
  dplyr::select(1:121, temperature)

###### for 84 dataset with physiochemical properties for 8 angstrom residues
origdat <- read_csv("data/residue_extraction/8_angstrom_84_OleA_physical_aa_features_extracted.csv") %>%
  full_join(read_csv("data/84_OleA_temps_noNAs.csv"))
rawdat <- origdat %>%
  dplyr::filter(!is.na(temperature_range)) %>%
  dplyr::select(1:166, temperature)

###### for 84 dataset with physiochemical properties for 10 angstrom residues
origdat <- read_csv("data/residue_extraction/10_angstrom_84_OleA_physical_aa_features_extracted.csv") %>%
  full_join(read_csv("data/84_OleA_temps_noNAs.csv"))
rawdat <- origdat %>%
  dplyr::filter(!is.na(temperature_range)) %>%
  dplyr::select(1:251, temperature)

###### for 84 dataset with physiochemical properties for 12 angstrom residues
origdat <- read_csv("data/residue_extraction/12_angstrom_84_OleA_physical_aa_features_extracted.csv") %>%
  full_join(read_csv("data/84_OleA_temps_noNAs.csv"))
rawdat <- origdat %>%
  dplyr::filter(!is.na(temperature_range)) %>%
  dplyr::select(1:421, temperature)

rf_list <- list()

# Only keep variables with nonzero variance
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 
which_rem

# Check for duplicates(safety check)
dat <- rawdat[!duplicated(rawdat),]


# Set random seed 
set.seed(1234)

# Split into test and training data
for(i in 1:1){
dat_split <- rsample::initial_split(rawdat, strata = "temperature", prop = 0.8)
dat_train <- rsample::training(dat_split)
dat_test  <- rsample::testing(dat_split)
nrow(dat_test)
nrow(dat_train)/nrow(dat)

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
                          # repeats = 3,
                           verboseIter = T, 
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")

rf_list[[i]] <- getTrainPerf(rf)
rf_list[[i]]$oob_error <- rf$finalModel$prediction.error # not positive what this is for
rf_list[[i]]$test_r2 <- rf$finalModel$r.squared # same with this 
rf_pred <- predict(rf, form_test)
rf_list[[i]]$test_rmse <- Metrics::rmse(rf_pred, y_test)
}
# Calculate performance
getTrainPerf(rf) # RMSE is 11.16

rf_imp <- varImp(rf)

ggplot(rf_imp, top = 32) + 
  xlab("") +
  theme_classic()

# Predict for new data
rf_pred <- predict(rf, newdata = form_test)

Metrics::rmse(y_test, rf_pred) 
rf_df <- data.frame(cbind(rf_pred, y_test))
rownames(rf_df) <- rownames(form_test)

my.formula <- y ~ x
summary(lm(rf_df$y_test ~ rf_pred))
summary(lm(rf_pred ~ rf_df$y_test))


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

# You can label the ones that are far from the line
ggplot(rf_df, aes(x = y_test, y = rf_pred)) +
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", method = "lm",
              lty = 2, lwd = 1, alpha = .5) +
  theme_classic() +
  xlab("Observed temperature") +
  ylab("Predicted temperature") +
  ggpmisc::stat_poly_eq(formula = my.formula,
                        aes(label = paste(..rr.label.., sep = "~~~")),
                        parse = TRUE) +
  geom_label_repel(label = rownames(rf_df),
                   force = 10, size = 1.5)

# Make a residual plot
rf_df_resid <- rf_df %>%
  mutate(resid = y_test - rf_pred)

ggplot(rf_df_resid, aes(x = rf_pred, y = resid)) +
  geom_point(alpha = .3) + 
  geom_line(col = "red", y = 0.0,
            lty = 2, lwd = 1, alpha = .5) +
  theme_classic() + 
  xlab("Predicted temperature") +
  ylab("Residuals")  +
  geom_label_repel(label = rownames(rf_df),
                   force = 10, size = 1.5)

# ?ggrepel
## Challenge 1. Play around with the 'ggrepel' package to display points on the residual plot.
# What is a residual plot/what are residuals? Which points have the highest residuals? 

# Challenge 2. Look up what RMSE, MAE, and R-squared mean

# Challenge 3. Try different random seeds and/or training/test split proportions and see how 
# these metrics and residual and observed v. predicted plots change.

# Super challenge 4. Try making residual and observed v. predicted plots for the training data?

# Challenge 5. How would you statistically test if temperature has a normal distribution.
# If it doesn't have one, what types of transformations could you do to normalize it? 

# Challenge 6. If you do any transformations to the dependent variable challenge 5, 
# do your results change at all?