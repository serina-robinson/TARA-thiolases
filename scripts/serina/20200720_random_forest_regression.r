# Install packages
pacman::p_load("tidyverse", "caret", "Metrics", "ggpmisc", "ggrepel")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/TARA-thiolases/")

# Read in the dataset
origdat <- read_csv("data/123_OleA_allstats.csv") 
prop <- read_csv("data/50_protein_props.csv")
head(dat)

# Look at the distribution...is it normal?
hist(origdat$temperature, col = 'skyblue3', breaks = 5) # Looks fairly normal with bin size of 5
hist(origdat$temperature, col = 'skyblue3', breaks = 40) # Less normal with bin size of 40


# Prepare data for machine learning
rawdat <- origdat %>%
  dplyr::filter(!is.na(temperature_range)) %>%
  dplyr::select(contains(colnames(prop)), 11:43, temperature,
                -temperature_range, -sqs, -acc)


# Only keep variables with nonzero variance
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] # didn't remove any!
which_rem

# Check for duplicates(safety check)
dat <- rawdat[!duplicated(rawdat),]


# Set random seed 
set.seed(5678)

# Split into test and training data
dat_split <- rsample::initial_split(rawdat, strata = "temperature", prop = 0.7)
dat_train <- rsample::training(dat_split)
dat_test  <- rsample::testing(dat_split)
nrow(dat_test)
nrow(dat_train)/nrow(dat) # 82.9% of the data is in dat_train

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
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, 
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")


# Calculate perfomrmance
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

# You can label the ones that are far frmo the line
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
  ylab("Residuals") 

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