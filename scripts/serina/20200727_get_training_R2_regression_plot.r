# Install packages
pacman::p_load("tidyverse", "caret", "Metrics", "ggpmisc", "ggrepel")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/TARA-thiolases/")

# Read in the dataset
# for 84 dataset with physiochemical properties for channels
origdat <- read_csv("data/residue_extraction/channelAB_84_OleA_physical_aa_features_extracted.csv") %>%
  full_join(read_csv("data/84_OleA_temps_noNAs.csv"))
rawdat <- origdat %>%
  dplyr::filter(!is.na(temperature_range)) %>%
  dplyr::select(1:121, temperature)

# Only keep variables with nonzero variance
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 
which_rem

# Check for duplicates(safety check)
dat <- rawdat[!duplicated(rawdat),]

# Set random seed 
set.seed(20200728) # today's date

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

train_plot <- data.frame(pred = rf$pred$pred,
                         obs = rf$pred$obs)
my.formula <- y ~ x
summ <- summary(lm(train_plot$obs ~ train_plot$pred))


ggplot(data = train_plot, aes(x = obs, y = pred)) +
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", method = "lm",
              lty = 2, lwd = 1, alpha = .5) +
  theme_classic() +
  xlab("Observed temperature") +
  ylab("Predicted temperature") +
  ggpmisc::stat_poly_eq(formula = my.formula,
                        aes(label = paste(..rr.label.., sep = "~~~")),
                        parse = TRUE)
