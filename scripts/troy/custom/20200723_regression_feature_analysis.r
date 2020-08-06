# Load packages
pacman::p_load("tidyverse")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases")
# Read data
dat <- read_csv("data/residue_extraction/OleA_machine_learning_results.csv")
colnames(dat)
dat_stats <- dat %>% 
  # select(feature_set, tt_split, training_rmse) %>%
  group_by(feature_set, tt_split) %>%
  summarize(tr_rmse_m = mean(training_rmse),
            tr_rmse_sd = sd(training_rmse),
            tr_r2_m = mean(training_r2),
            tr_r2_sd = sd(training_r2),
            oob_m = mean(oob_error),
            oob_sd = sd(oob_error),
            te_rmse_m = mean(testing_rmse),
            te_rmse_sd = sd(testing_rmse),
            te_r2_m = mean(testing_r2),
            te_r2_sd = sd(testing_r2),)

dat_stats
write_csv(dat_stats, "output/residue_extraction/regression_summary_stats.csv")
channels <- sort(c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246, 253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353))
view(data.frame(channels))



