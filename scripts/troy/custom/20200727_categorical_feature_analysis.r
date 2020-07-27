# Load packages
pacman::p_load("tidyverse")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases")
# Read data
dat <- read_csv("output/residue_extraction/OleA_machine_learning_results_categorical.csv")
colnames(dat)
dat_stats <- dat %>% 
  group_by(groups, threshold, tt_split) %>%
  summarize(tr_acc_m = mean(training_accuracy),
            tr_acc_sd = sd(training_accuracy),
            oob_m = mean(oob_error),
            oob_sd = sd(oob_error),
            te_acc_m = mean(testing_accuracy),
            te_acc_sd = sd(testing_accuracy))

dat_stats
write_csv(dat_stats, "output/residue_extraction/categorical_summary_stats.csv")
channels <- sort(c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246, 253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353))
view(data.frame(channels))



