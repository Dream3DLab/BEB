# This script processes the area under the curve data derived from the 
# FITC_Quantification python script and summarizes it into a boxplot.

# Load packages and define function to read csv files -------------------------
library(plyr)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)

read_plus <- function(flnm) {
  read_csv(flnm, skip=0) %>% 
    mutate(filename = flnm)
}



#  Set working directory and load data. Replace "YourFolderPath"  --------
working_directory <- "/YourFolderPath"
setwd(working_directory)  

# import Area data and make dataframe
pat = ".csv"
files <- list.files(path = working_directory, pattern = pat, full.names = TRUE, recursive = TRUE)
Area_Under_The_Curve <- ldply(files, read_plus)

colnames(Area_Under_The_Curve) <- c("Experiment", "Condition", "Area", "Filename")


# Plot data per experiment ---------------------------------------------------------------
plot_1 <- ggplot(Area_Under_The_Curve, aes(x = Condition, y = Area, fill = Condition)) +
  geom_bar(stat = "identity") +  # Removed position="fill"
  labs(title = "Increase of area under the curve",
       x = "Condition",
       y = "Area",
       fill = "Group") +
  facet_grid(~ Experiment) +
  theme_minimal()

print(plot_1)

# Plot data accross experiments --------------------------------------------
# Calculate average area under the curve per condition across all experiments
average_area_per_condition <- Area_Under_The_Curve %>%
  group_by(Condition) %>%
  summarise(Average_Area = mean(Area, na.rm = TRUE))  # na.rm = TRUE to remove NA values if any

# Box plot to show the distribution of 'Area' for each 'Condition'
plot_2 <- ggplot(Area_Under_The_Curve, aes(x = Condition, y = Area, fill = Condition)) +
  geom_boxplot() +  # Using geom_boxplot to display distributions
  labs(title = "Distribution of Area Under the Curve Across Experiments",
       x = "Condition",
       y = "Area",
       fill = "Group") +
  theme_minimal()

print(plot_2)


# Calculate statistics ----------------------------------------------------
summary_Area_Under_The_Curve <- Area_Under_The_Curve %>%
  group_by(Condition) %>%
  summarize(
    Average_Area = mean(Area),
    sd_Area = sd(Area),
    n = n()
  )

# Perform ANOVA
anova_result_Area_Under_The_Curve <- aov(Area ~ Condition, data = Area_Under_The_Curve)
summary(anova_result_Area_Under_The_Curve)

# Perform Tukey's HSD post-hoc test
tukey_result <- TukeyHSD(anova_result_Area_Under_The_Curve)

print(tukey_result)



