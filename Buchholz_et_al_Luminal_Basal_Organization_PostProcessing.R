### This script post processes the data obtained from Luminal_Basal_Organization.py and collected in a .csv file.

# Load packages, define file reading function and load data ---------------
library(plyr)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)

#Function to read csv file and import filename
read_plus <- function(flnm) {
  read_csv(flnm, skip=0) %>% 
    mutate(filename = flnm)
}

#Set working directory and replace "YourFolderPath" by the folder path your .csv file is stored in
working_directory <- "/YourFolderPath"
setwd(working_directory)  

# import Luminal_Basal data from csv and make dataframe
pat = ".csv"
files <- list.files(path = working_directory, pattern = pat, full.names = TRUE, recursive = TRUE)
Luminal_Basal <- ldply(files, read_plus)

colnames(Luminal_Basal) <- c("Experiment", "Condition", "Luminal_First", "Correct_Organization", "filename")

# Ensure Casted_Printed is a factor and set the levels in the desired order
Luminal_Basal$Condition <- factor(Luminal_Basal$Condition, 
                                  levels = c("Static", "Perfused"))


# Plot data ---------------------------------------------------------------
# Plot Luminal_First data
plot_1 <- ggplot(Luminal_Basal, aes(x = Condition, y = Luminal_First, fill = Condition)) +
  geom_boxplot() +  # Use geom_boxplot() instead of boxplot()
  labs(title = "Percentage of inner edge of the duct covered by K8/18 positive cells",
       x = "Condition",
       y = "Percentage of first cell layer that is K8/18 positive",
       fill = "Group") +
  theme_minimal()

print(plot_1)

# Plot Correct organization data
plot_2 <- ggplot(Luminal_Basal, aes(x = Condition, y = Correct_Organization, fill = Condition)) +
  geom_boxplot() +  # Use geom_boxplot() instead of boxplot()
  labs(title = "Percentage of duct surface covered first by luminal cells and then by basal cells",
       x = "Condition",
       y = "Percentage of duct surface covered first by luminal cells and then by basal cells",
       fill = "Group") +
  theme_minimal()

print(plot_2)


# Calculate statistics --------------------------------------------------------------
# Calculate summary statistics
summary_Luminal_Basal <- Luminal_Basal %>%
  group_by(Condition) %>%
  summarize(
    Average_Luminal_First = mean(Luminal_First),
    sd_Luminal_First = sd(Luminal_First),
    n = n(),
    Average_Correct_Organization = mean(Correct_Organization),
    sd_Correct_Organization = sd(Correct_Organization)
  )

# Perform ANOVA for Luminal first data
anova_result_Luminal_Basal_Luminal_First <- aov(Luminal_First ~ Condition, data = Luminal_Basal)
summary(anova_result_Luminal_Basal_Luminal_First)

# Perform Tukey's HSD post-hoc test for luminal first data
tukey_result_luminal_first <- TukeyHSD(anova_result_Luminal_Basal_Luminal_First)

print(tukey_result_luminal_first)

# Perform ANOVA for correct organization data
anova_result_Luminal_Basal_Correct_Organization <- aov(Correct_Organization ~ Condition, data = Luminal_Basal)
summary(anova_result_Luminal_Basal_Correct_Organization)

# Perform Tukey's HSD post-hoc test for correct organization data
tukey_result_Correct_Organization <- TukeyHSD(anova_result_Luminal_Basal_Correct_Organization)

print(tukey_result_Correct_Organization)

