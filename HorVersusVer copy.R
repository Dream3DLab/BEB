# This script will quantify cell alignment in imaged printed chips based on segmentation data from Imaris. Export Bounding box dimensions
# in the axes of interest as csv and put them into one merged folder. Be sure that the filenaming is consistent and clear.

#Load all necessary packages
library(plyr)
library(readr)
library(dplyr)
library(ggplot2)

#Function to read csv files in the directory and skip the first 3 rows (Row 4 contains column headers). Adjust this number if necessary
read_plus <- function(flnm, filename) {
  data <- read_csv(flnm, skip = 3)
  data$filename <- filename  # Add filename as a new column
  return(data)
}

working_directory <- "/Users/majbrittbuchholz/surfdrive/Shared/PhD/Thesis/Chapter_5_Breast_on_a_Chip/Analysis/Cell orientation/Stats_horvsver"
setwd(working_directory)

# import sum_red
pat <- "Length_X"
files <- list.files(path = working_directory, pattern = pat, full.names = TRUE, recursive = TRUE)
Circumferential_Length <- ldply(Map(read_plus, files, basename(files)))

# import sum_red
pat <- "Length_Z"
files <- list.files(path = working_directory, pattern = pat, full.names = TRUE, recursive = TRUE)
Longitudinal_Length <- ldply(Map(read_plus, files, basename(files)))

# Combine data frames using column names
Axis_Lengths <- cbind(Circumferential_Length[, c("BoundingBoxAA Length X", "ID", "filename")], 
                      Longitudinal_Length[, c("BoundingBoxAA Length Z", "ID", "filename")])

# Rename columns for clarity
colnames(Axis_Lengths) <- c("Circumferential_Length", "Circumferential_Length_ID", "Circumferential_Filename",
                            "Longitudinal_Length", "Longitudinal_Length_ID", "Longitudinal_Filename")

# Delete the second column
Axis_Lengths <- subset(Axis_Lengths, select = -c(2,3))

# Combine Filename and ID into Unique_ID
Axis_Lengths$Unique_ID <- paste(Axis_Lengths$Longitudinal_Filename, Axis_Lengths$Longitudinal_Length_ID, sep = "_")

# Delete the second column
Axis_Lengths <- subset(Axis_Lengths, select = -c(3))

Axis_Lengths$Celltype <- Axis_Lengths$Unique_ID
Axis_Lengths$Orientation <- Axis_Lengths$Unique_ID

Axis_Lengths$Celltype <- ifelse(grepl("HU", Axis_Lengths$Celltype), "HUVECs", Axis_Lengths$Celltype)
valid_orgs <- c("HUVECs")
Axis_Lengths$Celltype <- ifelse(Axis_Lengths$Celltype %in% valid_orgs, Axis_Lengths$Celltype, "MCF10A")

# Create the Alignment column based on the condition
Axis_Lengths$Alignment <- ifelse(Axis_Lengths$Longitudinal_Length > Axis_Lengths$Circumferential_Length,
                                 "Longitudinal", "Circumferential")

Axis_Lengths$Orientation <- ifelse(grepl("Exp064", Axis_Lengths$Orientation), "Horizontal", Axis_Lengths$Orientation)
Axis_Lengths$Orientation <- ifelse(grepl("Exp063", Axis_Lengths$Orientation), "Horizontal", Axis_Lengths$Orientation)
valid_orgs <- c("Horizontal")
Axis_Lengths$Orientation <- ifelse(Axis_Lengths$Orientation %in% valid_orgs, Axis_Lengths$Orientation, "Vertical")

# Calculate the percentage of Unique_ID for which Alignment is Longitudinal for each Longitudinal_Filename,
# including Celltype and Orientation
percentage_longitudinal <- Axis_Lengths %>%
  group_by(Longitudinal_Filename, Celltype, Orientation) %>%
  summarise(percentage_longitudinal = mean(Alignment == "Longitudinal") * 100)


# Create box plot
ggplot(percentage_longitudinal, aes(x = paste(Celltype, Orientation), y = percentage_longitudinal, fill = Orientation)) +
  geom_boxplot() +
  labs(title = "Percentage of Unique_ID with Longitudinal Alignment",
       x = "Celltype & Orientation",
       y = "Percentage Longitudinal") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############ Statistics
# Subset the data for MCF10A Celltype
subset_MCF10A <- subset(percentage_longitudinal, Celltype == "MCF10A")
# Perform statistical test separately for MCF10A Celltype
p_value_MCF10A <- t.test(percentage_longitudinal ~ Orientation, data = subset_MCF10A)$p.value

# Subset the data for HUVEC Celltype
subset_HUVEC <- subset(percentage_longitudinal, Celltype == "HUVECs")
# Perform statistical test separately for HUVEC Celltype
p_value_HUVEC <- t.test(percentage_longitudinal ~ Orientation, data = subset_HUVEC)$p.value

# Create box plot with color based on Orientation
p <- ggplot(percentage_longitudinal, aes(x = paste(Celltype, Orientation), y = percentage_longitudinal, fill = Orientation)) +
  geom_boxplot() +
  geom_point(aes(shape = Celltype), position = position_jitterdodge(dodge.width = 0.75), size = 3, fill = "black") +  # Map shape to Celltype
  labs(title = "Percentage of Unique_ID with Longitudinal Alignment",
       x = "Celltype & Orientation",
       y = "Percentage Longitudinal") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add lines for each Celltype
#p <- p +
#  geom_segment(data = subset_HUVEC, aes(x = 3, y = max(subset_HUVEC$percentage_longitudinal) + 0.05 * max(subset_HUVEC$percentage_longitudinal), 
#                                        xend = 4, yend = max(subset_HUVEC$percentage_longitudinal) + 0.05 * max(subset_HUVEC$percentage_longitudinal)),
#               linetype = "solid", color = "black") +
#  geom_segment(data = subset_MCF10A, aes(x = 1, y = max(subset_MCF10A$percentage_longitudinal) + 0.05 * max(subset_MCF10A$percentage_longitudinal), 
#                                         xend = 2, yend = max(subset_MCF10A$percentage_longitudinal) + 0.05 * max(subset_MCF10A$percentage_longitudinal)),
 #              linetype = "solid", color = "black") 

p <- p +
  geom_segment(data = subset_HUVEC, aes(x = 3, y = 102), 
                                        xend = 4, yend = 102,
               linetype = "solid", color = "black") +
  geom_segment(data = subset_MCF10A, aes(x = 1, y = 102), 
                                         xend = 2, yend = 102,
               linetype = "solid", color = "black") 
p <- p +
  annotate("text", x = c(1.5, 3.5), 
           y = 107, 
           label = c(
             ifelse(p_value_HUVEC < 0.05, ifelse(p_value_HUVEC < 0.01, ifelse(p_value_HUVEC < 0.001, ifelse(p_value_HUVEC < 0.0001, "****", "***"), "**"), "*"), "NS"), 
             ifelse(p_value_MCF10A < 0.05, ifelse(p_value_MCF10A < 0.01, ifelse(p_value_MCF10A < 0.001, ifelse(p_value_MCF10A < 0.0001, "****", "***"), "**"), "*"), "NS")
           ), 
           color = c(
             ifelse(p_value_HUVEC < 0.05, ifelse(p_value_HUVEC < 0.01, ifelse(p_value_HUVEC < 0.001, ifelse(p_value_HUVEC < 0.0001, "black", "black"), "black"), "black"), "black"), 
             ifelse(p_value_MCF10A < 0.05, ifelse(p_value_MCF10A < 0.01, ifelse(p_value_MCF10A < 0.001, ifelse(p_value_MCF10A < 0.0001, "black", "black"), "black"), "black"), "black")
           ), 
           size = 6)


p
