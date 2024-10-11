# This script uses nominal vectors exported from CloudCompare after nominal actual comparison of CAD models and 3D models of 
# 3D prints derived via light sheet imaging. Export this data from CloudCompare as .csv file. The script will take the nominal
# vectors to calculate the relative deviation of printed object points to the original CAD file and plot percentage deviations
# in all three axes. It further exports these percentage deviations as .csv file. You can import it back into CloudCompare as
#a Scalar Field to visualize percentage deviations in 3D

# Load all packages -------------------------------------------------------
library(plyr)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)


# Define file reading function and working directories and import data --------------------
## function to import nominal vectors from csv files extracted from CloudCompare
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}
## directory where the files are located; replace "YourInputFolderPath" with the path of the folder where your .csv file is located
## and "YourOutputFolderPath" with the folder path where you want to save the resulting .csv file
working_directory <- "/YourInputFolderPath"  ## Insert here the direction to the example dataset directory
setwd(working_directory)
result_directory <- "/YourOutputFolderPath"  ## Insert here the direction to the example dataset directory

# import nominal vectors
pat = ".csv"
files <- list.files(path = working_directory, pattern = pat,  recursive = TRUE, ignore.case = TRUE)
Scalar_Field <- ldply(files, read_plus)
colnames(Scalar_Field) <- c("X", "Y", "Z", "C2M_signed_distances", "Nx", "Ny", "Nz", "Filename")


# Calculate percentage deviations and export as .csv file to import into CloudCompare as Scalar field ---------------------------------------------------------------
#Calculate deviation vector for each axis separately, Vx, Vy, Vz. V_reverse = (-d * Nx, -d * Ny, -d * Nz)
Scalar_Field$Vx <- -Scalar_Field$C2M_signed_distances*Scalar_Field$Nx
Scalar_Field$Vy <- -Scalar_Field$C2M_signed_distances*Scalar_Field$Ny
Scalar_Field$Vz <- -Scalar_Field$C2M_signed_distances*Scalar_Field$Nz

# Compute original reference coordinates Rx, Ry, Rz
Scalar_Field$Rx <- Scalar_Field$X + Scalar_Field$Vx
Scalar_Field$Ry <- Scalar_Field$Y + Scalar_Field$Vy
Scalar_Field$Rz <- Scalar_Field$Z + Scalar_Field$Vz

# Calculate absolute deviation ADx, ADy, ADz for each coordinate
Scalar_Field$ADx <- abs(Scalar_Field$X - Scalar_Field$Rx)
Scalar_Field$ADy <- abs(Scalar_Field$Y - Scalar_Field$Ry)
Scalar_Field$ADz <- abs(Scalar_Field$Z - Scalar_Field$Rz)

# Calculate percentage deviation for each coordinate: Px, Py, Pz
Scalar_Field$Px <- (Scalar_Field$ADx / abs(Scalar_Field$Rx)) * 100
Scalar_Field$Py <- (Scalar_Field$ADy / abs(Scalar_Field$Ry)) * 100
Scalar_Field$Pz <- (Scalar_Field$ADz / abs(Scalar_Field$Rz)) * 100

# Assign sign based on whether the deviation is positive or negative and generate Percentage Deviations in all axes: PDx, PDy, PDz
Scalar_Field$PDx <- ifelse(Scalar_Field$X > Scalar_Field$Rx, Scalar_Field$Px, -Scalar_Field$Px)
Scalar_Field$PDy <- ifelse(Scalar_Field$Y > Scalar_Field$Ry, Scalar_Field$Py, -Scalar_Field$Py)
Scalar_Field$PDz <- ifelse(Scalar_Field$Z > Scalar_Field$Rz, Scalar_Field$Pz, -Scalar_Field$Pz)

Scalar_Field_subset <- Scalar_Field[, c(1:4, 12:14, 18:20)]

# Export results to CSV
# Select the required columns
export_data <- dplyr::select(Scalar_Field, X, Y, Z, PDx, PDy, PDz)
# Export to CSV
write.csv(export_data, "percentage_deviation_point_cloud.csv", row.names = FALSE)


# Prepare data for plotting by axis ------------------------------------------

# Select the required columns and pivot longer
Scalar_Field_long <- Scalar_Field %>%
  dplyr::select(PDx, PDy, PDz) %>%
  pivot_longer(cols = c(PDx, PDy, PDz),
               names_to = "Axis",
               values_to = "Percentage_Deviation")
Scalar_Field_long$Axis <- ifelse(grepl("PDx", Scalar_Field_long$Axis), "X", Scalar_Field_long$Axis)
Scalar_Field_long$Axis <- ifelse(grepl("PDy", Scalar_Field_long$Axis), "Y", Scalar_Field_long$Axis)
Scalar_Field_long$Axis <- ifelse(grepl("PDz", Scalar_Field_long$Axis), "Z", Scalar_Field_long$Axis)

# Data for X axis
data_X <- Scalar_Field_long %>%
  filter(Axis == "X")

# Data for Y axis
data_Y <- Scalar_Field_long %>%
  filter(Axis == "Y")

# Data for Z axis
data_Z <- Scalar_Field_long %>%
  filter(Axis == "Z")


# Define and filter outliers ----------------------------------------------
# Calculate the quartiles and IQR for the X axis
Q1_X <- quantile(data_X$Percentage_Deviation, 0.25)
Q3_X <- quantile(data_X$Percentage_Deviation, 0.75)
IQR_X <- Q3_X - Q1_X

# Define outlier boundaries
lower_bound_X <- Q1_X - 1.5 * IQR_X
upper_bound_X <- Q3_X + 1.5 * IQR_X

# Filter data to exclude outliers
filtered_data_X <- data_X %>%
  filter(Percentage_Deviation >= lower_bound_X & Percentage_Deviation <= upper_bound_X)
# Calculate the quartiles and IQR for the Y axis
Q1_Y <- quantile(data_Y$Percentage_Deviation, 0.25)
Q3_Y <- quantile(data_Y$Percentage_Deviation, 0.75)
IQR_Y <- Q3_Y - Q1_Y

# Define outlier boundaries
lower_bound_Y <- Q1_Y - 1.5 * IQR_Y
upper_bound_Y <- Q3_Y + 1.5 * IQR_Y

# Filter data to exclude outliers
filtered_data_Y <- data_Y %>%
  filter(Percentage_Deviation >= lower_bound_Y & Percentage_Deviation <= upper_bound_Y)

# Calculate the quartiles and IQR for the Z axis
Q1_Z <- quantile(data_Z$Percentage_Deviation, 0.25)
Q3_Z <- quantile(data_Z$Percentage_Deviation, 0.75)
IQR_Z <- Q3_Z - Q1_Z

# Define outlier boundaries
lower_bound_Z <- Q1_Z - 1.5 * IQR_Z
upper_bound_Z <- Q3_Z + 1.5 * IQR_Z

# Filter data to exclude outliers
filtered_data_Z <- data_Z %>%
  filter(Percentage_Deviation >= lower_bound_Z & Percentage_Deviation <= upper_bound_Z)


# Plot filtered and unfiltered data ---------------------------------------
ggplot(data_X, aes(x = Percentage_Deviation)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.7, color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Percentage Deviations for X Axis",
       x = "Percentage Deviation",
       y = "Frequency")
ggplot(data_Y, aes(x = Percentage_Deviation)) +
  geom_histogram(binwidth = 1, fill = "green", alpha = 0.7, color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Percentage Deviations for Y Axis",
       x = "Percentage Deviation",
       y = "Frequency")

ggplot(data_Z, aes(x = Percentage_Deviation)) +
  geom_histogram(binwidth = 1, fill = "red", alpha = 0.7, color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Percentage Deviations for Z Axis",
       x = "Percentage Deviation",
       y = "Frequency")
####
ggplot(filtered_data_X, aes(x = Percentage_Deviation)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.7, color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Percentage Deviations for X Axis",
       x = "Percentage Deviation",
       y = "Frequency")
ggplot(filtered_data_Y, aes(x = Percentage_Deviation)) +
  geom_histogram(binwidth = 1, fill = "green", alpha = 0.7, color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Percentage Deviations for Y Axis",
       x = "Percentage Deviation",
       y = "Frequency")

ggplot(filtered_data_Z, aes(x = Percentage_Deviation)) +
  geom_histogram(binwidth = 1, fill = "red", alpha = 0.7, color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Percentage Deviations for Z Axis",
       x = "Percentage Deviation",
       y = "Frequency")
####
ggplot(data_X, aes(x = Percentage_Deviation, fill = Axis)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Axis, scales = "free") +
  theme_minimal() +
  labs(title = "Density of Percentage Deviations by Axis",
       x = "Percentage Deviation",
       y = "Density")

ggplot(data_Y, aes(x = Percentage_Deviation, fill = Axis)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Axis, scales = "free") +
  theme_minimal() +
  labs(title = "Density of Percentage Deviations by Axis",
       x = "Percentage Deviation",
       y = "Density")

ggplot(data_Z, aes(x = Percentage_Deviation, fill = Axis)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Axis, scales = "free") +
  theme_minimal() +
  labs(title = "Density of Percentage Deviations by Axis",
       x = "Percentage Deviation",
       y = "Density")
###
ggplot(filtered_data_X, aes(x = Percentage_Deviation, fill = Axis)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Axis, scales = "free") +
  theme_minimal() +
  labs(title = "Density of Percentage Deviations by Axis",
       x = "Percentage Deviation",
       y = "Density")

ggplot(filtered_data_Y, aes(x = Percentage_Deviation, fill = Axis)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Axis, scales = "free") +
  theme_minimal() +
  labs(title = "Density of Percentage Deviations by Axis",
       x = "Percentage Deviation",
       y = "Density")

ggplot(filtered_data_Z, aes(x = Percentage_Deviation, fill = Axis)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Axis, scales = "free") +
  theme_minimal() +
  labs(title = "Density of Percentage Deviations by Axis",
       x = "Percentage Deviation",
       y = "Density")
###

# Average deviations and standard deviations per axis -----------------------------------------
setwd(result_directory)

# import Percentage_Deviation results
pat = ".csv"
files <- list.files(path = result_directory, pattern = pat,  recursive = TRUE, ignore.case = TRUE)
Deviations <- ldply(files, read_plus)
colnames(Deviations) <- c("Percentage_Deviation_X", "Percentage_Deviation_Y", "Percentage_Deviation_Z")

Deviations <- Deviations %>%
  dplyr::select(Percentage_Deviation_X, Percentage_Deviation_Y, Percentage_Deviation_Z)

# Reshape the data into long format
Deviations_long <- Deviations %>%
  pivot_longer(cols = everything(),
               names_to = "Axis",
               values_to = "Percentage_Deviation")
  
# Create violin plot
ggplot(Deviations_long, aes(x = Axis, y = Percentage_Deviation)) +
  geom_violin(fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Percentage Deviations by Axis",
       x = "Axis",
       y = "Percentage Deviation")
# Create boxplot
ggplot(Deviations_long, aes(x = Axis, y = Percentage_Deviation)) +
  geom_boxplot(fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Percentage Deviations by Axis",
       x = "Axis",
       y = "Percentage Deviation")

# Calculate average negative deviation and standard deviation
summary_negative_deviation <- Scalar_Field %>%
  summarise(
    avg_neg_dev_X = mean(pmin(PDx, 0), na.rm = TRUE),
    sd_neg_dev_X = sd(pmin(PDx, 0), na.rm = TRUE),
    avg_neg_dev_Y = mean(pmin(PDy, 0), na.rm = TRUE),
    sd_neg_dev_Y = sd(pmin(PDy, 0), na.rm = TRUE),
    avg_neg_dev_Z = mean(pmin(PDz, 0), na.rm = TRUE),
    sd_neg_dev_Z = sd(pmin(PDz, 0), na.rm = TRUE)
  )

# Calculate average positive deviation and standard deviation
summary_positive_deviation <- Scalar_Field %>%
  summarise(
    avg_pos_dev_X = mean(pmax(PDx, 0), na.rm = TRUE),
    sd_pos_dev_X = sd(pmax(PDx, 0), na.rm = TRUE),
    avg_pos_dev_Y = mean(pmax(PDy, 0), na.rm = TRUE),
    sd_pos_dev_Y = sd(pmax(PDy, 0), na.rm = TRUE),
    avg_pos_dev_Z = mean(pmax(PDz, 0), na.rm = TRUE),
    sd_pos_dev_Z = sd(pmax(PDz, 0), na.rm = TRUE)
  )






