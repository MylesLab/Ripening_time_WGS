# This script will be used to calculate sliding windows for read depth across samples. 

# For this script, I will apply a non-overlapping sliding window alogrythm. I will use 100bp windows. 

# The output from this script is the mean value from each of the 100bp windows across the genome of the apple. 

# Get the sample name from the command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

#set the basepath for the file 
base_path = "/home/tdavies/scratch/wgs_gwas_backup/wg_rd_analysis/all_rd_files/"

just_samp = strsplit(sample_name, split = ".txt")[[1]][1]

# Set the file name and read in the data
data <- read.table(file = paste(base_path, sample_name, sep = ""), header = FALSE, sep = "\t")

# Define the window size and step size
window_size <- 100
step_size <- 100

# Create a function to calculate the mean of a window
mean_window <- function(data, start_index, end_index) {
  mean(data[start_index:end_index, 3])
}

# Initialize variables
start_index <- 1
end_index <- window_size
window_num <- 1
means <- data.frame(Window_ID = numeric(), Window_loc = character(), Mean_Values = numeric())

# Loop through the data and calculate the means for each window
while (end_index <= nrow(data)) {
  mean_value <- mean_window(data, start_index, end_index)
  means[window_num, ] <- c(window_num, paste(start_index, "-",end_index, sep = ""), mean_value)
  #Print progress tracker 
  print(start_index)
  # Move the window
  start_index <- start_index + step_size
  end_index <- end_index + step_size
  window_num <- window_num + 1
}

colnames(means) <- c("Window_ID", "Window_loc", sample_name)

# Write the means to a new file
output_dir = "/home/tdavies/scratch/wgs_gwas_backup/wg_rd_analysis/all_rd_files/100bp_win_means"
write.table(means, file = paste(output_dir,just_samp,"_rd_100bp_win",".txt", sep = ""), row.names = FALSE)
