# This script will be used to calculate sliding windows for read depth across samples. 

# This is a second attempt at making this calculation. The other script used was too slow - it was an inefficient while loop written by hand. Upon further research, another group had made a package for this exact algorithm. 


#install.packages("TTR")
library(TTR)
# Get the sample name from the command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

#set the basepath for the file 
base_path = "/home/tdavies/scratch/wgs_gwas_backup/wg_rd_analysis/all_rd_files/"

just_samp = strsplit(sample_name, split = ".txt")[[1]][1]
print(just_samp)
# Set the file name and read in the data
data <- read.table(file = paste(base_path, sample_name, sep = ""), header = FALSE, sep = "\t")

# Calculate sliding window averages at EVERY POSITION 
a = SMA(data$V3, 100)
# Take only the 100th element, and then every 100th element thereafter. 
b = a[seq.int(100L, length(a), 100L)]
# name each window 
win_id = c(1:length(b))

# create table to write 
mean_tab = data.frame(win_id)
mean_tab$plotting_bp = ((mean_tab$win_id)*100)-50 # names the window as the centre of the window
mean_tab$win_mean = b

# Write the means to a new file
output_dir = "/home/tdavies/scratch/wgs_gwas_backup/wg_rd_analysis/all_rd_files/100bp_win_means/"
write.table(mean_tab, file = paste(output_dir,just_samp,"_rd_100bp_win",".txt", sep = ""), row.names = FALSE)