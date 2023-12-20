# This script will standardize windows per sample. The mean and standard deviation of all windows needs to be found. Then, the mean will be subtracted from each window and each window will be divided by the SD. 

# Get the sample name from the command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

#set the basepath for the file 
base_path = "/home/tdavies/scratch/wgs_gwas_backup/wg_rd_analysis/all_rd_files/100bp_win_means/"

just_samp = strsplit(sample_name, split = ".txt")[[1]][1]
print(just_samp)
# Set the file name and read in the data
data <- read.table(file = paste(base_path, sample_name, sep = ""), header = FALSE)
# add some column names 
colnames(data) <- c("window", "plotting_bp", "win_mean")

# calculate the genome wide mean read depth across windows 
rd_mean <- mean(data$win_mean)
# calculater the genome wide SD of the read depth across windows
rd_SD <- sd(data$win_mean)

# Now normalize the read depth data 

data$rd_norm = ((data$win_mean)-rd_mean)/rd_SD

# add back in the chromosome numbers 
# number of reps has been hand adjusted based on intermediate files 
chroms = c(rep("Chr00", 527283), rep("Chr01", 326254), rep("Chr02", 375777), rep("Chr03",375240), rep("Chr04", 323018), rep("Chr05", 479524), rep("Chr06", 371372), rep("Chr07",366911), rep("Chr08", 316092), rep("Chr09", 376049), rep("Chr10", 417624), rep("Chr11", 430598), rep("Chr12", 330500), rep("Chr13", 443395), rep("Chr14", 325134), rep("Chr15", 549454), rep("Chr16", 413894), rep("Chr17", 347487))

data = cbind(chroms, data)


# Write the means to a new file
output_dir = "/home/tdavies/scratch/wgs_gwas_backup/wg_rd_analysis/all_rd_files/100bp_win_means/normalized/"
write.table(data, file = paste(output_dir,just_samp,"_rd_100bp_win_norm",".txt", sep = ""), row.names = FALSE)
