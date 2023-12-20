# This script will be used to run a correlation test between the mean read depts of 100bp windows  all samples in the WGS-GWAS. 

# Here, a pearson correlation test is used to correlate each 100bp window across the genome to the ripening time data. In this script, I correlate the data from each 100bp window (from all samples) to the phenotype data from all samples. The correlation being made is between the rd values at each window to the phenotype data. 

# Read in the read depth data 
data = read.table("/home/tdavies/scratch/wgs_gwas_backup/wg_rd_analysis/all_rd_files/100bp_win_means/normalized/all_samps_full_genome_normalized.txt", header = T)
#data = read.table("~/Downloads/norm_test.txt", header = T)
proper_header = read.table("/home/tdavies/scratch/wgs_gwas_backup/wg_rd_analysis/all_rd_files/100bp_win_means/normalized/sample_list.txt", header = F)
proper_header = proper_header$V1[1:97]
# rename the header of the file
proper_header = c("chroms", "window", "plotting_bp", "win_mean", proper_header)
colnames(data) <- proper_header

# Read in the phenotype data
# NOTE: THE PHENOTYPE DATA HERE WILL BE MATCHED TO GENCOVE IDs. This is the opposite of the GWAS analysis. This is for ease of analysis. 

harv_dat_df = read.table("/home/tdavies/scratch/wgs_gwas_backup/wg_rd_analysis/all_rd_files/100bp_win_means/normalized/harvest_date_by_gencoveID.txt", header = T)
# columns of the real data need to be ordered first, then I need to run the apple not on file head but the real data frame 
# Create a vector with the column names in the right order 
header_copy <- c("chroms", "window", "plotting_bp", "win_mean", harv_dat_df$gencove_id)
# make sure the harv data and the rd data is ordered correctly 
header_copy == colnames(data)
#yes
colnames(data)[5:101] == harv_dat_df$gencove_id
#yes 
# They already match 
print("this is line 25")
# now run a corr.test for each row of the table with the phenotype data

huge_list = apply(data[,5:101], MARGIN = 1, function(x) {cor.test(x, harv_dat_df$date_17_harv, method = "pearson")})
print("this is line 29")
# extract required data from the huge list 
pvalues <- unlist(lapply(huge_list, function(x) x[3]))
r_vals <- unlist(lapply(huge_list, function(x) x[4]))

print("this is line 34")
# take the first two columns from the RD data table, and append this data so that we can plot it later. 
correlation_df = data[,1:3]
correlation_df$rvals = r_vals
correlation_df$pvals = pvalues
print("this is line 39")
# finally, create a column that represents a BP that is roughly in the middle of the window

print("this is line 42")
#write this table to file. 

write.table(correlation_df, "/home/tdavies/scratch/wgs_gwas_backup/wg_rd_analysis/all_rd_files/100bp_win_means/normalized/corr_results/100bp_win_whole_genome_correlation_results.txt", sep = '\t', quote = F, row.names = F)



