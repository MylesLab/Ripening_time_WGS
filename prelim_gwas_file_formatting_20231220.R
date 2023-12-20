# This file was written to create the appropriate files for running a GWAS with a subset of samples from the entire ABC. The samples were chosen using SVcollector. 

#First, make a text file with two columns that are the apple IDs (both columns are the same)

apple_IDs = read.csv("~/Desktop/DALMSc/kmer-gwas/prelim_gwas/20220624_gencove_20x_apple_ids.csv")

apple_id_tab = as.data.frame(apple_IDs)
apple_id_tab[,2] = apple_IDs
write.table(apple_id_tab, file = "~/Desktop/DALMSc/kmer-gwas/prelim_gwas/kmer-prelim-apple-ids.txt", sep = "\t", row.names = F,col.names = F, quote = F)

#Second, make a table where the name of the first column is <Trait> and the second is harvest_date (AKA Ripening time)

# get harv date date 

pheno_data = read.csv("~/Desktop/DALMSc/kmer-gwas/prelim_gwas/Watts_et_al_2021_supp_data copy.csv")
IDs_to_keep = as.vector(apple_IDs$x)

harv_dat = pheno_data[which(pheno_data$apple_id %in% IDs_to_keep),]

harv_dat = harv_dat[c("apple_id", "date_17_harv")]

colnames(harv_dat) <- c("<Trait>", "harvest_date")
write.table(harv_dat, file = "~/Desktop/DALMSc/kmer-gwas/prelim_gwas/kmer-prelim-harv_dates.txt", sep = "\t", row.names = F,col.names = T, quote = F)
#remove NA's
harv_no_nas = harv_dat[-which(is.na(harv_dat$harvest_date)),]
write.table(harv_no_nas, file = "~/Desktop/DALMSc/kmer-gwas/prelim_gwas/kmer-prelim-harv_dates_no_NA.txt", sep = "\t", row.names = F,col.names = T, quote = F)

#Look at the HD distribution 
# the entire ABC dist: 
hist(pheno_data$date_17_harv, breaks = 50)
max(pheno_data$date_17_harv, na.rm = T)
#290.2631
min(pheno_data$date_17_harv, na.rm = T)
#224.9375

# subset for WGS
hist(harv_no_nas$harvest_date, breaks = 50)
max(harv_no_nas$harvest_date, na.rm = T)
#285.156
min(harv_no_nas$harvest_date, na.rm = T)
#225.0926

pdf('~/Desktop/DALMSc/kmer-gwas/prelim_gwas/harvest_dates_comparison.pdf')
par(mfrow=c(2,1))
hist(pheno_data$date_17_harv, breaks = 50, main = "Harvest Dates: Entire ABC", xlab = "Harvest date (julian day)", xlim = c(220,300))
hist(harv_no_nas$harvest_date, breaks = 50, main = "Harvest Dates: WGS panel", xlab = "Harvest date (julian day)", xlim = c(220,300))
dev.off()



write.table(harv_no_nas, file = "~/Desktop/DALMSc/kmer-gwas/prelim_gwas/kmer-prelim-harv_dates_no_NA_spc_sep.txt", sep = " ", row.names = F,col.names = T, quote = F)
# Therefor remake the ID file
ids_no_nas = harv_no_nas$`<Trait>`
apple_ids_no_nas = as.data.frame(ids_no_nas)
apple_ids_no_nas[,2] = apple_ids_no_nas
write.table(apple_ids_no_nas, file = "~/Desktop/DALMSc/kmer-gwas/prelim_gwas/kmer-prelim-apple-ids_no_NAs.txt", sep = "\t", row.names = F,col.names = F, quote = F)

