# This script is mean to plot the genome wide correlation results 

# Manhattan and QQ plot creation 
# This script will be used to generate manhattan and QQ plots for GWAS results from the WGS reference panel. 

# The manhattan function needed slight adjusting in order to plot correlation data rather than GWAS result data.

install.packages("qqman")
library("qqman")

# Input data 
## test data 
#test_dat = read.table("~/Desktop/DALMSc/kmer-gwas/classic-gwas/read_depth_analysis/sample_norm_corr.txt", header = F, sep = '\t')




#real data 
test_dat = read.table("/home/tdavies/projects/def-smyles/tdavies/wgs-gwas-pipeline/extract_RD_from_bams/rd_scripts_and_results/100bp_win_whole_genome_correlation_results.txt", header = T, sep = '\t')


colnames(test_dat)<- c("chroms",	"window",	"plotting_bp",	"rvals",	"pvals")

test_dat$chroms = gsub("Chr", "", test_dat$chroms)
test_dat$chroms = as.numeric(test_dat$chroms)
test_dat[which(test_dat$chroms == 00), "chroms"] <- 18
test_dat$rvals = (test_dat$rvals)+1
#Create manhattan plot 

#Signifiance line = 0.05/7095607
sig_line = 0.05/7095607
#change chr0 to R


png(filename = "harvest_date_norm_rd_full_genome_manhattans.png", width = 1440, height = 480, units = "px")
par(mfrow=c(2,1),mai = c(0.5, 0.5, 0.1, 0.2))
manhattan(
  test_dat, main = "Harvest Date",
  chr = "chroms",
  bp = "plotting_bp",
  p = "pvals",
  snp = "window",
  col = c("gray10", "gray60"),
  suggestiveline = -log10(sig_line),
  genomewideline = -log10(sig_line),
  highlight = NULL,
  logp = TRUE,
  annotatePval = NULL,
  annotateTop = TRUE,
  chrlabs = c(c(1:17), "R")
)

manhattan(
  test_dat, main = NULL,
  chr = "chroms",
  bp = "plotting_bp",
  p = "rvals",
  snp = "window",
  col = c("gray10", "gray60"),
  suggestiveline = -log10(sig_line),
  genomewideline = -log10(sig_line),
  highlight = NULL,
  logp = F,
  annotatePval = NULL,
  annotateTop = TRUE,
  chrlabs = c(c(1:17), "R")
)
dev.off()


