# Manhattan and QQ plot creation 
# This script will be used to generate manhattan and QQ plots for GWAS results from the 76 WGS samples. 

# This will take the output file from GEMMA ("*.assoc.txt") and run it through the qqman package to generate a genome wide manhattan plot and a qqplot for the test. 

install.packages("qqman")
library(qqman)

test_dat = read.table("/home/ubuntu/Output/Linear_Mixed_Model/date_17_harv/date_17_harv_20230529_1704591/date_17_harv_mod_sub_wgs_gwas_harv_date_table_ref_panel_vars_filtered.assoc.txt", header = T, sep = '\t')

#Create manhattan plot 

#Signifiance line = 0.05/31242575
sig_line = 0.05/25769077
#change chr0 to R
test_dat[which(test_dat$chr == 0), "chr"] <- 18

png(filename = "Harv_date_manhattan_001.png", width = 1440, height = 480, units = "px")
manhattan(
  test_dat, main = "Harvest Date",
  chr = "chr",
  bp = "ps",
  p = "p_lrt",
  snp = "rs",
  col = c("gray10", "gray60"),
  suggestiveline = -log10(sig_line),
  genomewideline = -log10(sig_line),
  highlight = NULL,
  logp = TRUE,
  annotatePval = NULL,
  annotateTop = TRUE,
  chrlabs = c(c(1:17), "R")
)
dev.off()

# Now create QQplot
png(filename = "Harv_date_QQ_001.png", width = 480, height = 480, units = "px")
p_vector = test_dat$p_lrt
qq(pvector = p_vector)
dev.off()


