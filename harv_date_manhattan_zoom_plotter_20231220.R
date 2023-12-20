# This file will zoom in on GWAS result data for the ripening time GWAS. 
# I will use ggplot to look at very minute regions of the genome, in order to generate zoomed-in manhattan plots of the results from the GWAS. 

# Zoom in manhattan plot for chromosome 3 ripening time GWAS 

install.packages("qqman")
library(qqman)


#400 kb window
test_data = read.table("~/Downloads/harv_date_chr3_peak_zoom.txt", header = F, sep = '\t')
import_head = c("chr",	"rs",	"ps",	"n_miss",	"allele1",	"allele0",	"af",	"beta",	"se",	"p_score")
colnames(test_data)<- import_head
sig_line = 0.05/25769077

plot(test_data$ps, -log10(test_data$p_score), cex = 0.1, main = 'Harvest Date', xlab = "BP (Chr3)", ylab = "-log10(p)")
abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = , col = "red")


#162kb window
test_data = read.table("~/Downloads/harv_date_chr3_162kb_zoom.txt", header = F, sep = '\t')
import_head = c("chr",	"rs",	"ps",	"n_miss",	"allele1",	"allele0",	"af",	"beta",	"se",	"p_score")
colnames(test_data)<- import_head
sig_line = 0.05/25769077

plot(test_data$ps, -log10(test_data$p_score), cex = 0.1, main = 'Harvest Date', xlab = "BP (Chr3)", ylab = "-log10(p)")
abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region

#REAL REFERENCE PANEL DATA
test_data = read.table("~/Desktop/DALMSc/kmer-gwas/classic-gwas/results/ref_panel_var_results/harvest_date/harv_date_ref_panel_chr3_77kb_zoom.txt", header = F, sep = '\t')
import_head = c("chr",	"rs",	"ps",	"n_miss",	"allele1",	"allele0",	"af",	"logl_H1",	"l_mle",	"p_lrt")
min(test_data$ps)
max(test_data$ps)
range = max(test_data$ps) -  min(test_data$ps)
range # 77570 77.57 kb
colnames(test_data)<- import_head
sig_line = 0.05/25769077

plot(test_data$ps, -log10(test_data$p_lrt), cex = 0.1, main = 'Harvest Date', xlab = "BP (Chr3)", ylab = "-log10(p)")
abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region
abline(v = 30738599, col = "orange") # Other NAC
abline(v = 30740479, col = "orange") # Other NAC

abline(v = 30709586, col = "yellow") # top hit 1
abline(v = 30709888, col = "orange") # top hit 2
abline(v = 30709593, col = "purple") # top hit 3

