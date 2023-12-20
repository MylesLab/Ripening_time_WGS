# Read depth correlation analysis. 
# The data used here was generated via a sliding window script written by Tommy Davies. The data is the mean read depth from 100bp windows on chromosome 3 - which was identified as containing a signal after examining the entire genome. The windows were then run through a correlation script to test the correlation between mean window read depth and harvest date. 

# This script will be used to look at the correlation results for the read depth data in the WGS-GWAS 97 sample reference panel. 
library(ggplot2)

cor_data = read.table("~/Desktop/DALMSc/kmer-gwas/classic-gwas/read_depth_analysis/100bp_win_Chr03_correlation_results.txt", header = F)
colnames(cor_data) <- c("chroms",	"window",	"plotting_bp",	"rvals",	"pvals")
#what does the dist of r-values look like ?
summary(cor_data$rvals)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.608008 -0.071224  0.002696  0.005988  0.081699  0.745267
summary(-log10(cor_data$pvals))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00000  0.07793  0.29061  0.48566  0.67787 13.93686

# lets look at the max value
max(cor_data$rvals, na.rm = T)
#0.7452668
which(cor_data$rvals == max(cor_data$rvals, na.rm = T))
#row 307217
cor_data[307217,]
# Window location : 30721701-30721800

# Lets look at the min value 
min(cor_data$rvals, na.rm = T)
which(cor_data$rvals == min(cor_data$rvals, na.rm = T))

#row 131036
cor_data[131036,]
# Window location : 13103600-13103700

#4.462747 -- this is the top 0.1% of all p values

par(mfrow=c(2,1),mai = c(0.5, 0.5, 0.1, 0.2))
plot(cor_data$plotting_bp, cor_data$rvals, cex = 0.15, xlab = "Chr03 BP (Mb)", ylab = "r value")
plot(cor_data$plotting_bp, -log10(cor_data$pvals), cex = 0.15, xlab = "Chr03 BP (Mb)", ylab = "-log10(p)")
abline(h=4.46, col="red")
# In order to line up the plots with the GWAS data, I need to extract data from the same window 

corr_zoom_dat = cor_data[306309:307930,]
plot(corr_zoom_dat$plotting_bp, corr_zoom_dat$rvals, cex = 0.2, xlim = c(30630955, 30793029))
abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region

# Zoom in more 

corr_zoom_dat = cor_data[306801:307577,]
corr_full_zoom = plot(corr_zoom_dat$plotting_bp, corr_zoom_dat$rvals, cex = 0.2, xlim = c(30680168, 30757738), xlab = "BP (Chr3)", ylab = "R-value")

abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region

# I would like to overlay the two plots, so I need to import the results from the GWAS data 
test_data = read.table("~/Desktop/DALMSc/kmer-gwas/classic-gwas/results/ref_panel_var_results/harvest_date/harv_date_ref_panel_chr3_77kb_zoom.txt", header = F, sep = '\t')
import_head = c("chr",	"rs",	"ps",	"n_miss",	"allele1",	"allele0",	"af",	"logl_H1",	"l_mle",	"p_lrt")
colnames(test_data)<- import_head
sig_line = 0.05/25769077

gwas_full_zoom = plot(test_data$ps, -log10(test_data$p_lrt), cex = 0.1, main = 'Harvest Date', xlab = "BP (Chr3)", ylab = "-log10(p)")
gwas_full_zoom
abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region

#Now stack
par(mfrow=c(2,1),mai = c(0.5, 0.5, 0.1, 0.2))
plot(test_data$ps, -log10(test_data$p_lrt), cex = 0.1, main = 'Harvest Date', xlab = "BP (Chr3)", ylab = "-log10(p)")
abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region
abline(v = 30738599, col = "orange") # Other NAC
abline(v = 30740479, col = "orange") # Other NAC
plot(corr_zoom_dat$plotting_bp, corr_zoom_dat$rvals, cex = 0.3, xlim = c(30680168, 30757738), xlab = "BP (Chr3)", ylab = "R-value", pch = 15)
abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region
abline(v = 30738599, col = "orange") # Other NAC
abline(v = 30740479, col = "orange") # Other NAC

# Now, lets look at that one strongest correlation window 
strong_win = read.table("~/Desktop/DALMSc/kmer-gwas/classic-gwas/read_depth_analysis/chr03_most_corr_normed_window.txt", header = F)
# Remove columns 1, 2, and 4 as they are metadata, metadata and a duplicate, respectively. 
strong_win = strong_win[,c(-1,-2,-3,-4)]
hist(as.numeric(strong_win[1,]), breaks = 40)
#get pheno data 
harv_dat = read.table("~/Desktop/DALMSc/kmer-gwas/classic-gwas/read_depth_analysis/harvest_date_by_gencoveID.txt", header = T)
dev.off()
plot(as.numeric(strong_win[1,]), harv_dat$date_17_harv, cex = 0.2)
harv_dat$rd = as.numeric(strong_win[1,])
ggplot(harv_dat, aes(x = rd, y = date_17_harv)) + geom_point() + geom_smooth(method="lm", col="red", se=F)+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("Normalized Read Depth")+ylab("Harvest Date")+ggtitle("Chr03: Window 30721701-30721800bp")


# Lets look at the mean depths 
tot_rd = read.delim('~/Desktop/DALMSc/kmer-gwas/classic-gwas/pipeline-stats/total_depth_all_samples.txt', header = F, sep= ' ')
# calculate mean depth by dividing total depth across all positions by the total number of positions in the genome
tot_rd$mean_rd = (tot_rd$V2)/709561391
#extract ony the sample name from the first column and replace the first column with just the sample name 
sample_id =c()
for (i in tot_rd$V1){
  proper_name = strsplit(i, split = '/')[[1]][11]
  sample_id = c(sample_id, proper_name)
}
#set the first column to the sample names
tot_rd$V1 = sample_id

# lets take the read depths at this position, and then divide them by each samples mean read depth and see what happens.  

harv_dat$rd_over_samp_mean = harv_dat$rd/tot_rd$mean_rd
cor.test(harv_dat$date_17_harv, harv_dat$rd)
cor.test(harv_dat$date_17_harv, harv_dat$rd_over_samp_mean)

ggplot(harv_dat, aes(x = rd_over_samp_mean, y = date_17_harv)) + geom_point() + geom_smooth(method="lm", col="red", se=F)+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("Read depth / sample mean read depth")+ylab("Harvest Date")+ggtitle("Window 30721701-30721800bp")

# Plot both single window plots side by side. 

a = ggplot(harv_dat, aes(x = rd, y = date_17_harv)) + geom_point() + geom_smooth(method="lm", col="red", se=F)+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("Read depth")+ylab("Harvest Date")+ggtitle("Window 30721701-30721800bp")
b = ggplot(harv_dat, aes(x = rd_over_samp_mean, y = date_17_harv)) + geom_point() + geom_smooth(method="lm", col="red", se=F)+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("Read depth / sample mean read depth")+ylab("Harvest Date")+ggtitle("Window 30721701-30721800bp")

plot_grid(a,b, labels = c('A', 'B'), label_size = 12)







