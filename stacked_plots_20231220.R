# Ripening time GWAS and Read depth analysis plotting script 

# The goal is to have the WGS gwas plot, the 100 bp sliding window corr plot, and the pool-seq read depth plot all stacked together for the ripening time NAC window. 

# This script is purely for plotting purposes, there are no formal analyses covered in this script


install.packages("qqman")
library(qqman)
install.packages("ggplot2")
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
install.packages("cowplot")
library(cowplot)

# PART A: WGS GWAS 
# Zoom in manhattan plot for chromosome 3 harvest date GWAS 

#REAL REFERENCE PANEL DATA
test_data = read.table("~/Desktop/DALMSc/kmer-gwas/classic-gwas/results/ref_panel_var_results/harvest_date/harv_date_ref_panel_chr3_77kb_zoom.txt", header = F, sep = '\t')
import_head = c("chr",	"rs",	"ps",	"n_miss",	"allele1",	"allele0",	"af",	"logl_H1",	"l_mle",	"p_lrt")
colnames(test_data)<- import_head
min(test_data$ps)
max(test_data$ps)
range = max(test_data$ps) -  min(test_data$ps)
range # 77570 77.57 kb
test_data$p_lrt = -log10(test_data$p_lrt)
sig_line = 0.05/25769077

plot(test_data$ps, (test_data$p_lrt), cex = 0.1, main = 'Harvest Date', xlab = "BP (Chr3)", ylab = "-log10(p)")
abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region
abline(v = 30738599, col = "orange") # Other NAC
abline(v = 30740479, col = "orange") # Other NAC
abline(v = 30723682, col = "pink") # SV1 DUP
abline(v = 30723947, col = "pink") # SV1 DUP
abline(v = 30744310, col = "purple") # SV2 DEL
abline(v = 30744410, col = "purple") # SV2 DEL
abline(v = 30752219, col = "green") # SV3 DEL
abline(v = 30753916, col = "green") # SV3 DEL

#PART B: 100 BP SLIDING WINDOW CORRELATION PLOT 

cor_data = read.table("~/Desktop/DALMSc/kmer-gwas/classic-gwas/read_depth_analysis/100bp_win_Chr03_correlation_results.txt", header = F)
colnames(cor_data) <- c("chroms",	"window",	"plotting_bp",	"rvals",	"pvals")
# find top 1% of pvals 
#what is the 1% pval cutoff?
all_pvals = sort(cor_data$pvals)
all_pvals = -log10(all_pvals)
quantile(all_pvals, 0.9999)


corr_zoom_dat = cor_data[306801:307577,]
corr_full_zoom = plot(corr_zoom_dat$plotting_bp, -log10(corr_zoom_dat$pvals), cex = 0.2, xlim = c(30680168, 30757738), xlab = "BP (Chr3)", ylab = "-log10(p)")

abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region
abline(v = 30738599, col = "orange") # Other NAC
abline(v = 30740479, col = "orange") # Other NAC

#PART C: POOL-SEQ READ DEPTH BETWEEN POOLS PLOT 

poolseq_nac_RD = read.table('~/Desktop/DALMSc/abc_xpgwas/pool_seq_zoom_data/NAC_region_RD_poolseq.txt')
early_harv = data.frame(poolseq_nac_RD[,c(2,3)])
colnames(early_harv) = c('pos', "rd")
late_harv = data.frame(poolseq_nac_RD[,c(2,4)])
colnames(late_harv) = c('pos', "rd")
sig_line = 0.05/25769077

ggplot(early_harv,aes(pos,rd)) +geom_point( alpha = 0.6)+theme_bw()+ labs( x = "Position on Chromosome 3 (Kb)",y = "Read depth")
pool_seq_rd_plot = ggplot(early_harv,aes(pos,rd)) +geom_point( alpha = 0.3,size = 0.2, colour = 'gray')+geom_point(data=late_harv,colour='black',size = 0.7, alpha = 0.3)+theme_bw()+ labs( x = "Position on Chromosome 3 (Kb)",y = "Read depth")+ theme(plot.margin = unit(c(0,1,0,1),"cm"),panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2))+ geom_hline(yintercept=148, linetype="dashed", color = "dark gray")+ geom_hline(yintercept=158, linetype="dashed", color = "black")+ annotate("rect", xmin = 30726727, xmax = 30727185, ymin = 0, ymax = 400, fill = "blue", alpha=0.35)+ annotate("rect", xmin = 30696191, xmax = 30698216, ymin = 0, ymax = 400, fill = "red", alpha=0.35)+ annotate("rect", xmin = 30719000, xmax = 30724000, ymin = 0, ymax = 400, fill = "orange", alpha=0.35)


pool_seq_rd_plot

#PART Cb - RD plot but standardized scores. 
SD_data = read.delim("~/Downloads/harvest_date_standardized_RD_nac_region.txt", sep = ' ', header = F) 
SD_data2 = SD_data[which(SD_data$V2 <= 30704000),]
SD_data2$V5 = (SD_data2$V5*(1-SD_data2$V5))
plot(SD_data2$V2, SD_data2$V5, cex = 0.1)
plot(SD_data$V2, SD_data$V5, cex = 0.1)
SD_data1 = SD_data[which(SD_data$V2 > 30704000),]
all_SD_data = rbind(SD_data2, SD_data1)
plot(all_SD_data$V2, all_SD_data$V5, cex = 0.2)

pool_seq_rd_sd_plot = ggplot(all_SD_data,aes(V2,V5)) +geom_point(colour='black',size = 0.7, alpha = 0.3)+theme_bw()+ labs( x = "Position on Chromosome 3 (Kb)",y = "Absolute difference in scaled RD")+ theme(plot.margin = unit(c(0,1,0,1),"cm"),panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2))+ geom_hline(yintercept=0.04, linetype="dashed", color = "dark gray")+ annotate("rect", xmin = 30726727, xmax = 30727185, ymin = 0, ymax = 0.3, fill = "blue", alpha=0.35)+ annotate("rect", xmin = 30696191, xmax = 30698216, ymin = 0, ymax = 0.3, fill = "red", alpha=0.35)+ annotate("rect", xmin = 30719000, xmax = 30724000, ymin = 0, ymax = 0.3, fill = "orange", alpha=0.35)+ geom_hline(yintercept=0.185, linetype="dashed", color = "red")
pool_seq_rd_sd_plot

 #PART D 
# Examine the GWAS pvalues for the kmerGWAS 

#read in the kmers with pvals
install.packages("qqman")
library(qqman)
kmer_pvals = read.table("~/Downloads/blast_scores_with_pvals.txt")
colnames(kmer_pvals) = c("kmer",  "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "pval")

#now subset for only CHR03
kmer_pvals = kmer_pvals[which(kmer_pvals$sseqid =="Chr03"),]
kmer_pvals$pval = -log10(kmer_pvals$pval)
plot(kmer_pvals$send, kmer_pvals$pval)
# window of examination for other gwas plots:

#30680168 .. 30757738

#subset based on $send for kmers in the window
kmer_pvals = kmer_pvals[which(kmer_pvals$send >= 30680168),]
kmer_pvals = kmer_pvals[which(kmer_pvals$send <= 30757738),]
plot(kmer_pvals$send, kmer_pvals$pval)
plot(kmer_pvals$send, (kmer_pvals$pval), cex = 0.3, pch = 15, main = 'Harvest Date', xlab = "BP (Chr3)", ylab = "-log10(p)", xlim = c(30680168, 30757738))
abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region
abline(v = 30738599, col = "orange") # Other NAC
abline(v = 30740479, col = "orange") # Other NAC



# Now stack all plots 
par(mfrow=c(3,1),mai = c(0.5, 0.5, 0.1, 0.2))

#A
part_A = ggplot(test_data,aes(ps,(p_lrt))) +geom_point( alpha = 1.0,size = 0.2, colour = 'black')+theme_bw()+ labs( x = "Position on Chromosome 3 (Kb)",y = "-log10(p)")+ theme(plot.margin = unit(c(0,1,0,1),"cm"),panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2))+ annotate("rect", xmin = 30726727, xmax = 30727185, ymin = 0, ymax = 11, fill = "blue", alpha=0.35)+ annotate("rect", xmin = 30696191, xmax = 30698216, ymin = 0, ymax = 11, fill = "red", alpha=0.35)+ annotate("rect", xmin = 30719000, xmax = 30724000, ymin = 0, ymax = 11, fill = "orange", alpha=0.35)+scale_y_continuous(breaks = c(0,2,4,6,8,10,12), limits = c(0,11))+annotate("rect", xmin = 30700000, xmax = 30720000, ymin = 8.70, ymax = 8.72, fill = "red")

#B

part_B = ggplot(corr_zoom_dat,aes(plotting_bp,-log10(pvals))) +geom_point( alpha = 1.0,size = 0.4, colour = 'black')+theme_bw()+ labs( x = "Position on Chromosome 3 (Kb)",y = "-log10(p)")+ theme(plot.margin = unit(c(0,1,0,1),"cm"),panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2))+ annotate("rect", xmin = 30726727, xmax = 30727185, ymin = 0, ymax = 14, fill = "blue", alpha=0.35)+ annotate("rect", xmin = 30696191, xmax = 30698216, ymin = 0, ymax = 14, fill = "red", alpha=0.35)+ annotate("rect", xmin = 30719000, xmax = 30724000, ymin = 0, ymax = 14, fill = "orange", alpha=0.35)+scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14), limits = c(0,14))+annotate("rect", xmin = 30700000, xmax = 30757650, ymin = 9.47, ymax = 9.49, fill = "red")
part_B

#C
pool_seq_rd_plot


#D
part_D = ggplot(kmer_pvals, aes(send, pval)) +geom_point( alpha = 1.0,size = 0.4, colour = 'black')+theme_bw()+ labs( x = "Position on Chromosome 3 (Kb)",y = "-log10(p)")+ theme(plot.margin = unit(c(0,1,0,1),"cm"),panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2)) +xlim(30680168, 30757738)+ annotate("rect", xmin = 30726727, xmax = 30727185, ymin = 9, ymax = 12, fill = "blue", alpha=0.35)+ annotate("rect", xmin = 30696191, xmax = 30698216, ymin = 9, ymax = 12, fill = "red", alpha=0.35)+ annotate("rect", xmin = 30738599, xmax = 30740479, ymin = 9, ymax = 12, fill = "orange", alpha=0.35)+scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14), limits = c(9,12))

part_D







plot_grid(part_A, part_B, pool_seq_rd_plot, part_D, nrow = 4, align = "v", labels = c("A", "B", "C", "D"))


plot_grid(part_A, part_B, pool_seq_rd_sd_plot, nrow = 3, align = "v", labels = c("A", "B", "C"))


# Some additional plots for supp materuials # 

# S1 - Normalized scores for pool-seq read depths at the NAC window
stand_rds_poolsed = read.table('~/Desktop/DALMSc/abc_xpgwas/test_data/harvest_date_read_depth_standardized_scores_nac_window.txt')
colnames(stand_rds_poolsed) = c("CHR", "BP", "P1", "P2", "DIFF", "SNP", "P1scaled", "P2scaled")

# red line indicates the top 1% of standardized difference values 
ggplot(stand_rds_poolsed,aes(BP,DIFF)) +geom_point( alpha = 1.0,size = 0.2, colour = 'black')+theme_bw()+ labs( x = "Position on Chromosome 3 (Kb)",y = "|Standardized difference|")+ theme(plot.margin = unit(c(0,1,0,1),"cm"),panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2))+ annotate("rect", xmin = 30726727, xmax = 30727185, ymin = 0, ymax = 0.4, fill = "blue", alpha=0.35)+ annotate("rect", xmin = 30696191, xmax = 30698216, ymin = 0, ymax = 0.4, fill = "red", alpha=0.35)+ annotate("rect", xmin = 30738599, xmax = 30740479, ymin = 0, ymax = 0.4, fill = "orange", alpha=0.35)+annotate("rect", xmin = 30700000, xmax = 30757650, ymin = 0.185, ymax = 0.187, fill = "red")          
       