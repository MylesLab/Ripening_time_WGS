# LD plots for the GWAS on 76 WGS individuals vs GBS data across entire GBS. 

# The objective of this script is to compare the LD between GBS data (from previous works) and from the WGS data generated for this project. 

# To do this, I will create LD plots for both groups, and compare. I will also place LD plots in context of their respective GWAS plots. 


#load packages
library("readr")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("gplots")
library("viridis")
library("scales")
library("tidyverse")
library("cowplot")
library("gridExtra")
library("grid")
library("genetics")
install.packages("LDheatmap")
library("LDheatmap")
library("snpStats")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("snpStats")
library("devtools")
library("emma")
library("mlmm")
library("data.table")
library("ape")
library("readxl")
# Install the latest release version from CRAN and the
# imported/suggested BioConductor packages with
install.packages("LDheatmap")
source("https://bioconductor.org/biocLite.R")
biocLite(c("snpStats","rtracklayer","GenomicRanges","GenomInfoDb","IRanges"))

# Install the latest development version from GitHub with
devtools::install_github("SFUStatgen/LDheatmap", force = TRUE)

#read hapmap of genetic data.
snp<- read_delim("~/Downloads/nac_region_only_vcf.hmp.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# now filter it down to only contain the SNPs that are in the GBS data, so that you can do an apples to apples comparison 
snp2<- read_delim("~/Downloads/abc_GBS_chr03.hmp.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

snp <- snp[which(snp$pos %in% snp2$pos), ]

# now remove the samples that didnt make the gwas! 
# load in the samples that are in the ref panel 
ref_samples = read.delim("~/Downloads/ref_panel_samples.txt", header = F)
#which have pheno data?
pheno_dat = read.delim("~/Desktop/DALMSc/kmer-gwas/prelim_gwas/Watts_et_al_2021_supp_data copy.csv",header = T, sep = ',')
pheno_dat <- pheno_dat[,c("apple_id", "date_17_harv")]
pheno_dat <- pheno_dat[which(pheno_dat$apple_id %in% ref_samples$V1),]
pheno_dat <- pheno_dat[-which(is.na(pheno_dat$date_17_harv)),]
complete_set <- paste(pheno_dat$apple_id, "_", pheno_dat$apple_id, sep = "")
and_this <- which(colnames(snp) %in% complete_set)
#subset the snp file to include all info cols plus only the samples in the GWAS 
final_subset_index <- c(1:11, and_this)
snp <- snp[,final_subset_index]

hapmap <- snp

hapmap[12:dim(hapmap)[2]]<-apply(hapmap[12:dim(hapmap)[2]],2,as.character)
snp_dis<-as.numeric(as.character(hapmap$pos))
hapmap[hapmap=="A"]<-"A/A"
hapmap[hapmap=="T"]<-"T/T"
hapmap[hapmap=="C"]<-"C/C"
hapmap[hapmap=="G"]<-"G/G"
hapmap[hapmap=="M"]<-"A/C"
hapmap[hapmap=="R"]<-"A/G"
hapmap[hapmap=="W"]<-"A/T"
hapmap[hapmap=="S"]<-"C/G"
hapmap[hapmap=="Y"]<-"C/T"
hapmap[hapmap=="K"]<-"G/T"
snp_infor<-hapmap[,1:11]
snp_data<-hapmap[,12:dim(hapmap)[2]]
snp_data<-t(snp_data)
snp_data<-as.data.frame(snp_data)
names(snp_data)<-as.character(snp_infor$'rs#')
num<-ncol(snp_data)

for(i in 1:num){
  alleles1 <- as.vector(snp_infor[i,2])
  snp_data[,i]<-as.genotype(snp_data[,i])
}

#make heatmap
rgb.palette <- colorRampPalette(rev(brewer.pal(9,"YlGnBu")))
pdf("~/Downloads/WGS_nac_region_heatmap.pdf",width=5,height=5)
WGSHeatmap<-LDheatmap(snp_data,snp_dis,flip=TRUE,color=rgb.palette(20),title="",)
dev.off()

# WGS GWAS LD STATS
mean(WGSHeatmap[["LDmatrix"]], na.rm = T)
#0.4235322
median(WGSHeatmap[["LDmatrix"]], na.rm = T)
#0.217083
#Quantiles 
# 0%        25%        50%        75%       100% 
# 0.01047856 0.07631449 0.21708297 0.83106807 0.99952229 
hist(WGSHeatmap[["LDmatrix"]], na.rm = T, breaks = 100)



#Now compare to the GBS data across the population 
snp2<- read_delim("~/Downloads/abc_GBS_chr03.hmp.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

hapmap <- snp2

hapmap[12:dim(hapmap)[2]]<-apply(hapmap[12:dim(hapmap)[2]],2,as.character)

hapmap[hapmap=="A"]<-"A/A"
hapmap[hapmap=="T"]<-"T/T"
hapmap[hapmap=="C"]<-"C/C"
hapmap[hapmap=="G"]<-"G/G"
hapmap[hapmap=="M"]<-"A/C"
hapmap[hapmap=="R"]<-"A/G"
hapmap[hapmap=="W"]<-"A/T"
hapmap[hapmap=="S"]<-"C/G"
hapmap[hapmap=="Y"]<-"C/T"
hapmap[hapmap=="K"]<-"G/T"
snp2_infor<-hapmap[,1:11]
snp2_dis<-as.numeric(as.character(snp2_infor$pos))
snp2_data<-hapmap[,12:dim(hapmap)[2]]
snp2_data<-t(snp2_data)
snp2_data<-as.data.frame(snp2_data)
names(snp2_data)<-as.character(snp2_infor$'rs#')
num<-ncol(snp2_data)

for(i in 1:num){
  print(i)
  alleles1 <- as.vector(snp2_infor[i,2])
  print(i+1)
  snp2_data[,i]<-as.genotype(snp2_data[,i])
}

rgb.palette <- colorRampPalette(rev(brewer.pal(9,"YlGnBu")))
pdf("~/Downloads/ya_right9.pdf",width=6,height=6)
MyHeatmap<-LDheatmap(snp2_data,snp2_dis,flip=TRUE,color=rgb.palette(20),title="",)
dev.off()
# GBS HEATMAP STATS: 
mean(MyHeatmap[["LDmatrix"]], na.rm = T)
#0.3759633
median(MyHeatmap[["LDmatrix"]], na.rm = T)
#0.08223491
quantile(MyHeatmap[["LDmatrix"]], na.rm = T)
# 0%         25%         50%         75%        100% 
# 0.001428179 0.029184229 0.082234912 0.796132454 0.999240414 
hist(MyHeatmap[["LDmatrix"]], na.rm = T, breaks = 100)



# load in data from GWAS done by Watts on entire population 
gbs_dat = read.delim("~/Downloads/gbs_harv_dat_gwas_res.txt")
#keep only positions in the window we are talking about
gbs_dat = gbs_dat[which(gbs_dat$Pos >= 30680000),]
gbs_dat = gbs_dat[which(gbs_dat$Pos <= 30720000),]
gbs_dat = gbs_dat[which(gbs_dat$Chr == 3),]

# make a ggplot with that data 
pdf("~/Downloads/ya_right3.pdf",width=6,height=6)
gbs_plot = ggplot(gbs_dat,aes(Pos,-log10(p))) +geom_point( alpha = 0.3,size = 2.0, colour = 'black')+theme_bw()+ labs(title = "GBS GWAS", x = "Position on Chromosome 3 (Mb)",y = "-log10(p)")+ theme(plot.margin = unit(c(0,1,0,1),"cm"),panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2))+ annotate("rect", xmin = 30726727, xmax = 30727185, ymin = 0, ymax = 40, fill = "blue", alpha=0.35)+ annotate("rect", xmin = 30696191, xmax = 30698216, ymin = 0, ymax = 40, fill = "red", alpha=0.35)+scale_y_continuous(breaks = c(0,10,20,30,40), limits = c(0,40))+annotate("rect", xmin = 30696191, xmax = 30720000, ymin = 8.00, ymax = 8.05, fill = "red")+scale_x_continuous(limits = c(30695000,30720000 ))
gbs_plot
dev.off()

# plot the WGS GWAS data with the same layout beside it so that a direct comparison LD plot can be made. 

#WGS GWAS plot 
wgs_data = read.table("~/Desktop/DALMSc/kmer-gwas/classic-gwas/results/ref_panel_var_results/harvest_date/harv_date_ref_panel_chr3_77kb_zoom.txt", header = F, sep = '\t')
import_head = c("chr",	"rs",	"ps",	"n_miss",	"allele1",	"allele0",	"af",	"logl_H1",	"l_mle",	"p_lrt")
colnames(wgs_data)<- import_head
min(wgs_data$ps)
max(wgs_data$ps)
range = max(wgs_data$ps) -  min(wgs_data$ps)
range # 77570 77.57 kb

sig_line = 0.05/25769077

# be sure to make the ggplot as wide as possible. To do this, I will zoom in further. The new window for this zoom in plot will be: 
#30680168 - 30727571
# 47.4 kb 

import_head = c("chr",	"rs",	"ps",	"n_miss",	"allele1",	"allele0",	"af",	"logl_H1",	"l_mle",	"p_lrt")
colnames(wgs_data)<- import_head
sig_line = 0.05/25769077
wgs_data$p_lrt = -log10(wgs_data$p_lrt)
#make test data only have snps from the GBS 
wgs_data = wgs_data[which(wgs_data$ps %in% gbs_dat$Pos),]

wgs_gwas_plot = ggplot(wgs_data,aes(ps,(p_lrt))) +geom_point( alpha = 0.3,size = 2, colour = 'black')+theme_bw()+ labs(title = "WGS GWAS", x = "Position on Chromosome 3 (Mb)",y = "-log10(p)")+ theme(plot.margin = unit(c(0,1,0,1),"cm"),panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2))+ annotate("rect", xmin = 30696191, xmax = 30698216, ymin = 0, ymax = 10.5, fill = "red", alpha=0.35)+scale_y_continuous( limits = c(0,11))+scale_x_continuous(limits = c(30695000,30720000 ))
wgs_gwas_plot


pdf("~/Downloads/ya_right4.pdf",width=12,height=6)
plot_grid(gbs_plot,wgs_gwas_plot, nrow = 1, ncol = 2, labels = c("A", "B"))
dev.off()
