# MLMM on ripening time - Prelim data
# This script was made to run a preliminary GWAS to ensure that the experiment will yeild enough statistical power to detect a large effect locus on Chromosome 3. 

# Description: This file runs the GWAS analysis ripening time data

# This file is to be run on Compute Canada Beluga Servers. The file will also produce manhattan plots. All input to this file were generated in other scripts or on the command line on the Beluga cluster. 

# load the command line arguments
args = commandArgs(trailingOnly = TRUE)
pheno_arg = args[1]

# the base directory where all the relevant files and directories are
BASE_DIR = "/home/tdavies/projects/def-smyles/myles_lab/wgs-gwas-prelim/"

# create the output directory for given phenotype
dir.create(file.path(BASE_DIR, "autoplot"), recursive=TRUE)


######################################
## TABLE OF CONTENTS 	            ##
######################################
## 1. LIBRARY IMPORTS
## 2. KINSHIP MATRIX FORMATTING
## 3. PHENOTYPE MATRIX FORMATTING
## 4. GENOTYPE MATRIX FORMATTING
## 5. PCA MATRIX FORMATTING
## 6. MAPPING MATRIX FORMATTING
## 7. RUNNING GWAS
## 8. WRITING GWAS OUTPUT
## 9. EXPORTING PLOTS & TABLES

########################
## 1. LIBRARY IMPORTS ##
#########################

# loading the required libraries
library(tidyverse)
library(emma)
library(mlmm)
library(qqman)
library(scales)
library(xlsx)

##################################
## 2. KINSHIP MATRIX FORMATTING ##
##################################

# read the kinship matrix
kinship <- read_delim(
  paste0(BASE_DIR, "kmer-prelim-kinship.txt"), 
  "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE, skip = 3
)
kinship <- as.data.frame(kinship)
order <- kinship$X1
rownames(kinship) <- kinship[,1]
kinship <- kinship[,2:ncol(kinship)]
colnames(kinship) <- rownames(kinship)
kinship <- as.matrix(kinship)

####################################
## 3. PHENOTYPE MATRIX FORMATTING ##
####################################

# read the phenotype matrix
phenotype <- read_delim(
  paste0(
  	paste0(BASE_DIR,"kmer-prelim-pheno-for-R-gwas"),
  	".txt"
  ), 
  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE
)
phenotype <- phenotype[match(order, phenotype$X1),]
pheno_vec <- as.numeric(phenotype$X2)
names(pheno_vec) <- phenotype$X1


###################################
## 4. GENOTYPE MATRIX FORMATTING ##
###################################

# read the genotype matrix
tictoc::tic()
geno_dat <- read.table(
  paste0(BASE_DIR,"prelim-kmer-gwas-snps-MAF-filtered.raw"),
  header=T
)
tictoc::toc()
geno_dat <- geno_dat[match(order, geno_dat$FID),]
rownames(geno_dat) <- geno_dat[,1]
geno_dat <- geno_dat[,7:ncol(geno_dat)]

# save column names from genotype table.
snp_names <- colnames(geno_dat)
geno_dat <- as.matrix(geno_dat)

# need to save the snp names from the raw file that have the allele appended to 
# the end and then rename the snps in the map file to match the ones in the raw 
# file. First arrange snp names into df.
snp_names_tab <- as.data.frame(snp_names)
snp_names_tab <- snp_names_tab %>%
  rename(new_snp = snp_names) %>%
  mutate_all(~as.character(.)) %>%
  mutate(trimmed = str_sub(new_snp, 2, -3))

##############################
## 5. PCA MATRIX FORMATTING ##
##############################


# read the PCA matrix
pcs <- read_delim(
  paste0(BASE_DIR,"prelim-kmer-gwas-snps-MAF-filtered-pca-tassel1.txt"), 
  "\t", escape_double = FALSE, trim_ws = TRUE, skip = 2)
pcs <- as.data.frame(pcs)
pcs <- pcs[match(order, pcs$Taxa),]
rownames(pcs) <- pcs[,1]
pcs <- pcs[,2:ncol(pcs)]
pcs <- as.matrix(pcs)

##################################
## 6. MAPPING MATRIX FORMATTING ##
##################################

# read the SNPs info
map <- read_delim(
  paste0(BASE_DIR,"prelim-kmer-gwas-snps-MAF-filtered.map"),
  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE
)
map <- map %>%
  select(X1, X2, X4) %>%
  rename(SNP=X2, Chr=X1, Pos=X4) %>%
  mutate(Chr = replace(Chr, Chr == 0, 18)) %>%
  select(SNP, Chr, Pos)

map <- left_join(snp_names_tab, map, by = c("trimmed"="SNP"))

map <- map %>%
  rename(SNP = new_snp) %>%
  select(SNP, Chr, Pos)

map <- as.data.frame(map)

#####################
## 7. RUNNING GWAS ##
#####################

gwas_output <- mlmm_cof(Y=pheno_vec, X=geno_dat, cofs=pcs, K=kinship,
                        nbchunks=10, maxsteps=10)


############################
## 8. WRITING GWAS OUTPUT ##
############################

#################################
## 9. EXPORTING PLOTS & TABLES ##
#################################

# exporting the standard Manhattan plot
jpeg(
  file = paste0(BASE_DIR,"prelim","_std_manhat.jpeg"), 
  width = 800, height = 400
)
plot_fwd_GWAS(
  gwas_output, 1, map,1,main="Standard MLM", 
  abline(h=(-log10(0.05/250579)), col="black", lty=2)
)
dev.off()

# exporting the standard p-values table
std_pvals <- gwas_output[["pval_step"]][[1]][["out"]]
write.table(
  std_pvals, 
  file = paste0(BASE_DIR,"prelim","_std_pvals.csv"), 
  quote = F, sep = ",", row.names = F
)

# exporting the standard QQ plot
jpeg(
  file = paste0(BASE_DIR,"prelim","_std_qq.jpeg"), 
  width = 500, height = 500
)
qqplot <- qq(std_pvals[,2], main = "Standard MLM")
dev.off()

# exporting the MLMM Manhattan plot
jpeg(
  file = paste0(BASE_DIR,"prelim","_mlmm_manhat.jpeg"), 
  width = 800, height = 400
)
plot_opt_GWAS(
  gwas_output, 'extBIC', map, 1, main="optimal MLMM", 
  abline(h=(-log10(0.05/250579)), col="black", lty=2)
)
dev.off()

# exporting the MLMM p-values file
mlmm_pvals <- gwas_output[["opt_extBIC"]][["out"]]
write.table(
  mlmm_pvals, 
  file = paste0(BASE_DIR,"prelim","_mlmm_pvals.csv"), 
  quote = F, sep = ",", row.names = F
)

# exporting the MLMM QQ plot
jpeg(
  file = paste0(BASE_DIR,"prelim","_mlmm_qq.jpeg"), 
  width = 500, height = 500
)
qqplot <- qq(mlmm_pvals[,2], main = "optimal MLMM")
dev.off()


# exporting the SNP variation
var <- as.data.frame(gwas_output[["RSSout"]])
write.table(
  var, 
  file = paste0(BASE_DIR,"prelim","_snp_variation.csv"), 
  quote = F, sep = ",", row.names = F
)

