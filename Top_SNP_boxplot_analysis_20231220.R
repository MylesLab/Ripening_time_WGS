# Examine the genotypes of the top three hits from the WGS GWAS (chapter 4)
# Top three hits 
#30709586 (T/C SNP), 30709593 (C/T SNP), and 30709888 (T/A SNP)

install.packages("dplyr")
library(dplyr)
library(ggplot2)

genotypes = read.delim("~/Downloads/top_three_hits_genotypes_clean.txt", header = T)

genotypes[genotypes=="0/0"]<-"0"
genotypes[genotypes=="0|0"]<-"0"
genotypes[genotypes=="0/1"]<-"1"
genotypes[genotypes=="0|1"]<-"1"
genotypes[genotypes=="1/0"]<-"1"
genotypes[genotypes=="1|0"]<-"1"
genotypes[genotypes=="1/1"]<-"2"
genotypes[genotypes=="1|1"]<-"2"

samples = as.vector(colnames(genotypes[,5:101]))
samples <- sub("^X", "", samples)

pheno_data = read.csv("~/Desktop/DALMSc/kmer-gwas/prelim_gwas/Watts_et_al_2021_supp_data copy.csv")
pheno_data = pheno_data[,c(2,37)]

#Subset phenodata only for samples we have genotyped. 
pheno_data = pheno_data[pheno_data$apple_id %in% samples, ]
genotypes = t(genotypes)
genotypes = as.data.frame(genotypes)
genotypes = cbind(genotypes, empty_column=NA)
genotypes$empty_column = as.vector(row.names(genotypes))
colnames(genotypes) = c("SNP586","SNP593","SNP888","samp")
genotypes = genotypes[-c(1:4),]
genotypes$samp <- sub("^X", "", genotypes$samp)

# Arrange df2 based on the order of df1$apple_ID
pheno_and_geno <- pheno_data[order(match(pheno_data$apple_id, genotypes$samp)), ]

#combine the tables so that we have both genotype and phenotype data in a single table 
pheno_and_geno <- cbind(pheno_and_geno, genotypes)

ggplot(pheno_and_geno, aes(x=SNP586, y=date_17_harv, fill=SNP586)) + 
  geom_boxplot(notch=F)+theme_bw()+ labs( x = "SNP Identity",y = "Ripening time (Julian Date)")+theme(panel.border = element_blank(), axis.text=element_text(colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2), legend.position = "none")+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+scale_fill_manual(values=c("#696868", "#9e9e9e", "#e6e3e3"))

# I need to know how many homozygotes, heterozygotes are in the groups
length(which(genotypes$SNP586 == 0)) # 48
length(which(genotypes$SNP586 == 1)) # 36
length(which(genotypes$SNP586 == 2)) # 13
length(which(is.na(pheno_and_geno$date_17_harv))) # 21 --- makes sense. 97-21 = 76 which is what the GWAS ran on due to phenotype data missingness. 

