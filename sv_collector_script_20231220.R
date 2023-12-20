# This is a script modified from code written by Sophie Watts (sophia.watts@dal.ca) which takes the SV collector data and plots all samples from the ABC in PC space, and then highlights the samples chosen for WGS. 

library(tidyverse)

#look at outputs from SV collector

greedy_output <- read_delim("/Users/sophie/Documents/myles_lab/sv_collector/tdavies/abc_combined_maf001_harv_maf001.recode.vcf.greedy", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

topn_output <- read_delim("/Users/sophie/Documents/myles_lab/sv_collector/abc_combined_maf001_harv_maf001.recode.vcf.topN", 
              delim = "\t", escape_double = FALSE, 
              trim_ws = TRUE)

geno_pheno_meta_data <- read_csv("/Users/sophie/Documents/myles_lab/abc-gwas/outputs/geno_pheno_meta_data.csv")

greedy_pheno <- left_join(greedy_output, geno_pheno_meta_data, by = c("Sample" = "apple_id"))
greedy_pheno <- greedy_pheno %>% select(1:12, 19, 24, 107, 33:49, 51)

topn_pheno <- left_join(topn_output, geno_pheno_meta_data, by = c("Sample" = "apple_id"))

table <- greedy_pheno %>% select(1:10)
table$pheno <- rownames(table)
table <- table %>% rename(count = table)
table<- table %>% select(pheno, count)

write_csv(table, file = "abc_sv_collector_greedy.csv", na = "NA", quote = "none")

write_csv(table, file = "phenotype_counts.csv", na = "NA", quote = "none")


###########################################################################################
#make pca plot with the 107 apple ids that sv collector picked.

#load output form sv collector
output_collector_107 <- read_table("Documents/myles_lab/sv_collector/abc_combined_maf001_harv_maf001_107.recode.vcf.greedy")

sample <- output_collector_107 %>% select("Sample")
sample <- sample %>% mutate(sv_samples = "YES")

#load tassel pca output
abc_pcs <- read_delim("Documents/myles_lab/abc-gwas/data/abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_pruned1.txt", 
              delim = "\t", escape_double = FALSE, 
              trim_ws = TRUE, skip = 2)

#load phenotype data
geno_pheno_meta_data <- read_csv("Documents/myles_lab/abc-gwas/outputs/geno_pheno_meta_data.csv")
#join
output_collector_107 <- left_join(output_collector_107, geno_pheno_meta_data, by = c("Sample" = "apple_id"))

sample_names <- output_collector_107 %>% select("Sample", "PLANTID", "#SVs", "species", "country")

write_csv(sample_names, file = "sv_collector_samples.csv", na = "NA", quote = "none")


abc_pcs <- abc_pcs %>% select(Taxa, PC1, PC2)

#make a column that says whether it was picked by the sv_collector
abc_pcs <- left_join(abc_pcs, sample, by = c("Taxa" = "Sample"))

#Make NAs a value
abc_pcs <- abc_pcs %>% mutate(sv_samples = replace_na(sv_samples, "NO"))


#remove rows that are NA for harvest date
yes_samples <- abc_pcs %>% filter(sv_samples == "YES") 

no_samples <- abc_pcs %>% filter(sv_samples == "NO")

data <- rbind(no_samples, yes_samples)

pca_sv_collector <- ggplot(data, aes(x=PC1, y=PC2))+
  geom_point(aes(colour=sv_samples),size=3, stroke=0, alpha=0.8)+
  theme_bw()+
  coord_fixed()+
  scale_color_manual(name = "ABC Samples Selected for GWAS", values=c("#A4A4A4", "#eba236"))+
  labs(x=paste0("PC1 (5.5%)"), y=paste0("PC2 (3.1%)"))+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))
ggsave("pca_sv_collector.pdf", plot = pca_sv_collector)

