# CRISPR Sites within the Harvest Date GWAS peak

# The objective of this script is to examine the peak region around NAC18.1 for CRISPR sites within the apple genome. 

# CRISPRscan was used in 5 batches to efficiently scan the entire peak region for CRISPR sites. 

# load in the 5 regions that make up this window individually. 
crispr_sites1 = read.delim("~/Downloads/crisprscan_undefined1.tsv", sep = "\t")
crispr_sites1$start = crispr_sites1$start+30680020
crispr_sites2 = read.delim("~/Downloads/crisprscan_undefined2.tsv", sep = "\t")
crispr_sites2$start = crispr_sites2$start+30690020
crispr_sites3 = read.delim("~/Downloads/crisprscan_undefined3.tsv", sep = "\t")
crispr_sites3$start = crispr_sites3$start+30700020
crispr_sites4 = read.delim("~/Downloads/crisprscan_undefined4.tsv", sep = "\t")
crispr_sites4$start = crispr_sites4$start+30710020
crispr_sites5 = read.delim("~/Downloads/crisprscan_undefined5.tsv", sep = "\t")
crispr_sites5$start = crispr_sites5$start+30720020
# combine all files 
crispr_sites <- rbind(crispr_sites1, crispr_sites2, crispr_sites3, crispr_sites4, crispr_sites5)
# There are 4054 total editing sites 
length(which(crispr_sites$score_crisprscan >= 40 ))
#2008 medium or higher
length(which(crispr_sites$score_crisprscan >= 60 ))
#407 high quality

plot(crispr_sites$start, crispr_sites$score_crisprscan, cex = 0.2)
abline(v = 30696191, col = "red")
abline(v = 30698216, col = "red")
abline(v = 30727185, col = "blue") # Peach SV homologos region 
abline(v = 30726727, col = "blue") # Peach SV homologos region

hist(crispr_sites$score_crisprscan, breaks = 100)
hist(crispr_sites$start, breaks = 500)

# How many sites are within the coding region of NAC18.1 which is from 30696191-30698216



length(which(crispr_sites$start >= 30696191 & crispr_sites$start <= 30698216 & crispr_sites$score_crisprscan >= 40))
#231 in total
# 139 if you only take those with CRISPR scans > 40 
# 48 if you take those with CRISPR scans > 60 

# create a dot/line plot showing how many potential target sites there are per kb 
max(crispr_sites$start)-min(crispr_sites$start)
per_kb = seq(min(crispr_sites$start), max(crispr_sites$start), 1000)
per_kb <- append(per_kb, 30730027)
counter_all = c()
counter_med =c()
counter_high=c()
for (i in per_kb){
  new_count_all <- length(which(crispr_sites$start <= i))
  counter_all <- append(counter_all, new_count_all)
  new_count_med <- length(which(crispr_sites$start <= i & crispr_sites$score_crisprscan >= 40))
  counter_med <- append(counter_med, new_count_med)
  new_count_high <- length(which(crispr_sites$start <= i & crispr_sites$score_crisprscan >= 60))
  counter_high <- append(counter_high, new_count_high)
}

crispr_sites_df = data.frame(per_kb, per_kb_ref, counter_all, counter_med, counter_high)

dot_plot = ggplot(data=crispr_sites_df, aes(x=(per_kb_ref), counter_all))+geom_point(color = "#e2e2e2", shape = "square")+geom_point(aes(y=counter_med),color = "#6a6a6a", shape = "square")+geom_point(y=counter_high, color = "#000000", shape = "square")+theme_bw()+ labs( title = "Cumulative number of CRISPR sites in GWAS peak", x = "Distance in GWAS peak (kb)",y = "CRISPR site count")+ theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2))+theme(legend.position="none") + scale_x_continuous(expand = c(0, 0), limits = c(0,51000)) + scale_y_continuous(expand = c(0, 0), limits = c(0,4100))


ggplot(data = crispr_sites, aes(x=start)) + geom_histogram(binwidth = 500, color = "light gray", fill = "light gray")+theme_bw()+ labs( title = "CRISPR sites in GWAS peak", x = "Position on Chromosome 3 (kb)",y = "CRISPR site count")+ theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2))+scale_x_continuous(expand = c(0, 0))+scale_y_continuous(expand = c(0, 0), limits = c(0,80))+ annotate("rect", xmin = 30726727, xmax = 30727185, ymin = 0, ymax = 80, fill = "blue", alpha=0.35)+ annotate("rect", xmin = 30696191, xmax = 30698216, ymin = 0, ymax = 80, fill = "red", alpha=0.35)

crispr_sites_med = crispr_sites[which(crispr_sites$score_crisprscan > 40),]
crispr_sites_high = crispr_sites[which(crispr_sites$score_crisprscan > 60),]

hist_plot = ggplot(data = crispr_sites, aes(x=start)) + geom_histogram(binwidth = 500, color = "#e2e2e2", fill = "#e2e2e2") + geom_histogram(data = crispr_sites_med, y = start, binwidth = 500, color = "#6a6a6a", fill = "#6a6a6a") + geom_histogram(data = crispr_sites_high, y = start, binwidth = 500, color = "#000000", fill = "#000000") +theme_bw()+ labs( title = "CRISPR sites in GWAS peak", x = "Position on Chromosome 3 (kb)",y = "CRISPR site count")+ theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.2))+scale_x_continuous(expand = c(0, 0))+scale_y_continuous(expand = c(0, 0), limits = c(0,80))+ annotate("rect", xmin = 30726727, xmax = 30727185, ymin = 0, ymax = 80, fill = "blue", alpha=0.35)+ annotate("rect", xmin = 30696191, xmax = 30698216, ymin = 0, ymax = 80, fill = "red", alpha=0.35)+theme(plot.margin = unit(c(0,0.2,0,0.5),"cm"))

plot_grid(hist_plot, dot_plot, nrow = 2, labels = c("A", "B")) 
