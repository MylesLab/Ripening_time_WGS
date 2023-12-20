# Script for GATK variant calling

# This script analyzes the WGS raw reads and read depth from the 97 samples returned from Gencove sequencing centre. 

setwd('~/Desktop/DALMSc/kmer-gwas/classic-gwas/pipeline-stats')

######
# RAW READ COUNT
#####
# Files called "forward_reads" and "reverse_reads" are the wc -l counts from the raw fastq files 
# each line in the file is the raw line count of the fastq. Divide by 4 to get total number of reads. 
# forward is R1. Reverse is R2. They arent really forward or reverse, they are just the first sequence from a particular molecule. Perhaps better to combine. 

R1= read.delim('forward_reads.txt', sep = ' ', header = F)
R2= read.delim('reverse_reads.txt', sep = ' ', header = F)
samp_names= read.delim('sample_names.txt', sep = ' ', header = F)

R1$R2 = R2
R1$sample = samp_names

#now get true total read count by dividing each of the numbers by 4 

R1$V1= R1$V1/4
R1$R2 = R1$R2/4
R1$V1 == R1$R2 # all TRUE -makes sense as every single read captured in one file was captured by its counterpart in the other file 
# Because they are all equal, they can be doubled to find the total number of raw reads for each sample  
R1$total_reads = R1$V1*2
summary(R1$total_reads)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#44906942  89950640 111315898 111319067 134432614 237662374
#Min = 44M - DAL6187
#Max = 237M - 6119
#Mean = 111M 
#***HOWEVER : All reads in this experiment were paired; only PE sequencing was done. Therefor, the number of reads should be shown in "read pairs"

hist(R1$total_reads, breaks = 100 )
hist(R1$V1, breaks = 100 )
R1= read.delim('forward_reads.txt', sep = ' ', header = F) # Re-read file as numeric for plotting
R1$V1 = R1$V1/4
summary(R1$V1)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 22453471  44975320  55657949  55659533  67216307 118831187 
ggplot(R1, aes(x=V1)) + geom_histogram(fill = 'black', alpha = 1 , bins = 50)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Total raw read pairs by sample")+xlab('Total raw read pairs')+ylab('Sample count')+geom_vline(xintercept = mean(R1$V1), linetype="dotted",color = "Red", size=0.5)



######
#DUPLICATION METRICS 
######
dup_metrics = read.delim('master_metrics.txt', sep = '\t', row.names = NULL)
#first, add the sample names 
dup_metrics$LIBRARY = as.vector(samp_names$V1)

# the master metrics file is from Picard MarkDuplicates duplication metrics. 
# This was run on each of the sam files after alignment and duplicate marking. 
# See the Picard docs for info on each col. 
# I took just the line that has all the info shown in the picard documantation. 

# Next, I will look through each of the categories, and for the ones that contain relevant info, I will produce distributions 
#UNPAIRED_READS_EXAMINED
#The number of mapped reads examined which did not have a mapped mate pair, either because the read is unpaired, or the read is paired to an unmapped mate. # I believe in my case, all reads are paired, so this is the number of read pairs in which one or more read is unmapped
summary(dup_metrics$UNPAIRED_READS_EXAMINED)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 854776 1157944 1370005 1,493,159 1668957 4043281 
hist(dup_metrics$UNPAIRED_READS_EXAMINED, breaks = 100)

#READ_PAIRS_EXAMINED
# The number of mapped read pairs examined. (Primary, non- supplemental)
hist(dup_metrics$READ_PAIRS_EXAMINED, breaks = 100)
#Mean = 51,245,798
ggplot(dup_metrics, aes(x=READ_PAIRS_EXAMINED)) + geom_histogram(fill = 'black', alpha = 1 , bins = 50)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Total mapped read pairs by sample")+xlab('Mapped read pairs')+ylab('Sample count')+geom_vline(xintercept = mean(dup_metrics$READ_PAIRS_EXAMINED), linetype="dotted",color = "Red", size=0.5)


sum(dup_metrics[1,2:3])/R1$V1[1]
#1.00131 - the duplicate metrics is slightly inflated because

#SECONDARY_OR_SUPPLEMENTARY_RDS
#The number of reads that were either secondary or supplementary
hist(dup_metrics$SECONDARY_OR_SUPPLEMENTARY_RDS, breaks = 100)

#Possibly the previous two on a single overlapping plot?

#UNMAPPED_READS
# The total number of unmapped reads examined. (Primary, non- supplemental)
hist(dup_metrics$UNMAPPED_READS, breaks = 100)

#UNPAIRED_READ_DUPLICATES
# The number of fragments that were marked as duplicates.
hist(dup_metrics$UNPAIRED_READ_DUPLICATES, breaks = 100)

#READ_PAIR_DUPLICATES
#The number of read pairs that were marked as duplicates.
hist(dup_metrics$READ_PAIR_DUPLICATES, breaks = 100)

#READ_PAIR-OPTICAL_DUPLICATES
#The number of read pairs duplicates that were caused by optical duplication. Value is always < READ_PAIR_DUPLICATES, which counts all duplicates regardless of source.
hist(dup_metrics$READ_PAIR_OPTICAL_DUPLICATES, breaks = 100)

#PERCENT_DUPLICATION
#The fraction of mapped sequence that is marked as duplicate.
hist(dup_metrics$PERCENT_DUPLICATION*100, breaks = 100 )
summary(dup_metrics$PERCENT_DUPLICATION*100)


ggplot(dup_metrics, aes(x=PERCENT_DUPLICATION)) + geom_histogram(fill = 'black', alpha = 1 , bins = 50)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Percent of mapped sequence marked duplicate ")+xlab('% marked duplicate')+ylab('Sample count')+geom_vline(xintercept = mean(dup_metrics$PERCENT_DUPLICATION), linetype="dotted",color = "Red", size=0.5)

##### Could do a % of reads found with an unpaired mate (at least one of two unpaired?)

#ESTIMATED_LIBRARY_SIZE
#The estimated number of unique molecules in the library based on PE duplication.
hist(dup_metrics$ESTIMATED_LIBRARY_SIZE, breaks = 100 )




dup_no_samp = dup_metrics[,2:ncol(dup_metrics)]
a = cor(dup_no_samp, y = NULL, use = "everything", method = c('pearson'))
heatmap(a, Rowv = NA, Colv = NA, labRow = c("this","is","a","test","b","c","D","e","f"))

cor.test(tot_rd$mean_rd, dup_metrics$ESTIMATED_LIBRARY_SIZE) #SIG 

plot(tot_rd$mean_rd, dup_metrics$ESTIMATED_LIBRARY_SIZE)
hist(dup_metrics$UNMAPPED_READS)

#####
#AVERAGE READ DEPTH
#####
# the total read depth file is the output of samtools depth for each file that has then been summed. This will need to be divided by 709561391 (the number of positions in the genome)

# read file
tot_rd = read.delim('total_depth_all_samples.txt', header = F, sep= ' ')
t# calculate mean depth by dividing total depth across all positions by the total number of positions in the genome
tot_rd$mean_rd = (tot_rd$V2)/709561391
#extract ony the sample name from the first column and replace the first column with just the sample name 
sample_id =c()
for (i in tot_rd$V1){
  proper_name = strsplit(i, split = '/')[[1]][11]
  sample_id = c(sample_id, proper_name)
}
#set the first column to the sample names
tot_rd$V1 = sample_id
summary(tot_rd$mean_rd)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.223  14.471  18.011  17.560  20.852  36.274 
#Min sample - 8165
#Max sample 12103
#hist mean read depth 
hist(tot_rd$mean_rd, breaks = 100)

#note that all samples with DAL in the name have much higher average read depths 
ggplot(tot_rd, aes(x=mean_rd)) + geom_histogram(fill = 'black', alpha = 1 , bins = 50)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Mean Read Depth Across Samples")+xlab('Mean read depth')+ylab('Sample count')+geom_vline(xintercept = 17, linetype="dotted",color = "Red", size=0.5)

# compare raw read count with total read depth 
R1$sample == tot_rd$V1 # proof that the samples are in the same order 

plot(R1$V1, tot_rd$mean_rd)
cor.test(R1$total_reads, tot_rd$mean_rd) # no correlation between total raw reads and average read depth
