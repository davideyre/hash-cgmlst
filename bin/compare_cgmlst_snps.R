#!/usr/bin/env Rscript
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggthemes)
theme_set(theme_bw() + theme(axis.title = element_text(size = 10.5)))

setwd("~/Drive/academic/infrastructure/hash-cgmlst/comparison_study_data/")

## REPLICATES ##
#read in replicate list
replicates = read.csv("replicate_list.csv", stringsAsFactors = F) %>% as_tibble()
replicates$id1 = substring(replicates$id1, 1, 8)
replicates$id2 = substring(replicates$id2, 1, 8)

#read in read qc data
read.data = read.table("replicates_reads_summary.txt", header = T, stringsAsFactors = F) %>% as_tibble()
read.data$id = substring(read.data$id, 0, 8)
rd.1 = read.data %>% select(id1="id", cov.1 = "coverage")
rd.2 = read.data %>% select(id2="id", cov.2 = "coverage")

#read in kraken2 data
spp = read.table("replicates_kraken_summary.txt", header=T, stringsAsFactors = F, sep="\t") %>% as_tibble()
spp$id = substring(spp$id, 0, 8)
spp$cdiff[which(spp$cdiff>1)]=1
spp.1 = spp %>% select(id1="id", cdiff.1="cdiff")
spp.2 = spp %>% select(id2="id", cdiff.2="cdiff")

# read in read length data
rl = read.table("replicates_read_length_summary.txt", header = T, stringsAsFactors = F, sep="\t") %>% as_tibble()
rl$id = substring(rl$id, 0, 8)
rl.1 = rl %>% select(id1="id", max_length.1="max_length")
rl.2 = rl %>% select(id2="id", max_length.2="max_length")

## SNPS ##
#read in snps
snps = read.table("merge_distances.txt", sep="\t", stringsAsFactors = F, header=T) %>% as_tibble()
#duplicate snps to enable matching with cgmlst
snps = rbind(snps, setNames(snps, c("id2", "id1", "pw_snps", "ml_snps", "cf_snps", "merge_snps"))) 

#read in pctACGT
acgt = read.table("pctACGT.txt", sep="\t", stringsAsFactors = F) %>% as_tibble()
#manually add 4 missing rows
acgt = rbind(acgt, c("43b90a59", NA, 0), c("43372d19" , NA, 0), c("8285610e", NA, 0), c("d65d4d8b", NA, 0))
colnames(acgt) = c("id", "file", "pct.acgt")
acgt.1 = acgt %>% select(c(id1 = "id", acgt.1 = "pct.acgt"))
acgt.2 = acgt %>% select(c(id2 = "id", acgt.2 = "pct.acgt"))

#check for missing snp data
id.list = unique(c(replicates$id1, replicates$id2)) 
missing = id.list[which(!(id.list %in% acgt$id))]
print("Misisng SNP data:")
print(missing)


## CGMLST ##
#read in cgmlst
cgmlst = read.table("replicates_compare.txt", sep="\t", header=T, stringsAsFactors = F)  %>% as_tibble()
cgmlst$id1 = substring(cgmlst$id1, 1, 8)
cgmlst$id2 = substring(cgmlst$id2, 1, 8)
#duplicate to enable matching based on reversed id1 and id2
cgmlst = rbind(cgmlst, setNames(cgmlst, c("id2", "id1", "loci_compared", "differences", "dist")))

# check for missing cgmlst data
cgmlst.check = unique(c(cgmlst$id1, cgmlst$id2))
missing = id.list[which(!(id.list %in% cgmlst.check))]
print("Misisng cgMLST data:")
print(missing)


#read in cgmlst with excluded sites
#read in replicate list with excluded genes
cgmlst.exclude = read.table("replicates_compare_exclude.txt", sep="\t", header=T, stringsAsFactors = F)  %>% as_tibble()
cgmlst.exclude$id1 = substring(cgmlst.exclude$id1, 1, 8)
cgmlst.exclude$id2 = substring(cgmlst.exclude$id2, 1, 8)
colnames(cgmlst.exclude) = c("id1", "id2", "exclude_loci_compared", "exclude_differences", "exclude_dist")
#duplicate to enable matching based on reversed id1 and id2
cgmlst.exclude = rbind(cgmlst.exclude, setNames(cgmlst.exclude, c("id2", "id1", "exclude_loci_compared", "exclude_differences", "exclude_dist")))

#read in assembly stats
assembly.stats = read.table("assembly_stats.txt", stringsAsFactors = F, header=T) %>% as_tibble()
assembly.stats$id = unlist(lapply(assembly.stats$filename, function(x) {substring(unlist(strsplit(x, "/"))[7], 1, 8)}))
stats.1 = assembly.stats %>% select(c(id1 = "id", gc.1="gc_avg", contig_bp.1 = "contig_bp", n50.1 = "ctg_N50", contigs.1="n_contigs"))
stats.2 = assembly.stats %>% select(c(id2 = "id", gc.2="gc_avg", contig_bp.2 = "contig_bp", n50.2 = "ctg_N50", contigs.2="n_contigs"))


#joint cgmlst and snp data
#result = inner_join(replicates, merge.dist, by=c("id1", "id2"))
result = left_join(replicates, snps, by=c("id1", "id2"))
result = left_join(result, cgmlst, by=c("id1", "id2"))
result = left_join(result, cgmlst.exclude, by=c("id1", "id2"))
#add pctACGT stats
result = left_join(result, acgt.1, by=c("id1"))
result = left_join(result, acgt.2, by=c("id2"))
#add assembly stats to results
result = left_join(result, stats.1, by=c("id1"))
result = left_join(result, stats.2, by=c("id2"))
#add read stats
result = left_join(result, rd.1, by=c("id1"))
result = left_join(result, rd.2, by=c("id2"))
#add kraken data
result = left_join(result, spp.1, by=c("id1"))
result = left_join(result, spp.2, by=c("id2"))
#add read length data
result = left_join(result, rl.1, by=c("id1"))
result = left_join(result, rl.2, by=c("id2"))

result.complete = result %>% filter(!is.na(pw_snps) & !is.na(differences))
result.incomplete = result %>% filter(is.na(pw_snps) | is.na(differences))

#filter out genomes with contig lengths not within +/- 10% of the median
med.bp = median(assembly.stats$contig_bp)
lower.bp = med.bp *.9
upper.bp = med.bp * 1.1

result.filtered = result.complete %>% 
  filter(contig_bp.1>lower.bp & contig_bp.1<upper.bp & 
           contig_bp.2>lower.bp & contig_bp.2<upper.bp &
           cov.1 >=50 & cov.2 >=50 &
           max_length.1!=501 & max_length.2!=501)

print(table(result.filtered$differences))

result.filtered$contig_bp.max = apply(cbind(abs(result.filtered$contig_bp.1-med.bp),
                                            abs(result.filtered$contig_bp.2-med.bp)),1,FUN=max)
result.filtered$n50.max = apply(cbind(result.filtered$n50.1, result.filtered$n50.2),1,FUN=max)
result.filtered$contigs.max = apply(cbind(result.filtered$contigs.1, result.filtered$contigs.2),1,FUN=max)
result.filtered$acgt.min = apply(cbind(result.filtered$acgt.1, result.filtered$acgt.2),1,FUN=min)
result.filtered$cov.min = apply(cbind(result.filtered$cov.1, result.filtered$cov.2),1,FUN=min)
result.filtered$cdiff.min = apply(cbind(result.filtered$cdiff.1, result.filtered$cdiff.2),1,FUN=min)
result.filtered$rl.min = apply(cbind(result.filtered$max_length.1, result.filtered$max_length.2),1,FUN=min)
result.filtered$rl.min[which(result.filtered$rl.min==108)] = 100



write.csv(result.filtered, "results_filtered.csv", row.names = F )


## get genes with most differences
g = read.table("genes_with_differences.txt", header=T, stringsAsFactors = F, sep="\t") %>% as_tibble()
gs = g %>% select(c("sample", "gene")) %>% distinct()
gsp = g %>% select(c("sample", "gene", "gene_length", "positions")) %>% distinct()
gene_count = gs %>% group_by(gene) %>% summarise(unique_samples = n())
gsp_n = inner_join(gsp, gene_count) %>% arrange(desc(unique_samples), gene)
write.csv(gsp_n, "gene_differences_summary.csv")

## Mapping based QC metrics - omit as don't add much

#OMIT (mapping based) - relationship between pctACGT and differences
ggplot(result.filtered, aes(x=as.numeric(acgt.min), y=differences, color=cov.min)) +
  geom_point()
kruskal.test(acgt.min ~ differences, data=result.filtered)

# OMIT - no real relationship between coverage and pct called
ggplot(result.filtered, aes(x=as.numeric(cov.min), y=as.numeric(acgt.min))) +
  geom_point()












#### FIGURES ####
setwd("/Users/davideyre/Drive/academic/infrastructure/cgmlst_archive/manuscript/figures")

#get samples with same pool identifier
result.filtered$pool = factor(ifelse(1:nrow(result.filtered) %in% grep('_p', result.filtered$samplename),1,0), 
                              label=c("Same isolate", "Same DNA pool"))

### FIGURE 1 - distributions of differences and SNPs
col = tableau_color_pal('Tableau 10')(10)[c(4,1)]

# group differences
result.filtered$diff_cut = cut(result.filtered$differences, 
                               breaks=c(-Inf, 0, 1, 2, 3, 4, 5, 9, 14, 19, Inf), 
                               labels=c("0", "1", "2", "3", "4", "5", "6-9", "10-14", "15-19", "20+"))
result.filtered$diff_cut_exclude = cut(result.filtered$exclude_differences, 
                               breaks=c(-Inf, 0, 1, 2, 3, 4, 5, 9, 14, 19, Inf), 
                               labels=c("0", "1", "2", "3", "4", "5", "6-9", "10-14", "15-19", "20+"))
result.filtered$snp_cut = cut(result.filtered$pw_snps, 
                              breaks=c(-Inf, 0, 1, 2, 3, 4, 5, 9, 14, 19, Inf), 
                              labels=c("0", "1", "2", "3", "4", "5", "6-9", "10-14", "15-19", "20+"))

p1a = ggplot(result.filtered, aes(x=diff_cut, fill=pool)) +
  geom_bar() +
  scale_fill_manual(values=col) +
  labs(y="Frequency", x="cgMLST gene differences\nbetween replicate sequences",
       fill="DNA sequenced") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,300) +
  theme(legend.position=c(0.77,0.84)) 


p1b = ggplot(result.filtered, aes(x=snp_cut, fill=pool)) +
  geom_bar() + 
  scale_fill_manual(values=col, drop=F) +
  scale_x_discrete(drop=FALSE) +
  labs(y="Frequency", x="SNP differences\nbetween replicate sequences",
       fill="DNA sequenced") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,300) +
  theme(legend.position=c(0.77,0.84)) 

p1 = grid.arrange(p1a, p1b, ncol = 2, nrow = 1)

ggsave("figure1.pdf", p1, width = 20, height = 10, units="cm", useDingbats=FALSE)


p1as = ggplot(result.filtered, aes(x=diff_cut, fill=pool)) +
  geom_bar() +
  scale_fill_manual(values=col) +
  labs(y="Frequency", x="cgMLST gene differences\nbetween replicate sequences\n",
       fill="DNA sequenced") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,250) +
theme(legend.position=c(0.77,0.84)) 

p1c = ggplot(result.filtered, aes(x=diff_cut_exclude, fill=pool)) +
  geom_bar() +
  scale_fill_manual(values=col, drop=FALSE) +
  scale_x_discrete(drop=FALSE) +
  labs(y="Frequency", x="cgMLST gene differences\nbetween replicate sequences\nExcluding 26 potentially mis-assembled genes",
       fill="DNA sequenced") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,250) +
  theme(legend.position=c(0.77,0.84)) 

ps1 = grid.arrange(p1as, p1c, ncol = 2, nrow = 1)
ggsave("figureS1.pdf", ps1, width = 20, height = 12, units="cm", useDingbats=FALSE)

#explore numbers of gene differences after filtering out potentially mis-assembled genes
table(result.filtered$diff_cut)
table(result.filtered$diff_cut_exclude)
length(which(result.filtered$differences>2))
length(which(result.filtered$exclude_differences>2))
length(which(result.filtered$differences>2 & result.filtered$pool=="Same DNA pool"))
length(which(result.filtered$exclude_differences>2 & result.filtered$pool=="Same DNA pool"))


### FIGURE 2 - relationship between cgMLST gene differences and coverage and read length
col = tableau_color_pal('Tableau 10')(10)[c(10,1,2)]
#coverage and differences
result.filtered$rl.min.gp = as.factor(round((result.filtered$rl.min)/50)*50) #covert minimum read lengths to factor

p2 = ggplot(result.filtered, aes(x=as.numeric(cov.min), y=differences, color=rl.min.gp)) +
  geom_jitter() +
  scale_color_manual(values=col) +
  labs(y="cgMLST gene differences\nbetween replicate sequence pairs", 
       x="Minimum average genome coverage in sequence pair",
       color="Minimum read length in sequence pair") +
  theme(legend.position="bottom")

kruskal.test(cov.min ~ differences, data=result.filtered)
cor.test(x=result.filtered$cov.min, y=result.filtered$differences, method='spearman')


kruskal.test(rl.min.gp ~ differences, data=result.filtered)

result.filtered$diff.above2 = ifelse(result.filtered$differences>2,1,0)

result.filtered %>% group_by(rl.min.gp) %>% summarise(n = length(differences), above2=sum(diff.above2))


ggsave("figure2.pdf", p2, width = 17, height = 10, units="cm", useDingbats=FALSE)

### FIGURE 3 - relationship between cgMLST gene differences and de novo assembly metrics and kraken2 classification
col = tableau_color_pal('Tableau 10')(2)



#assembly size and differences
p3a = ggplot(result.filtered, aes(x=contig_bp.max/med.bp*100, y=differences, color=pool)) +
  geom_jitter() +
  scale_color_manual(values = col) +
  theme(legend.position="bottom") +
  labs(y="cgMLST gene differences\nbetween replicate sequence pairs", 
       x=paste("Maximum absolute percentage deviation from\noverall median assemby size in sequence pair", sep=""),
       color="Isolate DNA")

kruskal.test(contig_bp.max ~ differences, data=result.filtered)
cor.test(x=result.filtered$contig_bp.max, y=result.filtered$differences, method='spearman')
result.filtered$contig_mean = (result.filtered$contig_bp.1 + result.filtered$contig_bp.2)/2
cor.test(x=result.filtered$contig_mean, y=result.filtered$differences, method='spearman')
p3b = ggplot(result.filtered, aes(x=as.numeric(n50.max), y=differences, color=pool)) +
  geom_jitter() +
  scale_color_manual(values = col) +
  theme(legend.position="bottom") +
  labs(y="cgMLST gene differences\nbetween replicate sequence pairs", 
       x="Maximum L50\nin sequence pair",
       color="Isolate DNA")
kruskal.test(n50.max ~ differences, data=result.filtered)
cor.test(x=result.filtered$n50.max, y=result.filtered$differences, method='spearman')

#false differences>2 by N50 <= 125
result.filtered$n50.above125 = ifelse(result.filtered$n50.max<=125,0,1)
result.filtered %>% group_by(n50.above125) %>% summarise(n = length(differences), above2=sum(diff.above2))
result.filtered %>% group_by(n50.above125, pool) %>% summarise(n = length(differences), above2=sum(diff.above2))

#contig count and differences
p3c = ggplot(result.filtered, aes(x=contigs.max, y=differences, color=pool)) +
  geom_jitter() +
  scale_color_manual(values = col) +
  theme(legend.position="bottom") +
  labs(y="cgMLST gene differences\nbetween replicate sequence pairs", 
       x="Maximum number of contigs in sequence pair",
       color="Isolate DNA")
kruskal.test(contigs.max ~ differences, data=result.filtered)
cor.test(x=result.filtered$contigs.max, y=result.filtered$differences, method='spearman')

#kraken and differences
p3d = ggplot(result.filtered, aes(x=as.numeric(cdiff.min), y=differences, color=pool)) +
  geom_jitter() +
  scale_color_manual(values = col) +
  theme(legend.position="bottom") +
  labs(y="cgMLST gene differences\nbetween replicate sequence pairs", 
       x="Minimum proportion of sequenced reads\nclassified as C. difficile in sequence pair",
       color="Isolate DNA") +
  xlim(0.85, 1)
kruskal.test(cdiff.min ~ differences, data=result.filtered)
cor.test(x=result.filtered$cdiff.min, y=result.filtered$differences, method='spearman')
#!! remember to add to legned have removed 2 points !!

p = grid.arrange(p3a, p3b, p3c, p3d, ncol = 2, nrow = 2)
ggsave("figure3.pdf", p, width = 20, height = 20, units="cm", useDingbats=FALSE)



#### SUMMARY NUMBERS ####
#number of unique isolates
length(unique(c(result.filtered$samplename)))

#number of unique sequences
length(unique(c(result.filtered$id1, result.filtered$id2)))

#number of comparisons
nrow(result.filtered)

#range of replicates
tally.by.samplename = result.filtered %>% group_by(samplename) %>% tally()
table(tally.by.samplename$n+1)
quantile(tally.by.samplename$n+1)

#read lenghts ~50
unique(result.filtered$id1[which(result.filtered$max_length.1<60)], 
       result.filtered$id2[which(result.filtered$max_length.2<60)])

#overall error rates
table(result.filtered$differences)
mean(result.filtered$differences)
1/mean(result.filtered$differences)

table(result.filtered$pw_snps)
mean(result.filtered$pw_snps)
1/mean(result.filtered$pw_snps)

#subset from same pool
result.subset = result.filtered[grep('_p', result.filtered$samplename),]
nrow(result.subset)
table(result.subset$differences)
mean(result.subset$differences)
1/mean(result.subset$differences)

table(result.subset$pw_snps)
mean(result.subset$pw_snps)
1/mean(result.subset$pw_snps)

#median assembly size
med.bp


#### BENCHMARKING ####
trace = read.table("trace.txt", header=TRUE, sep='\t', stringsAsFactors = F) %>% as_tibble()
trace = trace %>% select(c("tag", "process", "realtime"))
trace.wide = trace %>% spread(process, realtime)
trace.wide$qc = trace.wide$bbDuk + trace.wide$rawFastQC + trace.wide$cleanFastQC + trace.wide$kraken2
trace.wide = trace.wide %>% select(c("qc", "spades", "cgmlst"))

round(quantile(trace.wide$qc/1000/60, na.rm=T),1)
round(quantile(trace.wide$spades/1000/60, na.rm=T),1)
round(quantile(trace.wide$cgmlst/1000, na.rm=T),1)

#time to compare - hash cgmlst
t = 92.73947334289551
n = 229504
t/n*1e5

# standard mlst
t.cgmlst = 88.82461214065552
n = 229504
t.cgmlst/n*1e5
