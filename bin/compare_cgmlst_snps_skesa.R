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
cgmlst = read.table("replicates_compare_skesa_spades.txt", sep="\t", header=T, stringsAsFactors = F)  %>% as_tibble()

cgmlst = cgmlst %>% 
  separate(id1, c("id1", "assembler1"), sep="_") %>%
  separate(id2, c("id2", "assembler2"), sep="_")
cgmlst$id1 = substring(cgmlst$id1, 1, 8)
cgmlst$id2 = substring(cgmlst$id2, 1, 8)

#duplicate to enable matching based on reversed id1 and id2
cgmlst = rbind(cgmlst, setNames(cgmlst, c("id2", "assembler2", "id1", "assembler1", "loci_compared", "differences", "dist")))

# check for missing cgmlst data
cgmlst.check = unique(c(cgmlst$id1, cgmlst$id2))
missing = id.list[which(!(id.list %in% cgmlst.check))]
print("Misisng cgMLST data:")
print(missing)

#read in assembly stats
assembly.stats = read.table("replicates_assembly_stats.txt", stringsAsFactors = F, header=T) %>% as_tibble()
assembly.stats$id = unlist(lapply(assembly.stats$filename, function(x) {substring(unlist(strsplit(x, "/"))[7], 1, 8)}))
stats.1 = assembly.stats %>% select(c(id1 = "id", gc.1="gc_avg", contig_bp.1 = "contig_bp", n50.1 = "ctg_N50", contigs.1="n_contigs"))
stats.2 = assembly.stats %>% select(c(id2 = "id", gc.2="gc_avg", contig_bp.2 = "contig_bp", n50.2 = "ctg_N50", contigs.2="n_contigs"))


#joint cgmlst and snp data
#result = inner_join(replicates, merge.dist, by=c("id1", "id2"))
result = left_join(replicates, snps, by=c("id1", "id2"))
result = left_join(result, cgmlst, by=c("id1", "id2"))

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






#### FIGURES ####
setwd("/Users/davideyre/Drive/academic/infrastructure/cgmlst_archive/manuscript/figures")

#get samples with same pool identifier
result.filtered$pool = factor(ifelse(1:nrow(result.filtered) %in% grep('_p', result.filtered$samplename),1,0), 
                              label=c("Same isolate", "Same DNA pool"))

### FIGURE 1 - distributions of differences and SNPs
col = tableau_color_pal('Tableau 10')(10)[c(4,1)]

# group differences
result.filtered$diff_cut = cut(result.filtered$differences, 
                               breaks=c(-Inf, 0, 1, 2, 3, 4, 5, 9, Inf), 
                               labels=c("0", "1", "2", "3", "4", "5", "6-9", "10+"))
result.filtered$snp_cut = cut(result.filtered$pw_snps, 
                              breaks=c(-Inf, 0, 1, 2, 3, 4, 5, 9, Inf), 
                              labels=c("0", "1", "2", "3", "4", "5", "6-9", "10+"))

r.spades = result.filtered %>% filter(assembler1=="spades", assembler2=="spades")
r.skesa = result.filtered %>% filter(assembler1=="skesa", assembler2=="skesa")

p1a = ggplot(r.spades, aes(x=diff_cut, fill=pool)) +
  geom_bar() +
  scale_fill_manual(values=col) +
  labs(y="Frequency", x="cgMLST gene differences\nbetween replicate sequences\n(Spades)",
       fill="DNA sequenced") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,300) +
  theme(legend.position=c(0.77,0.84)) 

p1b = ggplot(r.skesa, aes(x=diff_cut, fill=pool)) +
  geom_bar() +
  scale_fill_manual(values=col, drop=F) + scale_x_discrete(drop=FALSE) +
  labs(y="Frequency", x="cgMLST gene differences\nbetween replicate sequences\n(Skesa)",
       fill="DNA sequenced") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,300) +
  theme(legend.position=c(0.77,0.84)) 

p1 = grid.arrange(p1a, p1b, ncol = 2, nrow = 1)

ggsave("figure1z.pdf", p1a, width = 20, height = 12, units="cm", useDingbats=FALSE)

