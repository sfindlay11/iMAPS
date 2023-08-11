# compare conservation of rep sites to position matched sites

library(data.table)
library(tidyverse)

path_data_rep <- ""
path_data_cons <- ""
path_data_fasta <- ""
path_data_plots <- ""
path_data_cds <- ""
path_data_out <- ""

names_BED_std <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

options(scipen = 999)

#### SET UP #### =================

## conservation data #### ====
setwd(path_data_cons)
cons <- fread("hg38.prot_cod.pADB_utr3_last_cds.gencode_utr3.phyloP100way.gz")

## rep sites #### ====
setwd(path_data_rep)
reps <- fread("eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.all_possible_snvs.txt")
# subset 3'UTRs
reps <- reps[top_motif_in_cust_utr3 == 1 & top_motif_in_cust_cds == 0 & top_motif_in_cust_utr3_intron == 0]

# add pos_id
reps[, c("chrom_extra", "chromStart", "chromEnd", "strand_extra") := tstrsplit(base_id, "_")]
reps[, pos_id := paste(chrom_extra, chromStart, chromEnd, sep = "_")]
reps[, c("chrom_extra", "strand_extra") := NULL]

# add peak id
reps[, peak_id := paste(chrom_peak, chromStart_peak, chromEnd_peak, strand_peak, RBP_peak, sep = "_")]

# add conservation annotations
setkey(reps, pos_id)
setkey(cons, pos_id)
reps <- cons[reps]

# add base change
reps[, base_change := paste(ref_base, alt_base, sep = ">")]
table(reps$base_change)

# for rep sites with peaks in both cell lines, randomly select one
reps[, rep_id := paste(chrom_peak, chromStart_top, chromEnd_top, strand_peak, RBP_peak, sep = "_")]
reps[, cell_count := length(unique(cell_line_peak)), by = .(rep_id)]
table(reps$cell_count)
reps[cell_count == 1, cell_sel := cell_line_peak]
set.seed(314883)
reps[cell_count == 2, cell_sel := sample(unique(cell_line_peak), 1, replace = FALSE), by = .(rep_id)]

# check for CDS overlap
setwd(path_data_cds)
fwrite(unique(reps[, .(chrom_peak, chromStart_top, chromEnd_top, name = rep_id, score = 0, strand_peak)]), "rep_sites_for_int.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
system("bedtools intersect -loj -a rep_sites_for_int.bed -b merged_nat_sorted_hg38_kg_gencode_v32_main_chr_prot_cod_CDS_exons.bed > rep_sites_cds_int.txt")
ints <- fread("rep_sites_cds_int.txt")
# are there any intersections?
ints[V9 != -1] # none

"merged_nat_sorted_hg38_kg_gencode_v32_main_chr_prot_cod_CDS_exons.bed"

## position matches #### ====
setwd(path_data_rep)
matches <- fread("eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.pos_matched_CONTROL_motifs.all_possible_snvs.pm_50.txt")

# process
matches[, c("RBP_peak", "cell_line_peak", "chrom_peak", "chromStart_peak", "chromEnd_peak", "strand_peak", "extra1", "extra2", "chromStart_five_end", "chromEnd_five_end") := tstrsplit(name, "_")]
matches[, peak_id := paste(chrom_peak, chromStart_peak, chromEnd_peak, strand_peak, RBP_peak, sep = "_")]
# subset for only matches of 3'UTR rep sites
matches <- matches[peak_id %in% reps_summ$peak_id]

# add pos_id
matches[, c("chrom_extra", "chromStart", "chromEnd", "strand_extra") := tstrsplit(base_id, "_")]
matches[, pos_id := paste(chrom_extra, chromStart, chromEnd, sep = "_")]
matches[, c("chrom_extra", "strand_extra") := NULL]
# add conservation annotations
setkey(cons, pos_id)
setkey(matches, pos_id)
matches <- cons[matches]

# add base change
matches[, base_change := paste(ref_base, alt_base, sep = ">")]
table(matches$base_change)

# check for CDS overlap
matches[, match_id := paste(chrom_peak, chromStart_top, chromEnd_top, strand_peak, RBP_peak, sep = "_")]
setwd(path_data_cds)
fwrite(unique(matches[, .(chrom_peak, chromStart_top, chromEnd_top, name = match_id, score = 0, strand_peak)]), "match_sites_for_int.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
system("bedtools intersect -loj -a match_sites_for_int.bed -b merged_nat_sorted_hg38_kg_gencode_v32_main_chr_prot_cod_CDS_exons.bed > match_sites_cds_int.txt")
ints <- fread("match_sites_cds_int.txt")
# are there any intersections?
ints[V9 != -1] # 15k

# mark these
matches[, in_cds := ifelse(match_id %in% ints[V9 != -1, V4], 1, 0)]
table(matches$in_cds)

#### SUMMARIZE and MERGE #### ====================

## summarize reps #### ====

reps_summ <- reps[cell_line_peak == cell_sel & base_change %in% c("A>C", "C>G", "G>T", "T>A"), # since ReP site table is fully expanded to include all possible SNVs, and phyloP is by position
                  .(phylop_100_max = max(phylop_100, na.rm = TRUE),
                    phylop_100_mean = mean(phylop_100, na.rm = TRUE),
                    phylop_100_min_2_count = sum(phylop_100 > 2),
                    phylop_100_min_3_count = sum(phylop_100 > 3),
                    phylop_100_min_4_count = sum(phylop_100 > 4),
                    phylop_100_min_5_count = sum(phylop_100 > 5),
                    phylop_100_min_6_count = sum(phylop_100 > 6),
                    phylop_100_min_7_count = sum(phylop_100 > 7),
                    phylop_100_min_8_count = sum(phylop_100 > 8),
                    phylop_100_min_2_prop = sum(phylop_100 > 2) / .N,
                    phylop_100_min_3_prop = sum(phylop_100 > 3) / .N,
                    phylop_100_min_4_prop = sum(phylop_100 > 4) / .N,
                    phylop_100_min_5_prop = sum(phylop_100 > 5) / .N,
                    phastcons_100_mean = mean(phastcons_100, na.rm = TRUE),
                    phastcons_100_min_0d99_prop = sum(phastcons_100 > 0.99) / .N,
                    peak_id_count = .N
                  ),
                  by = .(peak_id, aff_top_ref)]
reps_summ

## summarize matches #### ====

match_summ <- matches[base_change %in% c("A>C", "C>G", "G>T", "T>A") & in_cds == 0,
                      .(phylop_100_max_MATCH = max(phylop_100, na.rm = TRUE),
                        phylop_100_mean_MATCH = mean(phylop_100, na.rm = TRUE),
                        phylop_100_min_2_count_MATCH = sum(phylop_100 > 2),
                        phylop_100_min_3_count_MATCH = sum(phylop_100 > 3),
                        phylop_100_min_4_count_MATCH = sum(phylop_100 > 4),
                        phylop_100_min_5_count_MATCH = sum(phylop_100 > 5),
                        phylop_100_min_6_count_MATCH = sum(phylop_100 > 6),
                        phylop_100_min_7_count_MATCH = sum(phylop_100 > 7),
                        phylop_100_min_8_count_MATCH = sum(phylop_100 > 8),
                        phylop_100_min_2_prop_MATCH = sum(phylop_100 > 2) / .N,
                        phylop_100_min_3_prop_MATCH = sum(phylop_100 > 3) / .N,
                        phylop_100_min_4_prop_MATCH = sum(phylop_100 > 4) / .N,
                        phylop_100_min_5_prop_MATCH = sum(phylop_100 > 5) / .N,
                        phastcons_100_mean_MATCH = mean(phastcons_100, na.rm = TRUE),
                        phastcons_100_min_0d99_prop_MATCH = sum(phastcons_100 > 0.99) / .N,
                        peak_id_count_MATCH = .N
                        ),
                      by = .(peak_id)]
match_summ

## merge #### ====

# first filter reps_summ so peaks with matches overlapping cds are excluded
reps_summ <- reps_summ[peak_id %in% match_summ$peak_id]

setkey(match_summ, peak_id)
setkey(reps_summ, peak_id)
summ_all <- match_summ[reps_summ]


## overall odds ratio approach
# combine all reps and matches

all_build <- NULL

for (aff_min_curr in c(0, 0.01, 0.05, 0.1, 0.2)) {
  
  for (phylop_min_curr in c(2, 3, 4, 5, 6)) {
    
    # testing
    #aff_min_curr <- 0.2
    #phylop_min_curr <- 5
    
    all_curr <- summ_all[aff_top_ref > aff_min_curr &
                           is.na(get(paste("phylop_100_min_", phylop_min_curr, "_count", sep = ""))) == FALSE & is.na(get(paste("phylop_100_min_", phylop_min_curr, "_count_MATCH", sep = ""))) == FALSE &
                           peak_id_count == peak_id_count_MATCH]
    
    all_curr_summ <- all_curr[, .(rep_cons_count = sum(get(paste("phylop_100_min_", phylop_min_curr, "_count", sep = ""))),
                                  rep_total_count = sum(peak_id_count),
                                  match_cons_count = sum(get(paste("phylop_100_min_", phylop_min_curr, "_count_MATCH", sep = ""))),
                                  match_total_count = sum(peak_id_count_MATCH))]
    
    all_curr_summ[, rep_non_count := rep_total_count - rep_cons_count]
    all_curr_summ[, match_non_count := match_total_count - match_cons_count]
    all_curr_summ[, fisher_or := fisher.test(x = matrix(data = c(rep_cons_count, rep_non_count, match_cons_count, match_non_count), ncol = 2))$estimate]
    all_curr_summ[, fisher_p := fisher.test(x = matrix(data = c(rep_cons_count, rep_non_count, match_cons_count, match_non_count), ncol = 2))$p.value]
    all_curr_summ[, fisher_ci_low := fisher.test(x = matrix(data = c(rep_cons_count, rep_non_count, match_cons_count, match_non_count), ncol = 2))$conf.int[1]]
    all_curr_summ[, fisher_ci_high := fisher.test(x = matrix(data = c(rep_cons_count, rep_non_count, match_cons_count, match_non_count), ncol = 2))$conf.int[2]]
    
    # add info
    all_curr_summ[, aff_min := aff_min_curr]
    all_curr_summ[, phylop_min := paste("phyloP >", phylop_min_curr)]
    
    # add to build
    all_build <- rbind(all_build, all_curr_summ)
    
  }
  
}

all_build[, .(aff_min, phylop_min, fisher_or, round(fisher_p, 10), fisher_ci_low, fisher_ci_high)]


plot_curr <- ggplot(data = all_build, aes(x = aff_min, y = log2(fisher_or), colour = phylop_min)) +
  geom_point(size = 2) +
  geom_line() +
  geom_ribbon(aes(ymin = log2(fisher_ci_low), ymax = log2(fisher_ci_high), fill = phylop_min, colour = NULL), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 1, colour = "gray30") +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_colour_discrete(labels = c("2", "3", "4", "5", "6"), name = "minimum\nPhyloP") +
  scale_fill_discrete(labels = c("2", "3", "4", "5", "6"), name = "minimum\nPhyloP") +
  theme_minimal(base_size = 16) +
  #theme(axis.text=element_text(size = 10)) +
  labs(x = "\nminimum affinity of ReP site", y = "odds ratio, log2\n(conserved position in ReP vs. match)\n")

plot_curr
ggsave("reps.pos_matches.cons_OR.by_min_phylop.by_min_aff.pdf", path = path_data_plots, width = 5, height = 5, units = "in")

all_build[phylop_min == "phyloP > 6" & aff_min == 0.2]

# save
setwd(path_data_out)
#fwrite(summ_all, "sum_all.rep_clip_cons.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#fwrite(all_build, "all_build.rep_clip_cons.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#all_build <- fread("all_build.rep_clip_cons.txt")





