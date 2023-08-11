# identify the motif / putative RBP binding site that "explains" each eCLIP peak

# considers motifs within 75 bases of 5' end of eCLIP peak
# takes the closest motif to the 5' end in the case of ties. 
# also considers the next highest affinity peak when the highest one is closer to another peak's 5' end (iteratively)
# this way one motif is not assigned to multiple peaks

## set up ------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(pbapply)
library(RColorBrewer)
library(parallel)

path_data_in <- ""
path_data_out <- ""
path_plot_out <- ""
path_data_utr3_space <- ""
path_data_cds_space <- ""
path_data_plots <- ""

options(scipen = 999)

names_BED_std <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

## -----------------------------------------------------------------------


## identify best motif for each eCLIP peak -----------------------------------------------------------------------

# read in (select sequence-based ("old") or sequence+structure ("new") motifs)
setwd(path_data_in)
#peaks <- readRDS("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_110.full_info.expanded.rds")
peaks <- readRDS("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_110.full_info.expanded.rds")
table(peaks$RBP_peak)

# add number of peaks per rbp cell line pair
peaks[, peak_count_by_pair := .N, by = .(RBP_peak, cell_line_peak)]

# add number of peaks per rbp (averaged across cell lines where applicable)
peaks[, avg_peak_count_by_rbp := mean(unique(peak_count_by_pair)), by = .(RBP_peak)]

# confirm extension distance used in data set
ext_dist <- as.integer(names(table(peaks[, nchar(ref_seq_ext) - ((k_size -1) *2) - 1]))) / 2 
ext_dist

# how long is each affinity list?
peaks[, aff_len := pbsapply(affinity_list, length)]
table(peaks$aff_len)
# this corresponds to affinities for the first (5') base of the peak, ext_dist + (k-1) (109/110) affinities upstream, and ext_dist (100) affinities downstream
# i.e. for k = 11 RBPs, affinity[106] represents the affinity of the motif centred around the 5' end of the peak

# create masked version for assignment (within specified ext_dist - trim_in distance of the 5' end)
mask_window <- function(aff_list_in, trim_in) {
  aff_list <- as.numeric(aff_list_in)
  aff_list[1:trim_in] <- -1
  aff_list[(length(aff_list) - (trim_in -1)):length(aff_list)] <- -1
  return(aff_list)
}
peaks[, affinity_list_masked := pbmapply(mask_window, aff_list_in = affinity_list, trim_in = 25)] # change trim dist here

# add max
peaks[, max_aff := sapply(affinity_list_masked, max)]

# which position is max 
peaks[, pos_max_aff := mapply(function(al, ma) which(al == ma), al = affinity_list_masked, ma = max_aff)]

# assess ties (select closest to 5' end in case of ties)
peaks[, count_pos_max_aff := sapply(pos_max_aff, length)]
table(peaks[, count_pos_max_aff, by = RBP_peak])

# record position of the first base of a motif that would be centred around the 5' end of the peak (using first base to match the index of the motif affinities)
# need to choose base to left or right for k = 10 RBPs (assign randomly for all bases, wont impact k = 11)
set.seed(009008)
peaks[, adjust_rand := sample(1:2, nrow(peaks), replace = T)]
peaks[adjust_rand == 1, pos_rel_five_base := pbsapply(k_size, function(k) ext_dist + k - floor(median(1:k) -1))]
peaks[adjust_rand == 2, pos_rel_five_base := pbsapply(k_size, function(k) ext_dist + k - ceiling(median(1:k) -1))]
table(peaks$k_size, peaks$pos_rel_five_base)
table(peaks$adjust_rand, peaks$pos_rel_five_base)

# for ties, select closest motif to 5' end of peak (apply to all data to avoid undesired interpretation of input by mapply)
peaks[, pos_max_aff_single := pbmapply(function(pos, pos_ref) pos[abs(pos - pos_ref) == min(abs(pos - pos_ref))], pos = pos_max_aff, pos_ref = pos_rel_five_base)]

# might still be ties
peaks[, count_pos_max_aff_single := sapply(pos_max_aff_single, length)]
table(peaks[, count_pos_max_aff_single])
# now can select randomly
set.seed(411311)
peaks[, pos_max_aff_single := sapply(pos_max_aff_single, function(x) ifelse(length(x) > 1, sample(x, 1), x))]

# tidy
peaks[, c("count_pos_max_aff", "count_pos_max_aff_single") := NULL]

# add position relative to 5' end of peak
peaks[, pos_max_aff_single_rel_five_end := pos_max_aff_single - pos_rel_five_base]
# expand to cover whole motif
peaks[, adjust := sapply(k_size, function(x) median(1:x))]
peaks[adjust_rand == 1, to_add_up := k_size - floor(adjust)]
peaks[adjust_rand == 1, to_add_down := k_size - ceiling(adjust)]
peaks[adjust_rand == 2, to_add_up := k_size - ceiling(adjust)]
peaks[adjust_rand == 2, to_add_down := k_size - floor(adjust)]
peaks[, pos_max_aff_single_rel_five_end_exp := pbmapply(function(pos, up, down) seq(pos - up, pos + down, by = 1), pos = pos_max_aff_single_rel_five_end, up = to_add_up, down = to_add_down, SIMPLIFY = F)]

## ensure unique assignment of motifs to peaks -----------------------------

# use pos_max_aff_single to extract coordinates of top motif
peaks[strand_peak == "+", chromStart_top := chromStart_ext + (pos_max_aff_single -1)]
peaks[strand_peak == "+", chromEnd_top := chromStart_top + k_size]

peaks[strand_peak == "-", chromEnd_top := chromEnd_ext - (pos_max_aff_single -1)]
peaks[strand_peak == "-", chromStart_top := chromEnd_top - k_size]

# check
peaks[RBP_peak == "RBFOX2" & max_aff == 1 & strand_peak == "+", ][1]
peaks[RBP_peak == "RBFOX2" & max_aff == 1 & strand_peak == "-", ][1]

# make id for top motif
peaks[, top_id := paste(RBP_peak, cell_line_peak, chrom_peak, chromStart_top, chromEnd_top, strand_peak, sep = "_")]

# how many motifs are assigned to multiple peaks?
peaks[, top_id_count := .N, by = .(top_id)]
table(peaks[, top_id_count])
length(unique(peaks[top_id_count == 1, top_id]))
length(unique(peaks[top_id_count == 1, top_id])) / length(unique(peaks$top_id)) # 89% unique
# flag these to keep
peaks[top_id_count == 1, keep := T]

# for those with more than one peak, keep the closest peak
peaks[top_id_count > 1, keep := ifelse(abs(pos_max_aff_single_rel_five_end) == min(abs(pos_max_aff_single_rel_five_end)), T, F), by = .(top_id)]
# randomly break ties
peaks[top_id_count > 1 & keep == T, keep_count := .N, by = top_id]
table(peaks[top_id_count > 1 & keep == T, keep_count])
peaks[keep_count == 2, ]
set.seed(455349)
peaks[keep_count == 2, keep_rand := rep(sample(1:2, 1), 2), by = top_id]
peaks[keep_count == 2 & keep_rand == 1 & pos_max_aff_single_rel_five_end > 0, keep := F]
peaks[keep_count == 2 & keep_rand == 2 & pos_max_aff_single_rel_five_end < 0, keep := F]

# what proportion of peaks are now "properly" assigned to a motif?
nrow(peaks[keep == T, ]) / nrow(peaks) # only 87%

# separate them
keeps <- peaks[keep == T, ]
cands_curr <- peaks[keep == F, ]

# are all keeps unique?
length(unique(keeps$top_id)) == nrow(keeps)

# mask unkept positions (replace affinity with -1 placeholder)
cands_curr[, affinity_list_masked := pbmapply(function(aff_list, to_mask) {aff_list[to_mask] <- -1; return(aff_list)}, aff_list = affinity_list_masked, to_mask = pos_max_aff_single)]
# reset
cands_curr[, pos_max_aff_single := NULL]

## iterate ---------------------------------------------------------------

round_curr <- 1

while (round_curr < Inf) {
  
  round_curr <- round_curr + 1
  
  print(paste("round", round_curr, sep = " "))
  
  # similar to code above
  
  table(cands_curr$RBP_peak)
  
  # update max (now that masking positions not kept in previous rounds has been applied)
  cands_curr[, max_aff := sapply(affinity_list_masked, max)]
  
  # which position is max (select closest to 5' end in case of ties)
  cands_curr[, pos_max_aff := mapply(function(al, ma) which(al == ma), al = affinity_list_masked, ma = max_aff, SIMPLIFY = F)]
  
  # assess ties
  cands_curr[, count_pos_max_aff := sapply(pos_max_aff, length)]
  table(cands_curr[, count_pos_max_aff, by = RBP_peak])
  
  # confirm extension distance used in data set
  ext_dist <- as.integer(names(table(cands_curr[, nchar(ref_seq_ext) - ((k_size -1) *2) - 1]))) / 2 
  
  # record position of the first base of a motif that would be centred around the 5' end of the peak (using first base to match the index of the motif affinities)
  # need to choose base to left or right for k = 10 RBPs (assign randomly for all bases, wont impact k = 11)
  # use adjust_rand already assigned
  cands_curr[adjust_rand == 1, pos_rel_five_base := sapply(k_size, function(k) ext_dist + k - floor(median(1:k) -1))]
  cands_curr[adjust_rand == 2, pos_rel_five_base := sapply(k_size, function(k) ext_dist + k - ceiling(median(1:k) -1))]
  table(cands_curr$k_size, cands_curr$pos_rel_five_base)
  
  # for ties, select closest motif to 5' end of peak (apply to all data to avoid problem with unintended simplification of input from mapply)
  cands_curr[, pos_max_aff_single := mapply(function(pos, pos_ref) pos[abs(pos - pos_ref) == min(abs(pos - pos_ref))], pos = pos_max_aff, pos_ref = pos_rel_five_base)]
  
  # might still be ties
  cands_curr[, count_pos_max_aff_single := sapply(pos_max_aff_single, length)]
  table(cands_curr[, count_pos_max_aff_single])
  # now can select randomly (again, apply to all rows, remember to avoid sampling single number because of an undesirable "convenience" feature of the sample function)
  set.seed(411311 * round_curr)
  cands_curr[, pos_max_aff_single := sapply(pos_max_aff_single, function(x) ifelse(length(x) > 1, sample(x, 1), x))]
  
  # tidy
  cands_curr[, c("count_pos_max_aff", "count_pos_max_aff_single") := NULL]
  
  # add position relative to 5' end of peak
  cands_curr[, pos_max_aff_single_rel_five_end := pos_max_aff_single - pos_rel_five_base]
  # expand to cover whole motif
  cands_curr[, pos_max_aff_single_rel_five_end_exp := NULL] # first reset previous
  cands_curr[, pos_max_aff_single_rel_five_end_exp := mapply(function(pos, up, down) seq(pos - up, pos + down, by = 1), pos = pos_max_aff_single_rel_five_end, up = to_add_up, down = to_add_down, SIMPLIFY = FALSE)]
  
  ## ensure unique assignment of motifs to cands_curr -----------------------------
  
  # use pos_max_aff_single to extract coordinates of top motif
  cands_curr[strand_peak == "+", chromStart_top := chromStart_ext + (pos_max_aff_single -1)]
  cands_curr[strand_peak == "+", chromEnd_top := chromStart_top + k_size]
  
  cands_curr[strand_peak == "-", chromEnd_top := chromEnd_ext - (pos_max_aff_single -1)]
  cands_curr[strand_peak == "-", chromStart_top := chromEnd_top - k_size]
  
  # make id for top motif
  cands_curr[, top_id := paste(RBP_peak, cell_line_peak, chrom_peak, chromStart_top, chromEnd_top, strand_peak, sep = "_")]
  
  # how many motifs are assigned to multiple cands_curr?
  cands_curr[, top_id_count := .N, by = .(top_id)]
  table(cands_curr[, top_id_count])
  length(unique(cands_curr[top_id_count == 1, top_id]))
  length(unique(cands_curr[top_id_count == 1, top_id])) / length(unique(cands_curr$top_id)) # round 2: 84% unique
  # flag these to keep
  cands_curr[top_id_count == 1, keep := T]
  
  # for those with more than one peak, keep the closest peak
  cands_curr[top_id_count > 1, keep := ifelse(abs(pos_max_aff_single_rel_five_end) == min(abs(pos_max_aff_single_rel_five_end)), T, F), by = .(top_id)]
  # randomly break ties
  cands_curr[, keep_count := NULL] # reset first
  cands_curr[top_id_count > 1 & keep == T, keep_count := .N, by = top_id]
  table(cands_curr[top_id_count > 1 & keep == T, keep_count])
  cands_curr[keep_count == 2, ]
  set.seed(328572 * round_curr)
  cands_curr[keep_count == 2, keep_rand := rep(sample(1:2, 1), 2), by = top_id]
  cands_curr[keep_count == 2 & keep_rand == 1 & pos_max_aff_single_rel_five_end > 0, keep := F]
  cands_curr[keep_count == 2 & keep_rand == 2 & pos_max_aff_single_rel_five_end < 0, keep := F]
  
  # among keeps, check for clashes with already kept motifs
  sum(unique(cands_curr[keep == T, top_id]) %in% keeps$top_id)
  sum(unique(cands_curr[keep == T, top_id]) %in% keeps$top_id) / length(unique(cands_curr[keep == T, top_id]))
  # flag them for removal
  cands_curr[keep == T & top_id %in% keeps$top_id, keep := F]
  
  # separate them
  print(paste(nrow(cands_curr[keep == T, ]), "more entries have been assigned"))
  keeps <- rbind(keeps, cands_curr[keep == T, ])
  
  cands_curr <- cands_curr[keep == F, ]
  print(paste(nrow(cands_curr), "entries remain unassigned"))
  
  if (nrow(cands_curr) == 0) {
    print("done!")
    break
  }
  
  # mask unkept positions (replace affinity with -1 placeholder)
  cands_curr[, affinity_list_masked := mapply(function(aff_in, to_mask) {aff_list <- unlist(aff_in); aff_list[to_mask] <- -1; return(aff_list)}, aff_in = affinity_list_masked, to_mask = pos_max_aff_single, SIMPLIFY = FALSE)]
  # reset
  cands_curr[, pos_max_aff_single := NULL]
  cands_curr[, keep := NULL]
  
}

# tidy
rm(cands_curr)

# quality check
# are all motifs unique?
length(unique(keeps$top_id)) == nrow(keeps)
keeps[, top_id_count_check := .N, by = .(top_id)]
keeps[top_id_count_check > 1, ]

# are all peaks represented?
sum(peaks$name_bed %in% keeps$name_bed) == nrow(peaks)

# do any affinities == -1 (i.e. ran out of non-masked sites)?
keeps[max_aff == -1, ]

# tidy
rm(peaks, round_curr)
keeps[, c("pos_max_aff", "adjust_rand", "top_id_count", "keep", "keep_count", "keep_rand", "top_id_count_check") := NULL]

## end of assigning top motifs to eCLIP peaks -----------------------------------------------------------------------

## save

# as bed file focused on top motifs (choose old or new motifs)
setwd(path_data_out)
#fwrite(keeps[, .(chrom_peak, chromStart_top, chromEnd_top, name = paste(name_bed, score_bed, sep = "_"), score = paste(max_aff, pos_max_aff_single_rel_five_end, sep = "_"), strand_peak)], "eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#fwrite(keeps[, .(chrom_peak, chromStart_top, chromEnd_top, name = paste(name_bed, score_bed, sep = "_"), score = paste(max_aff, pos_max_aff_single_rel_five_end, sep = "_"), strand_peak)], "eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# split
#library(strex)
#keeps[, name_bed := str_before_nth(name, "[_]", 10)]
#keeps[, score_bed := str_after_nth(name, "[_]", 10)]
#keeps[, c("gene_id_prot_cod_string", "gene_count", "prot_cod_prop", "region_count", "region_top") := tstrsplit(score_bed, "[_]")]
#keeps[, c("RBP_peak", "cell_line_peak", "chrom_peak", "chromStart_peak", "chromEnd_peak", "strand_peak", "enrichment_peak", "neg_log_10_p_value_peak", "chromStart_five_end", "chromEnd_five_end") := tstrsplit(name_bed, "_", fixed = TRUE)]

## update region classification-----------------------------------------------------------------------

# copy 3'UTR space to dir **3'UTRs downstream of stop codons only**
system(paste("cp ", path_data_utr3_space, "merged_by_strand_nat_sorted_hg38_pADB_v3_matched_to_ensembl_101_LAST_EXONS_STRANDED_plus_gencode_V32_UTR3s_minus_CDS_exons_for_space_filter.bed ", path_data_in, sep = ""))
system("head merged_by_strand_nat_sorted_hg38_pADB_v3_matched_to_ensembl_101_LAST_EXONS_STRANDED_plus_gencode_V32_UTR3s_minus_CDS_exons_for_space_filter.bed")
# intersect top motifs with 3' UTR space (use -u option to just report entries with any overlap)
#system("bedtools intersect -u -s -a eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.bed -b merged_by_strand_nat_sorted_hg38_pADB_v3_matched_to_ensembl_101_LAST_EXONS_STRANDED_plus_gencode_V32_UTR3s_minus_CDS_exons_for_space_filter.bed > eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_utr3s.bed")
#system("bedtools intersect -u -s -a eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.bed -b merged_by_strand_nat_sorted_hg38_pADB_v3_matched_to_ensembl_101_LAST_EXONS_STRANDED_plus_gencode_V32_UTR3s_minus_CDS_exons_for_space_filter.bed > eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3s.bed")
system("wc -l eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.bed")
system("wc -l eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3s.bed") # 20% of peaks
# mark these
setwd(path_data_in)
#tops_utr3 <- fread("eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_utr3s.bed")
#tops_utr3 <- fread("eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3s.bed")
names(tops_utr3) <- names_BED_std
keeps[, match_id := paste(name_bed, score_bed, sep = "_")]
keeps[, top_motif_in_cust_utr3 := ifelse(match_id %in% tops_utr3[, name], 1, 0)]
table(keeps[, top_motif_in_cust_utr3])
table(keeps[, top_motif_in_cust_utr3], keeps[, region_top], useNA = "ifany") # remember original region was based on peak, not top motif

# do same for CDS 
system(paste("cp ", path_data_cds_space, "merged_nat_sorted_hg38_kg_gencode_v32_main_chr_prot_cod_CDS_exons.bed ", path_data_in, sep = ""))
system("head merged_nat_sorted_hg38_kg_gencode_v32_main_chr_prot_cod_CDS_exons.bed")
# intersect top motifs with 3' UTR space (use -u option to just report entries with any overlap; any strand here)
setwd(path_data_in)
#system("bedtools intersect -u -a eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.bed -b merged_nat_sorted_hg38_kg_gencode_v32_main_chr_prot_cod_CDS_exons.bed > eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_cds.bed")
#system("bedtools intersect -u -a eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.bed -b merged_nat_sorted_hg38_kg_gencode_v32_main_chr_prot_cod_CDS_exons.bed > eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_cds.bed")
system("wc -l eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.bed")
system("wc -l eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_cds.bed") # 10% of peaks
# mark these
setwd(path_data_in)
#tops_cds <- fread("eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_cds.bed")
#tops_cds <- fread("eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_cds.bed")
names(tops_cds) <- names_BED_std
keeps[, top_motif_in_cust_cds := ifelse(match_id %in% tops_cds[, name], 1, 0)]
table(keeps[, top_motif_in_cust_cds])
table(keeps[, top_motif_in_cust_cds], keeps[, region_top], useNA = "ifany") # remember original region was based on peak, not top motif

# add 3'UTR intron category
system(paste("cp ", path_data_utr3_space, "gencode_v32_UTR3s_in_prot_cod_genes_main_chr_basic_INTRON_SEGMENTS.bed ", path_data_in, sep = ""))
# intersect top motifs (use -u option to just report entries with any overlap)
#system("bedtools intersect -u -s -a eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.bed -b gencode_v32_UTR3s_in_prot_cod_genes_main_chr_basic_INTRON_SEGMENTS.bed > eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_utr3_introns.bed")
#system("bedtools intersect -u -s -a eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.bed -b gencode_v32_UTR3s_in_prot_cod_genes_main_chr_basic_INTRON_SEGMENTS.bed > eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3_introns.bed")
system("wc -l eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.bed")
system("wc -l eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3_introns.bed") # 1.5% of peaks
# mark these
setwd(path_data_in)
#tops_utr3_intron <- fread("eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_utr3_introns.bed")
#tops_utr3_intron <- fread("eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3_introns.bed")
names(tops_utr3_intron) <- names_BED_std
keeps[, top_motif_in_cust_utr3_intron := ifelse(match_id %in% tops_utr3_intron[, name], 1, 0)]
table(keeps[, top_motif_in_cust_utr3_intron])
table(keeps[, top_motif_in_cust_utr3_intron], keeps[, region_top], useNA = "ifany") # remember original region was based on peak, not top motif

# original breakdown
table(keeps$region_top)
# new breakdown (based on any overlap with top motif position)
table(utr3 = keeps$top_motif_in_cust_utr3, cds = keeps$top_motif_in_cust_cds)
table(utr3 = keeps$top_motif_in_cust_utr3, utr3_intron = keeps$top_motif_in_cust_utr3_intron)
# therefore get an over 10% increase / boost in number of top motifs in 3'UTRs relative to original classification based on clip region files

# tidy
rm(tops_cds, tops_utr3, tops_utr3_intron)

## add meta-position matched coordinates for mock top motifs if positional information was all that mattered for top motif relative to eCLIP peak. 
# don't allow any overlap with top motif
find_top_pos_match <- function(rbp_in, cell_in, pos_exp_in) {
  
  cands <- keeps[RBP_peak == rbp_in & cell_line_peak == cell_in, .(top_id, pos_max_aff_single_rel_five_end_exp, pos_max_aff_single_rel_five_end)]
  
  # try picking one entry to calc overlap at a time for more efficient processing...
  for(test_round in 1:1000) {
    
    winner <- cands[][sample(.N, 1)]
    winner[, overlap := sapply(pos_max_aff_single_rel_five_end_exp, function(x) ifelse(sum(x %in% pos_exp_in) > 0, T, F))]
    
    if(winner[, overlap] == F) {
      return(winner[, pos_max_aff_single_rel_five_end]) # can be used to generate coordinates outside of function
    }
    
    if(test_round == 1000) {
      return(stop("custom error: no non-overlapping match after 1000 tries"))
    }
    
  }
  
  rm(cands, elig, winner)  
  
}
  
set.seed(231990)
keeps[, pos_MATCH_rel_five_end := mcmapply(find_top_pos_match, rbp_in = RBP_peak, cell_in = cell_line_peak, pos_exp_in = pos_max_aff_single_rel_five_end_exp, mc.cores = 8L)]
hist(keeps$pos_max_aff_single_rel_five_end, breaks = 100)
hist(keeps$pos_MATCH_rel_five_end, breaks = 100) # shift because of overlap requirement
hist(keeps[, pos_MATCH_rel_five_end - pos_max_aff_single_rel_five_end], breaks = 150)

# add position-matched random motif coords
keeps[strand_peak == "+", chromStart_pos_MATCH := chromStart_top + (pos_MATCH_rel_five_end - pos_max_aff_single_rel_five_end)] # pos_MATCH_rel_five_end - pos_max_aff_single_rel_five_end > 0 is downstream of top
keeps[strand_peak == "+", chromEnd_pos_MATCH := chromEnd_top + (pos_MATCH_rel_five_end - pos_max_aff_single_rel_five_end)]
keeps[strand_peak == "-", chromStart_pos_MATCH := chromStart_top - (pos_MATCH_rel_five_end - pos_max_aff_single_rel_five_end)]
keeps[strand_peak == "-", chromEnd_pos_MATCH := chromEnd_top - (pos_MATCH_rel_five_end - pos_max_aff_single_rel_five_end)]
keeps[, .(strand_peak, pos_max_aff_single_rel_five_end, chromStart_top, chromEnd_top, pos_MATCH_rel_five_end, chromStart_pos_MATCH, chromEnd_pos_MATCH)][sample(.N, 5)]
keeps[, match_id := paste(RBP_peak, cell_line_peak, chrom_peak, chromStart_pos_MATCH, chromEnd_pos_MATCH, strand_peak, sep = "_")]

## save ("OLD" and "NEW" motifs in consecutive chunks of code)

# OLD: =========================================================
setwd(path_data_out)
# all
saveRDS(keeps, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.rds")
# update region to reflect top motif
# as bed file focused on top motifs
## all
fwrite(keeps[, .(chrom_peak, chromStart_top, chromEnd_top, name = paste(name_bed, score_bed, top_motif_in_cust_utr3, top_motif_in_cust_cds, top_motif_in_cust_utr3_intron, sep = "_"), score = paste(max_aff, pos_max_aff_single_rel_five_end, sep = "_"), strand_peak)],
       "eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.top_region.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
## in 3'UTRs
fwrite(keeps[top_motif_in_cust_utr3 == 1 & top_motif_in_cust_cds == 0, .(chrom_peak, chromStart_top, chromEnd_top, name = paste(name_bed, score_bed, top_motif_in_cust_utr3, top_motif_in_cust_cds, top_motif_in_cust_utr3_intron, sep = "_"), score = paste(max_aff, pos_max_aff_single_rel_five_end, sep = "_"), strand_peak)],
       "eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
## merge for space filtering
system("sort -V eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed > nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed")
system("bedtools merge -i nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed > merged.nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed")
system("rm nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed")

# as bed file focused on position-matched CONTROL motifs
## all
fwrite(keeps[, .(chrom_peak, chromStart_pos_MATCH, chromEnd_pos_MATCH, name = name_bed, score = paste(pos_max_aff_single_rel_five_end, pos_MATCH_rel_five_end, sep = "_"), strand_peak)],
       "eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.pos_matched_CONTROL_motifs.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
## in 3'UTRs (based on top motif)
fwrite(keeps[top_motif_in_cust_utr3 == 1 & top_motif_in_cust_cds == 0, .(chrom_peak, chromStart_pos_MATCH, chromEnd_pos_MATCH, name = name_bed, score = paste(pos_max_aff_single_rel_five_end, pos_MATCH_rel_five_end, sep = "_"), strand_peak)],
       "eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
## merge for space filtering
system("sort -V eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed > nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed")
system("bedtools merge -i nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed > merged.nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed")
system("rm nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed")

# NEW: =========================================================
setwd(path_data_out)
# all
saveRDS(keeps, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.rds")
# update region to reflect top motif
# as bed file focused on top motifs
## all
fwrite(keeps[, .(chrom_peak, chromStart_top, chromEnd_top, name = paste(name_bed, score_bed, top_motif_in_cust_utr3, top_motif_in_cust_cds, top_motif_in_cust_utr3_intron, sep = "_"), score = paste(max_aff, pos_max_aff_single_rel_five_end, sep = "_"), strand_peak)],
       "eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.top_region.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
## in 3'UTRs
fwrite(keeps[top_motif_in_cust_utr3 == 1 & top_motif_in_cust_cds == 0, .(chrom_peak, chromStart_top, chromEnd_top, name = paste(name_bed, score_bed, top_motif_in_cust_utr3, top_motif_in_cust_cds, top_motif_in_cust_utr3_intron, sep = "_"), score = paste(max_aff, pos_max_aff_single_rel_five_end, sep = "_"), strand_peak)],
       "eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
## merge for space filtering
system("sort -V eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed > nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed")
system("bedtools merge -i nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed > merged.nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed")
system("rm nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.in_utr3_ONLY.bed")
# as bed file focused on position-matched CONTROL motifs
## all
fwrite(keeps[, .(chrom_peak, chromStart_pos_MATCH, chromEnd_pos_MATCH, name = name_bed, score = paste(pos_max_aff_single_rel_five_end, pos_MATCH_rel_five_end, sep = "_"), strand_peak)],
       "eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.pos_matched_CONTROL_motifs.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
## in 3'UTRs (based on top motif)
fwrite(keeps[top_motif_in_cust_utr3 == 1 & top_motif_in_cust_cds == 0, .(chrom_peak, chromStart_pos_MATCH, chromEnd_pos_MATCH, name = name_bed, score = paste(pos_max_aff_single_rel_five_end, pos_MATCH_rel_five_end, sep = "_"), strand_peak)],
       "eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
## merge for space filtering
system("sort -V eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed > nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed")
system("bedtools merge -i nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed > merged.nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed")
system("rm nat_sorted.eCLIP_peak.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.pos_matched_CONTROL_motifs.where_top_in_utr3_ONLY.bed")

## -----------------------------------------------------------------------


## make "super summary" table of affinity quantiles for each RBP -----------------------------------------------------------------------
# separate entries for each cell line (and for union) and 3'UTR only vs. all
# two versions for affinity values: 1) maxv = max version 2) sumv = sum version, average of all positions in top motif (to match approach used when assessing selection for variants in motifs)

## generate sum data:
# based on how the affinity list is structured (corresponds to first base of motif),
# first base of motif would be: (max_pos - (k-1)):max_pos
# second base of motif would be: (max_pos - (k-1) + 1):(maxpos + 1)
# etc. 
# function
get_aff_sum <- function(aff_list_in, k_in, pos_max_in) {
  
  aff_list_in <- as.numeric(aff_list_in)
  pos <- rep(NA, k_in)
  for(i in 0:(k_in -1)) {
    index_curr <- (pos_max_in + i - (k_in - 1)):(pos_max_in + i)
    pos[i+1] <- sum(aff_list_in[index_curr])
  }
  return(mean(pos))
}
# apply
keeps[, sum_aff_mean := pbmapply(get_aff_sum, aff_list_in = affinity_list, k_in = k_size, pos_max_in = pos_max_aff_single)]

options(scipen = 0)

## A) cell lines individually
rbp_quant_summ_pair_all <- NULL

# **change top line for different regions**
rbp_quant_summ_pair <- keeps[,
#rbp_quant_summ_pair <- keeps[top_motif_in_cust_utr3 == T & top_motif_in_cust_cds == F,
                             .(peak_count = .N,
                               min_aff_maxv = min(max_aff),
                               min_aff_sumv = min(sum_aff_mean),
                               d1_maxv = quantile(max_aff, 0.1),
                               d2_maxv = quantile(max_aff, 0.2),
                               d3_maxv = quantile(max_aff, 0.3),
                               d4_maxv = quantile(max_aff, 0.4),
                               d5_maxv = quantile(max_aff, 0.5),
                               d6_maxv = quantile(max_aff, 0.6),
                               d7_maxv = quantile(max_aff, 0.7),
                               d8_maxv = quantile(max_aff, 0.8),
                               d9_maxv = quantile(max_aff, 0.9),
                               q1_maxv = quantile(max_aff, 0.25),
                               q2_maxv = quantile(max_aff, 0.5),
                               q3_maxv = quantile(max_aff, 0.75),
                               d1_sumv = quantile(sum_aff_mean, 0.1),
                               d2_sumv = quantile(sum_aff_mean, 0.2),
                               d3_sumv = quantile(sum_aff_mean, 0.3),
                               d4_sumv = quantile(sum_aff_mean, 0.4),
                               d5_sumv = quantile(sum_aff_mean, 0.5),
                               d6_sumv = quantile(sum_aff_mean, 0.6),
                               d7_sumv = quantile(sum_aff_mean, 0.7),
                               d8_sumv = quantile(sum_aff_mean, 0.8),
                               d9_sumv = quantile(sum_aff_mean, 0.9),
                               q1_sumv = quantile(sum_aff_mean, 0.25),
                               q2_sumv = quantile(sum_aff_mean, 0.5),
                               q3_sumv = quantile(sum_aff_mean, 0.75),
                               max_aff_maxv = max(max_aff),
                               max_aff_sumv = max(sum_aff_mean)
                               ),
                             by = .(RBP_peak, cell_line_peak)]

rbp_quant_summ_pair[cell_line_peak == "HepG2", cell_line_peak := "H"]
rbp_quant_summ_pair[cell_line_peak == "K562", cell_line_peak := "K"]

cols_to_round <- c("min_aff_maxv", "min_aff_sumv", "d1_maxv", "d2_maxv", "d3_maxv", "d4_maxv", "d5_maxv", "d6_maxv", "d7_maxv", "d8_maxv", "d9_maxv", "q1_maxv", "q2_maxv", "q3_maxv", "d1_sumv", "d2_sumv", "d3_sumv", "d4_sumv", "d5_sumv", "d6_sumv", "d7_sumv", "d8_sumv", "d9_sumv", "q1_sumv", "q2_sumv", "q3_sumv", "max_aff_maxv", "max_aff_sumv")
rbp_quant_summ_pair[ , (cols_to_round) := pblapply(.SD, signif, digits = 4), .SDcols = cols_to_round]

# add region
rbp_quant_summ_pair[, regions := "all"]
#rbp_quant_summ_pair[, regions := "utr3_no_cds"]

# join
rbp_quant_summ_pair_all <- rbind(rbp_quant_summ_pair_all, rbp_quant_summ_pair)
rm(rbp_quant_summ_pair)


## B) union (i.e. don't double count top motifs bound in both cell lines)
rbp_quant_summ_rbp_all <- NULL

# can just use one cell line arbitrarily in DT[i,]. It doesn't matter which since all affinity data is the same. 
keeps[, cell_count := .N, by = .(RBP_peak, chrom_peak, chromStart_top, chromEnd_top, strand_peak)]
table(keeps[, cell_count])

# **change first line for different regions**
#rbp_quant_summ_rbp <- keeps[cell_count == 1 | (cell_count == 2 & cell_line_peak == "HepG2"),
rbp_quant_summ_rbp <- keeps[(cell_count == 1 | (cell_count == 2 & cell_line_peak == "HepG2")) & top_motif_in_cust_utr3 == T & top_motif_in_cust_cds == F,
                            .(peak_count = .N,
                              min_aff_maxv = min(max_aff),
                              min_aff_sumv = min(sum_aff_mean),
                              d1_maxv = quantile(max_aff, 0.1),
                              d2_maxv = quantile(max_aff, 0.2),
                              d3_maxv = quantile(max_aff, 0.3),
                              d4_maxv = quantile(max_aff, 0.4),
                              d5_maxv = quantile(max_aff, 0.5),
                              d6_maxv = quantile(max_aff, 0.6),
                              d7_maxv = quantile(max_aff, 0.7),
                              d8_maxv = quantile(max_aff, 0.8),
                              d9_maxv = quantile(max_aff, 0.9),
                              q1_maxv = quantile(max_aff, 0.25),
                              q2_maxv = quantile(max_aff, 0.5),
                              q3_maxv = quantile(max_aff, 0.75),
                              d1_sumv = quantile(sum_aff_mean, 0.1),
                              d2_sumv = quantile(sum_aff_mean, 0.2),
                              d3_sumv = quantile(sum_aff_mean, 0.3),
                              d4_sumv = quantile(sum_aff_mean, 0.4),
                              d5_sumv = quantile(sum_aff_mean, 0.5),
                              d6_sumv = quantile(sum_aff_mean, 0.6),
                              d7_sumv = quantile(sum_aff_mean, 0.7),
                              d8_sumv = quantile(sum_aff_mean, 0.8),
                              d9_sumv = quantile(sum_aff_mean, 0.9),
                              q1_sumv = quantile(sum_aff_mean, 0.25),
                              q2_sumv = quantile(sum_aff_mean, 0.5),
                              q3_sumv = quantile(sum_aff_mean, 0.75),
                              max_aff_maxv = max(max_aff),
                              max_aff_sumv = max(sum_aff_mean)
                              ),
                            by = .(RBP_peak)]

rbp_quant_summ_rbp[, cell_line_peak := "any"]

cols_to_round <- c("min_aff_maxv", "min_aff_sumv", "d1_maxv", "d2_maxv", "d3_maxv", "d4_maxv", "d5_maxv", "d6_maxv", "d7_maxv", "d8_maxv", "d9_maxv", "q1_maxv", "q2_maxv", "q3_maxv", "d1_sumv", "d2_sumv", "d3_sumv", "d4_sumv", "d5_sumv", "d6_sumv", "d7_sumv", "d8_sumv", "d9_sumv", "q1_sumv", "q2_sumv", "q3_sumv", "max_aff_maxv", "max_aff_sumv")
rbp_quant_summ_rbp[ , (cols_to_round) := pblapply(.SD, signif, digits = 4), .SDcols = cols_to_round]

# add region
#rbp_quant_summ_rbp[, regions := "all"]
rbp_quant_summ_rbp[, regions := "utr3_no_cds"]

# join
rbp_quant_summ_rbp_all <- rbind(rbp_quant_summ_rbp_all, rbp_quant_summ_rbp)
rm(rbp_quant_summ_rbp)

## combine A and B
rbp_quant_summ <- rbind(rbp_quant_summ_rbp_all, rbp_quant_summ_pair_all)
rm(rbp_quant_summ_rbp_all, rbp_quant_summ_pair_all)

# save 
setkey(rbp_quant_summ, RBP_peak, cell_line_peak)
setcolorder(rbp_quant_summ, c("RBP_peak", "cell_line_peak", "regions", "peak_count", "min_aff_maxv", "min_aff_sumv", "d1_maxv", "d2_maxv", "d3_maxv", "d4_maxv", "d5_maxv", "d6_maxv", "d7_maxv", "d8_maxv", "d9_maxv", "q1_maxv", "q2_maxv", "q3_maxv", "d1_sumv", "d2_sumv", "d3_sumv", "d4_sumv", "d5_sumv", "d6_sumv", "d7_sumv", "d8_sumv", "d9_sumv", "q1_sumv", "q2_sumv", "q3_sumv", "max_aff_maxv", "max_aff_sumv"))
setwd(path_data_out)
fwrite(rbp_quant_summ, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.aff_quantile_summ.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# read back
setwd(path_data_out)
#rbp_quant_summ <- fread("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.aff_quantile_summ.txt")
## -----------------------------------------------------------------------


## assess relationship between position of top motifs and 5' end of peak -----------------------------------------------------------------------

## **subset peaks for faster processing (repeat for each desired space)**

## 1) all 
#peaks_sub <- copy(keeps)

## OR

## 2) 3'UTR only (excluding those with any CDS overlap)
peaks_sub <- keeps[top_motif_in_cust_utr3 == T & top_motif_in_cust_cds == F, ] 

sort(table(peaks_sub[, RBP_peak]), decreasing = T)
peaks_sub

# define quantile boundaries for each RBP cell line pair (for desired region)

# choose A, B, or C below

## A) (for deciles):
#peaks_summ_int <- peaks_sub[, .("0" = min(max_aff),
#                                "0.1" = quantile(max_aff, 0.1),
#                                "0.2" = quantile(max_aff, 0.2),
#                                "0.3" = quantile(max_aff, 0.3),
#                                "0.4" = quantile(max_aff, 0.4),
#                                "0.5" = quantile(max_aff, 0.5),
#                                "0.6" = quantile(max_aff, 0.6),
#                                "0.7" = quantile(max_aff, 0.7),
#                                "0.8" = quantile(max_aff, 0.8),
#                                "0.9" = quantile(max_aff, 0.9),
#                                "1.0" = max(max_aff)
#                                ),
#                            by = .(RBP_peak, cell_line_peak)] # by cell line / RBP pair

# B) (for quartiles):
#peaks_summ_int <- peaks_sub[, .("0" = min(max_aff),
#                                "0.25" = quantile(max_aff, 0.25),
#                                "0.5" = quantile(max_aff, 0.5),
#                                "0.75" = quantile(max_aff, 0.75),
#                                "1.0" = max(max_aff)
#                                ),
#                            by = .(RBP_peak, cell_line_peak)] # by cell line / RBP pair

# OR C) (for all data):
peaks_summ_int <- peaks_sub[, .("0" = min(max_aff),
                                "1.0" = max(max_aff)
                                ),
                            by = .(RBP_peak, cell_line_peak)] # by cell line / RBP pair

# reformat
peaks_summ <- as.data.table(gather(peaks_summ_int, "quantile", "max_aff", "0":"1.0"))
peaks_summ[, quantile := as.numeric(quantile)]
rm(peaks_summ_int)

# update min so that > min and <= max can be used later (i.e. include lowest values in lowest quantile)
peaks_summ[quantile == 0, max_aff := max_aff - 10^-10]

# make new DT with RBPs and corresponding max aff positions
# function
get_pos_by_rbp <- function(pos_in, rbp_in, cell_line_in = "any") {
  return(data.table(rbp = rep(rbp_in, length(pos_in)), pos = pos_in, cell_line = cell_line_in))
}

# set up
plot_min <- -75
plot_max <- 75

samples_to_sim <- 10000

#set.seed(939824) # currently used for all quantiles cell each and cell union, all regions
set.seed(888721) # all quantiles, 3'UTR only, cell each
#set.seed(545691) # all quantiles, 3'UTR only, cell union
#set.seed(314589) # quartiles, all regions, cell each
#set.seed(979879) # quartiles, all regions, cell union
#set.seed(890001) # quartiles, 3'UTRs, cell each
#set.seed(778214) # quartiles, 3'UTRs, cell union

# confirm desired regions present in sub
table(peaks_sub$top_motif_in_cust_utr3, useNA = "ifany")

# confirm number of quantiles being used here (e.g. 10 = deciles)
length(unique(peaks_summ$quantile)) - 1

all_quants <- round(seq(from = 0, to = 1, by = 1 / (length(unique(peaks_summ$quantile)) - 1)), 1) # need round to avoid some kind of floating point problem 
#all_quants <- round(seq(from = 0, to = 1, by = 1 / (length(unique(peaks_summ$quantile)) - 1)), 2) # use 2 instead of 1 for quantiles with more than one digit (i.e. Q1 = 0.25)


## for EACH cell line: -----------------------------------------------------------------------

# reset
pos_all <- NULL

# A) loop through RBPs
for (curr_rbp in unique(peaks_summ[, RBP_peak])) {
  
  #curr_rbp <- "RBFOX2"
  print(paste("processing ", curr_rbp, "at", Sys.time()))
  print(paste("processing RBP", which(unique(peaks_summ[, RBP_peak]) == curr_rbp), "of", length(unique(peaks_summ[, RBP_peak]))))
  
  # B) loop through cell lines
  for (curr_cell_line in unique(peaks_summ[RBP_peak == curr_rbp, cell_line_peak])) {
    
    #curr_cell_line <- "HepG2"
    print(paste("processing ", curr_cell_line))
    
    # C) loop through quantiles
    for (p in 1:(length(all_quants) - 1)) {
      
      #p <- 1
      print(paste("processing quant", p))
      
      ## determine current min and max affs
      curr_maxaff_min <- peaks_summ[RBP_peak == curr_rbp & cell_line_peak == curr_cell_line & quantile == all_quants[p], max_aff]
      
      curr_maxaff_max <- peaks_summ[RBP_peak == curr_rbp & cell_line_peak == curr_cell_line & quantile == all_quants[p + 1], max_aff]
      
      # expand position list
      curr_positions <- mapply(get_pos_by_rbp,
                               pos_in = peaks_sub[RBP_peak == curr_rbp & cell_line_peak == curr_cell_line & max_aff > curr_maxaff_min & max_aff <= curr_maxaff_max, pos_max_aff_single_rel_five_end_exp],
                               rbp_in = peaks_sub[RBP_peak == curr_rbp & cell_line_peak == curr_cell_line & max_aff > curr_maxaff_min & max_aff <= curr_maxaff_max, RBP_peak],
                               cell_line_in = peaks_sub[RBP_peak == curr_rbp & cell_line_peak == curr_cell_line & max_aff > curr_maxaff_min & max_aff <= curr_maxaff_max, cell_line_peak],
                               SIMPLIFY = F)
      
      if (length(curr_positions) > 0) {
        
        curr_positions <- rbindlist(curr_positions)
        
        # trim "extra" extended regions where there will be "decay"
        curr_positions <- curr_positions[pos >= plot_min & pos <= plot_max]
        
        # store
        curr_positions[, quantile := p]
        curr_positions[, max_aff_min := curr_maxaff_min]
        curr_positions[, max_aff_max := curr_maxaff_max]
        
        # sample the same number of positions samples_to_sim times to simulate a p value for observing a peak equal to or higher than the highest peak
        
        all_samples <- data.table(sample = 1:samples_to_sim, any_pos_above_obs = NA)
        
        curr_max <- max(table(curr_positions[, pos]) / sum(table(curr_positions[, pos])))
        curr_k <- unique(peaks_sub[RBP_peak == curr_rbp, k_size])
        curr_motif_count <- nrow(peaks_sub[RBP_peak == curr_rbp & cell_line_peak == curr_cell_line & max_aff > curr_maxaff_min & max_aff <= curr_maxaff_max, ])
        
        samples_to_add_up <- max(peaks_sub$to_add_up)
        samples_to_add_down <- max(peaks_sub$to_add_up)
        
        for (curr_sample in all_samples$sample) {
          # sample from mock eligible pos_max_aff_single_rel_five_end positions
          curr_sample_mids <- sample((plot_min - (ceiling(median(1:curr_k)) -1)):(plot_max + (ceiling(median(1:curr_k)) -1)), curr_motif_count, replace = T)
          # expand **(assume k = 11 for simplicity. since curr_max is based on proportion and not raw count above, this simplification is appropriate)**
          curr_sample_pos <- rep(curr_sample_mids, each = (samples_to_add_up + samples_to_add_down + 1)) + (samples_to_add_down*-1):(samples_to_add_up)
          # trim
          curr_sample_pos <- curr_sample_pos[curr_sample_pos >= plot_min & curr_sample_pos <= plot_max]
          # evaluate
          all_samples[sample == curr_sample, any_pos_above_obs := ifelse(max(table(curr_sample_pos) / sum(table(curr_sample_pos))) >= curr_max, T, F)]
        }
        
        curr_positions[, p_val := sum(all_samples[, any_pos_above_obs]) / nrow(all_samples)]
        
        pos_all <- rbind(pos_all, curr_positions)
        
        # tidy
        rm(curr_positions, all_samples, curr_max, curr_k, curr_motif_count)
        
      }
      
    }
    
  }
}

# add adjusted p value (use BH correction considering number of RBPs x number of quantiles tests)
pos_all_temp <- pos_all[, .N, by = .(rbp, cell_line, quantile, p_val)]
pos_all_temp[, p_val_adj := p.adjust(p = p_val, method = "BH")]
mock <- copy(pos_all_temp)
mock[which(mock[, p_val == 0, ])[1], p_val := 1/samples_to_sim]
mock[, p_val_adj := p.adjust(p = p_val, method = "BH")]
mock_p_min <- mock[p_val == 1/samples_to_sim, min(p_val_adj)]
rm(mock)
pos_all_temp[, c("p_val", "N") := NULL]
setkey(pos_all_temp, rbp, cell_line, quantile)
setkey(pos_all, rbp, cell_line, quantile)
pos_all <- pos_all_temp[pos_all]

# add peak counts (check since counting across quantiles)
pos_all[, count_by_pair := .N, by = .(rbp, cell_line)]
pos_all[, count_by_pair_quant := .N, by = .(rbp, cell_line, quantile)]
pos_all[, count_by_rbp := mean(unique(count_by_pair)), by = rbp]
pos_all[, count_by_rbp_quant := mean(unique(count_by_pair_quant)), by = .(rbp, quantile)]

# save 
setwd(path_data_out)
#fwrite(pos_all, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.all_regions.cell_each.all_quantiles.max_aff_pos_list.p_sim_10k_bh_adj.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
fwrite(pos_all, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.UTR3.cell_each.all_quantiles.max_aff_pos_list.p_sim_10k_bh_adj.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#fwrite(pos_all, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.all_regions.cell_each.quartiles.max_aff_pos_list.p_sim_10k_bh_adj.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#fwrite(pos_all, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.UTR3.cell_each.quartiles.max_aff_pos_list.p_sim_10k_bh_adj.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
## -----------------------------------------------------------------------


## for UNION of both cell lines: -----------------------------------------------------------------------

# reset first
pos_all <- NULL

# updates
peaks_summ[, max_aff_by_rbp := max(max_aff), by = .(RBP_peak, quantile)]
peaks_sub[, pos_max_aff_single_rel_five_end_for_union := round(mean(pos_max_aff_single_rel_five_end)), by = .(RBP_peak, chrom_peak, chromStart_top, chromEnd_top, strand_peak)] # round to avoid .5s
peaks_sub[, pos_max_aff_single_rel_five_end_exp_for_union := mapply(function(pos, up, down) seq(pos - up, pos + down, by = 1), pos = pos_max_aff_single_rel_five_end_for_union, up = to_add_up, down = to_add_down, SIMPLIFY = F)]

# A) loop through RBPs
for (curr_rbp in unique(peaks_summ[, RBP_peak])) {

  print(paste("processing ", curr_rbp, "at", Sys.time()))
  print(paste("processing RBP", which(unique(peaks_summ[, RBP_peak]) == curr_rbp), "of", length(unique(peaks_summ[, RBP_peak]))))
  
  # B) loop through quantiles
  for (p in 1:(length(all_quants) - 1)) {
    
    #p <- 1
    print(paste("processing quant", p))
    
    ## determine current min and max affs
    curr_maxaff_min <- unique(peaks_summ[RBP_peak == curr_rbp & quantile == all_quants[p], max_aff_by_rbp])
    
    curr_maxaff_max <- unique(peaks_summ[RBP_peak == curr_rbp & quantile == all_quants[p + 1], max_aff_by_rbp])
    
    # expand position list (**arbitrarily use HepG2 data for cell_count == 2 now that same averaged position lists have been generated for each top motif**)
    curr_positions <- mapply(get_pos_by_rbp,
                             pos_in = peaks_sub[RBP_peak == curr_rbp & max_aff > curr_maxaff_min & max_aff <= curr_maxaff_max & (cell_count == 1 | (cell_count == 2 & cell_line_peak == "HepG2")), pos_max_aff_single_rel_five_end_exp_for_union],
                             rbp_in = peaks_sub[RBP_peak == curr_rbp & max_aff > curr_maxaff_min & max_aff <= curr_maxaff_max & (cell_count == 1 | (cell_count == 2 & cell_line_peak == "HepG2")), RBP_peak],
                             SIMPLIFY = F)
    
    if (length(curr_positions) > 0) {
      
      curr_positions <- rbindlist(curr_positions)
      
      # trim "extra" extended regions where there will be "decay"
      curr_positions <- curr_positions[pos >= plot_min & pos <= plot_max]
      
      # store
      curr_positions[, quantile := p]
      curr_positions[, max_aff_min := curr_maxaff_min]
      curr_positions[, max_aff_max := curr_maxaff_max]
      
      # sample the same number of positions samples_to_sim times to simulate a p value for observing a peak equal to or higher than the highest peak
      
      all_samples <- data.table(sample = 1:samples_to_sim, any_pos_above_obs = NA)
      
      curr_max <- max(table(curr_positions[, pos]) / sum(table(curr_positions[, pos])))
      curr_k <- unique(peaks_sub[RBP_peak == curr_rbp, k_size])
      curr_motif_count <- nrow(peaks_sub[RBP_peak == curr_rbp & max_aff > curr_maxaff_min & max_aff <= curr_maxaff_max & (cell_count == 1 | (cell_count == 2 & cell_line_peak == "HepG2")), ])
      
      samples_to_add_up <- max(peaks_sub$to_add_up)
      samples_to_add_down <- max(peaks_sub$to_add_up)
      
      for (curr_sample in all_samples$sample) {
        # sample from mock eligible pos_max_aff_single_rel_five_end positions
        curr_sample_mids <- sample((plot_min - (ceiling(median(1:curr_k)) -1)):(plot_max + (ceiling(median(1:curr_k)) -1)), curr_motif_count, replace = T)
        # expand **(assume k = 11 for simplicity. since curr_max is based on proportion and not raw count above, this simplification is appropriate)**
        curr_sample_pos <- rep(curr_sample_mids, each = (samples_to_add_up + samples_to_add_down + 1)) + (samples_to_add_down*-1):(samples_to_add_up)
        # trim
        curr_sample_pos <- curr_sample_pos[curr_sample_pos >= plot_min & curr_sample_pos <= plot_max]
        # evaluate
        all_samples[sample == curr_sample, any_pos_above_obs := ifelse(max(table(curr_sample_pos) / sum(table(curr_sample_pos))) >= curr_max, T, F)]
      }
      
      curr_positions[, p_val := sum(all_samples[, any_pos_above_obs]) / nrow(all_samples)]
      
      pos_all <- rbind(pos_all, curr_positions)
      
      # tidy
      rm(curr_positions, all_samples, curr_max, curr_k, curr_motif_count)
      
    }
    
  }
}

# add adjusted p value (use BH correction considering number of RBPs x number of quantiles tests)
pos_all_temp <- pos_all[, .N, by = .(rbp, cell_line, quantile, p_val)]
pos_all_temp[, p_val_adj := p.adjust(p = p_val, method = "BH")]
mock <- copy(pos_all_temp)
mock[which(mock[, p_val == 0, ])[1], p_val := 1/samples_to_sim]
mock[, p_val_adj := p.adjust(p = p_val, method = "BH")]
mock_p_min <- mock[p_val == 1/samples_to_sim, min(p_val_adj)]
rm(mock)
pos_all_temp[, c("p_val", "N") := NULL]
setkey(pos_all_temp, rbp, cell_line, quantile)
setkey(pos_all, rbp, cell_line, quantile)
pos_all <- pos_all_temp[pos_all]

# add peak counts (check since counting across quantiles)
pos_all[, count_by_pair := .N, by = .(rbp, cell_line)]
pos_all[, count_by_pair_quant := .N, by = .(rbp, cell_line, quantile)]
pos_all[, count_by_rbp := mean(unique(count_by_pair)), by = rbp]
pos_all[, count_by_rbp_quant := mean(unique(count_by_pair_quant)), by = .(rbp, quantile)]

# save 
setwd(path_data_out)
#fwrite(pos_all, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.all_regions.cell_union.all_quantiles.max_aff_pos_list.p_sim_10k_bh_adj.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
fwrite(pos_all, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.UTR3.cell_union.all_quantiles.max_aff_pos_list.p_sim_10k_bh_adj.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#fwrite(pos_all, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.all_regions.cell_union.quartiles.max_aff_pos_list.p_sim_10k_bh_adj.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#fwrite(pos_all, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.UTR3.cell_union.quartiles.max_aff_pos_list.p_sim_10k_bh_adj.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
## -----------------------------------------------------------------------

## save summary for congruent filtering
compat_summ <- pos_all[, .(quantile_for_p = "all"), by = .(RBP_peak = rbp, cell_line_peak = cell_line, count_by_pair, count_by_rbp, max_aff_min, max_aff_max, p_val_adj)]
compat_summ[, regions := "utr3_no_cds"]
setwd(path_data_out)
fwrite(compat_summ, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.UTR3.cell_union.all_quantiles.p_sim_10k_bh_adj.compat_summ.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

## start of plotting section:

## read in data to plot...
# read back
setwd(path_data_out)
pos_all <- fread("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.UTR3.cell_union.all_quantiles.max_aff_pos_list.p_sim_10k_bh_adj.txt")
#pos_all <- fread("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_75.top_motifs_assigned.UTR3.cell_each.all_quantiles.max_aff_pos_list.p_sim_10k_bh_adj.txt")


# summarize by rbp / cell line and calc prop
pos_all_summ <- pos_all[, .(count = .N), by = .(rbp, pos, count_by_pair, count_by_pair_quant, count_by_rbp, count_by_rbp_quant, cell_line, max_aff_min, max_aff_max, quantile, p_val_adj)]

# add missing (zero) data
for (curr_rbp in unique(pos_all_summ$rbp)) {
  print(curr_rbp)
  
  for (curr_cell_line in unique(pos_all_summ[rbp == curr_rbp, cell_line])){
    print(curr_cell_line)
    
    for (curr_quant in 1:(length(unique(peaks_summ$quantile)) - 1)) {
      #print(curr_quant)
      
      for (curr_pos in plot_min:plot_max) {
        #print(curr_pos)
        
        if(nrow(pos_all_summ[rbp == curr_rbp & cell_line == curr_cell_line & quantile == curr_quant & pos == curr_pos, ]) == 0) {
          to_add <- unique(pos_all_summ[rbp == curr_rbp & cell_line == curr_cell_line & quantile == curr_quant, .(rbp, pos = curr_pos, count_by_pair, count_by_pair_quant, count_by_rbp, count_by_rbp_quant, cell_line, max_aff_min, max_aff_max, quantile, p_val_adj, count = 0)])
          pos_all_summ <- rbind(pos_all_summ, to_add)
          print("row added!")
        }
        
      }
    }
  }
}
# add prop
pos_all_summ[, prop := count / count_by_pair_quant]
# all props sum to 1?
pos_all_summ[, sum(prop), by = .(rbp, cell_line)]
pos_all_summ[, cell_line := factor(as.character(cell_line), levels = sort(unique(cell_line)))]

# factor RBP to order by number of peaks (based on desired data)
temp_sub <- pos_all_summ
temp_levels <- unique(temp_sub[order(-rank(count_by_rbp), rbp), rbp])
pos_all_summ[, rbp := factor(as.character(rbp), levels = temp_levels)]

# add p_val_adj string for plot label
pos_all_summ[, p_val_adj_LABEL := as.character(signif(p_val_adj, 3))]
pos_all_summ[p_val_adj == 0, p_val_adj_LABEL := paste("P <", signif(mock_p_min, 2))]
pos_all_summ[p_val_adj != 0, p_val_adj_LABEL := paste("P =", signif(p_val_adj, 2))]
unique(pos_all_summ[, .(rbp, cell_line, p_val_adj, p_val_adj_LABEL)])

## PLOT ALL DATA
plot_colours <- c("#1C1C1BFF", "#CE4A7EFF")

pos_all_summ[, cell_line_count := length(unique(cell_line)), by = rbp]

pos_all_summ[, facet_label := NULL]

pos_all_summ[cell_line_count == 1, facet_label := paste0(rbp,
                                                         "\n", p_val_adj_LABEL), by = rbp]

pos_all_summ[cell_line_count == 2, facet_label := paste0(unique(rbp),
                                                         "\n", unique(p_val_adj_LABEL[cell_line == "HepG2"]),
                                                         "\n", unique(p_val_adj_LABEL[cell_line == "K562"])), by = rbp]

# additional adjust for union version
#pos_all_summ[, facet_label := paste(facet_label, "\n ")]

rbp_temp <- as.data.table(table(peaks_sub[, .(RBP_peak, cell_line_peak)]))
setorder(rbp_temp, -N)
rbps_to_sel <- unique(rbp_temp[N >= 10, RBP_peak]) # at least 10 peaks in at least one cell line

# order by average count
rbp_order <- rbp_temp[, .(N = mean(N)), by = RBP_peak]
setorder(rbp_order, -N)
rbp_order <- rbp_order$RBP_peak[rbp_order$RBP_peak %in% rbps_to_sel]
rm(rbp_temp)

pos_all_summ[, sel_rbps := factor(rbp, levels = rbp_order)]
pos_all_summ[, facet_labels_sel_rbps := factor(facet_label, levels = unique(pos_all_summ[order(sel_rbps), facet_label]))]

pos_plot <- ggplot(pos_all_summ[sel_rbps %in% rbps_to_sel, ], aes(x = pos, y = prop, colour = cell_line, fill = cell_line)) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank()) +
  geom_bar(stat = "identity", position = "dodge", size = 0, alpha = 0.75) +
  geom_line(size = 0.5, alpha = 0.75) +
  geom_vline(xintercept = 0, size = 0.5, linetype = "solid", color = "gray20") +
  geom_hline(yintercept = 1 / (plot_max - plot_min), size = 0.5, linetype = "dotted", color = "gray20") +
  coord_cartesian(ylim = c(0, 0.035)) +
  #scale_fill_manual(values = plot_colours) + 
  scale_fill_manual(values = plot_colours, labels = "HepG2\nOR\nK562") + 
  #scale_colour_manual(values = plot_colours) +
  scale_colour_manual(values = plot_colours, labels = "HepG2\nOR\nK562") + 
  #scale_x_continuous(breaks = seq(-100, 100, by = 100), minor_breaks = seq(-100, 100, by = 25)) +
  scale_x_continuous(breaks = seq(plot_min, plot_max, by = 50), minor_breaks = seq(plot_min, plot_max, by = 25)) +
  labs(x = "\nPosition relative to 5' end of peak", y = "Proportion of top affinity sites at position\n") +
  ggtitle("Location of highest affinity RBPamp site at eCLIP peak loci in 3' UTRs\n")

pos_plot_facet <- pos_plot + facet_wrap(~ facet_labels_sel_rbps, ncol = 6, scales = "fixed")
pos_plot_facet

# save plot
# top_aff_unique_assign.pos_prop.all_RBPs_w_new_RBPamp.ALL_DATA.UTR3.pdf # 8.5 x 14
# top_aff_unique_assign.pos_prop.all_RBPs_w_new_RBPamp.CELL_EACH.ALL_DATA.UTR3.SORTED.pdf # 8.5 x 14
# top_aff_unique_assign.pos_prop.all_RBPs_w_new_RBPamp.CELL_UNION.ALL_DATA.UTR3.SORTED.pdf # 8.5 x 14

# save RBP cell line pairs to include in overall analysis
#stat_summ <- NULL
stat_summ <- rbind(stat_summ, unique(pos_all_summ[, .(rbp, cell_line, num_quants = length(unique(peaks_summ$quantile)) - 1, quantile, p_val_adj)]))
setwd(path_data_out)
fwrite(stat_summ, "", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

## PLOT by RBP
for (curr_rbp in unique(pos_all_summ[, rbp])) {
  #for (curr_rbp in "RBFOX2") {
  
  print(curr_rbp)
  
  # add title / faceting label
  if (length(unique(pos_all_summ[rbp == curr_rbp, cell_line])) == 2) {
    
    pos_all_summ[rbp == curr_rbp,
                 facet_label := paste0("D", quantile,
                                       "\n", "affinity: ", signif(max_aff_min, 2), " to ", signif(max_aff_max, 2),
                                       "\n", p_val_adj_LABEL[cell_line == "HepG2"],
                                       "\n", p_val_adj_LABEL[cell_line == "K562"] 
                 ),
                 by = .(quantile)]
    
    # factor for order
    pos_all_summ[, facet_label := factor(facet_label, levels = unique(pos_all_summ[order(quantile),facet_label])[is.na(unique(pos_all_summ[order(quantile),facet_label])) == F])]
    
  } else {
    
    pos_all_summ[rbp == curr_rbp,
                 facet_label := paste0("D", quantile,
                                       "\n", "affinity: ", signif(max_aff_min, 2), " to ", signif(max_aff_max, 2),
                                       "\n", p_val_adj_LABEL 
                 ),
                 by = .(quantile)]
    
    # factor for order
    pos_all_summ[, facet_label := factor(facet_label, levels = unique(pos_all_summ[order(quantile),facet_label])[is.na(unique(pos_all_summ[order(quantile),facet_label])) == F])]
    
  }
  
  # set colours
  if(length(unique(pos_all_summ[rbp == curr_rbp, cell_line])) == 2) {
    plot_colours <- c("#1C1C1BFF", "#CE4A7EFF")
  } else if (length(unique(pos_all_summ[rbp == curr_rbp, cell_line])) == 1) {
    if (unique(pos_all_summ[rbp == curr_rbp, cell_line]) == "HepG2") {
      plot_colours <- c("#1C1C1BFF")
    } else if (unique(pos_all_summ[rbp == curr_rbp, cell_line]) == "K562") {
      plot_colours <- c("#CE4A7EFF")
    }
  }
  
  pos_plot <- ggplot(pos_all_summ[rbp == curr_rbp, ], aes(x = pos, y = prop, colour = cell_line, fill = cell_line)) +
    theme_minimal(base_size = 11) +
    theme(legend.title = element_blank()) +
    geom_bar(stat = "identity", position = "dodge", size = 0, alpha = 0.75) +
    geom_line(size = 0.5, alpha = 0.75) +
    geom_vline(xintercept = 0, size = 0.5, linetype = "solid", color = "gray20") +
    geom_hline(yintercept = 0.005, size = 0.5, linetype = "dotted", color = "gray20") +
    coord_cartesian(ylim =c(0, NA)) +
    scale_fill_manual(values = plot_colours) +
    scale_colour_manual(values = plot_colours) +
    #scale_x_continuous(breaks = seq(-100, 100, by = 100), minor_breaks = seq(-100, 100, by = 25)) +
    scale_x_continuous(breaks = seq(plot_min, plot_max, by = 50), minor_breaks = seq(plot_min, plot_max, by = 25)) +
    labs(x = "\nPosition relative to 5' end of peak", y = "Proportion of top affinity sites at position\n") +
    ggtitle(paste("Location of highest affinity RBPamp site at", curr_rbp, "eCLIP peak loci\n"))
  
  
  pos_plot_facet <- pos_plot + facet_wrap(~ facet_label, ncol = 5, scales = "fixed")
  #pos_plot + facet_wrap(~ rbp, ncol = 6, scales = "free")
  #pos_plot + facet_wrap(~ rbp, ncol = 6, scales = "free_x")
  
  # for min and max at first and last quantile
  #ggsave(paste("pos_prop_by_pos_by_decile_", curr_rbp, "_", unique(peaks[RBP_peak == curr_rbp, avg_peak_count_by_rbp]), "_avg_peaks.pdf", sep = ""), plot = pos_plot_facet, path = path_plot_out, width = 14, height = 8.5, units = "in")
  # for average across all quantiles
  ggsave(paste("pos_prop_by_pos_by_full_avg_decile_", curr_rbp, "_", unique(peaks[RBP_peak == curr_rbp, avg_peak_count_by_rbp]), "_avg_peaks.pdf", sep = ""), plot = pos_plot_facet, path = path_plot_out, width = 14, height = 8.5, units = "in")
  
  # tidy
  pos_all_summ[, facet_label := NULL]
  
}


## PLOT all RBPs


# pick thresholds
# average (middle) of lowest affinity quantile that has p < 0.01 (in both cell lines if assayed in both)
summ_for_thresh <- pos_all_summ[, .(.N), by = .(rbp, cell_line, p_val_adj, quantile, max_aff_min, max_aff_max)]
summ_for_thresh[, all_p_under := ifelse(sum(p_val_adj < 0.01) == .N, T, F), by = .(rbp, quantile)]
summ_for_thresh[all_p_under == T, quant_min_sig := ifelse(quantile == min(quantile), T, F), by = .(rbp)]
summ_for_thresh <- summ_for_thresh[quant_min_sig == T, ]
summ_for_thresh[, c("p_val_adj", "cell_line", "N", "all_p_under", "quant_min_sig") := NULL]
summ_for_thresh <- unique(summ_for_thresh)
summ_for_thresh[, max_aff_mid := round((max_aff_min + max_aff_max) / 2, 3)]
summ_for_thresh[, c("max_aff_min", "max_aff_max") := NULL]
# missing RBPs:
to_add <- data.table(rbp = as.character(unique(pos_all_summ$rbp)[!(unique(pos_all_summ$rbp) %in% summ_for_thresh$rbp)]), quantile = NA, max_aff_mid = NA)
summ_for_thresh <- rbind(summ_for_thresh, to_add)
rm(to_add)
# add aff_quant_min_sig to peaks
sub_to_merge <- summ_for_thresh[, .(rbp, aff_quant_min_sig = max_aff_mid)]
setkey(sub_to_merge, rbp)
setkey(peaks, RBP_peak)
peaks <- sub_to_merge[peaks]
rm(sub_to_merge)
# summarize peak data
peaks[, peaks_above_thresh_by_pair := sum(max_aff >= aff_quant_min_sig), by = .(rbp, cell_line_peak)]
peaks[, prop_above_thresh_by_pair := peaks_above_thresh_by_pair / peak_count_by_pair]
temp <- unique(peaks[, .(rbp, cell_line_peak, prop_above_thresh_by_pair)])
temp_summ <- temp[, .(prop_above_thresh_mean = round(mean(prop_above_thresh_by_pair), 3)), by = rbp]
rm(temp)
# add back
setkey(summ_for_thresh, rbp)
setkey(temp_summ, rbp)
summ_for_thresh <- temp_summ[summ_for_thresh]
setcolorder(summ_for_thresh, c("rbp", "quantile", "max_aff_mid", "prop_above_thresh_mean"))
summ_for_thresh
View(summ_for_thresh)

# save
setwd(path_data_out)
fwrite(summ_for_thresh, "SUMM.eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_100.full_avg_deciles.max_aff_pos.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


# note: missing RBPs have no significant quantile
# see if these RBPs have sig result overall... 

temp_levels <- summ_for_thresh[order(max_aff_mid, decreasing = T), rbp]
summ_for_thresh[, rbp := factor(as.character(rbp), levels = temp_levels)]

ggplot(summ_for_thresh, aes(x = rbp, y = max_aff_mid)) +
  geom_point(size = 4, colour = "dodger blue") +
  theme_minimal(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "\nRBP", y = "Minimum relative affinity\n") +
  ggtitle("Motif position analysis")

temp_levels <- summ_for_thresh[order(prop_above_thresh_mean, decreasing = T), rbp]
summ_for_thresh[, rbp := factor(as.character(rbp), levels = temp_levels)]

ggplot(summ_for_thresh, aes(x = rbp, y = prop_above_thresh_mean)) +
  geom_bar(fill = "dodger blue", colour = "dodger blue", stat = "identity") +
  theme_minimal(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "\nRBP", y = "Proportion of eCLIP peaks above threshold\n") +
  ggtitle("Motif position analysis")

# save
# pos_based_thresholds.min_aff.geom_point.pdf
# pos_based_thresholds.prop_peaks.geom_bar.pdf

