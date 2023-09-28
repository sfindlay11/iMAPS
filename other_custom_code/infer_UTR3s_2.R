
# set up ------------------------------------------------------------------

library(data.table)

names_BED_std <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

path_data_in <- ""
path_data_out <- ""

# -----------------------------------------------------------------------

# read in data
setwd(path_data_in)
features <- fread("hg38_liftover_inferred_prot_cod_UTR3s_from_polyA_DB_v3_sites_matched_to_closest_upstream_ensembl_101_prot_cod_trans_merged_exon_no_mismatch_intronic_stop_updated_sel_cols.txt", header = TRUE)
names(features) <- c("chrom", "chromStart", "chromEnd", "exon_gene_id", "exon_gene_symbol", "strand", "id_site", "mean_RPM", "PSE", "PAS_Signal", "conserved", "extended_prop", "last_cds_exon_ever", "intronic_no_stop")
features

# total stop codons / five prime ends per gene
features[strand == "+", five_end := paste(chrom, chromStart, sep = "_")]
features[strand == "-", five_end := paste(chrom, chromEnd, sep = "_")]
features[, five_end_count_by_gene := length(unique(five_end)), by = exon_gene_id]
temp <- features[, .(five_end_count_by_gene = length(unique(five_end))), by = exon_gene_id]
table(temp[, five_end_count_by_gene])
table(temp[, five_end_count_by_gene]) / nrow(temp) 
rm(temp)

# rank stop codons / five prime ends (1 == most proximal to five_end_count_by_gene = most distal)
features[strand == "+", five_end_rank_by_gene := frank(chromStart, ties.method = "dense"), by = exon_gene_id]
features[strand == "-", five_end_rank_by_gene := frank(-chromEnd, ties.method = "dense"), by = exon_gene_id] # negative (-) for descending order

# annotated stop codons / five prime ends per gene (with matched pA sites; remember unmatched don't show up here; see comment a few lines down)
features[last_cds_exon_ever == 1 & strand == "+", annot_five_end := paste(chrom, chromStart, sep = "_")]
features[last_cds_exon_ever == 1 & strand == "-", annot_five_end := paste(chrom, chromEnd, sep = "_")]
features[, annot_five_end_count_by_gene := length(unique(annot_five_end[last_cds_exon_ever == 1])), by = exon_gene_id] # adds entry for each row, so better than placing last_cds_exon_ever == 1 in i
temp <- features[last_cds_exon_ever == 1, .(annot_five_end_count_by_gene = length(unique(annot_five_end))), by = exon_gene_id]
table(temp[, annot_five_end_count_by_gene])
table(temp[, annot_five_end_count_by_gene]) / nrow(temp)
# 1 = 78%, 2 = 17%, 3 = 4%, 4 = 1 % (max is 11)
# 1 = 84%, 2 = 14%, 3 = 2%, 4 = 0.4%, (max is 6) for matching to 75
# for reference: annotated stop codons per gene in ensembl 75 is: 1 = 55%, 2 = 24%, 3 = 11%, 4 = 5% (see assign_pA_sites...)
# this suggests that maybe many stop codons are very close together?

# rank annotated stop codons / five prime ends (1 == most proximal to annot_five_end_count_by_gene = most distal)
features[last_cds_exon_ever == 1 & strand == "+", annot_five_end_rank_by_gene := frank(chromStart, ties.method = "dense"), by = exon_gene_id]
features[last_cds_exon_ever == 1 & strand == "-", annot_five_end_rank_by_gene := frank(-chromEnd, ties.method = "dense"), by = exon_gene_id] # negative (-) for descending order

# NOVEL internal / intronic / un-annotated stop codons / five prime ends per gene
features[last_cds_exon_ever == 0 & strand == "+", internal_five_end := paste(chrom, chromStart, sep = "_")]
features[last_cds_exon_ever == 0 & strand == "-", internal_five_end := paste(chrom, chromEnd, sep = "_")]
features[, internal_five_end_count_by_gene := length(unique(internal_five_end[last_cds_exon_ever == 0])), by = exon_gene_id] # adds entry for each row, so better than placing last_cds_exon_ever == 0 in i
temp <- features[last_cds_exon_ever == 0, .(internal_five_end_count_by_gene = length(unique(internal_five_end))), by = exon_gene_id]
table(temp[, internal_five_end_count_by_gene])
table(temp[, internal_five_end_count_by_gene]) / nrow(temp)
# 1 = 38%, 2 = 25%, 3 = 15%, 4 = 8%, (max is 21)
# 1 = 37%, 2 = 24%, 3 = 15%, 4 = 9%, (max is 21) (mathcing to ens 75)
rm(temp)

# rank NOVEL internal five prime ends (1 == most proximal to internal_five_end_count_by_gene = most distal)
features[last_cds_exon_ever == 0 & strand == "+", internal_five_end_rank_by_gene := frank(chromStart, ties.method = "dense"), by = exon_gene_id]
features[last_cds_exon_ever == 0 & strand == "-", internal_five_end_rank_by_gene := frank(-chromEnd, ties.method = "dense"), by = exon_gene_id] # negative (-) for descending order


### number of APAs / 3'UTRs per 5' end
## all 
# count
features[, APA_count_by_five_end := .N, by = five_end]
temp <- features[, .(APA_count_by_five_end = .N), by = five_end]
table(temp[, APA_count_by_five_end])
table(temp[, APA_count_by_five_end]) / nrow(temp) 
# 1 = 50%, 2 = 16%, 3 = 8%, 4 = 5%, 5= 4%
# 1 = 50%, 2 = 16%, 3 = 8%, 4 = 5%, 5 = 4% (matching to 75)
# rank (1 == most proximal, n == most distal)
features[strand == "+", APA_rank_by_five_end := frank(chromEnd, ties.method = "dense"), by = five_end]
features[strand == "-", APA_rank_by_five_end := frank(-chromStart, ties.method = "dense"), by = five_end] # negative (-) for descending order
# mark each entry as single, proximal, middle, distal, etc.
features[, APA_category := ifelse(APA_count_by_five_end == 1, "single", ifelse(APA_count_by_five_end == APA_rank_by_five_end, "distal", ifelse(APA_rank_by_five_end == 1, "proximal", "middle")))]
table(features[, last_cds_exon_ever], features[, APA_category])

## excluding UTRs above and below the cut-off lengths (to reduce potential mismatches / annotation errors and other problems)
features[, length_UTR3 := chromEnd - chromStart]
features[, len_over_12kb := ifelse(length_UTR3 > 12000, 1, 0)] # about 97.5%ile of annotated stop pADB 3'UTRs
features[, len_under_50 := ifelse(length_UTR3 < 50, 1, 0)] # about 2.5%ile
# count
features[len_over_12kb == 0 & len_under_50 == 0, APA_count_by_five_end_len_trim := .N, by = five_end]
temp <- features[len_over_12kb == 0 & len_under_50 == 0, .(APA_count_by_five_end_len_trim = .N), by = five_end]
table(temp[, APA_count_by_five_end_len_trim])
table(temp[, APA_count_by_five_end_len_trim]) / nrow(temp) 
# 1 = 49%, 2 = 16%, 3 = 8%, 4 = 5%, 5= 4% (basically the same)
# rank (1 == most proximal, n == most distal)
features[strand == "+" & len_over_12kb == 0 & len_under_50 == 0, APA_rank_by_five_end_len_trim := frank(chromEnd, ties.method = "dense"), by = five_end]
features[strand == "-" & len_over_12kb == 0 & len_under_50 == 0, APA_rank_by_five_end_len_trim := frank(-chromStart, ties.method = "dense"), by = five_end] # negative (-) for descending order
# mark each entry as single, proximal, middle, distal, etc.
features[, APA_category_len_trim := ifelse(APA_count_by_five_end_len_trim == 1, "single", ifelse(APA_count_by_five_end_len_trim == APA_rank_by_five_end_len_trim, "distal", ifelse(APA_rank_by_five_end_len_trim == 1, "proximal", "middle")))]
# replace NA with "filtered"
features[is.na(APA_category_len_trim), APA_category_len_trim := "filtered"]
table(features[, last_cds_exon_ever], features[, APA_category_len_trim])
table(features[, APA_category], features[, APA_category_len_trim], useNA = "ifany")

## annotated last
features[is.na(annot_five_end) == F, APA_count_by_annot_five_end := .N, by = annot_five_end]
temp <- features[is.na(annot_five_end) == F, .(APA_count_by_annot_five_end = .N), by = annot_five_end]
table(temp[, APA_count_by_annot_five_end])
table(temp[, APA_count_by_annot_five_end]) / nrow(temp)
# 1 = 25%, 2 = 15%, 3 = 11%, 4 = 8%, 5 = 7%
# 1 = 23%, 2 = 14%, 3 = 11%, 4 = 9%, 5 = 7% (matching to 75)
# rank
features[is.na(annot_five_end) == F & strand == "+", APA_rank_by_annot_five_end := frank(chromEnd, ties.method = "dense"), by = annot_five_end]
features[is.na(annot_five_end) == F & strand == "-", APA_rank_by_annot_five_end := frank(-chromStart, ties.method = "dense"), by = annot_five_end] # negative (-) for descending order

## intronic / unannotated
features[is.na(internal_five_end) == F, APA_count_by_internal_five_end := .N, by = internal_five_end]
temp <- features[is.na(internal_five_end) == F, .(APA_count_by_internal_five_end = .N), by = internal_five_end]
table(temp[, APA_count_by_internal_five_end])
table(temp[, APA_count_by_internal_five_end]) / nrow(temp) # 1 = 70%, 2 = 17%, 3 = 6%, 4 = 3%, 5 = 1.5%
# rank
features[is.na(internal_five_end) == F & strand == "+", APA_rank_by_internal_five_end := frank(chromEnd, ties.method = "dense"), by = internal_five_end]
features[is.na(internal_five_end) == F & strand == "-", APA_rank_by_internal_five_end := frank(-chromStart, ties.method = "dense"), by = internal_five_end] # negative (-) for descending order

#! mark highest RPM 3'UTR per stop codon and per gene in pADB samples 
features[, high_by_gene := ifelse(mean_RPM == max(mean_RPM), 1, 0), by = exon_gene_id]
features[, high_by_five_end := ifelse(mean_RPM == max(mean_RPM), 1, 0), by = five_end]
features[len_over_12kb == 0 & len_under_50 == 0, high_by_gene_len_trim := ifelse(mean_RPM == max(mean_RPM), 1, 0), by = exon_gene_id]
features[len_over_12kb == 0 & len_under_50 == 0, high_by_five_end_len_trim := ifelse(mean_RPM == max(mean_RPM), 1, 0), by = five_end]


## flag genes where upstream 3'UTRs are very close to the stop codon. It is difficult to assign pA sites confidently to one stop codon or another for these genes.
# for each five_end, report most downstream pA site
features[strand == "+", most_distal_coord := max(chromEnd), by = five_end]
features[strand == "-", most_distal_coord := min(chromStart), by = five_end]
# list all by gene
features[, most_distal_coord_list := list(list(unique(most_distal_coord))), by = exon_gene_id]
# add col for coordinate of five end
features[strand == "+", five_end_coord := chromStart]
features[strand == "-", five_end_coord := chromEnd]
# list all by gene
features[, five_end_coord_list := list(list(unique(five_end_coord))), by = exon_gene_id]
# then calculate distance between closest five end and upstream pA site (or vice versa) by five_end
get_closest_up_d <- function(distal_entry, distals, five_end_entry, strand) {
  if (strand == "+") {
    if (length(distals[distals < distal_entry]) == 0) {
      return(NA)
    } else {
      return(min(five_end_entry - distals[distals < distal_entry]))
    }
  } else if (strand == "-") {
    if (length(distals[distals > distal_entry]) == 0) {
      return(NA)
    } else {
      return(min(distals[distals > distal_entry] - five_end_entry))
    }
  }
}
features[, up_dist := mapply(get_closest_up_d, distal_entry = most_distal_coord, distals = most_distal_coord_list, five_end_entry = five_end_coord, strand = strand)]

get_closest_down_d <- function(distal_entry, five_ends, five_end_entry, strand) {
  if (strand == "+") {
    if (length(five_ends[five_ends > five_end_entry]) == 0) {
      return(NA)
    } else {
      return(min(five_ends[five_ends > five_end_entry] - distal_entry))
    }
  } else if (strand == "-") {
    if (length(five_ends[five_ends < five_end_entry]) == 0) {
      return(NA)
    } else {
      return(min(distal_entry - five_ends[five_ends < five_end_entry]))
    }
  }
}
features[, down_dist := mapply(get_closest_down_d, distal_entry = most_distal_coord, five_ends = five_end_coord_list, five_end_entry = five_end_coord, strand = strand)]

# set distance threshold from pA site to adjacent stop / five end
#! run for each cut dist
cut_dist <- 100
cut_dist <- 250
cut_dist <- 500
cut_dist <- 1000
cut_dist <- 2000

# flag by gene
genes_to_flag <- unique(features[up_dist < cut_dist | down_dist < cut_dist, exon_gene_id])
# number of genes
length(genes_to_flag)
# prop of genes
length(genes_to_flag) / length(unique(features[, exon_gene_id]))

# 100 filters 11% of genes
# 250 filters 16% of genes
# 500 filters 20% of genes
# 1000 filters 25% of genes
# 2000 filters 33% of genes

# add flags
features[, eval(paste("within", cut_dist, "filter", sep = "_")) := ifelse(exon_gene_id %in% genes_to_flag, 1, 0)]

## write to file for visualization in browser
setwd(path_data_out)
fwrite(features[, .(chrom, chromStart, chromEnd, strand)], "for_vis_101.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

## write full file (removing a few columns)
setwd(path_data_out)
fwrite(features[, -c("most_distal_coord", "most_distal_coord_list", "five_end_coord", "five_end_coord_list")], "hg38_liftover_inferred_prot_cod_UTR3s_from_polyA_DB_v3_sites_matched_to_closest_upstream_ensembl_101_prot_cod_trans_merged_exon_no_mismatch_intronic_stop_updated_sel_cols_plus_UTR3_info.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## read back
setwd(path_data_out)
features <- fread("hg38_liftover_inferred_prot_cod_UTR3s_from_polyA_DB_v3_sites_matched_to_closest_upstream_ensembl_101_prot_cod_trans_merged_exon_no_mismatch_intronic_stop_updated_sel_cols_plus_UTR3_info.txt", header = TRUE)
# explore
View(table(features[, PSE]))
summary(features[, PSE])
summary(features[, mean_RPM])

# explore
features[exon_gene_symbol == "NMT1", ]
View(features[exon_gene_symbol == "TIMP2", ])
features[exon_gene_symbol == "DNAJB1", ]


# examine relative expression among 3' UTRs that share a 5' end
summary(features[, mean_RPM])
features[, mean_RPM_rel := mean_RPM / max(mean_RPM), by = .(exon_gene_symbol, five_end)]
summary(features[, mean_RPM_rel])
plot(log10(features$mean_RPM), features$PSE)

# from pADB website
#Sample list
#Species	Sample	No. of replicates	No. of PAS reads
# Human	A549 cells	1	1,904,382
# Human	U87 cells	4	11,605,312
# Human	HEK293 cells	7	17,520,283
# Human	HeLa cells	25	78,904,647
# Human	Embroyonic stem cells	7	20,992,021
# Human	Tissue mix	3	11,224,423
# Human	MCF7 cells	4	5,836,425
# Human	293T cells	7	18,322,377
# Human	HTR8 cells	5	22,440,198
# Human	Neural stem cells	17	32,565,958
# Human	Cerebellum	4	28,862,575
# Human	JEG3 cells	5	20,985,499
# Human	LNCaP cells	6	18,700,102
# Human	MCF10A cells	1	8,133,793
# Human	MDA-MB-231 cells	1	10,349,391
# Human	MDA-MB-468 cells	1	6,670,321
# Human	Mammary epithelial cells	1	5,095,222
# Human	MeWo cells	6	21,647,307
# Human	SKBR3 cells	1	10,525,835
# Human	T474 cells	1	13,315,686

# total of 20 different samples
# total of 107 replicates
# based on the tabled results, it looks like PSE is calculated as a percent of the total replicates used (i.e. over-weighting HeLa and Neural stem cells). Tremendous source of bias as 3' UTRs with specific values are more likely to come from specific tissues.


