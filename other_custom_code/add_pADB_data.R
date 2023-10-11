# add data relevant to 3' UTRs / pA sites (from polyA database)

library(data.table)
library(stringr)
library(pbapply)

path_data_in <- ""
path_data_out <- ""
path_data_gen <- ""
path_data_pA <- ""

names_BED_std <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

# read in variant data
setwd(path_data_in)
variants <- fread("gnomad_genomes_r3.0_vars.in_main_chr_prot_cod.pADB_v3_ensembl_101_matched.plus_gencode_v32.UTR3.filter_pass.no_indels.info_exp.intersected.expanded.RNA_strand.ancest_deriv.k_11_seq.cov.methyl.dnase.cpgi.h3k9me3.cds.bed.gz", header = FALSE)
names(variants) <- names_BED_std
variants

# feature analysis to intersect with transcripts -----------------------------------------------------------------------

# extract features for calculation of transcript features on smaller non-redundant data
features <- data.table(feature_string = unique(variants[, score]))
nrow(features) # 151,286
head(features)

# give each feature string a lookup id
setkey(features, feature_string)
features[, feature_string_lookup_id := 1:nrow(features)]
# convert to letters
features[, feature_string_lookup_id := chartr("1234567890", "ABCDEFGHIJ", feature_string_lookup_id)]
features
# all unique?
length(unique(features[, feature_string_lookup_id])) == nrow(features)

# write to file to serve as reference 
setwd(path_data_out)
#lookup_table_file_name <- "hg38_pADB_inferred_UTR3s_plus_some_gencode_v32_UTR3s_var_feature_string_lookup_table_01.24.2022.txt"

# replace feature_string (score) with feature_string_lookup_id in variants
setkey(variants, score)
setkey(features, feature_string)
variants <- features[variants]
variants[, feature_string := NULL]
gc()

# expand name column
variants[, c("qual", "use_as_ancest", "use_as_deriv", "complex_ancest", "ancest_seq_strand", "ref_strand", "alt_strand", "locus_allele_n", "lcr", "use_as_deriv_count", "AN", "deriv_count_afr", "deriv_count_amr", "deriv_count_eas", "deriv_count_nfe", "deriv_count_sas", "AN_afr", "AN_amr", "AN_eas", "AN_nfe", "AN_sas", "ancest_k_seq", "deriv_k_seq", "gnomAD_v3_cov_mean", "gnomAD_v3_cov_med", "methyl_mean", "DNase", "CpG_island", "H3K9me3", "CDS_overlap") := tstrsplit(name, "_", fixed = TRUE)]
variants[, name := NULL]
variants
gc()

# reorder columns
setcolorder(variants, c("chrom", "chromStart", "chromEnd", "strand", "qual", "use_as_ancest", "use_as_deriv", "complex_ancest", "ancest_seq_strand", "ref_strand", "alt_strand", "locus_allele_n", "lcr", "use_as_deriv_count", "AN", "deriv_count_afr", "deriv_count_amr", "deriv_count_eas", "deriv_count_nfe", "deriv_count_sas", "AN_afr", "AN_amr", "AN_eas", "AN_nfe", "AN_sas", "ancest_k_seq", "deriv_k_seq", "gnomAD_v3_cov_mean", "gnomAD_v3_cov_med", "methyl_mean", "DNase", "CpG_island", "H3K9me3", "CDS_overlap", "feature_string_lookup_id"))

# write
setwd(path_data_out)
fwrite(variants, "gnomad_genomes_r3.0_vars.in_main_chr_prot_cod.pADB_v3_ensembl_101_matched.plus_gencode_v32.UTR3.filter_pass.no_indels.all_info.feature_string_lookup_id.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# expand features ---------------------------------------------------------

features
# expand features into list form where each entry is a data table of the original feature string (repeated) and its individual features
expand_features <- function(input_string, input_id) {
  return(data.table(feature_string_lookup_id = input_id, ind_feature = unique(str_split(input_string, fixed(","))[[1]])))
}
# apply
exp_features <- pbmapply(expand_features, input_string = features$feature_string, input_id = features$feature_string_lookup_id, SIMPLIFY = FALSE)
exp_features <- rbindlist(exp_features)
exp_features

# extract components of individual features
exp_features[, c("coord", "class", "gene_id", "gene_symbol", "info") := tstrsplit(ind_feature, "__", fixed = TRUE)]

# extract additional sub-components
exp_features[, c("chrom", "chromStart", "chromEnd", "strand") := tstrsplit(coord, "_", fixed = TRUE)]
exp_features[, c("last_cds_exon_ever", "protein_coding_exon_ever", "intronic_no_stop", "source") := tstrsplit(info, "_", fixed = TRUE)]

# calc length
exp_features[, length_UTR3 := as.integer(chromEnd) - as.integer(chromStart)]

# remove un-needed cols
exp_features[, c("info", "chrom", "chromStart", "chromEnd", "strand") := NULL]

# any feature strings with multiple genes?
exp_features[, gene_id_count_by_feature_string := length(unique(gene_id)), by = .(feature_string_lookup_id)]
table(exp_features$gene_id_count_by_feature_string)
exp_features[, gene_symbol_count_by_feature_string := length(unique(gene_symbol)), by = .(feature_string_lookup_id)]
table(exp_features$gene_symbol_count_by_feature_string)

# flag entries where there are both gencode and pADB data within a feature string
exp_features[, conflict_by_feature_string := ifelse(sum(source == "pADB") > 0 & sum(source == "gencode") > 0, 1, 0), by = feature_string_lookup_id]
table(exp_features[, conflict_by_feature_string])
table(exp_features[, conflict_by_feature_string]) / nrow(exp_features) # about 0.6% (was 1% for matching to ens 75)

# flag entries where there are both gencode and pADB data for a gene
exp_features[, conflict_by_gene_id := ifelse(sum(source == "pADB") > 0 & sum(source == "gencode") > 0, 1, 0), by = gene_id]
table(exp_features[, conflict_by_gene_id]) # none
exp_features[, conflict_by_gene_symbol := ifelse(sum(source == "pADB") > 0 & sum(source == "gencode") > 0, 1, 0), by = gene_symbol]
table(exp_features[, conflict_by_gene_symbol]) # 969
table(exp_features[, conflict_by_gene_symbol]) / nrow(exp_features) # about 0.1%

# flag entries where gene_id / gene_symbol pairs are not always consistent
exp_features[, gene_pair := paste(gene_id, gene_symbol, sep = "_")]
temp_summ <- exp_features[, length(unique(gene_pair)), by = gene_id]
temp_summ[V1 > 1, gene_id] # none
temp_summ <- exp_features[, length(unique(gene_pair)), by = gene_symbol]
length(unique(temp_summ[V1 > 1, gene_symbol])) # 11
exp_features[, conflict_by_gene_pair := ifelse(gene_symbol %in% temp_summ[V1 > 1, gene_symbol], 1, 0)]
table(exp_features$conflict_by_gene_pair) # 1029
rm(temp_summ)

# add filter column (keep in table for future filtering)
exp_features[, feature_conflict_filter_pass := ifelse(gene_id_count_by_feature_string > 1 | gene_symbol_count_by_feature_string > 1 | conflict_by_feature_string == 1 | conflict_by_gene_id == 1 | conflict_by_gene_symbol == 1 | conflict_by_gene_pair == 1, 0, 1)]
table(exp_features$feature_conflict_filter_pass) # 5530
table(exp_features$feature_conflict_filter_pass) / sum(table(exp_features$feature_conflict_filter_pass)) # 99.2% pass

# tidy
exp_features[, c("gene_id_count_by_feature_string", "gene_symbol_count_by_feature_string", "conflict_by_feature_string", "conflict_by_gene_id", "conflict_by_gene_symbol", "gene_pair", "conflict_by_gene_pair") := NULL]
setkey(exp_features, feature_string_lookup_id)
exp_features

# explore:
# what prop of data (by gene) is from pADB vs. gencode? (rough estimate, no filters)
exp_features[source == "pADB", length(unique(gene_symbol))] # 15,968 pADB
exp_features[source == "gencode", length(unique(gene_symbol))] # 3,098 gencode

# -----------------------------------------------------------------------


# add pADB data by ind_feature --------------------------------------------

# read in polyadenylation data
setwd(path_data_pA)
poly_info <- fread("hg38_liftover_inferred_prot_cod_UTR3s_from_polyA_DB_v3_sites_matched_to_closest_upstream_ensembl_101_prot_cod_trans_merged_exon_no_mismatch_intronic_stop_updated_sel_cols_plus_UTR3_info.txt", header = TRUE)
poly_info

# explore
# do distal pA tend to be most highly used? 
summ <- as.data.table(table(poly_info[last_cds_exon_ever == 1 & within_100_filter == 0, APA_category], poly_info[last_cds_exon_ever == 1 & within_100_filter == 0, high_by_gene]))
summ2 <- summ[, .(prim_count = N[V2 == 1], sec_count = N[V2 == 0]), by = V1]
summ2[, total_count := prim_count + sec_count]
summ2[, prim_ratio := prim_count / total_count]
summ2
rm(summ, summ2)

# make coord col for matching
poly_info[, id_site := NULL] # remove old (was actual pA site and not 3'UTR coords)
poly_info[, coord := paste(chrom, chromStart, chromEnd, strand, sep = "_")]

# add coordinate for pA site (downstream)
poly_info[strand == "+", down_pa_site_pADB := chromEnd]
poly_info[strand == "-", down_pa_site_pADB := chromStart]

# list all pA site coordinates by 5' end
temp_summ_plus <- poly_info[strand == "+", .(pa_site_list = (list(sort(as.integer(down_pa_site_pADB), decreasing = F))), temp_count = .N), by = .(five_end, strand)] # sort from most prox to most distal
temp_summ_minus <- poly_info[strand == "-", .(pa_site_list = (list(sort(as.integer(down_pa_site_pADB), decreasing = T))), temp_count = .N), by = .(five_end, strand)] # sort from most prox to most distal
temp_summ <- rbind(temp_summ_plus, temp_summ_minus)
rm(temp_summ_plus, temp_summ_minus)
setkey(temp_summ, five_end, strand)
setkey(poly_info, five_end, strand)
poly_info <- temp_summ[poly_info]
sum(poly_info[, temp_count == APA_count_by_five_end]) == nrow(poly_info)
poly_info[, temp_count := NULL]
rm(temp_summ)

# add upstream pA site coord for each 3' UTR
table(poly_info$APA_category)
poly_info[APA_category == "middle" | APA_category == "distal", up_pa_site_pADB := pbmapply(function(rank_in, list_in) unlist(list_in)[rank_in - 1], rank_in = APA_rank_by_five_end, list_in = pa_site_list)]
poly_info[APA_category == "single" | APA_category == "proximal", up_pa_site_pADB := NA] # line order is important
# confirm
poly_info[, .(five_end, strand, APA_category, APA_rank_by_five_end, pa_site_list, down_pa_site_pADB, up_pa_site_pADB)]

# select cols to merge (also have up_dist and down_dist to next pA sites if desired)
poly_info_sub <- poly_info[, .(coord,
                               mean_RPM_pADB = mean_RPM, 
                               PSE_pADB = PSE, # from polyA DB: "Percentage of samples with detected expression of all samples" (maybe not as straightforward as it sounds, so be careful using this)
                               PAS_hex_pADB = PAS_Signal, # polyA signal hexamer from polyA DB
                               conserved_pADB = conserved, # from polA DB
                               extended_prop_pADB = extended_prop, # from polyA DB (was the site in an extended region past the end of the annotation)
                               five_end_count_by_gene_annot = annot_five_end_count_by_gene,
                               five_end_rank_by_gene_annot = annot_five_end_rank_by_gene,
                               five_end_count_by_gene_pADB = five_end_count_by_gene, # marks the number of 3' UTR groups in the gene (a group of 3' UTRs share the same 5' end / stop codon)
                               five_end_rank_by_gene_pADB = five_end_rank_by_gene, # numbers these groups from 1 to n (1 = most proximal)
                               APA_count_by_five_end_annot = APA_count_by_annot_five_end,
                               APA_rank_by_five_end_annot = APA_rank_by_annot_five_end,
                               APA_count_by_five_end_pADB = APA_count_by_five_end, # marks the number of 3' UTRs / APAs in the group
                               APA_rank_by_five_end_pADB = APA_rank_by_five_end, # ranks from 1 to n (1 = most proximal)
                               high_by_gene_pADB = high_by_gene, # marks the isoform / APA with the highest expression in the gene
                               high_by_five_end_pADB = high_by_five_end, # marks the isoform / APA with the highest expression in the group of 3' UTRs that share the same 5' end / stop codon
                               APA_category_pADB = APA_category, 
                               five_end,
                               down_pa_site_pADB,
                               up_pa_site_pADB,
                               within_100_filter_pADB = within_100_filter, # does the group that this 3' UTR belongs to have another group within 100 bases (upstream or downstream)
                               within_250_filter_pADB = within_250_filter,
                               within_500_filter_pADB = within_500_filter,
                               within_1000_filter_pADB = within_1000_filter,
                               within_2000_filter_pADB = within_2000_filter)]

# tidy
#rm(poly_info)

# merge
setkey(poly_info_sub, coord)
setkey(exp_features, coord)
exp_features_pa <- poly_info_sub[exp_features]
rm(poly_info, poly_info_sub)
# save
setwd(path_data_out)
#fwrite(exp_features_pa, "gnomad_genomes_r3.0_vars.in_main_chr_prot_cod.pADB_v3_ensembl_101_matched.plus_gencode_v32.UTR3.exp_features_pa.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#read back
#exp_features_pa <- fread("gnomad_genomes_r3.0_vars.in_main_chr_prot_cod.pADB_v3_ensembl_101_matched.plus_gencode_v32.UTR3.exp_features_pa.txt")

# note: data not available for entries from gencode only.
# -----------------------------------------------------------------------

# summarize by feature_string_lookup_id ------------------------------------------------------------

#! add more columns here as desired

# leverage the fact that the "segment" region always corresponds to the shortest 3' UTR in the group
# use only pADB data for now, since gencode data does not have details needed

features_summ <- exp_features_pa[source == "pADB",
                                 .(RPM_mean_down_pADB = mean_RPM_pADB[length_UTR3 == min(length_UTR3)][1],
                                   RPM_mean_down_any_max_pADB = max(mean_RPM_pADB),
                                   PSE_down_pADB = PSE_pADB[length_UTR3 == min(length_UTR3)][1],
                                   PSE_down_any_max_pADB = max(PSE_pADB),
                                   PAS_hex_down_pADB = PAS_hex_pADB[length_UTR3 == min(length_UTR3)][1],
                                   PAS_hex_any_down_is_AWUAAA_pADB = ifelse(sum(PAS_hex_pADB == "AAUAAA" | PAS_hex_pADB == "AUUAAA") > 0, 1, 0),
                                   down_conserved_pADB = conserved_pADB[length_UTR3 == min(length_UTR3)][1],
                                   any_down_conserved_pADB = ifelse(sum(conserved_pADB) > 0, 1, 0),
                                   APA_location_pADB = ifelse(sum(APA_category_pADB == "single") > 0, "single", ifelse(sum(APA_category_pADB == "proximal") > 0, "common_all", ifelse(sum(APA_category_pADB == "distal") == .N, "distal_unique", "common_partial"))),
                                   down_pa_site_pADB = down_pa_site_pADB[length_UTR3 == min(length_UTR3)][1],
                                   up_pa_site_pADB = up_pa_site_pADB[length_UTR3 == min(length_UTR3)][1],
                                   down_pa_is_high_by_gene_pADB = ifelse(high_by_gene_pADB[length_UTR3 == min(length_UTR3)][1] == 1, 1, 0),
                                   down_pa_is_high_by_five_end_pADB = ifelse(high_by_five_end_pADB[length_UTR3 == min(length_UTR3)][1] == 1, 1, 0),
                                   in_high_by_gene_pADB = ifelse(sum(high_by_gene_pADB) > 0, 1, 0),
                                   in_high_by_five_end_pADB = ifelse(sum(high_by_five_end_pADB) > 0, 1, 0),
                                   down_APA_rank_by_five_end_annot = min(APA_rank_by_five_end_annot),
                                   down_APA_rank_by_five_end_pADB = min(APA_rank_by_five_end_pADB),
                                   len_utr3_short_pADB = min(length_UTR3),
                                   len_utr3_longest_pADB = max(length_UTR3),
                                   in_extended_pADB = ifelse(sum(extended_prop_pADB[length_UTR3 == min(length_UTR3)][1]) > 0, 1, 0),
                                   last_cds_exon_ever = ifelse(sum(last_cds_exon_ever == 1) > 0, 1, 0),
                                   protein_coding_exon_ever = ifelse(sum(protein_coding_exon_ever == 1) > 0, 1, 0),
                                   len_filt_pass = ifelse(max(length_UTR3) > 50000, 0, 1)), # allowing for potential filtering of very long UTRs (corresponds to top 0.5% of genes). see below for quantiles
                                 by = .(gene_id, gene_symbol, feature_string_lookup_id, feature_conflict_filter_pass, five_end_pADB = five_end, five_end_count_by_gene_annot, five_end_rank_by_gene_annot, five_end_count_by_gene_pADB, five_end_rank_by_gene_pADB, APA_count_by_five_end_annot, APA_count_by_five_end_pADB, within_100_filter_pADB, within_250_filter_pADB, within_500_filter_pADB, within_1000_filter_pADB, within_2000_filter_pADB)] # includes carry-over variables
features_summ # 145,867

# abbreviate APA_location_pADB
table(features_summ$APA_location_pADB)
features_summ[APA_location_pADB == "single", APA_location_pADB:= "s"]
features_summ[APA_location_pADB == "common_all", APA_location_pADB:= "ca"]
features_summ[APA_location_pADB == "common_partial", APA_location_pADB:= "cp"]
features_summ[APA_location_pADB == "distal_unique", APA_location_pADB:= "u"]

features_summ[, dup_count := .N, by = .(feature_string_lookup_id)]
features_summ[dup_count > 1, ]
features_summ[dup_count > 1 & feature_conflict_filter_pass == 1, ] # all are already filtered
features_summ[, dup_count := NULL]

# explore quantiles for above 
temp_summ <- exp_features[last_cds_exon_ever == 1 & protein_coding_exon_ever == 1, max(length_UTR3), by = gene_id]
round(quantile(temp_summ$V1, seq(0,1, by = 0.001)), digits = 0)
rm(temp_summ)

# save
setwd(path_data_out)
fwrite(features_summ, "ALL_FEATURE_INFO.lookup_table_v01.24.2022.pADB_infered_UTR3.variants.main_chr_prot_cod.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# merge with variants ------------------------------------------------------------

# select basic features for now. can use above file to add more detailed polyA data downstream as desired
features_summ_BASIC <- unique(features_summ[, .(feature_string_lookup_id, feature_conflict_filter_pass, len_filt_pass, gene_id, gene_symbol, APA_location_pADB, down_pa_site_pADB, down_pa_is_high_by_gene_pADB, down_pa_is_high_by_five_end_pADB, up_pa_site_pADB, within_100_filter_pADB, in_high_by_gene_pADB, in_high_by_five_end_pADB, last_cds_exon_ever, protein_coding_exon_ever)])

# merge features and gene info to variants
setkey(variants, feature_string_lookup_id)
setkey(features_summ_BASIC, feature_string_lookup_id)
variants <- features_summ_BASIC[variants]
variants

setwd(path_data_out)
fwrite(variants[last_cds_exon_ever == 1, ], "pADB_infered_UTR3.variants.main_chr_prot_cod.info.last_cds_ex_only.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#system("gzip -1 -f pADB_infered_UTR3.variants.main_chr_prot_cod.info.last_cds_ex_only.txt")


