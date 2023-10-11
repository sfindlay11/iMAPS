# set up ------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(pbapply)
library(ggrepel)
library(RColorBrewer)

path_data_in <- ""
path_data_out <- ""
path_data_cds <- ""
path_data_pa_info <- ""
path_data_plots <- ""

names_BED_std <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

#options(scipen = 999)

options(max.print = 10000)

# -----------------------------------------------------------------------

# read in variants 
setwd(path_data_in)
variants <- fread("pADB_infered_UTR3_plus_some_gencode.variants.main_chr_prot_cod.info.allele_freq_cats.context.calib_sel_intergenic_match_OR_vova_mut_rate_by_di_con_pair.extra_info.LAST_CDS_EXON_VARS.txt", header = TRUE) 
variants

# remove cols for this analysis
variants[, c("calib_count", "i.gene_symbol") := NULL] 
variants[, c("qual", "trinuc", "quadranuc", "hexanuc", "pentanuc", "AN_afr", "AN_amr", "AN_eas", "AN_nfe", "AN_sas", "gnomAD_v3_cov_mean", "deriv_count_afr", "deriv_count_amr", "deriv_count_eas", "deriv_count_nfe", "deriv_count_sas") := NULL]
gc()

# -----------------------------------------------------------------------

# add more detailed polyA data
setwd(path_data_pa_info)
features_summ <- fread("ALL_FEATURE_INFO.lookup_table_v01.24.2022.pADB_infered_UTR3.variants.main_chr_prot_cod.txt")
# add col
features_summ[, down_is_high_by_gene_pADB := ifelse(RPM_mean_down_pADB == max(RPM_mean_down_pADB), 1, 0), by = .(gene_id)]
names_to_sub <- c("feature_string_lookup_id", names(features_summ)[!(names(features_summ) %in% names(variants))])
features_summ_SUB <- unique(features_summ[, ..names_to_sub])

# merge
setkey(variants, feature_string_lookup_id)
setkey(features_summ_SUB, feature_string_lookup_id)
variants <- features_summ_SUB[variants] # some expansions but all with feature_conflict_filter_pass == 0
rm(features_summ, features_summ_SUB)

# calc distances to nearest pA sites
variants[strand == "+", dist_to_up_pa := chromEnd - up_pa_site_pADB]
variants[strand == "-", dist_to_up_pa := up_pa_site_pADB - chromStart]
variants[strand == "+", dist_to_down_pa := down_pa_site_pADB - chromEnd]
variants[strand == "-", dist_to_down_pa := chromStart - down_pa_site_pADB]

# group single and common_all regions together 
variants[, super_APA_location_pADB := APA_location_pADB]
variants[APA_location_pADB %in% c("s", "ca"), super_APA_location_pADB := "s|ca"]
table(variants$super_APA_location_pADB)

# motif analysis ----------------------------------------------------------

# set k of interest
k_mer_size <- 6

# confirm that k-mer size of ancest_k_seq and deriv_k_seq
default_k <- ( ( unique(nchar(variants[, ancest_k_seq])) - 1) / 2) + 1
default_k

# trim k-seq to desired window size
variants[, ancest_k_seq_curr := substr(ancest_k_seq, default_k - (k_mer_size -1), default_k + (k_mer_size - 1))]
variants[, deriv_k_seq_curr := substr(deriv_k_seq, default_k - (k_mer_size -1), default_k + (k_mer_size - 1))]
variants[, .(ancest_k_seq, ancest_k_seq_curr, deriv_k_seq, deriv_k_seq_curr)]

# motif roller function 
motif_roll <- function(input_seq, k) {
  
  motif_list <- NULL
  
  for (i in 1:(nchar(input_seq) - (k - 1))) {
    motif_list[i] <- substr(input_seq, i, i + (k-1))
  }
  
  return(motif_list)
}


# apply
variants[, ancest_motifs := pblapply(ancest_k_seq_curr, motif_roll, k = 6)]
variants[, deriv_motifs := pblapply(deriv_k_seq_curr, motif_roll, k = 6)]
variants

# tidy
variants[, c("ancest_k_seq", "deriv_k_seq") := NULL]
gc()

# -----------------------------------------------------------------------


# are top two PAS 6-mers present? (AAUAAA and AUUAAA; see Ni et al BMC Genomics 2013)
variants[, ancest_top_PAS_count := pbsapply(ancest_motifs, function(x) sum(x %in% c("AATAAA", "ATTAAA")))] # this format counts each occurrence of each 6-mer
variants[, deriv_top_PAS_count := pbsapply(deriv_motifs, function(x) sum(x %in% c("AATAAA", "ATTAAA")))]
# net top PAS change
variants[, top_PAS_net := deriv_top_PAS_count - ancest_top_PAS_count]

# are any enriched PAS 6-mers present? (see Ni et al BMC Genomics 2013)
variants[, ancest_any_PAS_count := pbsapply(ancest_motifs, function(x) sum(x %in% c("AATAAA", "ATTAAA", "AAAAAA", "ATAAAA", "AAATAA", "AAAATA", "ATAAAT", "ATAAAG", "TAAAAA")))]
variants[, deriv_any_PAS_count := pbsapply(deriv_motifs, function(x) sum(x %in% c("AATAAA", "ATTAAA", "AAAAAA", "ATAAAA", "AAATAA", "AAAATA", "ATAAAT", "ATAAAG", "TAAAAA")))]
# net any PAS change
variants[, any_PAS_net := deriv_any_PAS_count - ancest_any_PAS_count]
variants

table(variants$ancest_top_PAS_count)
table(variants$deriv_top_PAS_count)
table(variants$top_PAS_net)

# in PAS region (-30 to -15)
variants[, in_pas_region := ifelse(dist_to_down_pa >= 15 & dist_to_down_pa <= 30, 1, 0)]

# save
setwd(path_data_out)
#fwrite(variants, "pADB_infered_UTR3_plus_some_gencode.variants.main_chr_prot_cod.info.allele_freq_cats.context.calib_sel_intergenic_match_OR_vova_mut_rate_by_di_con_pair.extra_info.LAST_CDS_EXON_VARS.PAS_analysis.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

## OE analysis 

# set up ------------------------------------------------------------------


OE_all <- NULL

# primary filter
variants[, prim_filter_pass := ifelse(feature_conflict_filter_pass == 1 & len_filt_pass == 1 & last_cds_exon_ever == 1 & ancest_seq_strand == ref_strand & (!chrom %in% c("chrX", "chrY", "chrM")) & gnomAD_v3_cov_med > 24 & gnomAD_v3_cov_med < 43 & AN > 129064 & ((cat == "CpG_ti" & methyl_mean_cat != "N") | cat != "CpG_ti") & complex_ancest == 0 & lcr == 0 & CpG_island == 0 & CDS_overlap == 0 & protein_coding_exon_ever == 1, TRUE, FALSE)]

# add secondary / final filters
variants[, sec_filter_pass := ifelse(prim_filter_pass == TRUE & within_100_filter_pADB == 0 & base_change_pair != "A>G", TRUE, FALSE)]


# all ---------------------------------------------------------------------

OE <- variants[sec_filter_pass == 1,
               .(
                 observed_singletons = sum(singleton),
                 expected_singletons = sum(highest_sig_k_prop_sing),
                 n_variants = .N,
                 gene_list = list(gene_id)
                 )
               ]

OE[, maps := (observed_singletons - expected_singletons) / n_variants]
OE[, descr := "all"]
OE

# specific categories ---------------------------------------------------------------------

# add descriptions (always reset first)
variants[, descr := NULL]

# TOP:
variants[sec_filter_pass == T &
           ancest_top_PAS_count > 0 & deriv_any_PAS_count == 0,
         descr := "PAS_lost"]

variants[sec_filter_pass == T &
           ancest_top_PAS_count > 0 & deriv_any_PAS_count > 0,
         descr := "PAS_pres"]

## OR ##
# ANY:
#variants[sec_filter_pass == T &
#           ancest_any_PAS_count > 0 & deriv_any_PAS_count == 0,
#         descr := "PAS_lost"]
           
#variants[sec_filter_pass == T &
#           ancest_any_PAS_count > 0 & deriv_any_PAS_count > 0,
#         descr := "PAS_pres"]

# selection

## 1) all top PAS vars split by consequence and primary (15 to 30 upstream of most used pA site, or secondary) #### ============================

OE <- variants[sec_filter_pass == 1,
#OE <- variants[sec_filter_pass == 1 & APA_location_pADB == "u" & down_is_high_by_gene_pADB == 1 & down_conserved_pADB == 1, # just to see what scale max should be
               .(
                 observed_singletons = sum(singleton),
                 expected_singletons = sum(highest_sig_k_prop_sing),
                 n_variants = .N,
                 gene_list = list(gene_id)
               ),
               by = .(descr,
                      in_pas_region = ifelse(in_pas_region == 1 & down_pa_is_high_by_gene_pADB == 1, 1, 0))
               ]

OE[, maps := (observed_singletons - expected_singletons) / n_variants]
OE

# factor
OE[, in_pas_region := factor(as.character(in_pas_region), levels = c("1", "0"))]

# plot
OE[, in_pas_region := factor(in_pas_region, levels = c(1, 0))]
ggplot(data = OE[is.na(descr) == FALSE, ], aes(x = in_pas_region, y = maps, colour = in_pas_region, shape = descr)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour = "black") + 
  geom_point(size = 2, stroke = 2) +
  theme_classic(base_size = 16) +
  guides(colour = "none") + 
  scale_colour_manual(values = c("dodger blue", "#666666")) +
  scale_shape_manual(values = c(4,1), name = "PAS impact", labels = c("lost", "preserved")) +
  #scale_shape_manual(values = c(17,1), name = "downstream pA site", labels = c("primary &\nconserved", "non-primary&\nnon-conserved")) +
  scale_x_discrete(labels = c("primary", "secondary")) +
  scale_y_continuous(breaks = seq(0, 0.07, by = 0.01), minor_breaks = NULL) +
  coord_cartesian(ylim =c(0,0.07)) +
  labs(x = "\nPAS region", y = "iMAPS\n") +
  ggtitle("\n")

fwrite(OE[is.na(descr) == FALSE, ], paste(path_data_plots, "imaps.top_pas.by_lost_pres.by_pas_region_15_30_up_of_primary.geom_point.DATA.txt", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
ggsave("imaps.top_pas.by_lost_pres.by_pas_region_15_30_up_of_primary.geom_point.PLOT.pdf", path = path_data_plots, width = 5, height = 5, units = "in")
#setwd(path_data_plots)

# stats
# subset (if applicable)
OE_curr <- OE[is.na(descr) == FALSE, ]
# add category / comparison column
OE_curr[, cat := paste(descr, in_pas_region, sep = "&")]
# add non-singleton numbers
OE_curr[, observed_others := n_variants - observed_singletons]
OE_curr[, expected_others := n_variants - expected_singletons]
# singleton to non-sing ratio for each group
OE_curr[, exp_ratio := expected_singletons / expected_others]
# assign column to compare
OE_curr[, to_comp := cat]
# change categories here as desired (use y for expected higher maps)
cat_x <- "PAS_pres&1"
#cat_x <- "PAS_lost&0"
cat_y <- "PAS_lost&1"
# calc expected OR for null testing
OR_exp_curr <- OE_curr[to_comp == cat_y, exp_ratio] / OE_curr[to_comp == cat_x, exp_ratio]
OR_exp_curr
# fisher exact test with expected OR != 1, but to the expected OR based on calibration
fisher_curr <- fisher.test(x = matrix(c(OE_curr[to_comp == cat_x, observed_others], OE_curr[to_comp == cat_y, observed_others], OE_curr[to_comp == cat_x, observed_singletons], OE_curr[to_comp == cat_y, observed_singletons]), nrow = 2, ncol = 2), or = OR_exp_curr, alternative = "greater")
fisher_curr

## 2) PAS lost vars by pA conservation status and primary / secondary #### ============================

OE <- variants[sec_filter_pass == 1 & descr == "PAS_lost",
               .(
                 observed_singletons = sum(singleton),
                 expected_singletons = sum(highest_sig_k_prop_sing),
                 n_variants = .N,
                 gene_list = list(gene_id)
               ),
               by = .(descr,
                      prim = down_pa_is_high_by_gene_pADB == 1 & in_pas_region == 1,
                      down_conserved_pADB)
]

OE[, maps := (observed_singletons - expected_singletons) / n_variants]
OE

# factor
OE[, down_conserved_pADB := factor(as.character(down_conserved_pADB), levels = c("1", "0"))]
OE[, prim := factor(prim, levels = c(TRUE, FALSE))]
OE

# plot
ggplot(data = OE, aes(x = down_conserved_pADB, y = maps, colour = prim)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour = "black") +
  #geom_vline(xintercept = 2, linetype = "dashed", size = 1, colour = "black") + 
  geom_point(size = 2, stroke = 2, shape = 4, alpha = 1, position = position_jitter(seed = 89552, width = 0.03)) +
  theme_classic(base_size = 16) +
  #theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12), axis.text.x = element_text(vjust = 0.5)) +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  scale_x_discrete(labels = c("conserved", "non\nconserved")) +
  scale_y_continuous(breaks = seq(0, 0.07, by = 0.01), minor_breaks = NULL) +
  scale_colour_manual(values = c("dodger blue", "#666666"), name = "PAS region", labels = c("primary", "secondary")) +
  coord_cartesian(ylim =c(0,0.07)) +
  labs(x = "\ndownstream poly(A) site usage", y = "iMAPS\n") + 
  ggtitle("\n")

#fwrite(OE, paste(path_data_plots, "imaps.top_pas.lost.15_30_up_of_primary.by_conservation.geom_point.DATA.txt", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
fwrite(OE, paste(path_data_plots, "imaps.top_pas.lost.by_conservation.by_primary.geom_point.DATA.txt", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#ggsave("imaps.top_pas.lost.15_30_up_of_primary.by_conservation.geom_point.PLOT.pdf", path = path_data_plots, width = 5, height = 5, units = "in")
ggsave("imaps.top_pas.lost.by_conservation.by_primary.geom_point.PLOT.pdf", path = path_data_plots, width = 5, height = 5, units = "in")
#OE[, -"gene_list"]

# stats
# subset (if applicable)
OE_curr <- OE
#OE_curr <- OE[prim == TRUE, ]
#OE_curr <- OE[prim == FALSE, ]
# add category / comparison column
#OE_curr[, cat := down_conserved_pADB]
OE_curr[, cat := paste(prim, down_conserved_pADB, sep = "&")]
# add non-singleton numbers
OE_curr[, observed_others := n_variants - observed_singletons]
OE_curr[, expected_others := n_variants - expected_singletons]
# singleton to non-sing ratio for each group
OE_curr[, exp_ratio := expected_singletons / expected_others]
# assign column to compare
OE_curr[, to_comp := cat]
# change categories here as desired (use y for expected higher maps)
#cat_x <- "0"
#cat_y <- "1"
cat_x <- "FALSE&1"
cat_y <- "TRUE&1"
# calc expected OR for null testing
OR_exp_curr <- OE_curr[to_comp == cat_y, exp_ratio] / OE_curr[to_comp == cat_x, exp_ratio]
OR_exp_curr
# fisher exact test with expected OR != 1, but to the expected OR based on calibration
fisher_curr <- fisher.test(x = matrix(c(OE_curr[to_comp == cat_x, observed_others], OE_curr[to_comp == cat_y, observed_others], OE_curr[to_comp == cat_x, observed_singletons], OE_curr[to_comp == cat_y, observed_singletons]), nrow = 2, ncol = 2), or = OR_exp_curr, alternative = "greater")
#fisher_curr <- fisher.test(x = matrix(c(OE_curr[to_comp == cat_x, observed_others], OE_curr[to_comp == cat_y, observed_others], OE_curr[to_comp == cat_x, observed_singletons], OE_curr[to_comp == cat_y, observed_singletons]), nrow = 2, ncol = 2), or = OR_exp_curr, alternative = "two.sided")
fisher_curr



