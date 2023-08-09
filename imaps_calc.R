#### ==================== set up ==================== ####

print(paste("running iMAPS_calc.R at", Sys.time()))

user_in <- commandArgs(trailingOnly = TRUE)

# test if there are three arguments provided: if not, return an error
if (length(user_in) != 3) {
  stop("please supply a single bed file, a single allele version, and a single path to required data", call. = FALSE)
}

# required packages
packages <- c("data.table", "stringi", "ggplot2", "RColorBrewer", "R.utils")

#install all required packages that are not already installed
install.packages(setdiff(packages, rownames(installed.packages())))

# load required packages
library(data.table)
library(stringi)
library(ggplot2)
library(RColorBrewer)
library(R.utils)

names_BED_std <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")


#### ============= END OF set up ==================== ####


#### ==================== intersect ==================== ####

setwd(user_in[3])
system(paste("./feature_intersection.sh", user_in[1], user_in[2], sep = " "))

setwd(user_in[3])
ints <- fread("features_intersected.filtered.sel_cols.txt")

if (user_in[2] %in% c("sense", "plus")) {
  names(ints) <- c(names_BED_std, "user_category", "user_alleles")
} else if (user_in[2] == "none") {
  names(ints) <- c(names_BED_std, "user_category")
}

# expand name
ints[, c("gene_symbol", "ancest_seq_strand", "ref_strand", "alt_strand", "use_as_ancest", "use_as_deriv", "ancest_k_seq") := tstrsplit(name, "_", fixed = TRUE)]
ints[, name := NULL]

# expand score
ints[, c("use_as_deriv_count", "AN", "lcr", "gnomAD_v3_cov_med", "methyl_mean", "DNase", "CpG_island", "H3K9me3", "CDS_overlap", "within_100_filter_pADB", "in_high_by_gene_pADB", "down_pa_is_high_by_gene_pADB", "dist_to_down_pa", "down_conserved_pADB") := tstrsplit(score, "_", fixed = TRUE)]

# update data types
cols_to_int <- c("use_as_deriv_count", "AN", "lcr", "gnomAD_v3_cov_med", "methyl_mean", "DNase", "CpG_island", "H3K9me3", "CDS_overlap", "within_100_filter_pADB", "in_high_by_gene_pADB", "down_pa_is_high_by_gene_pADB", "dist_to_down_pa", "down_conserved_pADB")
ints[, (cols_to_int) := lapply(.SD, as.integer), .SDcols = cols_to_int]
rm(cols_to_int)

# match and filter variants
if (user_in[2] == "sense") {
  
  # extract alleles
  ints[, c("ref_user", "alt_user") := tstrsplit(user_alleles, "[|]")]
  ints[, ref_user := toupper(ref_user)]
  ints[, alt_user := toupper(alt_user)]
  
  # bedtools loj matches all gnomAD variants to all user variant alleles for a given position
  # results in expansion
  # remove non-matching entries
  ints <- ints[ref_user == ref_strand & alt_user == alt_strand]
  
  # tidy
  ints[, c("ref_user", "alt_user", "user_alleles") := NULL]
  
} else if (user_in[2] == "plus") {
  
  # extract alleles
  ints[, c("ref_user", "alt_user") := tstrsplit(user_alleles, "[|]")]
  ints[, ref_user := toupper(ref_user)]
  ints[, alt_user := toupper(alt_user)]
  
  # update user alleles to match 3'UTR strand
  ints[strand == "-", ref_user := chartr("ACGT", "TGCA", ref_user)]
  ints[strand == "-", alt_user := chartr("ACGT", "TGCA", alt_user)]
  
  # bedtools loj matches all gnomAD variants to all user variant alleles for a given position
  # results in expansion
  # remove non-matching entries
  ints <- ints[ref_user == ref_strand & alt_user == alt_strand]
  
  # tidy
  ints[, c("ref_user", "alt_user", "user_alleles") := NULL]
  
} else if (user_in[2] == "none") {
  
} else {
  stop("invalid input allele version specified", call. = FALSE)
}

#### ============= END OF intersect ==================== ####


#### ==================== calibrate ==================== ####

print("calibrating variants for iMAPS calculation")

# add base change
ints[, base_change := paste(use_as_ancest, use_as_deriv, sep = ">")]
ints[use_as_ancest %in% c("A", "C"), base_change_pair := base_change]
ints[use_as_ancest %in% c("G", "T"), base_change_pair := chartr("ACGT", "TGCA", base_change)] # not "reversed" since doing both bases together

## DInucleotide context (including methylation level for CpG transitions)

# confirm k_mer size
k_curr <- 3
# extract dinucleotide context
ints[, dinuc := substr(ancest_k_seq, ((nchar(ancest_k_seq) - k_curr) / 2) + 1, nchar(ancest_k_seq) - ((nchar(ancest_k_seq) - k_curr) / 2))]

## add category and super category 

# add context category (CpG transition, non-CpG transition, transversion)
ints[dinuc %in% c("ACG", "CCG", "GCG", "TCG", "CGT", "CGG", "CGC", "CGA") & base_change %in% c("C>T", "G>A"), cat := "CpG-ti"]
ints[!dinuc %in% c("ACG", "CCG", "GCG", "TCG", "CGT", "CGG", "CGC", "CGA") & base_change %in% c("C>T", "G>A", "T>C", "A>G"), cat := "other-ti"]
ints[base_change %in% c("A>C", "A>T", "C>A", "C>G", "G>C", "G>T", "T>A", "T>G"), cat := "tv"]
#table(ints$cat)

# add methylation bin (for CpG transitions only)
# from Karczewski et al. 2020: low methylation, missing or < 0.2; medium, 0.2-0.6; and high, > 0.6
# here, don't use data with missing methylation status (filtered later)
ints[cat == "CpG-ti", methyl_mean_cat := ifelse(is.na(methyl_mean) == TRUE, "N", ifelse(methyl_mean > 60, "H", ifelse(methyl_mean > 20, "M", "L")))] # N for none / missing data
ints[cat != "CpG-ti", methyl_mean_cat := "N"] # assign N for "none" to all non CpG ti variants (not considering impact of any potential methylation for this category)
#table(ints[ , methyl_mean_cat])

## TRInucleotide context

# confirm k_mer size
k_curr <- 5
# extract trinucleotide context
ints[, trinuc := substr(ancest_k_seq, ((nchar(ancest_k_seq) - k_curr) / 2) + 1, nchar(ancest_k_seq) - ((nchar(ancest_k_seq) - k_curr) / 2))]

## TETRAnucleotide context

# confirm k_mer size
k_curr <- 7
# extract quadranucleotide context
ints[, tetranuc := substr(ancest_k_seq, ((nchar(ancest_k_seq) - k_curr) / 2) + 1, nchar(ancest_k_seq) - ((nchar(ancest_k_seq) - k_curr) / 2))]

## PENTAnucleotide context

# confirm k_mer size used for k_seq
k_curr <- 9
# extract pentanucleotide context
ints[, pentanuc := substr(ancest_k_seq, ((nchar(ancest_k_seq) - k_curr) / 2) + 1, nchar(ancest_k_seq) - ((nchar(ancest_k_seq) - k_curr) / 2))]

## HEXAnucleotide context 

# confirm k_mer size used for k_seq
k_curr <- 11
# extract hexanucleotide context
ints[, hexanuc := substr(ancest_k_seq, ((nchar(ancest_k_seq) - k_curr) / 2) + 1, nchar(ancest_k_seq) - ((nchar(ancest_k_seq) - k_curr) / 2))]
ints[, hexanuc_rev_comp := stri_reverse(chartr("ACGT", "TGCA", hexanuc))]
# make new column of showing context and base change (for strand with ancestral base as A or C)
ints[use_as_ancest %in% c("A", "C"), hexa_con_pair := paste(substr(hexanuc, 1, (k_curr -1) / 2), "[", base_change_pair, "]", substr(hexanuc, (((k_curr -1) / 2) + 2), nchar(hexanuc)), sep = "")]
ints[use_as_ancest %in% c("G", "T"), hexa_con_pair := paste(substr(hexanuc_rev_comp, 1, (k_curr -1) / 2), "[", base_change_pair, "]", substr(hexanuc_rev_comp, (((k_curr -1) / 2) + 2), nchar(hexanuc_rev_comp)), sep = "")]
# update context to include methylation level
ints[, hexa_con_pair := paste(hexa_con_pair, methyl_mean_cat, sep = "-")]

# tidy
ints[, c("hexanuc_rev_comp", "trinuc", "tetranuc", "pentanuc", "hexanuc") := NULL]

## make id for matching to calibration variants

# convert DNase score to binary
ints[, DNase_binary := ifelse(DNase > 0, 1, 0)]
ints[, hexa_id_pair := paste(hexa_con_pair, DNase_binary, H3K9me3, sep = "_")]
ints[, hexa_con_pair := NULL]

setwd(user_in[3])
calib <- fread("intergenic_calibration_data.txt.gz", header = TRUE)

# add expected proportion singleton based on calibration data to variants
setkey(ints, hexa_id_pair)
setkey(calib, hexa_id_pair)
ints <- calib[ints]

# tidy
rm(calib)

# how many variants are calibrated at each level of extended nucleotide context?
#table(ints[, highest_sig_k])
#round(table(ints[, highest_sig_k]) / nrow(ints) * 100, 1)

#### ============= END OF calibrate ==================== ####


#### ==================== filter and calculate iMAPS ==================== ####

print("calculating iMAPS")

# mark singletons
ints[, singleton := ifelse(use_as_deriv_count == 1, 1, 0)]

# add complex ancestry (cases where ancestral allele is neither reference or alternative allele)
ints[, complex_ancest := ifelse(ancest_seq_strand == ref_strand | ancest_seq_strand == alt_strand, 0, 1)]

# set up filters
ints[, prim_filter_pass := ifelse(ancest_seq_strand == ref_strand & (!chrom %in% c("chrX", "chrY", "chrM")) & gnomAD_v3_cov_med > 24 & gnomAD_v3_cov_med < 43 & AN > 129064 & ((cat == "CpG-ti" & methyl_mean_cat != "N") | cat != "CpG-ti") & complex_ancest == 0 & lcr == 0 & CpG_island == 0 & CDS_overlap == 0 & base_change_pair != "A>G", 1, 0)]

# calculate iMAPS
OE <- ints[prim_filter_pass == 1,
               .(observed_singletons = sum(singleton),
                 expected_singletons = sum(highest_sig_k_prop_sing),
                 n_variants = .N),
           by = .(category = user_category)]
OE[, imaps := (observed_singletons - expected_singletons) / n_variants]

#### ============= END OF filter and calculate iMAPS ==================== ####


#### ==================== statistical testing ==================== ####

# ** manually change for desired comparisons **
# un-comment to run for demo variants

#  # calculate non-singleton or "other" variant counts
#  OE[, observed_others := n_variants - observed_singletons]
#  OE[, expected_others := n_variants - expected_singletons]
#  
#  # calculate singleton to non-singleton ratios
#  OE[, exp_ratio := expected_singletons / expected_others]
#  
#  # set desired "category" values to compare 
#  cat_x <- "pres"
#  cat_y <- "lost"
#  
#  # calc expected OR for null testing
#  OE_exp <- OE[category == cat_y, exp_ratio] / OE[category == cat_x, exp_ratio]
#  
#  # fisher exact test against expected OR
#  alt_hyp <- "greater" # change as desired
#  fisher_res <- fisher.test(x = matrix(c(OE[category == cat_x, observed_others], OE[category == cat_y, observed_others],
#                                         OE[category == cat_x, observed_singletons], OE[category == cat_y, observed_singletons]),
#                                       nrow = 2, ncol = 2),
#                            or = OE_exp, alternative = alt_hyp)
#  fisher_res$p.value
#  fisher_res$estimate
#  
#  # tidy
#  OE[, c("observed_others", "expected_others", "exp_ratio") := NULL]

#### ============= END OF statistical testing ==================== ####


#### ==================== plot and save ==================== ####

print("generating plot and saving data")

# add averages for 3' UTRs, synonymous, missense, and stop-gain variants for reference
OE <- rbind(OE, data.table(category = c("synonymous", "missense", "stop-gain", "all 3' UTR"),
                           observed_singletons = c(514995, 1215034, 54509, 3663120),
                           expected_singletons = c(485121.1, 1051709.1, 40162.9, 3533800),
                           n_variants = c(1122454, 2344700, 90188, 7167824)), fill = TRUE)
OE[, imaps := (observed_singletons - expected_singletons) / n_variants]

# set plot order
OE[, category := factor(as.character(category), levels = c("synonymous", "missense", "stop-gain", "all 3' UTR", sort(unique(ints[, user_category]))))]

# select plot colours
plot_cols <- c(brewer.pal(n = 8, name = "Dark2")[c(1,2,4)], brewer.pal(n = 9, name = "Blues")[c(8, rep(5, length(unique(ints[, user_category]))))])

# plot
ggplot(data = OE, aes(x = category, y = imaps, fill = category, colour = category, alpha = category)) +
  geom_bar(position = "dodge", stat = "identity", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.75, colour = "black") + 
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.2, 0.2), legend.key.height = unit(1.5,"line"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = plot_cols, guide = "none", na.value = NA) +
  scale_colour_manual(values = plot_cols, guide = "none") +
  scale_alpha_manual(values = c(rep(0.025, 3), rep(1, 4)), guide = "none") +
  scale_y_continuous(breaks = c(seq(0, 0.16, by = 0.04))) +
  coord_cartesian(ylim = c(0,0.16)) +
  labs(x = "\n", y = "iMAPS\n")

# save
fwrite(OE, paste(user_in[3], "imaps.by_category.DATA.txt", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
ggsave("imaps.by_category.PLOT.pdf", path = user_in[3], width = 5, height = 5, units = "in")

#### ============= END OF plot and save ==================== ####




