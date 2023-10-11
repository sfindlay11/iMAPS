#################### notes ####################

# uses IDR peaks resulting from Yeo Lab's merge_peaks function
# sequence-based ("old") or sequence+structure ("new") motif models used for RBPamp

#################### END of notes ####################


#################### set up ####################
library(data.table)
library(pbapply)

path_data_in <- ""
path_data_out <- ""
path_data_fasta <- ""
path_data_rbp_amp_old_motifs <- ""
path_data_rbp_amp_new_motifs <- ""

names_BED_std <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

options(scipen = 999)

# set desired extension window (plus max(k)-1 to avoid decay at edges)
aff_ext_len <- 100 + (11 -1) # for max aff location plot (set for k = 11 here and RBPs with k < 11 are adjusted later)

#################### END of set up ####################


#################### process peaks ####################

# read in peak 5' ends
setwd(path_data_in)
peaks <- fread("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_name.gene_region_score.bed", header = FALSE)
names(peaks) <- names_BED_std
peaks

# add bases to each side of peak 5' end (to allow for different windows relative to peak)
peaks[, chromStart_ext := chromStart - aff_ext_len]
peaks[, chromEnd_ext := chromEnd + aff_ext_len]
summary(peaks[, chromEnd_ext - chromStart_ext])

# write to file
setwd(path_data_out)
setkey(peaks, chrom, chromStart_ext, name, strand)
fwrite(peaks[, .(chrom, chromStart_ext, chromEnd_ext, name, score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_name.gene_region_score.pm_ext.for_aff.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
rm(peaks)

# fetch reference seq for each extended peak
# copy eCLIP data to fasta dir
setwd(path_data_out)
system("head eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_name.gene_region_score.pm_ext.for_aff.bed")
system(paste("cp eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_name.gene_region_score.pm_ext.for_aff.bed", path_data_fasta, sep = " "))
setwd(path_data_fasta)
system("bedtools getfasta -bedOut -s -fi hg38.fa -bed eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_name.gene_region_score.pm_ext.for_aff.bed > eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_name.gene_region_score.pm_ext.for_aff.PLUS_ref_seq.txt")
system("rm eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_name.gene_region_score.pm_ext.for_aff.bed") # removes copy
# move back result
system(paste("mv eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_name.gene_region_score.pm_ext.for_aff.PLUS_ref_seq.txt", path_data_in))

setwd(path_data_in)
peaks <- fread("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_name.gene_region_score.pm_ext.for_aff.PLUS_ref_seq.txt", header = FALSE)
names(peaks) <- c("chrom", "chromStart_ext", "chromEnd_ext", "name", "score", "strand", "ref_seq_ext")
peaks

# convert all seq to upper case
peaks[, ref_seq_ext := toupper(ref_seq_ext)]

# check lengths
table(peaks[, nchar(ref_seq_ext)])

# extract name info
peaks[, c("RBP_peak", "cell_line_peak", "chrom_peak", "chromStart_peak", "chromEnd_peak", "strand_peak", "enrichment_peak", "neg_log_10_p_value_peak", "chromStart_five_end", "chromEnd_five_end") := tstrsplit(name, "_", fixed = TRUE)]

# extract score info
peaks[, c("gene_id_prot_cod_string", "gene_count", "prot_cod_prop", "region_count", "region_top") := tstrsplit(score, "_", fixed = TRUE)]

# explore
sort(table(peaks[region_top == "utr3" & cell_line_peak == "HepG2", RBP_peak]), decreasing = TRUE)
sort(table(peaks[region_top == "utr3" & cell_line_peak == "K562", RBP_peak]), decreasing = TRUE)
temp <- peaks[, .(cell_count = length(unique(cell_line_peak))), by = RBP_peak]
duals <- temp[cell_count == 2, RBP_peak]
sort(table(peaks[region_top == "utr3" & RBP_peak %in% duals, RBP_peak]), decreasing = TRUE)
#summary(peaks[RBP_peak == "RBFOX2", as.numeric(enrichment_peak)])
#quantile(peaks[RBP_peak == "RBFOX2" & grepl("GCATG|CGACG", ref_seq_ext) == T, as.numeric(enrichment_peak)], seq(0,1, by = 0.1))
#quantile(peaks[RBP_peak == "RBFOX2" & grepl("GCATG|GCACG", ref_seq_ext) == F, as.numeric(enrichment_peak)], seq(0,1, by = 0.1))

# flag RBPs for which we have RBPamp models
setwd(path_data_rbp_amp_old_motifs)
old_motifs <- system("ls *.txt | cut -f 1 -d '.'", intern = TRUE)
setwd(path_data_rbp_amp_new_motifs)
new_motifs <- system("ls *.txt | cut -f 1 -d '.'", intern = TRUE)
peaks[, RBPamp_seq_only := ifelse(RBP_peak %in% old_motifs, T, F)]
table(peaks$RBPamp_seq_only) / nrow(peaks) # have RBPamp old motif for 17% of peaks 
peaks[, RBPamp_seq_and_struct := ifelse(RBP_peak %in% new_motifs, T, F)]
table(peaks$RBPamp_seq_and_struct) / nrow(peaks) # have RBPamp new motif for 18% of peaks 
peaks

# explore
# number of peaks per RBP and whether or not we have RBPamp for it
#summ <- as.data.table((sort(table(peaks[cell_line_peak == "HepG2" & region_top == "utr3", RBP_peak]), decreasing = T)))
summ <- as.data.table((sort(table(peaks[cell_line_peak == "K562" & region_top == "utr3", RBP_peak]), decreasing = T)))
names(summ) <- c("RBP_peak", "count_peak")
#summ[, cell_line_peak := "HepG2"]
summ[, cell_line_peak := "K562"]
summ[, RBPamp_seq_only := ifelse(RBP_peak %in% old_motifs, T, F)]
summ[, RBPamp_seq_and_struct := ifelse(RBP_peak %in% new_motifs, T, F)]
# proportion of RBPs (with more than 10 peaks) with RBPamp data
length(unique(summ[count_peak >= 10 & (RBPamp_seq_only == T | RBPamp_seq_and_struct == T), RBP_peak])) / length(unique(summ[count_peak >= 10, RBP_peak])) # 18% / 19% of RBPs
length(unique(summ[count_peak >= 1000 & (RBPamp_seq_only == T | RBPamp_seq_and_struct == T), RBP_peak])) / length(unique(summ[count_peak >= 1000, RBP_peak])) # 31% / 42% of "top" RBPs
# proportion of 3'UTR peaks with RBPamp data
sum(summ[RBPamp_seq_only == T | RBPamp_seq_and_struct == T, count_peak]) / sum(summ[, count_peak]) # 21% / 30% of peaks 

# write them to files (make sure key is set for sorting)
setwd(path_data_out)
setkey(peaks, RBP_peak, cell_line_peak, chrom, chromStart_ext, chromEnd_ext)
fwrite(peaks[RBPamp_seq_only == 1, .(name = paste(name, score, sep = "&"), RBP_peak, seq = ref_seq_ext)], "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.name_score.RBP.ref_seq_ext.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
fwrite(peaks[RBPamp_seq_and_struct == 1, .(name = paste(name, score, sep = "&"), RBP_peak, seq = ref_seq_ext)], "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.name_score.RBP.ref_seq_ext.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# read back
#peaks <- fread("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.name_score.RBP.ref_seq_ext.txt")
peaks[, c("name_short", "score") := tstrsplit(name, "&", fixed = TRUE)]
peaks[, c("RBP_peak", "cell_line_peak", "chrom_peak", "chromStart_peak", "chromEnd_peak", "strand_peak", "enrichment_peak", "neg_log_10_p_value_peak", "chromStart_five_end", "chromEnd_five_end") := tstrsplit(name_short, "_", fixed = TRUE)]
peaks[, c("gene_id_prot_cod_string", "gene_count", "prot_cod_prop", "region_count", "region_top") := tstrsplit(score, "_", fixed = TRUE)]
peaks[region_top == "utr3", ]
sort(table(peaks[region_top == "utr3", RBP_peak]))
length(unique(peaks[region_top == "utr3", RBP_peak]))
length(unique(peaks[region_top == "utr3" & RBPamp_seq_only == 1, RBP_peak]))
length(unique(peaks[region_top == "utr3" & RBPamp_seq_and_struct == 1, RBP_peak]))
View(as.data.table(table(peaks[region_top == "utr3", RBP_peak])))
#################### END of process peaks ####################


#################### fetch RBPamp affinities ####################

# ensured motif PSAMs were (manually) SORTED (descending A0) so as to avoid scaling problems by program

# rbp_amp_eclip_fetch_general.py
# python rbp_amp_eclip_fetch_general.py eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.name_score.RBP.ref_seq_ext.txt /net/eofe-data010/data001/burgelab/nevermind/data/nm/sfindlay/rbpamp/old_motifs_from_HJ_PRIOR_Nov_25_2020_SORTED/
# mv rbp_amp_eclip_fetch_general_OUT.txt eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.name_score.RBP.ref_seq_ext.k_size.affinity_by_pos.txt
# python rbp_amp_eclip_fetch_general.py eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.name_score.RBP.ref_seq_ext.txt /net/eofe-data010/data001/burgelab/nevermind/data/nm/sfindlay/rbpamp/new_motifs_from_HJ_AFTER_Nov_25_2020_SORTED/
# mv rbp_amp_eclip_fetch_general_OUT.txt eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.name_score.RBP.ref_seq_ext.k_size.affinity_by_pos.txt

#################### END of fetch RBPamp affinities ####################


#################### process RBPamp output ####################

# **(REPEAT THIS SECTION FOR "OLD" AND "NEW" MOTIFS SEPARATELY)** -----------------------------------------------------------

# read in
setwd(path_data_in)
#amp_peaks <- fread("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.name_score.RBP.ref_seq_ext.k_size.affinity_by_pos.txt", header = TRUE)
amp_peaks <- fread("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.name_score.RBP.ref_seq_ext.k_size.affinity_by_pos.txt", header = TRUE)
names(amp_peaks) <- c("name", "RBP_peak", "ref_seq_ext", "aff_string", "k_size")
amp_peaks

# convert affinity_array string into list
amp_peaks[, affinity_list := pblapply(aff_string, function(x) as.numeric(strsplit(x, split = ",", fixed = TRUE)[[1]]))]
amp_peaks[, aff_string := NULL]
amp_peaks

# length check
table(amp_peaks[, nchar(ref_seq_ext)]) # aff_ext_len * 2 + 1
amp_peaks[, len_check := sapply(affinity_list, length)]
table(amp_peaks[, len_check]) # should be L - (k -1)
amp_peaks[, len_check := NULL]

# expand name
amp_peaks[, c("name_bed", "score_bed") := tstrsplit(name, "&", fixed = TRUE)]
amp_peaks[, c("RBP_peak_extra", "cell_line_peak", "chrom_peak", "chromStart_peak", "chromEnd_peak", "strand_peak", "enrichment_peak", "neg_log_10_p_value_peak", "chromStart_ref_point", "chromEnd_ref_point") := tstrsplit(name_bed, "_", fixed = TRUE)]
amp_peaks[, c("gene_id_prot_cod_string", "gene_count", "prot_cod_prop", "region_count", "region_top") := tstrsplit(score_bed, "_", fixed = TRUE)]
amp_peaks[, chromStart_peak := as.integer(chromStart_peak)]
amp_peaks[, chromEnd_peak := as.integer(chromEnd_peak)]
amp_peaks[, chromStart_ref_point := as.integer(chromStart_ref_point)]
amp_peaks[, chromEnd_ref_point := as.integer(chromEnd_ref_point)]
amp_peaks[, gene_count := as.integer(gene_count)]
amp_peaks[, prot_cod_prop := as.numeric(prot_cod_prop)]
amp_peaks[, region_count := as.integer(region_count)]
amp_peaks[, c("RBP_peak_extra", "name") := NULL]

# update start and end to represent extension
amp_peaks[, chromStart_ext := chromStart_ref_point - aff_ext_len]
amp_peaks[, chromEnd_ext := chromEnd_ref_point + aff_ext_len]
table(amp_peaks[, chromEnd_ext - chromStart_ext])

# set new ext for each entry (based on k)
amp_peaks[k_size ==  10, ind_aff_ext_len := (aff_ext_len * 2) -2] 
amp_peaks[k_size ==  11, ind_aff_ext_len := aff_ext_len * 2] 
table(amp_peaks[, ind_aff_ext_len])
# trim k = 10 entries (first and last)
amp_peaks[k_size == 10, ref_seq_ext := substr(ref_seq_ext, 2, nchar(ref_seq_ext) -1)]
amp_peaks[k_size == 10, affinity_list := pblapply(affinity_list, function(x) x[2:(length(x) - 1)])]
# make sure to update ext coordinates as well
amp_peaks[k_size == 10, chromStart_ext := chromStart_ext + 1]
amp_peaks[k_size == 10, chromEnd_ext := chromEnd_ext - 1]

# write to file for separate analysis..
# manually change aff_ext_len in file name
setwd(path_data_out)
#saveRDS(amp_peaks, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_OLD_motifs.pm_110.full_info.expanded.rds")
#saveRDS(amp_peaks, "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.RBPs_with_amp_NEW_motifs.pm_110.full_info.expanded.rds")


#################### END of process RBPamp output ####################

