
# set up ------------------------------------------------------------------

library(data.table)
library(stringr)
library(tidyverse)

names_BED_std <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")
names_chr_main <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")

path_data_in <- ""
path_data_out <- ""
path_data_ensembl <- ""

options(scipen = 999)

# -----------------------------------------------------------------------


# read in polyA data base v3 site info ------------------------------------

setwd(path_data_in)
polys <- fread("polyA_DB_v3_hg19.txt", header = TRUE)
nrow(polys) # 311,594
names(polys)
names(polys) <- c("PAS_id", "chrom", "pos", "strand_site", "mean_RPM", "Intron_exon_location", "Ensemble_ID", "RefSeq_Gene_ID", "Gene_Symbol", "Gene_Name", "FAMTOM_ID", "FANTOM_Category", "Extension", "PAS_type", "PSE", "PAS_Signal", "Conservation", "intergenic_TE")
sum(polys[, Gene_Name == ""])
sum(polys[, RefSeq_Gene_ID == ""])
sum(polys[, Ensemble_ID == ""])
sum(polys[, Intron_exon_location == "intergenic" & Ensemble_ID == ""])
polys[, c("Gene_Name") := NULL]
polys[Ensemble_ID == "", Ensemble_ID := NA]
polys[Ensemble_ID == "na", Ensemble_ID := NA]
polys[RefSeq_Gene_ID == "", RefSeq_Gene_ID := NA]
polys[RefSeq_Gene_ID == "na", RefSeq_Gene_ID := NA]
polys[Gene_Symbol == "", Gene_Symbol := NA]
polys[Gene_Symbol == "na", Gene_Symbol := NA]
polys[, fantom_id := FAMTOM_ID] # correct typo
polys[fantom_id == "na", fantom_id := NA]
polys[, FAMTOM_ID := NULL]
polys[, fantom_category := FANTOM_Category] 
polys[fantom_category == "na", fantom_category := NA] 
polys[, FANTOM_Category := NULL]
polys[RefSeq_Gene_ID == "na", RefSeq_Gene_ID := NA] 

polys[, conserved := ifelse(Conservation == "Yes", 1, ifelse(Conservation == "No", 0, NA))]
polys[, Conservation := NULL]
polys[, extended := ifelse(Extension == "YES", 1, ifelse(Extension == "NO", 0, NA))]
polys[, Extension := NULL]
polys[, mean_RPM := round(mean_RPM, 2)]
polys[, PSE := round(as.numeric(substr(PSE, 1, str_locate(PSE, "%") -1)) / 100, 3)]
polys[, chromStart_site := pos - 1]
polys[, chromEnd_site := pos]
polys[, name_site := paste(chrom, chromStart_site, chromEnd_site, strand_site, sep = "_")]
polys[, score_site := "."]
polys

# -----------------------------------------------------------------------


# convert polyA site data to hg38 ------------------------------------

# write to bed file for liftover to hg38
polys[, chromStart_site := pos - 1]
polys[, chromEnd_site := pos]
polys[, name_site := paste(chrom, chromStart_site, chromEnd_site, strand_site, sep = "_")]
polys[, score_site := "."]
setwd(path_data_in)
fwrite(polys[, .(chrom, chromStart_site, chromEnd_site, name_site, score_site, strand_site)], "polyA_DB_v3_hg19_for_liftover.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# liftover at http://genome.ucsc.edu/cgi-bin/hgLiftOver
# successfully converted 311500 records
# failed on 94 records
# read in results
setwd(path_data_in)
lift_over_results <- fread("polyA_DB_v3_liftover_hg38.bed", header = FALSE)
names(lift_over_results) <- c("chrom_site", "chromStart_site", "chromEnd_site", "name_site", "score_site", "strand_site")
lift_over_results

# any duplicates?
sum(duplicated(lift_over_results)) # yes
#remove
lift_over_results <- lift_over_results[duplicated(lift_over_results) == FALSE, ]

# filter non main chromosomes
lift_over_results <- lift_over_results[chrom_site %in% names_chr_main, ]

# update to hg38 coordinates by name_site
polys[, c("chrom", "chromStart_site", "chromEnd_site", "score_site", "strand_site") := NULL]
setkey(polys, name_site)
setkey(lift_over_results, name_site)
polys <- lift_over_results[polys, nomatch = 0] # inner join
rm(lift_over_results)
polys

#! IMPORTANT
# update previous name_site to hg38 (used to be hg19 coordinates)
polys[, name_site := NULL]
polys[, name_site := paste(chrom_site, chromStart_site, chromEnd_site, strand_site, sep = "_")]

# -----------------------------------------------------------------------

# write unique pA sites that have been matched to protein coding genes in polys to bed file for bedtools closest
table(polys[, Intron_exon_location], polys[, PAS_type])
prot_cod_polys_bed <- polys[Intron_exon_location != "intergenic" & PAS_type != "LncRNA(FANTOM5)" & PAS_type != "Pseudogene", .(name = paste(Ensemble_ID, collapse = "_"), score = paste(Gene_Symbol, collapse = "_")), by = .(chrom = chrom_site, chromStart_site, chromEnd_site, strand_site)]
setcolorder(prot_cod_polys_bed, c("chrom", "chromStart_site", "chromEnd_site", "name", "score", "strand_site"))
prot_cod_polys_bed
setwd(path_data_out)
fwrite(prot_cod_polys_bed, "unique_prot_cod_polyA_DB_v3_hg38_liftover_sites.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# link each polyA site to its upstream exon / stop codon, so records can be converted to 3' UTRs

# process Ensembl data ------------------------------------

# read in Ensembl data (Wang et al 2018 used release 75) ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
setwd(path_data_ensembl)
#ens <- fread("Homo_sapiens.GRCh37.75.gtf", header = FALSE, skip = 5)
ens <- fread("Homo_sapiens.GRCh38.101.gtf.gz", header = FALSE, skip = 5) # current as of Oct 2020
names(ens) <- c("chrom", "gene_source", "component", "chromStart", "chromEnd", "score", "strand", "frame", "info")
ens

#! ens .gtf file is 1-based. convert to 0-based.
ens$chromStart <- ens$chromStart - 1

# update chrom names
ens[, chrom := paste("chr", chrom, sep = "")]
ens[chrom == "chrMT", chrom := "chrM"]

# filter non main chromosomes
ens <- ens[chrom %in% names_chr_main, ]

# extract gene id
info_tag_string <- "gene_id"
ens$gene_id <- str_match(ens$info, paste(info_tag_string, ' ".*?"?(;|$)', sep = "")) [ ,1] # ? stops at first instance () is for ; or end of line (so last info field can be captured) # only need first column, second column is capture group
head(ens)
# inspect for uniformity
table(substr(ens$gene_id, 1, nchar(info_tag_string) + 2))
table(substr(ens$gene_id, nchar(ens$gene_id) -1, nchar(ens$gene_id)))
ens$gene_id <- substr(ens$gene_id, nchar(info_tag_string) + 3, nchar(ens$gene_id) -2)
head(ens)

# extract transcript_id
info_tag_string <- "transcript_id"
ens$transcript_id <- str_match(ens$info, paste(info_tag_string, ' ".*?"?(;|$)', sep = "")) [ ,1] # ? stops at first instance () is for ; or end of line (so last info field can be captured) # only need first column, second column is capture group
head(ens)
# inspect for uniformity
table(substr(ens$transcript_id, 1, nchar(info_tag_string) + 2))
table(substr(ens$transcript_id, nchar(ens$transcript_id) -1, nchar(ens$transcript_id)))
ens$transcript_id <- substr(ens$transcript_id, nchar(info_tag_string) + 3, nchar(ens$transcript_id) -2)
head(ens)

# extract gene biotype
info_tag_string <- "gene_biotype"
ens$gene_biotype <- str_match(ens$info, paste(info_tag_string, ' ".*?"?(;|$)', sep = "")) [ ,1] # ? stops at first instance () is for ; or end of line (so last info field can be captured) # only need first column, second column is capture group
head(ens)
# inspect for uniformity
table(substr(ens$gene_biotype, 1, nchar(info_tag_string) + 2))
table(substr(ens$gene_biotype, nchar(ens$gene_biotype) -1, nchar(ens$gene_biotype)))
ens$gene_biotype <- substr(ens$gene_biotype, nchar(info_tag_string) + 3, nchar(ens$gene_biotype) -2)
head(ens)

# extract gene name
info_tag_string <- "gene_name"
ens$gene_name <- str_match(ens$info, paste(info_tag_string, ' ".*?"?(;|$)', sep = "")) [ ,1] # ? stops at first instance () is for ; or end of line (so last info field can be captured) # only need first column, second column is capture group
head(ens)
# inspect for uniformity
table(substr(ens$gene_name, 1, nchar(info_tag_string) + 2))
table(substr(ens$gene_name, nchar(ens$gene_name) -1, nchar(ens$gene_name)))
ens$gene_name <- substr(ens$gene_name, nchar(info_tag_string) + 3, nchar(ens$gene_name) -2)
head(ens)

# extract transcript_biotype
info_tag_string <- "transcript_biotype"
ens$transcript_biotype <- str_match(ens$info, paste(info_tag_string, ' ".*?"?(;|$)', sep = "")) [ ,1] # ? stops at first instance () is for ; or end of line (so last info field can be captured) # only need first column, second column is capture group
head(ens)
# inspect for uniformity
table(substr(ens$transcript_biotype, 1, nchar(info_tag_string) + 2))
table(substr(ens$transcript_biotype, nchar(ens$transcript_biotype) -1, nchar(ens$transcript_biotype)))
ens$transcript_biotype <- substr(ens$transcript_biotype, nchar(info_tag_string) + 3, nchar(ens$transcript_biotype) -2)
table(ens$transcript_biotype)
head(ens)

table(ens$gene_biotype, ens$transcript_biotype)

# update
ens[, c("info") := NULL]
ens

# -----------------------------------------------------------------------


# Extract and process ensembl subset data -------------------------------------

# note: 101 is different from 75 in that it uses gene biotype and transcript biotype, but not gene_cat

# make coding stop codon DT
stops <- ens[component == "stop_codon" & gene_biotype == "protein_coding", ] 
table(stops$transcript_biotype)
names(stops) <- c("chrom_stop", "source_stop", "component", "chromStart_stop", "chromEnd_stop", "score_stop", "strand_stop", "frame_stop", "ens_gene_id_stop", "ens_trans_id_stop", "gene_biotype_stop", "gene_name_stop", "transcript_biotype_stop")
stops[, c("component", "source_stop", "score_stop", "gene_biotype_stop") := NULL]
stops
# sometimes there is more than one stop codon entry per transcript. these are spliced / split stop codons
stops[, transcript_count := .N, by = ens_trans_id_stop]
# keep the most distal portion as the stop codon
stops[transcript_count == 1, keep := 1]
stops[transcript_count > 1 & strand_stop == "+", keep := ifelse(chromEnd_stop == max(chromEnd_stop), 1, 0), by = ens_trans_id_stop]
stops[transcript_count > 1 & strand_stop == "-", keep := ifelse(chromStart_stop == min(chromStart_stop), 1, 0), by = ens_trans_id_stop]
stops[transcript_count > 1, ]
stops <- stops[keep == 1, ]
stops # 80,525 of 80,973 kept
stops[, c("keep", "transcript_count") := NULL]
stops
# how many unique stops per gene?
temp <- stops[, .(count = length(unique(paste(chrom_stop, chromStart_stop, chromEnd_stop, strand_stop, sep = "_")))), by = ens_gene_id_stop]
table(temp$count)
table(temp$count) / sum(table(temp$count))
# 1 = 46%, 2 = 23%, 3 = 14%, 4 = 7%, 5+ = 10% (for ens 101)
# 1 = 55%, 2 = 24%, 3 = 11%, 4 = 5%, 5+ =  5% (for ens 75)
temp[count > 1, ]

# make exon DT (very important to use both gene_cat and gene_biotype appropriately. take "CDS" for anything, and "exons" for protein coding genes (includes non-coding transcripts within them))
exons <- ens[component == "CDS" | (gene_biotype %in% unique(ens[component == "CDS", gene_biotype]) & transcript_biotype != "protein_coding" & component == "exon"), ]
names(exons) <- c("chrom_exon", "source_exon", "component", "chromStart_exon", "chromEnd_exon", "score_exon", "strand_exon", "frame_exon", "ens_gene_id_exon", "ens_trans_id_exon", "gene_biotype_exon", "gene_name_exon", "transcript_biotype_exon")
exons

# add stop codon annotations to each cds exon entry by transcript id
length(unique(stops[, ens_trans_id_stop])) # 80,525 (in ens 75: 66,304) transcripts
length(unique(exons[, ens_trans_id_exon])) # 155,148 (146,110) transcripts
stops_sub <- stops[, .(chromStart_stop, chromEnd_stop, ens_trans_id_stop, transcript_biotype_stop)]
setkey(stops_sub, ens_trans_id_stop)
setkey(exons, ens_trans_id_exon)
# merge is NOT limited by stop codon data. Limiting by stop codon data caused problems. need to consider non-coding transcripts in protein coding genes and filter these later
exons_stops <- stops_sub[exons] 
exons_stops
rm(stops_sub)
# do categories always agree between stop and exon? (except where there is no matching stop codon)
sum(exons_stops[, transcript_biotype_stop == transcript_biotype_exon], na.rm = TRUE) == nrow(exons_stops[is.na(transcript_biotype_stop) == FALSE, ])
# yes, so use cat_exon only
exons_stops[, transcript_biotype_stop := NULL]
# what proportion of unique exons does this include? (only used if limiting by stop codon)
#exons[, temp_id := paste(chrom_exon, chromStart_exon, chromEnd_exon, strand_exon, sep = "_")]
#exons_stops[, temp_id := paste(chrom_exon, chromStart_exon, chromEnd_exon, strand_exon, sep = "_")]
#length(unique(exons_stops[, temp_id])) / length(unique(exons[, temp_id])) # 84%
#exons_stops[, temp_id := NULL]
#exons[, temp_id := NULL]

# label last exons (i.e. adjacent to or overlapping stop codon; overlapping only includes another few hundred cases)
exons_stops[component == "CDS" & strand_exon == "+", last_cds_exon := ifelse(chromEnd_exon >= chromStart_stop, 1, 0)]
exons_stops[component == "CDS" & strand_exon == "-", last_cds_exon := ifelse(chromStart_exon <= chromEnd_stop, 1, 0)]
table(exons_stops[, last_cds_exon]) # 79,735 (66,084 for 75)
sum(exons_stops[, is.na(last_cds_exon)])
table(exons_stops[, last_cds_exon], exons_stops[, frame_exon], useNA = "ifany")
# update last exons to end at end of stop codon (so UTR can be defined starting after stop)
exons_stops[last_cds_exon == 1 & strand_exon == "+", chromEnd_exon := chromEnd_stop]
exons_stops[last_cds_exon == 1 & strand_exon == "-", chromStart_exon := chromStart_stop]

# also need to truncate exons from non coding transcripts that overlap with stop codon
# for this, compare BY GENE. 
exons_stops[, gene_chromStart_stops := list(list(unique(chromStart_stop[is.na(chromStart_stop) == F]))), by = ens_gene_id_exon]
exons_stops[, gene_chromEnd_stops := list(list(unique(chromEnd_stop[is.na(chromEnd_stop) == F]))), by = ens_gene_id_exon]
# mark entries
exons_stops[component == "exon", nc_exon_stop_overlap := mapply(function(exStart,exEnd,stpStart,stpEnd) ifelse(sum(exStart <= unlist(stpStart) & exEnd >= unlist(stpEnd)) > 0, 1, 0), exStart = chromStart_exon, exEnd = chromEnd_exon, stpStart = gene_chromStart_stops, stpEnd = gene_chromEnd_stops)]
table(exons_stops[, nc_exon_stop_overlap])
exons_stops[nc_exon_stop_overlap == 1, ]
# could also trim these exons to end at stop codon, but it's a bit tricky. just mask them and don't match to them
# remove stop columns now that last exons are marked
exons_stops[, c("chromStart_stop", "chromEnd_stop") := NULL]
exons_stops

# write to file for bedtools closest (exclude chrom M since not present in pA data)
# nc_exon_stop_overlap exons are excluded here
exons_stops_sel <- exons_stops[chrom_exon != "chrM" & nc_exon_stop_overlap %in% c(NA, 0), .(gene_biotype_exon, transcript_biotype_exon, component, chrom_exon, chromStart_exon, chromEnd_exon, ens_gene_id_exon, gene_name_exon, strand_exon, frame_exon, last_cds_exon)]

# for the purpose of this script, need to update frame to be in relation to the 3' end rather than the default 5' end
exons_stops_sel[, frame_three_exon := (chromEnd_exon - chromStart_exon - as.integer(frame_exon)) %% 3]

# summarize by gene (not all frames and last exon statuses are the same)
exons_stops_sel_summ <- exons_stops_sel[, .(frames_three_exon = list(sort(unique(frame_three_exon[is.na(frame_three_exon) == F]))), last_cds_exon_prop = sum(last_cds_exon, na.rm = TRUE) / .N, protein_coding_exon_prop = sum(transcript_biotype_exon == "protein_coding") / .N), by = .(gene_biotype_exon, chrom_exon, chromStart_exon, chromEnd_exon, ens_gene_id_exon, gene_name_exon, strand_exon)]
# 399,698
#View(exons_stops_sel_sum[!last_cds_exon_prop %in% c(0,1),])

# take longest exon for exons sharing a 3' end (by gene, don't merge across genes)
exons_stops_sel_summ[strand_exon == "+", three_prime_end := chromEnd_exon]
exons_stops_sel_summ[strand_exon == "-", three_prime_end := chromStart_exon]
exons_stops_longest <- exons_stops_sel_summ[, .(chromStart_exon_min = min(chromStart_exon), chromEnd_exon_max = max(chromEnd_exon), frames_three_exon = list(sort(unique(unlist(frames_three_exon)[unlist(is.na(frames_three_exon) == F)]))), last_cds_exon_ever = ifelse(sum(last_cds_exon_prop) > 0, 1, 0), protein_coding_exon_ever = ifelse(sum(protein_coding_exon_prop) > 0, 1, 0)), by = .(gene_biotype_exon, chrom_exon, ens_gene_id_exon, gene_name_exon, three_prime_end, strand_exon)]
# 324,457
exons_stops_longest[, frames_three_exon := sapply(frames_three_exon, function(x) paste(unlist(x), collapse = "|"))]
exons_stops_longest[frames_three_exon == "", frames_three_exon := "none"]
table(exons_stops_longest[, frames_three_exon])
exons_stops_longest[strand_exon == "+", chromStart_exon := chromStart_exon_min]
exons_stops_longest[strand_exon == "+", chromEnd_exon := three_prime_end]
exons_stops_longest[strand_exon == "-", chromStart_exon := three_prime_end]
exons_stops_longest[strand_exon == "-", chromEnd_exon := chromEnd_exon_max]
exons_stops_longest[, c("three_prime_end", "chromStart_exon_min", "chromEnd_exon_max") := NULL]
# careful with delims here. "-" is present in fusion genes, "_" present in "protein_coding" etc, "|" used in frames
exons_stops_longest_bed <- exons_stops_longest[, .(chrom_exon, chromStart_exon, chromEnd_exon, name = paste(gene_biotype_exon, ens_gene_id_exon, gene_name_exon, sep = "__"), score = paste(frames_three_exon, last_cds_exon_ever, protein_coding_exon_ever, sep = "_"), strand_exon)]
exons_stops_longest_bed # 324,457
setwd(path_data_out)
fwrite(exons_stops_longest_bed, "ensembl_101_prot_cod_merged_cds_exons.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) # note "protein_coding" name isn't entirely technically accurate here

# -----------------------------------------------------------------------


# Match polyA sites to their closest upstream exon -------------------------------------

# match unique pA sites to closest upstream exon using bedtools closest
# if the gene that matches is a match in the pADB data, keep it. if not, remove it.
setwd(path_data_out)
system("sort -V ensembl_101_prot_cod_merged_cds_exons.bed > nat_sorted_ensembl_101_prot_cod_merged_cds_exons.bed")
system("sort -V unique_prot_cod_polyA_DB_v3_hg38_liftover_sites.bed > nat_sorted_unique_prot_cod_polyA_DB_v3_hg38_liftover_sites.bed")
system("bedtools closest -a nat_sorted_unique_prot_cod_polyA_DB_v3_hg38_liftover_sites.bed -b nat_sorted_ensembl_101_prot_cod_merged_cds_exons.bed -s -id -io -D a > unique_prot_cod_polyA_DB_v3_hg38_liftover_sites_closest_upstream_prot_cod_cds_exon_ens_101.bed")

# read results
setwd(path_data_out)
site_exon_pairs <- fread("unique_prot_cod_polyA_DB_v3_hg38_liftover_sites_closest_upstream_prot_cod_cds_exon_ens_101.bed", header = FALSE)
names(site_exon_pairs) <- c(paste(names_BED_std, "site", sep = "_"), paste(names_BED_std, "exon", sep = "_"), "distance")

# check for sites that have no match
site_exon_pairs[chromStart_exon == -1, ] # 6
# remove them
site_exon_pairs <- site_exon_pairs[chromStart_exon != -1, ]

# extract gene ids and symbols
site_exon_pairs[, site_gene_ids := lapply(name_site, function(x) strsplit(x,  split = "_", fixed = TRUE)[[1]])]
site_exon_pairs[, site_gene_ids := lapply(site_gene_ids, function(x) x[is.na(x) == FALSE & x != "NA" & x != "NO"])]
site_exon_pairs[, site_gene_symbols := lapply(score_site, function(x) strsplit(x,  split = "_", fixed = TRUE)[[1]])]
site_exon_pairs[, site_gene_symbols := lapply(site_gene_symbols, function(x) x[is.na(x) == FALSE & x != "NA" & x != "NO"])]
site_exon_pairs[, c("name_site", "score_site") := NULL]
site_exon_pairs[, exon_gene_biotype := sapply(strsplit(name_exon, "__"), "[", 1)] # different delim
site_exon_pairs[, exon_gene_id := sapply(strsplit(name_exon, "__"), "[", 2)]
site_exon_pairs[, exon_gene_symbol := sapply(strsplit(name_exon, "__"), "[", 3)]
site_exon_pairs[, frames_three_exon := sapply(strsplit(score_exon, "_"), "[", 1)]
table(site_exon_pairs$frames_three_exon)
site_exon_pairs[, frames_three_exon := lapply(frames_three_exon, function(x) strsplit(x, split = "|", fixed = TRUE)[[1]])]
site_exon_pairs[, last_cds_exon_ever := sapply(strsplit(score_exon, "_"), "[", 2)]
site_exon_pairs[, protein_coding_exon_ever := sapply(strsplit(score_exon, "_"), "[", 3)]
site_exon_pairs[, c("name_exon", "score_exon") := NULL]
site_exon_pairs

# how many times is each site id in the table?
site_exon_pairs[, id_site := paste(chrom_site, chromStart_site, chromEnd_site, strand_site, sep = "_")]
site_exon_pairs[, id_site_count := .N, by = id_site]
table(site_exon_pairs[, id_site_count])
site_exon_pairs

# mark matches (any gene id or symbol from site in exon)
site_exon_pairs[, matching_id := mapply(function(s,e) ifelse(sum(s %in% e) > 0, 1, 0), s = site_gene_ids, e = exon_gene_id)]
site_exon_pairs[, matching_sym := mapply(function(s,e) ifelse(sum(s %in% e) > 0, 1, 0), s = site_gene_symbols, e = exon_gene_symbol)]
site_exon_pairs[, match_both := ifelse(matching_id == 1 & matching_sym == 1, 1, 0)]
table(site_exon_pairs[, match_both]) 
table(site_exon_pairs[, match_both]) / nrow(site_exon_pairs) # 93% matching

# filter mismatches and assess
site_exon_pairs_matched <- site_exon_pairs[match_both == 1, ]
site_exon_pairs_matched[, id_site_count := .N, by = id_site]
table(site_exon_pairs_matched[, id_site_count])
site_exon_pairs_matched[id_site_count > 1, ] # 1793 entries
# many of these are fusions. remove the fusion genes
site_exon_pairs_matched <- site_exon_pairs_matched[!grepl("-", site_exon_pairs_matched[, exon_gene_symbol]), ]
site_exon_pairs_matched[, id_site_count := .N, by = id_site]
table(site_exon_pairs_matched[, id_site_count])
site_exon_pairs_matched[id_site_count > 1, ] # 1079 entries
# many of these are genes with related names
# make note
site_exon_pairs_matched[id_site_count > 1, rep_gene := "yes"]
# sort by gene name, number them, keep first one
setkey(site_exon_pairs_matched, exon_gene_symbol)
site_exon_pairs_matched[, temp_num := 1:.N, by = id_site]
site_exon_pairs_matched[id_site_count > 1, ]
site_exon_pairs_matched <- site_exon_pairs_matched[temp_num == 1, ]
site_exon_pairs_matched[, id_site_count := .N, by = id_site]
table(site_exon_pairs_matched[, id_site_count])
site_exon_pairs_matched[id_site_count > 1, ] # 0 entries

# what prop of sites remain after filtering?
length(unique(site_exon_pairs_matched[, id_site])) # 173,350
prot_cod_polys_bed[, id_site := paste(chrom, chromStart_site, chromEnd_site, strand_site, sep = "_")]
length(unique(prot_cod_polys_bed[, id_site])) # 185,053
prot_cod_polys_bed[, id_site := NULL]
# = 94% of unique sites still used

# prep for merge back to original data (note some gene IDs / symbols are representative of multiple genes)
site_exon_pairs_matched_sel <- site_exon_pairs_matched[, .(chrom = chrom_site, chromStart_site, chromEnd_site, strand = strand_site, id_site, chromStart_exon, chromEnd_exon, exon_gene_id, exon_gene_symbol, frames_three_exon, last_cds_exon_ever, protein_coding_exon_ever, gene_biotype = exon_gene_biotype)]
polys_sel <- polys[, .(extended_prop = sum(extended) / .N, mean_RPM = mean(mean_RPM), PSE = mean(PSE)), by = .(id_site = name_site, PAS_Signal, conserved)]
polys_sel[, id_site_count := .N, by = id_site]
table(polys_sel[, id_site_count]) # should all be 1
polys_sel[id_site_count > 1, ]
polys_sel <-polys_sel[id_site_count == 1, ]
polys_sel[, id_site_count := NULL]
# merge
setkey(site_exon_pairs_matched_sel, id_site)
setkey(polys_sel, id_site)
matches <- polys_sel[site_exon_pairs_matched_sel]
matches

# write to file
setwd(path_data_out)
saveRDS(matches, "unique_prot_cod_polyA_DB_v3_hg38_liftover_sites_closest_upstream_prot_cod_cds_exon_ens_101_no_mismatch.rds")





