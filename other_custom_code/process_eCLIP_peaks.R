# see file descriptions at end of script for broad overview 

# load packages
library(data.table)
library(stringr)
library(pbapply)
library(parallel)
library(EnvStats)

# set paths
path_data_in <- ""
path_data_out <- ""
path_data_regions <- ""
path_data_genes <- ""

# set main values
names_BED_std <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")
names_chr_main <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")

# download all bed narrowPeak files from ENCODE (https://www.encodeproject.org/cart-view/)
# returns files.txt 
# rename encode_eCLIP_files_all_BED_narow_peak.txt and move to ENCODE_eCLIP_dta sub dir in proj dir

# a lot of lines are first time only and are currently commented 

setwd(path_data_in)

#! use curl to get download files (first time only)
#system("xargs -L 1 curl -O -L < encode_eCLIP_files_all_BED_narow_peak.txt")

## read in associated metadata and select desired files to keep
metadata <- fread("metadata.tsv", head = TRUE)
metadata

# first time only...

## select desired files
nrow(metadata)

table(metadata[`File assembly` == "GRCh38", `Biological replicate(s)`])
table(metadata[`File assembly` == "GRCh38", `Biosample term name`])

## hg38 and biological replicates 1 and 2 ("1,2") (IDR-passing peaks) =======================================
files_to_select <- metadata[`File assembly` == "GRCh38" & `Biological replicate(s)` == "1, 2" & `Biosample term name` %in% c("K562", "HepG2"), `File accession`]
length(files_to_select)

other_files <- metadata[!`File accession` %in% files_to_select, `File accession`]

# format
files_to_select <- paste(files_to_select, ".bed.gz", sep = "")
files_to_select <- paste(files_to_select, collapse = " ")

other_files <- paste(other_files, ".bed.gz", sep = "")
length(other_files)
other_files <- paste(other_files, collapse = " ")

# copy selected files to new sub dir (first time only)
setwd(path_data_in)
#system(paste("cp", files_to_select, "eCLIP_GRCh38_IDR_of_rep1_rep2", sep = " "))
# note: had to manually download ENCFF230QOU (was missing)

# cat (first time only)
setwd(paste(path_data_in, "eCLIP_GRCh38_IDR_of_rep1_rep2/", sep = ""))
#system(paste("cat", files_to_select, "> eCLIP_GRCh38_IDR_of_rep1_rep2_narrow_peak_all.bed.gz", sep = " "))

# read in
setwd(paste(path_data_in, "eCLIP_GRCh38_IDR_of_rep1_rep2/", sep = ""))
clip_peaks <- fread(file = "eCLIP_GRCh38_IDR_of_rep1_rep2_narrow_peak_all.bed.gz", header = FALSE)


## hg38 and biological replicates 1 OR 2 (peaks called from CLIPer BEFORE IDR) =======================================
files_to_select <- metadata[`File assembly` == "GRCh38" & `Biological replicate(s)` %in% c("1","2") & `Biosample term name` %in% c("K562", "HepG2"), `File accession`]
length(files_to_select)

other_files <- metadata[!`File accession` %in% files_to_select, `File accession`]

# format
files_to_select <- paste(files_to_select, ".bed.gz", sep = "")
files_to_select <- paste(files_to_select, collapse = " ")

other_files <- paste(other_files, ".bed.gz", sep = "")
length(other_files)
other_files <- paste(other_files, collapse = " ")

# copy selected files to new sub dir (first time only)
setwd(path_data_in)
#system(paste("cp", files_to_select, "eCLIP_GRCh38_rep1_rep2", sep = " "))
# note: had to manually download ENCFF230QOU (was missing)

# cat (first time only)
setwd(paste(path_data_in, "eCLIP_GRCh38_rep1_rep2/", sep = ""))
#system(paste("cat", files_to_select, "> eCLIP_GRCh38_IDR_of_rep1_rep2_narrow_peak_all.bed.gz", sep = " "))

# read in
setwd(paste(path_data_in, "eCLIP_GRCh38_rep1_rep2/", sep = ""))
clip_peaks <- fread(file = "eCLIP_GRCh38_IDR_of_rep1_rep2_narrow_peak_all.bed.gz", header = FALSE)

# column names described in: https://genome.ucsc.edu/goldenPath/help/bigNarrowPeak.html

## string chrom;        "Reference sequence chromosome or scaffold"
## uint   chromStart;   "Start position in chromosome"
## uint   chromEnd;     "End position in chromosome"
## string name;	 "Name given to a region (preferably unique). Use . if no name is assigned"
## uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "
## char[1]  strand;     "+ or - or . for unknown"
## float  signalValue;  "Measurement of average enrichment for the region"
## float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
## float  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
## int   peak;         "Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called."

## start of common portion =======================================

# process
names(clip_peaks) <- c(names_BED_std, "enrichment", "neg_log_10_p_value", "neg_log_10_q_value", "peak_point_source")
clip_peaks[, enrichment := signif(enrichment, digits = 4)]
clip_peaks[, neg_log_10_p_value := signif(neg_log_10_p_value, digits = 4)]
nrow(clip_peaks) # 848,077 peaks for IDR; 58,470,810
clip_peaks

# score feature 
table(clip_peaks[, score])

# extract name components
# for IDR
clip_peaks[, c("RBP", "cell_line", "IDR_string") := tstrsplit(name, "_", fixed = TRUE)]
clip_peaks[, c("name", "score", "IDR_string", "neg_log_10_q_value", "peak_point_source") := NULL]
clip_peaks[, score := 0]
# for rep1 rep2
#clip_peaks[, c("RBP", "cell_line", "rep_string") := tstrsplit(name, "_", fixed = TRUE)]
#clip_peaks[, replicate := gsub("rep0", "", rep_string, fixed = T)]
#clip_peaks[, c("name", "rep_string") := NULL]
#clip_peaks

# filter non-main chromosomes
clip_peaks <- clip_peaks[chrom %in% names_chr_main, ]
nrow(clip_peaks) # 847,805 or 58,450,084
clip_peaks

# explore
View(as.data.table(table(clip_peaks[, RBP])))

# store peak info in name
# for IDR
clip_peaks[, name := paste(RBP, cell_line, chrom, chromStart, chromEnd, strand, enrichment, neg_log_10_p_value, sep = "_")]
clip_peaks
# for rep1 rep2
#clip_peaks[, name := paste(RBP, cell_line, chrom, chromStart, chromEnd, strand, enrichment, neg_log_10_p_value, sep = "_")]
#clip_peaks[, name := paste(name, replicate, sep = "|")]
#clip_peaks

# any duplicate names?
table(table(clip_peaks$name))

# write
## rep 1 rep 2
### peaks
setwd(path_data_out)
fwrite(clip_peaks[, ..names_BED_std], "eCLIP_peaks.GRCh38.rep1_rep2.main_chr.full_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[cell_line == "HepG2", ..names_BED_std], "eCLIP_peaks.GRCh38.rep1_rep2.main_chr.HepG2.full_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[cell_line == "K562", ..names_BED_std], "eCLIP_peaks.GRCh38.rep1_rep2.main_chr.K562.full_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[replicate == "1", ..names_BED_std], "eCLIP_peaks.GRCh38.rep1.main_chr.full_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[replicate == "2", ..names_BED_std], "eCLIP_peaks.GRCh38.rep2.main_chr.full_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

## IDR pass 
### peaks
setwd(path_data_out)
fwrite(clip_peaks[, ..names_BED_std], "eCLIP_peaks.GRCh38.IDR_pass.main_chr.full_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[cell_line == "HepG2", ..names_BED_std], "eCLIP_peaks.GRCh38.IDR_pass.main_chr.HepG2.full_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[cell_line == "K562", ..names_BED_std], "eCLIP_peaks.GRCh38.IDR_pass.main_chr.K562.full_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
### RBPs
setwd(path_data_out)
fwrite(data.table(RBP = sort(unique(clip_peaks[, RBP]))), "RBPs_with_GRCh38_IDR_pass_eCLIP_peaks_in_main_chr.bed", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
fwrite(data.table(RBP = sort(unique(clip_peaks[cell_line == "HepG2", RBP]))), "RBPs_with_GRCh38_IDR_pass_eCLIP_peaks_in_main_chr_HepG2.bed", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
fwrite(data.table(RBP = sort(unique(clip_peaks[cell_line == "K562", RBP]))), "RBPs_with_GRCh38_IDR_pass_eCLIP_peaks_in_main_chr_K562.bed", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# explore
setwd(path_data_out)
clip_peaks <- fread("eCLIP_peaks.GRCh38.IDR_pass.main_chr.full_name.bed")
names(clip_peaks) <- names_BED_std
summary(clip_peaks[, chromEnd - chromStart])

#################### re-format relative to 5' end of peaks ####################

# write bed files for 5' ENDS of eCLIP peaks
clip_peaks[strand == "+", chromStart_five_end := chromStart]
clip_peaks[strand == "+", chromEnd_five_end := chromStart_five_end + 1]
clip_peaks[strand == "-", chromEnd_five_end := chromEnd]
clip_peaks[strand == "-", chromStart_five_end := chromEnd_five_end - 1]

## for rep 1 rep 2
setwd(path_data_out)
fwrite(clip_peaks[, .(chrom, chromStart_five_end, chromEnd_five_end, paste(name, chromStart_five_end, chromEnd_five_end, sep = "_"), score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.full_peak_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[cell_line == "HepG2", .(chrom, chromStart_five_end, chromEnd_five_end, paste(name, chromStart_five_end, chromEnd_five_end, sep = "_"), score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.HepG2.full_peak_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[cell_line == "K562", .(chrom, chromStart_five_end, chromEnd_five_end, paste(name, chromStart_five_end, chromEnd_five_end, sep = "_"), score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.K562.full_peak_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[replicate == "1", .(chrom, chromStart_five_end, chromEnd_five_end, paste(name, chromStart_five_end, chromEnd_five_end, sep = "_"), score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.rep1.main_chr.full_peak_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[replicate == "2", .(chrom, chromStart_five_end, chromEnd_five_end, paste(name, chromStart_five_end, chromEnd_five_end, sep = "_"), score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.rep2.main_chr.full_peak_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

## for IDR
setwd(path_data_out)
fwrite(clip_peaks[, .(chrom, chromStart_five_end, chromEnd_five_end, paste(name, chromStart_five_end, chromEnd_five_end, sep = "_"), score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_peak_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[cell_line == "HepG2", .(chrom, chromStart_five_end, chromEnd_five_end, paste(name, chromStart_five_end, chromEnd_five_end, sep = "_"), score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.HepG2.full_peak_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(clip_peaks[cell_line == "K562", .(chrom, chromStart_five_end, chromEnd_five_end, paste(name, chromStart_five_end, chromEnd_five_end, sep = "_"), score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.K562.full_peak_name.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# add gene and gene name and region (for 5' end base) based on CLIPPer definitions / genome segmentation
# see clipper_GRCh38_genes.bed, clipper_GRCh38_cds.bed, etc.

# join all regions
# (in terminal)
# cd ~/_/R_projects/ENCODE_eCLIP/data_files
# awk -v OFS='\t' '{print $1, $2, $3, $4, "utr5", $6}' clipper_GRCh38_five_prime_utrs.bed > clipper_GRCh38_five_prime_utrs_mod.bed
# awk -v OFS='\t' '{print $1, $2, $3, $4, "cds", $6}' clipper_GRCh38_cds.bed > clipper_GRCh38_cds_mod.bed
# awk -v OFS='\t' '{print $1, $2, $3, $4, "intron", $6}' clipper_GRCh38_introns.bed > clipper_GRCh38_introns_mod.bed
# awk -v OFS='\t' '{print $1, $2, $3, $4, "utr3", $6}' clipper_GRCh38_three_prime_utrs.bed > clipper_GRCh38_three_prime_utrs_mod.bed
# cat clipper_GRCh38_five_prime_utrs_mod.bed clipper_GRCh38_cds_mod.bed clipper_GRCh38_introns_mod.bed clipper_GRCh38_three_prime_utrs_mod.bed > clipper_GRCh38_main_annots_joined.bed

# intersect eclip 5' ends (not peaks!) with regions
# note: using loj to return peaks without matches to protein coding genes. will be recoded as NA below. presumably from non-coding genes
setwd(path_data_out)
system("bedtools intersect -s -loj -a eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_peak_name.bed -b clipper_GRCh38_main_annots_joined.bed > eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_peak_name.region.txt")
system("bedtools intersect -s -loj -a eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.full_peak_name.bed -b clipper_GRCh38_main_annots_joined.bed > eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.full_peak_name.region.txt")
# awk -v OFS='\t' '{print $1, $2, $3, $4, $10"_"$11, $6}' eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_peak_name.region.txt > eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_peak_name.region.select_cols.bed
# awk -v OFS='\t' '{print $1, $2, $3, $4, $10"_"$11, $6}' eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.full_peak_name.region.txt > eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.full_peak_name.region.select_cols.bed
system("rm eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_peak_name.region.txt")
system("rm eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.full_peak_name.region.txt")

# tidy (if running rep 1, rep 2)
rm(clip_peaks)
gc()

# read in
setwd(path_data_out)
fives <- fread("eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_peak_name.region.select_cols.bed")
#fives <- fread("eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.full_peak_name.region.select_cols.bed")
names(fives) <- c("chrom", "chromStart_five_end", "chromEnd_five_end", "name", "score", "strand")

# reformat entries with no matching region in b
View(as.data.table(table(fives[, score])))
fives[score == "._-1", score := NA]

# extract name and score info
fives[, c("RBP", "cell_line", "chrom_extra", "chromStart_peak", "chromEnd_peak", "strand_extra", "enrichment", "neg_log_10_p_value", "chromStart_five_end", "chromEnd_five_end") := tstrsplit(name, "_", fixed = TRUE)]
fives[, c("gene_id", "region") := tstrsplit(score, "_", fixed = TRUE)]
fives[, c("chrom_extra", "strand_extra") := NULL]
fives[, gene_id_less := gsub("[.]\\d*", "", gene_id)]

# add gene type
setwd(path_data_genes)
genes <- fread("hg38_kg_gencode_v32_main_chr_prot_cod_all_info.txt")

fives[, gene_type := ifelse(is.na(gene_id), NA, ifelse(gene_id_less %in% genes[gene_type == "protein_coding", gene_id_less], "protein_coding", "other"))]
table(fives[, gene_type], useNA = "ifany")

# rank regions
fives[region == "cds", region_rank := 1]
fives[region == "utr3", region_rank := 2]
fives[region == "utr5", region_rank := 3]
fives[region == "intron", region_rank := 4]

# summarize by peak (ie. some peaks might have multiple intersections)
fives_summ <- fives[, .(gene_count = length(unique(gene_id_less)),
                        prot_cod_prop = signif(sum(gene_type == "protein_coding") / .N, 3),
                        region_count = length(unique(region)),
                        region_top = unique(region[region_rank == min(region_rank)]),
                        gene_id_prot_cod_string = paste(sort(unique(gene_id_less[gene_type == "protein_coding"])), collapse = "|")
                        ),
                    by = .(chrom, chromStart_five_end, chromEnd_five_end, name, strand, cell_line)]

# correct string NAs for genes with no region match
fives_summ[name %in% fives[is.na(gene_id), name], gene_count := NA]
fives_summ[name %in% fives[is.na(gene_id), name], region_count := NA]
fives_summ[name %in% fives[is.na(gene_id), name], prot_cod_prop := NA]
fives_summ[name %in% fives[is.na(gene_id), name], region_top := NA]
fives_summ[name %in% fives[is.na(gene_id), name], gene_id_prot_cod_string := NA]
fives_summ[gene_id_prot_cod_string == "", gene_id_prot_cod_string := NA] # for peaks with intersections but not to any protein coding genes

# explore
table(fives_summ[, gene_count], useNA = "ifany") # 
table(fives_summ[, prot_cod_prop], useNA = "ifany") # 
table(fives_summ[, gene_count], fives_summ[, prot_cod_prop], useNA = "ifany")
table(fives_summ[, region_count], useNA = "ifany")

# repack score
fives_summ[, score := paste(gene_id_prot_cod_string, gene_count, prot_cod_prop, region_count, region_top, sep = "_")]

# write to file
setwd(path_data_out)

## rep1, rep 2
fwrite(fives_summ[, .(chrom, chromStart_five_end, chromEnd_five_end, name, score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.full_name.gene_region_score.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(fives_summ[cell_line == "HepG2", .(chrom, chromStart_five_end, chromEnd_five_end, name, score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.HepG2.full_name.gene_region_score.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(fives_summ[cell_line == "K562", .(chrom, chromStart_five_end, chromEnd_five_end, name, score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.rep1_rep2.main_chr.K562.full_name.gene_region_score.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#fwrite(fives_summ[replicate == "1", .(chrom, chromStart_five_end, chromEnd_five_end, name, score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.rep1.main_chr.K562.full_name.gene_region_score.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#fwrite(fives_summ[replicate == "2", .(chrom, chromStart_five_end, chromEnd_five_end, name, score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.rep2.main_chr.K562.full_name.gene_region_score.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# (need to expand replicate upstream for these last two versions)

## IDR
fwrite(fives_summ[, .(chrom, chromStart_five_end, chromEnd_five_end, name, score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.full_name.gene_region_score.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(fives_summ[cell_line == "HepG2", .(chrom, chromStart_five_end, chromEnd_five_end, name, score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.HepG2.full_name.gene_region_score.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(fives_summ[cell_line == "K562", .(chrom, chromStart_five_end, chromEnd_five_end, name, score, strand)], "eCLIP_peak_five_prime_ends.GRCh38.IDR_pass.main_chr.K562.full_name.gene_region_score.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#################### END of 5' end of peaks ####################

