# this script is to compare observed and expected proportion singleton (or maps) values for k-mers using different calibration approaches
# 5-mers with variant at middle position

#### initial set up ------------------------------------------------------------------

library(data.table)
library(tidyverse)

path_data_in <- ""
path_data_proj <- ""
path_data_out <- ""
path_data_plots <- ""
path_data_calibs <- ""

names_BED_std <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

options(max.print = 5000)

stats_all <- NULL

#### END OF initial set up -----------------------------------------------------------------------


#### process variants ------------------------------------------------------------------

#setwd(path_data_in)
setwd(path_data_proj) # other than inlcuding the 3'UTR version of MAPS, only difference to above commented out file is addition of optional columns. can be added if desired.
variants <- fread("FOR_REVISION.pADB_infered_UTR3_plus_some_gencode.variants.main_chr_prot_cod.info.last_cds_ex_only.allele_freq_cats.context.calib_sel_intergenic.calib_utr3.txt.gz", header = TRUE)

variants[, di_con_pair := paste(substr(hexa_con_pair, 5, 11), substr(hexa_con_pair, 16, 17), sep = "")]

#### END OF process variants ------------------------------------------------------------------


#### add filters to motif variants ------------------------------------------------------------------

### set up shared filter (for all versions; additional filters added as appropriate for each version)
variants[, prim_filter_pass := ifelse(feature_conflict_filter_pass == 1 & len_filt_pass == 1 & last_cds_exon_ever == 1 & (!chrom %in% c("chrX", "chrY", "chrM")) & gnomAD_v3_cov_med > 24 & gnomAD_v3_cov_med < 43 & AN > 129064 & ((cat == "CpG_ti" & methyl_mean_cat != "N") | cat != "CpG_ti") & lcr == 0 & CDS_overlap == 0 & protein_coding_exon_ever == 1, TRUE, FALSE)]

# mark CpGs
variants[, is_cpg := grepl(".CG|CG.", dinuc)]

#### END OF add filters to motif variants ------------------------------------------------------------------


# add maps calibration method values
setwd(path_data_calibs)
calibs <- fread("maps_related_calibrations_by_dinuc_base_change.txt")
setkey(calibs, di_con_pair)
setkey(variants, di_con_pair)
variants <- calibs[, .(di_con_pair, prop_sing_for_calib_maps_all, prop_sing_for_calib_maps_split, prop_sing_for_calib_maps_raw)][variants]
rm(calibs)


#################### all 3'UTR ####################


# selection
OE <- variants[prim_filter_pass == TRUE
               & ancest_seq_strand == ref_strand
               
               & CpG_island == 0
               & base_change_pair != "A>G"
               
               #& is.na(prop_sing_for_calib_maps_all) == FALSE
               
               ,
               .(observed_singletons = sum(singleton),
                 
                 expected_singletons = sum(highest_sig_k_prop_sing_INT), # iMAPS
                 #expected_singletons = sum(highest_sig_k_prop_sing_UTR3), # 3'UTR MAPS
                 #expected_singletons = sum(prop_sing_for_calib_maps_all), # MAPS
                 #expected_singletons = sum(prop_sing_for_calib_maps_split), # MAPS split regressions
                 #expected_singletons = sum(prop_sing_for_calib_maps_raw), # just take syn prop singletons
                 
                 n_variants = .N),
               
               by = .(trinuc, is_cpg)]
               #by = .(trinuc, base_change, is_cpg)]

OE[, expected_prop_sing := expected_singletons / n_variants]
OE[, observed_prop_sing := observed_singletons / n_variants]
OE[, obs_div_exp := observed_prop_sing / expected_prop_sing]

# add linear models
# https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R
lin_mod <- lm(observed_prop_sing ~ expected_prop_sing, data = OE)
summary(lin_mod)
str(summary(lin_mod))

inter_curr <- summary(lin_mod)$coefficients[1, 1]
slope_curr <- summary(lin_mod)$coefficients[2, 1]
residuals_curr <- unname(summary(lin_mod)$residuals)
abs_resid_mean_curr <- mean(abs(residuals_curr))
sigma_curr <- summary(lin_mod)$sigma # residual standard error
r_sq_curr <- summary(lin_mod)$r.squared

# manually select version
#version_curr <- "imaps"
#version_curr <- "maps"
#version_curr <- "syn-raw"

stats_all <- rbind(stats_all, data.table(version = version_curr,
                                           intercept = inter_curr,
                                           slope = slope_curr,
                                           resid_list = list(residuals_curr),
                                           resid_std_err = sigma_curr,
                                           r_sq = r_sq_curr))
stats_all[, -"resid_list"]

# after all versions
#stats_all[, abs_resid_mean := sapply(resid_list, function(x) mean(abs(x)))]
# save
setwd(path_data_out)
#fwrite(stats_all, "diff_maps_calib_by_trinuc_regression_stats.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# save all versions of this plot

#ggplot(data = OE, aes(x = expected_singletons, y = observed_singletons)) +
ggplot(data = OE, aes(x = expected_prop_sing, y = observed_prop_sing)) +
  geom_point(alpha = 0.1) +
  #geom_point(size = 2, alpha = 0.25) +
  geom_abline(slope = slope_curr, intercept = inter_curr) +
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray70", linewidth = 0.5) +
  coord_cartesian(xlim = c(0.1, 0.65), ylim = c(0.1, 0.65)) + # all
  #coord_cartesian(xlim = c(0.1, 0.7), ylim = c(0.1, 0.7)) + # all (trinuc + base change)
  #coord_cartesian(xlim = c(0.4, 0.65), ylim = c(0.3, 0.675)) + # excluding CpG ti
  theme_minimal(base_size = 14) +
  labs(x = "\nExpected proportion singleton", y = "Observed proportion singleton\n")
setwd(path_data_plots)
OE <- fread("obs_exp_prop_sing.maps.DATA.txt")
cor(OE$observed_prop_sing, OE$expected_prop_sing)^2
# save
#fwrite(OE, paste(path_data_plots, "obs_exp_prop_sing.imaps.DATA.txt", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#ggsave("obs_exp_prop_sing.imaps.WIDE.pdf", path = path_data_plots, width = 10, height = 5, units = "in")


## correlations
# all together
cor(x = OE$expected_prop_sing, y = OE$observed_prop_sing)^2 
# iMAPS: 0.9865
# MAPS: 0.9485
# MAPS split: 0.9563
# syn raw: 0.9555

# CpGs
cor(x = OE[is_cpg == TRUE, expected_prop_sing], y = OE[is_cpg == TRUE, observed_prop_sing])^2 
# iMAPS: 0.9316
# 3'UTR MAPS: 0.9780
# MAPS: 0.5752
# MAPS split: 0.5911
# syn raw: 0.5931

# non-CpGs
cor(x = OE[is_cpg == FALSE, expected_prop_sing], y = OE[is_cpg == FALSE, observed_prop_sing])^2 
# iMAPS: 0.6040
# 3'UTR MAPS: 0.6792
# MAPS: 0.3147
# MAPS split: 0.3147
# syn raw: 0.3382

## residual analysis
OE[, residual := observed_prop_sing - ((slope_curr * expected_prop_sing) + inter_curr)]

#OE_to_compare_1 <- OE[, .(trinuc, residual_1 = residual)] 
#OE_to_compare_2 <- OE[, .(trinuc, residual_2 = residual)] 
setkey(OE_to_compare_1, trinuc)
setkey(OE_to_compare_2, trinuc)
OE_to_compare <- OE_to_compare_2[OE_to_compare_1]
rm(OE_to_compare_1, OE_to_compare_2)

# compare
sum(OE_to_compare[, abs(residual_2) == abs(residual_1)])
sum(OE_to_compare[, abs(residual_2) < abs(residual_1)])
sum(OE_to_compare[, abs(residual_2) < abs(residual_1)]) / nrow(OE_to_compare)
sum(OE_to_compare[, abs(residual_2) > abs(residual_1)])
sum(OE_to_compare[, abs(residual_2) > abs(residual_1)]) / nrow(OE_to_compare)
odds <- sum(OE_to_compare[, abs(residual_2) < abs(residual_1)]) / sum(OE_to_compare[, abs(residual_2) > abs(residual_1)])
OE_to_compare[, mean_abs_improvement := abs(residual_1) - abs(residual_2)]
mean(OE_to_compare$mean_abs_improvement)
mean(abs(OE_to_compare$residual_1))
mean(OE_to_compare$mean_abs_improvement) / mean(abs(OE_to_compare$residual_1))
median(OE_to_compare$mean_abs_improvement)
median(OE_to_compare$mean_abs_improvement) / mean(abs(OE_to_compare$residual_1))
# imaps vs maps: 650 improved; 374 worse; odds = 1.737968
# mean = 0.005770643 (as percent of average residual = 0.005770643 / 0.01497162 = 0.3854386); median = 0.003002026 (as percent of average residual = 0.003002026 / 0.01497162 = 0.2005144)


# plot correlation
#ggplot(data = OE_to_compare, aes(x = residual_1, y = residual_2)) +
ggplot(data = OE_to_compare, aes(x = abs(residual_1), y = abs(residual_2))) +
  geom_point(alpha = 0.1) +
  #geom_point(size = 2, alpha = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.25) +
  geom_vline(xintercept = 0, linetype = "solid", size = 0.25) +
  #coord_cartesian(xlim = c(-0.11, 0.11), ylim = c(-0.11, 0.11)) +
  #coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) +
  coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 0.2)) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  #coord_cartesian(xlim = c(-0, 0.11), ylim = c(0, 0.11)) +
  #scale_x_continuous(breaks = c(0, 0.05, 0.1)) +
  #scale_y_continuous(breaks = c(0, 0.05, 0.1)) +
  theme_minimal(base_size = 14) +
  #labs(x = "\niMAPS residual", y = "3' UTR MAPS residual\n")
  #labs(x = "\niMAPS abs(residual)", y = "3' UTR MAPS abs(residual)\n")
  #labs(x = "\nMAPS residual", y = "iMAPS residual\n")
  labs(x = "\nMAPS abs(residual)", y = "iMAPS abs(residual)\n")

cor(x = OE_to_compare$residual_1, y = OE_to_compare$residual_2)^2 # 3'UTR maps vs imaps: 0.4786

#fwrite(OE_to_compare, paste(path_data_plots, "residuals.imaps_vs_maps.DATA.txt", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#ggsave("residuals.imaps_vs_maps.pdf", path = path_data_plots, width = 5, height = 5, units = "in")

# plot cdf of delta residual
ggplot(data = OE_to_compare, aes(x = residual_1, y = residual_2)) +
  stat_ecdf() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #coord_cartesian(xlim = c(0.1, 0.65)) +
  theme_minimal(base_size = 14) +
  labs(x = "\nExpected proportion singleton", y = "Residual\n")

# waterfall
OE_to_compare <- OE_to_compare[order(mean_abs_improvement, decreasing = TRUE),]
OE_to_compare[, order := 1:nrow(OE_to_compare)]
OE_to_compare[, plot_colour := ifelse(mean_abs_improvement > 0, "blue", "gray")]
ggplot(data = OE_to_compare, aes(x = order, y = mean_abs_improvement, colour = plot_colour, fill = plot_colour)) +
  geom_bar(stat = "identity", position = "identity", show.legend = FALSE) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = c(1, 1024), minor_breaks = NULL, labels = c(1, 1024)) +
  scale_colour_manual(values = c("dodger blue", "gray70")) +
  scale_fill_manual(values = c("dodger blue", "gray70")) +
  labs(x = "\n5-mer rank", y = "improvement of residual using iMAPS vs. MAPS\n") +
  coord_cartesian(ylim = c(-0.2, 0.2))
# calc areas
abs(sum(OE_to_compare[mean_abs_improvement > 0, mean_abs_improvement])) # 8.89
abs(sum(OE_to_compare[mean_abs_improvement < 0, mean_abs_improvement])) # 2.98

#ggsave("residual_improvements.utr3maps_vs_imaps.waterfall.pdf", path = path_data_plots, width = 5, height = 5, units = "in")

# merge
setkey(resid_imaps, trinuc)
setkey(resid_maps, trinuc)
setkey(resid_maps_split, trinuc)
setkey(resid_syn_raw, trinuc)
resid_all <- Reduce(merge,list(resid_imaps, resid_maps, resid_maps_split, resid_syn_raw))
# save
setwd(path_data_out)
#fwrite(resid_all, "diff_maps_calib_by_trinuc_regression_residuals.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


#################### end of all 3'UTR ####################



