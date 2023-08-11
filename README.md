# iMAPS calculator

Quantify negative selection for 3' UTR variants or regions of interest using the iMAPS approach described in Findlay et al. 2023: "Quantifying negative selection in human 3ʹ UTRs uncovers constrained targets of RNA-binding proteins"

## System requirements
  
R v4.2.0  
bedtools v2.29.1  
gsutil v 5.25  

Tested on:  
Mac OS 10.15.7  
Mac OS 13.4  
Ubuntu 20.04  

## Installation guide

```bash

# clone the repo
git clone https://github.com/sfindlay11/iMAPS.git

# download all required files into a single directory:

# from Google Cloud:
gsutil cp gs://imaps/gnomad_vars_annotated.bed .
gsutil cp gs://imaps/intergenic_calibration_data.txt.gz .
gsutil cp gs://imaps/demo_vars.bed .
# from this repository:
feature_intersection.sh

# file permissions for feature_intersection.sh may need to be updated using
chmod 777 feature_intersection.sh

```

## Demo

```bash

# run a demo calculating negative selection (iMAPS) for variants disrupting or preserving ReP sites with relative affinities >= 0.1
Rscript /path/to/imaps_calc.R /path/to/demo_vars.bed sense /path/to/required_files/

# expected output (approximate run time of 3 minutes):
imaps.by_category.DATA.txt

```

## Usage

```bash

Rscript /path/to/imaps_calc.R /path/to/demo_vars.bed [allele format] /path/to/required_files/ 

# allele fromat: "none", "sense", or "plus"; see below

# .bed format requirements:
# tab-separated values
# use 0-based hg38 coordinates
# provide single nucleotide variants or genomic ranges corresponding to features of interest
# no header

# columns:
# 1) chrom
# 2) chromStart
# 3) chromEnd
# 4) name = categorical variable on which to calculate iMAPS (e.g. "test" & "control"). If no comparison is desired, provide any string as a placeholder
# 5) score (optional; see below) = "reference_allele|alternative_allele"
# 6) strand (optional; see below)


# if providing features / genomic ranges of interest:
# use allele format "none" and provide a four column .bed file if genomic ranges are provided with no variant alleles

# if providing variants:
# use allele format "sense" and provide a six column .bed file with reference and alternative alleles corresponding to the sense strand in the 5th "score" column, separated by "|" (e.g. "A|C"). Provide the strand in column six.
# OR
# use allele format "plus" and provide a five column .bed file  with reference and alternative alleles corresponding to the plus strand in the 5th "score" column, separated by "|" (e.g. "A|C").

# ** CAUTION should be exercised when analyzing variants where the allele frequency spectrum may influence variant ascertainment or classification. **
# ** As one example, significant GWAS and QTL variants generally have relatively high allele frequencies as these analyses are not typically powered to detect effects from lower freuency variants. **
# ** Such sets of variants may not produce interpretable iMAPS results. **

```


## License

[MIT](https://choosealicense.com/licenses/mit/)