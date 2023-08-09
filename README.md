# iMAPS

Quantify negative selection for 3' UTR variants or regions of interest using our iMAPS approach described in Findlay et al. 2023: "Quantifying negative selection in human 3ʹ UTRs uncovers constrained targets of RNA-binding proteins"

## System requirements

Tested on:
R version 4.2.0 
bedtools v2.29.1 
gsutil v 5.25 

## Installation guide

```bash

# download required files
gsutil cp gs://imaps/gnomad_vars_annotated.bed .
gsutil cp gs://imaps/intergenic_calibration_data.txt.gz .
gsutil cp gs://imaps/demo_vars.bed .

```

## Demo

```bash

# run a demo calculating negative selection (iMAPS) for variants disrupting or preserving ReP sites with relative affinities >= 0.1
Rscript path/to/imaps_calc.R path/to/demo_vars.bed sense path/to/required_files/

# expected output (approximate run time of 3 minutes):
imaps.by_category.DATA.txt

```

## Usage

```bash

# provide imaps_calc.R with:
# 1) variant.bed file
# 2) the allele format used in for .bed file ("sense", "plus", or "none")
# 3) a directory containing data (will also be used for output)

Rscript path/to/imaps_calc.R path/to/demo_vars.bed [allele format] path/to/required_files/ 

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
# 5) score (optional) = "reference_allele|alternative_allele"
# 6) strand (optional)

# use allele format "none" and provide a four column .bed file if genomic ranges are provided with no variant alleles
# use allele format "sense" and provide a six column .bed file with reference and alternative alleles corresponding to the sense strand in the 5th "score" column, separated by "|" (e.g. "A|C"). Provide the strand in column six.
# use allele format "plus" and provide a five column .bed file  with reference and alternative alleles corresponding to the plus strand in the 5th "score" column, separated by "|" (e.g. "A|C").

```


## License

[MIT](https://choosealicense.com/licenses/mit/)