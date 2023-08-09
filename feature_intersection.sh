# intersects 3'UTR gnomAD variants with demo / user provided variants or regions

echo "running intersection for $1"

if [ "$2" = "sense" ];
then

  # strands must match for intersection
  bedtools intersect -s -loj -a gnomad_vars_annotated.bed -b "$1" > features_intersected.txt
  awk '{ if($7 != ".") { print }}' features_intersected.txt > features_intersected.filtered.txt
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $10, $11}' features_intersected.filtered.txt > features_intersected.filtered.sel_cols.txt

elif [ "$2" = "plus" ];
then

  # no requirement for strand matching
  bedtools intersect -loj -a gnomad_vars_annotated.bed -b "$1" > features_intersected.txt
  awk '{ if($7 != ".") { print }}' features_intersected.txt > features_intersected.filtered.txt
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $10, $11}' features_intersected.filtered.txt > features_intersected.filtered.sel_cols.txt

elif [ "$2" = "none" ];
then

  # no requirement for strand matching
  bedtools intersect -loj -a gnomad_vars_annotated.bed -b "$1" > features_intersected.txt
  awk '{ if($7 != ".") { print }}' features_intersected.txt > features_intersected.filtered.txt
  # no allels to output
  awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $10}' features_intersected.filtered.txt > features_intersected.filtered.sel_cols.txt

else
  echo "invalid allele version provided"  
fi

n_entries=$(wc -l $1 | awk '{print $1}')
n_int=$(wc -l features_intersected.filtered.txt | awk '{print $1}')
echo "there are ${n_int} 3'UTR gnomAD variants overlapping your ${n_entries} entries provided" # expansion is common and is assessed downstream

rm features_intersected.txt
rm features_intersected.filtered.txt