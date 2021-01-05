declare -a chrom_list=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
segmentation_fn=$(readlink -f $1)
output_fn=$(readlink -f $2)
rm -f $output_fn # so that we can create a new one
echo "Make sure: segmentation_fn is gzipped"
echo "output_fn has its folder created"

for chrom_index in "${chrom_list[@]}"
do
	this_chrom_length=$(zcat $segmentation_fn | awk -F'\t' -v chrom_symbol="chr${chrom_index}" '{if ($1 == chrom_symbol) print $0}' | sort -k3n | tail -1 | awk -F'\t' '{print $3}')
	echo -e "chr$chrom_index\t0\t$this_chrom_length" >> $output_fn
	echo "Done processing chromosome: $chrom_index"
done