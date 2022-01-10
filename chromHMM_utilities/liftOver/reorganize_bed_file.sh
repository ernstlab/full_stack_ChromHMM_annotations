echo "This file will reorganize your bed file, such that lines correspond to a chromosome are grouped together. And the starting coordinate in lines are ordered ascendingly."
echo "reorganize_bed_file.sh <input bed file> <output bed file>"
echo "input and output files should be in gzip format. You should always do that!"
echo "this file assume that there is a header line in the input file, which should NOT be copied to the output file"
echo "this code also assumes that the chrom_symbol in the input bed file is chr<chrom_index>. Please go into this file to fix the code if you want to adjust it. "
input_file=$(readlink -f $1)
output_fn=$(readlink -f $2)
rm -f $output_fn # so that we can write new one
rm -f ${output_fn}.gz
declare -a chrom_list=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

# copy the header line from the input file into the output file
# zcat $input_file | sed -n 1p > $output_fn 

# the following line of code is used for filtering snp gwax data into chromosome, ordering snps based on their genomic positions
for chrom_index in "${chrom_list[@]}"
do
	if [[ $input_file =~ \.t?gz$ ]]; # check if the input file has extension .gz or .tgz
	then
		zcat $input_file | sed -n '1,$p' | awk -F'\t' -v chrom_symbol="chr${chrom_index}" 'BEGIN {OFS="\t"} { if ($1 == chrom_symbol) print $0; }' | sort -k2n >> $output_fn # fix the chrom_symbol fomat if needed. It's mentioned in the instruction of this code
		echo "Done processing for chromosome $chrom_index"
	else
		cat $input_file | sed -n '1,$p' | awk -F'\t' -v chrom_symbol="chr${chrom_index}" 'BEGIN {OFS="\t"}{ if ($1 == chrom_symbol) print $0; }' | sort -k2n >> $output_fn # fix the chrom_symbol fomat if needed. It's mentioned in the instruction of this code
		echo "Done processing for chromosome $chrom_index"
	fi
done

echo "We zipped your output bed file. Therefore, the bed file name is now: ${output_fn}.gz"
