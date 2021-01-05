segmentation_fn=$(readlink -f $1)
output_folder=$(readlink -f $2)
prefix_for_state_annotation=$3 # states are annotated as E1 or U1 or something like that
num_state=$4
shDir=$(readlink -f $5)
shPrefix=divide_segmentation
job_index=$6

# echo ${segmentation_fn}
# echo ${output_folder}
# echo ${prefix_for_state_annotation}
# echo ${num_state}
# echo ${shDir}
# echo ${job_index}
declare -a chrom_list=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
function divide_segmentation_bed_by_state_chromosome_one_bed () {
	segmentation_fn=$1
	out_dir=$2
	mkdir -p $out_dir
	rm -f ${shDir}/${shPrefix}_${job_index}.sh
	# the following line of code is used for filtering snp gwax data into chromosome, ordering snps based on their genomic positions
	for chrom_index in "${chrom_list[@]}"
	do
		chrom_output_fn=$out_dir/chrom_${chrom_index}_all_states
		echo " 
		rm -f ${chrom_output_fn}.gz
		zcat $segmentation_fn | awk -v chrom_symbol='chr${chrom_index}' '{ if (\$1 == chrom_symbol) print \$0; }' | sort -n -k2 > $chrom_output_fn
		gzip $chrom_output_fn" >> ${shDir}/${shPrefix}_${job_index}.sh
	done

	# for state_index in `seq 1 $num_state`
	# do
	# 	state_folder=$out_dir/state_${state_index}/
	# 	mkdir -p $state_folder
	# 	for chrom_index in "${chrom_list[@]}"
	# 	do
	# 		this_chrom_segmentation_fn=$out_dir/chrom_${chrom_index}_all_states.gz
	# 		this_chrom_output_fn=$state_folder/segmentation_state_${state_index}_chr${chrom_index}.bed
	# 		echo "
	# 		rm $this_chrom_output_fn.gz
	# 		zcat $this_chrom_segmentation_fn | awk -v state_symbol='${prefix_for_state_annotation}${state_index}' '{ if (\$4 == state_symbol) print \$0; }' > $this_chrom_output_fn
	# 		gzip $this_chrom_output_fn" >> ${shDir}/${shPrefix}_${job_index}.sh
	# 	done
	# done

	# # remove the files that contain all segmentations in one chromosome because we already divide segmentation into chromosome and then state
	# for chrom_index in "${chrom_list[@]}"
	# do
	# 	echo "rm -f $out_dir/chrom_${chrom_index}_all_states.gz" >> ${shDir}/${shPrefix}_${job_index}.sh
	# done
	chmod +x ${shDir}/${shPrefix}_${job_index}.sh
	job_index=$(($job_index+1))
}

function divide_segmentation_by_state () {
	segmentation_fn=$1
	output_folder=$2
	mkdir -p $output_folder
	num_state=$3
	for state_index in `seq 1 $num_state`
	do
		this_state_output_fn=${output_folder}/state_${state_index}.bed
		echo "
		zcat $segmentation_fn | awk -F'\t' -v state_symbol='${prefix_for_state_annotation}${state_index}' '{ if (\$4 == state_symbol) print \$0; }' > $this_state_output_fn
		gzip $this_state_output_fn 
		" >> ${shDir}/${shPrefix}_${job_index}.sh
		chmod +x ${shDir}/${shPrefix}_${job_index}.sh
	done
}

divide_segmentation_bed_by_state_chromosome_one_bed ${segmentation_fn} $output_folder

# divide_segmentation_by_state $segmentation_fn $output_folder $num_state