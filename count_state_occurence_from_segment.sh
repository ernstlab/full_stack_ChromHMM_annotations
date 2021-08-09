input_segment_fn=$(readlink -f $1)
output_fn=$(readlink -f $2)
num_state=$3


###### this file will count the number of occurences of each of the chromHMM states into output file

# rm -f $output_fn # so that I can create a new output file
# for state_index in `seq 1 $num_state`
# do 
# 	num_bin_in_this_state=$(zcat $input_segment_fn | awk -F'\t' -v state_symbol="E${state_index}" -v x="chrX" -v y="chrY" '{if ($4 == state_symbol && $1!=x && $1 !=y) print $0}' | wc -l | awk '{print $1}')
# 	echo "$state_index,$num_bin_in_this_state" >> $output_fn
# 	echo "Done with processing state $state_index"
# done

# echo "Done counting the number of occurences of states"

# bedtools intersect -a ~/project-ernst/data/cosmic_genome/Cosmic_NCV_38_noCommon_wgs_bed_files/mutation_oc
# c1.bed.gz  -b $input_segment_fn -wb > ~/project
# -ernst/data/cosmic_genome/Cosmic_NCV_38_noCommon_wgs_bed_files/full_model/draft_occ1

for state_index in `seq 1 $num_state`
do 
	num_mutation_in_this_state=$(cat ~/project-ernst/data/cosmic_genome/Cosmic_NCV_38_noCommon_wgs_bed_files/full_model/draft_occ1 | awk -F'\t'  -v state_symbol="E${state_index}" '{if ($7 == state_symbol) print $0}' | wc -l | awk '{print $1}')
	echo "$state_index,$num_bin_in_this_state" >> ~/project-ernst/data/cosmic_genome/Cosmic_NCV_38_noCommon_wgs_bed_files/full_model/count_state_occ1
done
