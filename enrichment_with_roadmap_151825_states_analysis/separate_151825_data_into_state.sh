separate_state_code='divide_segmentation_bed_by_chrom_state.sh'
shDir=/u/home/h/havu73/project-ernst/source/job_jobs # TO BE UPDATED: change this folder to the relative path where you want to store all the .sh files that this code will produce. These .sh files will then be run on computing cluster or on your local computer. If you submit these jobs using the computing cluster, you can run them in parallel and that will save a lot fo time. 

# all suffixes are written with .gz
download_link_suffix_15='_15_coreMarks_segments.bed.gz'
download_link_suffix_18='_18_core_K27ac_segments.bed.gz'
download_link_suffix_25="_25_imputed12marks_segments.bed.gz"

segment_folder_15_state='/u/home/h/havu73/project-ernst/data/roadmap_epigenome/15_core_model_downloaded' # TO BE UPDATED: change this folder to the relative path where you want to store the 15-state concatenated annotations for 127 reference epigenomes. 
segment_folder_18_state=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded # TO BE UPDATED: change this folder to the relative path where you want to store the 18-state concatenated annotations for 98 reference epigenomes. 
segment_folder_25_state='/u/home/h/havu73/project-ernst/data/roadmap_epigenome/25_imputed12marks_model_downloaded' # TO BE UPDATED: change this folder to the relative path where you want to store the 25-state concatenated annotations for 127 reference epigenomes. 

job_index=1
for ct_folder in $segment_folder_25_state/*/
do
	ct=$(echo $ct_folder | awk -F'/' '{print $(NF-1) }')
	this_ct_state_segment_folder=${ct_folder}/state_segments/
	mkdir -p ${this_ct_state_segment_folder}
	this_ct_segment_fn=${ct_folder}/${ct}${download_link_suffix_25}
	if [ -f $this_ct_segment_fn ];
	then
		${separate_state_code}  $this_ct_segment_fn $this_ct_state_segment_folder 'E' 25 $shDir $job_index
		job_index=$(($job_index + 1))
	fi
done
