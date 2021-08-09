raw_segmentation_folder='/u/home/h/havu73/project-ernst/data/roadmap_epigenome/25_imputed12marks_model_downloaded/hg19/'
state_annot_fn=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/25_imputed12marks_model_downloaded/state_annotation.txt
online_data_dir=/u/project/ernst/PUBLIC_SHARED/havu73/full_stack/igv_source/
job_dir=/u/home/h/havu73/project-ernst/source/job_jobs
########################
##### QUERY REGION #####
NUM_BP_PER_WINDOW=10000000
region_chrom='chr10'
region_start=30000000
region_end=40000000
region_index=$(echo $(( ${region_start}/${NUM_BP_PER_WINDOW} )) ) # this is already rounded down in bash
region_out_dir=${online_data_dir}/${region_chrom}_${region_index}
mkdir -p ${region_out_dir}
########################

########################
## FUNCTION TO CREATE THE SEGMENTATION DATA FOR DIFFERENT CELL GROUPS AND PUT IT UP ONLINE
########################
function create_igv_segment_one_cell_group(){
	cell_group=$1
	cg_dir=${region_out_dir}/${cell_group}
	mkdir -p ${cg_dir}
	sample_fn=/u/home/h/havu73/project-ernst/diff_pete/roadmap/${cell_group}/sample.list
	reported_start=$(( $NUM_BP_PER_WINDOW * $region_index ))
	reported_end=$(( $NUM_BP_PER_WINDOW * ($region_index+1) ))
	while IFS= read -r ct 
	do
		raw_fn=${raw_segmentation_folder}/${ct}/${ct}_25_imputed12marks_segments.bed.gz
		output_fn=${cg_dir}/${ct}_${region_chrom}_${region_index}.bed
		code=/u/home/h/havu73/project-ernst/source_pete/draw_igv_plots/get_segmentation_file_with_color.py
		if [ ! -f ${output_fn}.gz ] # if the file does not exist already then write code to produce it
		then 
			echo "python ${code} ${raw_fn} ${state_annot_fn} ${output_fn} ${cell_group}_${ct} ${region_chrom} ${reported_start} ${reported_end} 
			gzip ${output_fn}" >> ${job_dir}/get_igv_plots
		fi
	done < ${sample_fn}
}
########################
########################
# second we create the segmentation for full stack state annotations
function create_igv_segment_full_stack() {
	code=/u/home/h/havu73/project-ernst/source/chromHMM_utilities/get_segmentation_file_with_color.py
	fS_segment_fn=/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/hg19_segmentations/genome_100_segments.bed.gz
	fS_state_annot_fn=/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/figures/supp_excel/state_annotations_processed.csv
	fS_out_dir=${region_out_dir}/full_stack/
	output_fn=${fS_out_dir}/full_stack.bed
	echo "python ${code} ${fS_segment_fn} ${fS_state_annot_fn} ${output_fn} full_stack ${region_chrom} ${region_start} ${region_end}" >> ${job_dir}/get_igv_plots
	echo "gzip ${output_fn}" >> ${job_dir}/get_igv_plots
}



cg_list='blood Brain Digestive Epithelial ESC Heart iPSC Mesench Muscle Myosat Neuron Skin Thymus'
rm -f ${job_dir}/get_igv_plots
for cg in ${cg_list}
do
	create_igv_segment_one_cell_group $cg
done	
create_igv_segment_full_stack 

########################
########################
# third we print out the link to the files for uploading onto UCSC genome browser
ha_url=https://public.hoffman2.idre.ucla.edu/ernst/2K9RS
upload_ucsc_fn=${job_dir}/upload_link_to_ucsc
rm -f ${upload_ucsc_fn}
for cg in ${cg_list}
do
	cg_dir=${region_out_dir}/${cg}/
	for f in ${cg_dir}/*.bed.gz
	do
		shortF=$(echo $f | awk -F'/' '{print $NF}')
		echo ${ha_url}/full_stack/igv_source/${region_chrom}_${region_index}/${cg}/${shortF} >> ${upload_ucsc_fn}
	done
done
echo ${ha_url}/full_stack/igv_source/${region_chrom}_${region_index}/full_stack/full_stack.bed.gz >> ${upload_ucsc_fn}
chmod +x ${job_dir}/get_igv_plots