download_folder_15_state='https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/'
download_folder_18_state='https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/'
download_folder_25_state="https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/"


download_link_suffix_15='_15_coreMarks_segments.bed'
download_link_suffix_18='_18_core_K27ac_segments.bed'
download_link_suffix_25="_25_imputed12marks_segments.bed.gz"

output_folder_15_state='/u/home/h/havu73/project-ernst/data/roadmap_epigenome/15_core_model_downloaded' # TO BE UPDATED: change this folder to the relative path where you want to store the 15-state concatenated annotations for 127 reference epigenomes. 
output_folder_18_state=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded # TO BE UPDATED: change this folder to the relative path where you want to store the 18-state concatenated annotations for 98 reference epigenomes. 
output_folder_25_state='/u/home/h/havu73/project-ernst/data/roadmap_epigenome/25_imputed12marks_model_downloaded' # TO BE UPDATED: change this folder to the relative path where you want to store the 25-state concatenated annotations for 127 reference epigenomes. 
cell_type_file="/u/home/h/havu73/project-ernst/data/roadmap_epigenome/all_127_celltype" # TO BE UPDATED: change this folder to the relative path of the file where we list the biosample code of the 127 reference epigenomes. This can be extracted from the Roadmap metadata (https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15). NOTE: we also provide this file in the github page associated with this project.
while read ct
do
	this_ct_link=${download_folder_25_state}/${ct}${download_link_suffix_25}
	this_ct_output_folder=${output_folder_25_state}/${ct}
	mkdir -p ${this_ct_output_folder}
	output_segment_Fn=${this_ct_output_folder}/${ct}${download_link_suffix_25}
	if [ ! -f  ${output_segment_Fn}.gz ];then
		wget -P ${this_ct_output_folder} $this_ct_link # download the data into output folder
		gzip $output_segment_Fn
		echo "Done downloading for cell type: $ct"
    fi 
done < $cell_type_file
