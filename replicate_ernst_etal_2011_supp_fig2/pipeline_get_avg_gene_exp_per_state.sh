input_data_folder=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/gene_exp/
raw_exp_fn=${input_data_folder}/57epigenomes.RPKM.pc.gz
gene_info_fn=${input_data_folder}/Ensembl_v65.Gencode_v10.ENSG.gene_info.gz
output_data_folder="/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/replicate_supp_fig2_ernst_etal_2011/one_avg_exp_per_state/"

# First, we have to clean the gene_info file. We want to keep only genes that are protein coding. And we create the file in the form of a bed file with the first three columnns being chromosome, start_bp, end_bp. 
# in the gene_info_fn: $6: gene type, $2: chromosome, $3: start_bp, $4: end_bp, $5: strand, $1: ensemble gene id
# zcat ${gene_info_fn} | awk -F'\t' -v pc="protein_coding" '{if ($6==pc) print "chr"$2"\t"$3"\t"$4"\t"$5"\t"$1}' > ${output_data_folder}/step1_output
echo "Done with getting protein coding genes"
echo ""

# Second, give the output file from step1, we will rearrange the file so that all genes are grouped based on the chromosome and arranged in ascending order of start_bp
reorganize_bed_file_code=/u/home/h/havu73/project-ernst/source/chromHMM_utilities/reorganize_bed_file.sh
command="$reorganize_bed_file_code ${output_data_folder}/step1_output ${output_data_folder}/step2_output" 
# $command
# rm -f ${output_data_folder}/step1_output 
# now the output of step2 would be ${output_data_folder}/step2_output.gz
echo "Done with reorganning the bed file for gene locations. Now things are in place"
echo ""

# thrid, we will will find the chromatin state that intersect with each genes
full_stack_segment_fn=/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/hg19_segmentations/genome_100_segments.bed.gz # chromsome, start_bp, end_bp, state (E1 --> E100)
# bedtools intersect -a $full_stack_segment_fn -b ${output_data_folder}/step2_output.gz -wa -wb > ${output_data_folder}/step3_output
# rm -f ${output_data_folder}/step2_output.gz
# gzip ${output_data_folder}/step3_output
echo "Done with intersecting the gene info data with chromatin state data. The output of this step is ${output_data_folder}/step3_output.gz"

# fourth, we will call on a python program that combine the raw gene expression data with the gene info- chromatin state data that we just obtained and calculate the average expression corresponding to each chromatin state
calculate_avg_exp_code='/u/home/h/havu73/project-ernst/source/replicate_ernst_etal_2011_supp_fig2/add_gene_info_to_gene_expression_data_roadmap.py'
num_states=100
python ${calculate_avg_exp_code} $raw_exp_fn ${output_data_folder}/step3_output.gz  $output_data_folder $num_states
echo "" 
echo "Done with calculating the average gene expression in all states in all cell types"
# rm -r ${output_data_folder}/step3_output.gz

# fifth, visualization:
echo "Please look into file /Users/vuthaiha/Desktop/window_hoff/source/replicate_ernst_etal_2011_supp_fig2/plot_avg_gene_exp_heatmap.Rmd for visualization code"