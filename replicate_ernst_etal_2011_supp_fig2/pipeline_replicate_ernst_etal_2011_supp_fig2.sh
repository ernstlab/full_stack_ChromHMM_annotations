roadmap_gene_exp_folder=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/gene_exp
# first, based on the gene info, we will extract the 50kb region around the TSS of each gene, along with the gene's expressions in different cell types.
# after this step, I confirm that all the genes in the 57RPKM file are protein-coding genes. Checked!
pypy get_50kb_around_TSS_roadmap_gene_info.py $roadmap_gene_exp_folder/Ensembl_v65.Gencode_v10.ENSG.gene_info.gz $roadmap_gene_exp_folder/gene_exp_TSS_region_Ensembl_v65.Gencode_v10.ENSG_57roadmap.gz $roadmap_gene_exp_folder/57epigenomes.RPKM.pc.gz ## wrote into output_fn a table with the following headers: ['chr', 'start_tss_region', 'end_tss_region', 'gene_ENS_id', 'strand', 'tss_bp'] + celltype_list. Strand is -1 --> tss is the max coordinate, strand is 1 --> tss is the min coordinate
echo 'Done extracting the TSS regions from gene_info file'
echo ''
echo ''
echo ''
# second, we will reorganize the resulting TSS regions, so that genes inside the same chromosome are put together, and ordered based on the starting coordinate, ascendingly
# /u/home/h/havu73/project-ernst/source/chromHMM_utilities/reorganize_bed_file.sh $roadmap_gene_exp_folder/gene_exp_TSS_region_Ensembl_v65.Gencode_v10.ENSG_57roadmap.gz $roadmap_gene_exp_folder/gene_exp_TSS_region_Ensembl_v65.Gencode_v10.ENSG_57roadmap_ordered
echo 'Done reorganizing the TSS region files so that all genes of the same chromosome are written one after another, ordered based on starting bp ascendingly'
echo ''
echo ''
echo ''
# third, we interesct the TSS region (with data of gene expression), with the segmentation data
hg19_fullstack_segmentation_fn=/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/hg19_segmentations/genome_100_segments.bed.gz
gene_exp_tss_folder=/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/replicate_supp_fig2_ernst_etal_2011/
mkdir -p $gene_exp_tss_folder
gene_exp_segmentation_output_fn=$gene_exp_tss_folder/gene_exp_TSS_region_Ensembl_v65.Gencode_v10.ENSG_57roadmap_segmentation
rm -f $gene_exp_segmentation_output_fn # so that we can create a new one
# bedtools intersect -a $hg19_fullstack_segmentation_fn -b $roadmap_gene_exp_folder/gene_exp_TSS_region_Ensembl_v65.Gencode_v10.ENSG_57roadmap_ordered.gz -wa -wb > $gene_exp_segmentation_output_fn
# gzip -f $gene_exp_segmentation_output_fn
echo 'Done combining the data of TSS regions (with gene expression data) and segmentation data from chromHMM_model'
echo ''
echo ''
echo ''

# forth, get rid of files that are temporary necessary to create the file fo gene expressions and segmentation that we want. We can always recreate these files if we need to
# rm -f ${roadmap_gene_exp_folder}/gene_exp_TSS_region_Ensembl_v65.Gencode_v10.ENSG_57roadmap.gz
# echo Removed $roadmap_gene_exp_folder/gene_exp_TSS_region_Ensembl_v65.Gencode_v10.ENSG_57roadmap.gz
# rm -f ${roadmap_gene_exp_folder}/gene_exp_TSS_region_Ensembl_v65.Gencode_v10.ENSG_57roadmap_ordered.gz
# echo Removed ${roadmap_gene_exp_folder}/gene_exp_TSS_region_Ensembl_v65.Gencode_v10.ENSG_57roadmap_ordered.gz
echo ''
echo ''
echo ''

# fifth, we will stratifiy the data from previous step into states. Previously, we had data of genes's 50kb regions, and the state segmentation segments overlapping those regions. Note: multiple genes can have their 50kb-surrounding-tss region overlapping. So the same segments of one state and overlaping with that region of multiple genes. We now put data of each state into its own file. Right now we comment this out because this step has been done adn takes a lot of time.  We do this stratification because it's beneficial for the next step processed in python

get_avg_gene_exp_code='/u/project/ernst/havu73/source/replicate_ernst_etal_2011_supp_fig2/get_gene_exp_by_state_by_distance_from_TSS.py'
gene_exp_tss_folder=/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/replicate_supp_fig2_ernst_etal_2011
output_folder=${gene_exp_tss_folder}/ct_avg_gExp_tss_relative/
ct_list_fn=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/gene_exp/header_57epigenomes_rpkm_protein_coding
num_processes=4
NUM_STATE=100

# for state_index in `seq 1 $NUM_STATE`
# do
# 	this_state_gene_exp_fn=${gene_exp_tss_folder}/state_${state_index}_gene_exp
# 	zcat ${gene_exp_segmentation_output_fn}.gz | awk -F'\t' -v state_symbol="E$state_index" '{if ($4==state_symbol) print $0;}' > $this_state_gene_exp_fn
# 	gzip -f $this_state_gene_exp_fn
# 	echo "Done filtering out state: $state_index"
# done
# python $get_avg_gene_exp_code ${gene_exp_tss_folder} ${output_folder} ${ct_list_fn} $NUM_STATE $num_processes

# sixth, calculate the average of processed data across different cell types
calculate_over_ct_gExp_tss_relative_code=/u/home/h/havu73/project-ernst/source/replicate_ernst_etal_2011_supp_fig2/average_over_ct_gExp_tss_relative.py
# python  $calculate_over_ct_gExp_tss_relative_code ${output_folder} ${output_folder}/avg_across_ct_avg_exp_tss_relative.gz