roadmap_gene_info_fn=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/gene_exp/gene_locations_reorganized.gz
tss_fn=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/gene_exp/tss.bed
# zcat $roadmap_gene_info_fn | awk -F'\t' '{if($5 == -1){print $1,$3,$4} else {print $1,
# $2,$4}}' | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$3}' > $tss_fn
# gzip $tss_fn
# in $roadmap_gene_info_fn
# $1: chromosome
# $2: smaller coordinate
# $3: bigger coodinate
# $4: gene name
# $5: strand 
# if strand =  1 --> tss is the smaller one
# if strand = -1 --> tss is the bigger one
# 
tss_fn=${tss_fn}.gz

outDir=/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/gene_ontology_analysis  
full_stack_segment_fn=/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/hg19_segmentations/genome_100_segments.bed.gz
tss_segmentation_fn=$outDir/tss_segment.bed
command="bedtools intersect -a $tss_fn -b $full_stack_segment_fn -wb | awk -F'\t' 'BEGIN {OFS=\"\t\"} {print \$1,\$2,\$3,\$4,\$8}' > $tss_segmentation_fn"
# echo $command 
# gzip $tss_segmentation_fn 
tss_segmentation_fn=${tss_segmentation_fn}.gz
# the output of the bedtools command is
# $1, $2, $3: locations of tss
# $4: gene names
# $5,$6,$7: locations of the segments
# $8: state 

most_popular_tss_state=$outDir/most_popular_tss_state
# zcat $tss_segmentation_fn | awk -F'\t' '{print $5}' | sort -k1 | uniq -c  | sort -k1n | awk 'BEGIN {OFS="\t"} {if($1>=100){print $1,$2}}' > $most_popular_tss_state
# the format of $tss_segmentation_fn is
# $1, $2, $3: locations of tss
# $4: gene names
# $5: state

while read line;  
do
	state=$(echo $line | awk '{print $2}')
	this_state_gene_list_fn=$outDir/${state}_gene_list
	zcat $tss_segmentation_fn | awk -F'\t' -v this_state=$state '{if($5==this_state){print $4}}' > $this_state_gene_list_fn
done < $most_popular_tss_state