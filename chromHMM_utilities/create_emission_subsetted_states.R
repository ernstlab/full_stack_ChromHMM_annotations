source("/Users/vuthaiha/Desktop/window_hoff/source/summary_analysis/common_functions.R")
library(heatmaply)
library(pheatmap)
CHROM_MARK_COLOR_CODE <- c('class1_acetyl' = '#E5EAE7', 'H2BK12ac'= '#E5EAE7', 'H2AK5ac'= '#E5EAE7', 'H2BK20ac'= '#E5EAE7', 'H3K27ac'= '#F7CB4D', 'H2BK5ac'= '#E5EAE7', 'H4K12ac'= '#E5EAE7', 'H4K5ac'= '#E5EAE7', 'H3K4me3'= '#F13712', 'H3K4me2'= '#FFA500', 'H3K4me1'= '#EDF732', 'H2A.Z'= '#C3D59C', 'H3K9me1'= '#F3AEAE', 'H3K23ac'= '#E5EAE7', 'H4K91ac'= '#E5EAE7', 'H3K4ac'= '#E5EAE7', 'H3T11ph'= '#F3AEAE', 'H3K23me2'= '#F3AEAE', 'H2BK120ac'= '#E5EAE7', 'H3K27me3'= '#A6A6A4', 'H3K79me1'= '#9DCEA8', 'H3K79me2'= '#377A45', 'H3K9ac'= '#F2A626', 'H3K18ac'= '#E5EAE7', 'K3BK15ac'= '#E5EAE7', 'H3K36me3'= '#49AC5E', 'H4K20me1'= '#6AD1C8', 'H3K56ac'= '#E5EAE7', 'DNase'= '#DBE680', 'H3K9me3'= '#677BF6', 'H2AK9ac'= '#E5EAE7', 'H3K14ac'= '#E5EAE7', 'H4K8ac'= '#E5EAE7', 'NA' = 'black', 'H2BK15ac' = '#E5EAE7', 'others' = '#F3AEAE')
ANATOMY_COLOR_CODE = c('BLOOD'= '#e34a33', 'ESC_DERIVED'= '#fdbb84', 'ESC'= '#fef0d9', 'BRAIN'= '#ffffd4', 'LUNG'= '#006837', 'SKIN'= '#78c679', 'MUSCLE'= '#c2e699', 'IPSC'= '#7a0177', 'GI_STOMACH' = '#f768a1', 'HEART'= '#fbb4b9', 'BREAST' = '#d7b5d8', 'GI_INTESTINE' = '#253494', 'GI_RECTUM' = "#2c7fb8", 'GI_COLON' = "#a1dab4", 'FAT' = "#54278f", 'LIVER' = '#9e9ac8', 'VASCULAR' = '#f2f0f7', 'PANCREAS'= '#1c9099', 'STROMAL_CONNECTIVE' = '#bdc9e1', 'THYMUS' = "#f6eff7", 'PLACENTA'= '#3182bd', 'GI_DUODENUM' = "#636363", 'CERVIX'= '#cccccc', 'BONE'= "#051C34", 'ADRENAL'= '#55A4F4', 'KIDNEY'= '#F4556C', 'MUSCLE_LEG' = '#9B7A7F', 'OVARY' = '#F2FA10', 'SPLEEN'= '#11FA10', 'GI_ESOPHAGUS' = '#10F7FA', 'NaN'= "#040404")
CELL_GROUP_COLOR_CODE = c('Neurosph'= '#FFD924', 'HSC & B-cell'= '#678C69', 'Mesench'= '#B65C73', 'Brain'= '#C5912B', 'Adipose'= '#AF5B39', 'Thymus'= '#DAB92E', 'Sm. Muscle'= '#F182BC', 'IMR90'= '#E41A1C', 'Myosat'= '#E67326', 'iPSC'= '#69608A', 'Muscle'= '#C2655D', 'Digestive'= '#C58DAA', 'ESC'= '#924965', 'Epithelial'= '#FF9D0C', 'Heart'= '#D56F80', 'ENCODE2012'= '#000000', 'ES-deriv'= '#4178AE', 'Other'= '#999999', 'Blood & T-cell'= '#55A354', 'NA' = 'black')
STATE_COLOR_DICT = c('quescient' = "#ffffff", "HET" = '#b19cd9', 'acetylations' = '#fffacd', 'enhancers'= '#FFA500', 'ct_spec enhancers' = '#FFA500', 'transcription' = '#006400', 'weak transcription' = '#228B22', 'weak enhancers' = '#ffff00', 'transcribed and enhancer' = '#ADFF2F', 'exon' = '#3cb371', 'promoters' = '#ff4500', 'TSS' = '#FF0000', 'weak promoters' = '#800080', 'polycomb repressed' = '#C0C0C0', 'others' = '#fff5ee', 'znf' = '#7fffd4', 'DNase' = '#fff44f')

get_rid_of_bedGZ_context_name <- function(enrichment_name){
	if (endsWith(enrichment_name, '.bed.gz')){
		result <- substring(enrichment_name, 1, nchar(enrichment_name) - 7)
		return(result)
	}
	return(enrichment_name)
}

read_state_annot_df <- function(state_annot_fn, states_to_plot){
	annot_df <- as.data.frame(read.csv(state_annot_fn, header = TRUE, stringsAsFactors = FALSE))
	tryCatch({
		annot_df <- annot_df %>% rename('state_type' = 'Group')
		}, error = function (e) {message("tried to change column names in annot_df and nothing worth worrying happened")})
	annot_df <- annot_df %>% slice(states_to_plot)
	# annot_df <- annot_df %>% arrange(state_order_by_group) # order rows based on the index that we get from state_order_by_group column
	return(annot_df)	
}

#############
### FUNCTIONS TO GET DATA FROM CHROMHMM EMISSION FORMAT INTO DATA FRAME
#############
# return emisison_df with attached metadata about the marks. Original data taken from ChromHMM emission matrix format
# mark_names, S<state_index>, ct, chrom_mark, GROUP, COLOR
get_emission_df_from_chromHMM_emission <- function(emission_fn){
	emission_df <- as.data.frame(read.table(emission_fn, header = FALSE, sep = "\t"))
	# first line of mark names are also part of the data, not the header
	emission_df <- as.data.frame(t(emission_df)) # transpos--> rows: marks, columns: states
	emission_colnames <- c("mark_names", sapply(emission_df[1,-1], paste_col_name, "S"))
	emission_df <- emission_df[-1,] # get rid of the first rows which are column names or somethings like that
	colnames(emission_df) <- emission_colnames # assign column names
	emission_df <- emission_df %>% mutate_all(.funs = as.character)%>% mutate_at(vars(starts_with("S")), funs(as.double)) %>% mutate_at(vars("mark_names"), funs(as.character)) 
	# have to change everything to character first because otherwise converting from factor to numeric will convert to indices of levels of factor
	emission_df["ct"] <- apply(emission_df['mark_names'], FUN = get_mark_ct, MARGIN = 1)
	emission_df["chrom_mark"] <- apply(emission_df['mark_names'], FUN = get_chrom_mark_name, MARGIN = 1)
	emission_df <- emission_df %>% dplyr::left_join(metadata, by = c("ct" = "CT_NAME"))
# combine data and group into cell groups and calculate mean of emission from states for each cell groups
	return(emission_df)
}

order_experiments_within_one_chrom_mark <- function(emission_df, chrM){
	# uniq_group is the list of unique cell groups that we will sort experiments within each chrom mark into
	this_chrM_df <- emission_df %>% filter(chrom_mark == chrM)
	uniq_group <- sort(unique(this_chrM_df$GROUP))
	final_chrM_df <- data.frame(matrix(ncol = ncol(emission_df), nrow = 0))# the final results of ordering experiments associated with this mark.
	for (group in uniq_group){
		this_group_df <- this_chrM_df %>% filter(GROUP == group)
		final_chrM_df <- bind_rows(final_chrM_df, this_group_df)
	}
	return(final_chrM_df)
}


pheatmap_emission_no_grouping <- function(emission_df, save_fn, state_annot_fn, states_to_plot){
	uniq_chrom_mark <- c('H3K9me3', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H4K20me1', 'H3K79me2', 'H3K79me1', 'H3K36me3', 'DNase', 'H3K27ac', 'H3K27me3', 'H2A.Z', 'H3K9ac', 'H2AK9ac', 'H4K12ac', 'H2BK20ac', 'H3K56ac', 'H4K5ac', 'K3BK15ac', 'H2BK12ac', 'H3K14ac', 'H4K91ac', 'H2AK5ac', 'H2BK120ac', 'H2BK5ac', 'H3K18ac', 'H3K23ac', 'H3K4ac', 'H2BK15ac', 'H4K8ac', 'class1_acetyl', 'H3T11ph', 'H3K23me2', 'H3K9me1', 'others') # put chromatin marks in an order that puts all the acetylation marks together --> we are putting the rarely studied acetylation marks together, and some other less well studied marks together #sort(unique(emission_df$chrom_mark)) # list of unique chromatin marks
	uniq_group <- sort(unique(emission_df$GROUP)) # list of unique tissue types
	# print(uniq_group)
	plot_df <- data.frame(matrix(ncol = ncol(emission_df), nrow = 0))
	colnames(plot_df) <- colnames(emission_df)
	# for (anaM in uniq_group){
	# 	this_anaM_df <- emission_df %>% filter(GROUP == anaM)
	# 	plot_df <- bind_rows(plot_df, this_anaM_df)
	# }
	# random_mark_order <- sample(seq(1, nrow(emission_df)))
	# plot_df <- emission_df %>% slice(random_mark_order)
	for (chrM in uniq_chrom_mark){
	  this_chrM_df <- order_experiments_within_one_chrom_mark(emission_df, chrM) # get experiments associated with this chromatin marks, and within this mark, arrange the experiments based on the associated cell groups
	  plot_df <- bind_rows(plot_df, this_chrM_df)
	}	
	cell_group_data <- plot_df$GROUP
	#plot_df <- plot_df %>% select(-ct, -chrom_mark, - GROUP, -COLOR, -TYPE, -Epig_name, -ANATOMY)
	state_annot_df <- read_state_annot_df(state_annot_fn, states_to_plot)
	emission_dist <- plot_df %>% select(paste0('S', state_annot_df$state))
	rownames(emission_dist) <- plot_df$mark_names # add rownames of marks, so that the pheatmap function can create plots with the chrom mark names on it
	colnames(emission_dist) <- state_annot_df$mneumonics
	chrom_mark_df <- data.frame(chrom_mark = apply(plot_df['mark_names'], FUN = get_chrom_mark_name, MARGIN = 1))
	chrom_mark_df['chrom_mark'] <- as.character(chrom_mark_df$chrom_mark)
	chrom_mark_df['GROUP'] <- as.character(cell_group_data)
	rownames(chrom_mark_df) <- as.character(plot_df$mark_names)
	# create the df for annotationof rows in our heatmap
	rownames(state_annot_df) <- state_annot_df$mneumonics
	state_annot_df <- state_annot_df %>% select(state_type)
	# cell_group_df <- data.frame(cell_group = plot_df$GROUP)
	emission_value_ranges <- seq(0,1, 0.01)
	hm_no_grouping <- pheatmap(t(emission_dist), breaks = emission_value_ranges, fontsize = 5, annotation_col = chrom_mark_df, annotation_row = state_annot_df, annotation_colors = list(chrom_mark = CHROM_MARK_COLOR_CODE, GROUP = CELL_GROUP_COLOR_CODE, state_type = STATE_COLOR_DICT), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, fontsize_col = 3, angle_col = 90 , cellheight = 5, filename = save_fn)
	return (hm_no_grouping)
}


args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3)
{
	stop("wrong command line argument format", call.=FALSE)
}
emission_fn <- args[1] # output of ChromHMM OverlapEnrichment, also input to this program
save_fn <- args[2] # where the figures should be stored
state_annot_fn <- args[3] # where annotations of states, the state group and state index by groups are stored
states_to_plot <- as.integer(args[4:length(args)])
print(states_to_plot)
emission_df <- get_emission_df_from_chromHMM_emission(emission_fn) 
pheatmap_emission_no_grouping(emission_df, save_fn, state_annot_fn, states_to_plot)