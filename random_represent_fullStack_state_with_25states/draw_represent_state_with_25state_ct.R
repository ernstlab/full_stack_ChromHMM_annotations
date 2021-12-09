source("/Users/vuthaiha/Desktop/window_hoff/source/summary_analysis/common_functions.R")
library("tidyverse")
library("tidyr")
library("dplyr")
library('ggplot2')
library(pheatmap)
CHROM_MARK_COLOR_CODE <- c('class1_acetyl' = '#E5EAE7', 'H3K4me3'= '#F13712', 'H3K4me2'= '#FFA500', 'H3K4me1'= '#EDF732', 'H2A.Z'= '#C3D59C', 'H3K27me3'= '#A6A6A4', 'H3K79me1'= '#9DCEA8', 'H3K79me2'= '#377A45', 'H3K9ac'= '#F2A626', 'H3K36me3'= '#49AC5E', 'H4K20me1'= '#6AD1C8', 'DNase'= '#DBE680', 'H3K9me3'= '#677BF6', 'H3K27ac'= '#F7CB4D')
ANATOMY_COLOR_CODE = c('BLOOD'= '#e34a33', 'ESC_DERIVED'= '#fdbb84', 'ESC'= '#fef0d9', 'BRAIN'= '#ffffd4', 'LUNG'= '#006837', 'SKIN'= '#78c679', 'MUSCLE'= '#c2e699', 'IPSC'= '#7a0177', 'GI_STOMACH' = '#f768a1', 'HEART'= '#fbb4b9', 'BREAST' = '#d7b5d8', 'GI_INTESTINE' = '#253494', 'GI_RECTUM' = "#2c7fb8", 'GI_COLON' = "#a1dab4", 'FAT' = "#54278f", 'LIVER' = '#9e9ac8', 'VASCULAR' = '#f2f0f7', 'PANCREAS'= '#1c9099', 'STROMAL_CONNECTIVE' = '#bdc9e1', 'THYMUS' = "#f6eff7", 'PLACENTA'= '#3182bd', 'GI_DUODENUM' = "#636363", 'CERVIX'= '#cccccc', 'BONE'= "#051C34", 'ADRENAL'= '#55A4F4', 'KIDNEY'= '#F4556C', 'MUSCLE_LEG' = '#9B7A7F', 'OVARY' = '#F2FA10', 'SPLEEN'= '#11FA10', 'GI_ESOPHAGUS' = '#10F7FA', 'NaN'= "#040404")
CELL_GROUP_COLOR_CODE = c('Neurosph'= '#FFD924', 'HSC & B-cell'= '#678C69', 'Mesench'= '#B65C73', 'Brain'= '#C5912B', 'Adipose'= '#AF5B39', 'Thymus'= '#DAB92E', 'Sm. Muscle'= '#F182BC', 'IMR90'= '#E41A1C', 'Myosat'= '#E67326', 'iPSC'= '#69608A', 'Muscle'= '#C2655D', 'Digestive'= '#C58DAA', 'ESC'= '#924965', 'Epithelial'= '#FF9D0C', 'Heart'= '#D56F80', 'ENCODE2012'= '#000000', 'ES-deriv'= '#4178AE', 'Other'= '#999999', 'Blood & T-cell'= '#55A354', 'NA' = 'black')
STATE_COLOR_DICT = c('quescient' = "#ffffff", "HET" = '#b19cd9', 'acetylations' = '#fffacd', 'enhancers'= '#FFA500', 'ct_spec enhancers' = '#FFA500', 'transcription' = '#006400', 'weak transcription' = '#228B22', 'weak enhancers' = '#ffff00', 'transcribed and enhancer' = '#ADFF2F', 'exon' = '#3cb371', 'promoters' = '#ff4500', 'TSS' = '#FF0000', 'weak promoters' = '#800080', 'polycomb repressed' = '#C0C0C0', 'others' = '#fff5ee', 'znf' = '#7fffd4', 'DNase' = '#fff44f')

convert_rgb_to_hex <- function(rgb_string){ # convert from 255,0,0 --> FF0000
	rgb_list <- strsplit(rgb_string, ',', fixed = TRUE)
	rgb_list <- as.integer(rgb_list[[1]])
	hex <- rgb(rgb_list[1],  rgb_list[2], rgb_list[3], maxColorValue = 255)
	return(hex)
}

get_25state_color_vector <- function(){ # process annotation of 25-state
	state25_annot_fn <- '~/Desktop/window_hoff/data/roadmap_epigenome/25_imputed12marks_model_downloaded/state_annot.txt'
	df <- as.data.frame(read.csv(state25_annot_fn, header = TRUE, stringsAsFactors = FALSE, sep = '\t'))
	df$hex_color <- sapply(df$rgb, convert_rgb_to_hex)
	return(df$hex_color) # a vector that will be passed into the pheatmap function in the color parameters.
}

read_full_stack_state_annot_df <- function(fS_state_annot_fn, states_to_plot){ # process annotation of fullstack states
	annot_df <- as.data.frame(read.csv(fS_state_annot_fn, header = TRUE, stringsAsFactors = FALSE))
	tryCatch({
		annot_df <- annot_df %>% rename('state_type' = 'Group')
		}, error = function (e) {message("tried to change column names in annot_df and nothing worth worrying happened")})
	annot_df <- annot_df %>% filter(state %in% states_to_plot)
	annot_df <- annot_df %>% arrange(state_order_by_group) # order rows based on the index that we get from state_order_by_group column
	return(annot_df)	
}

read_count_df <- function(count_fn){
	df <- as.data.frame(read.table(count_fn, header = TRUE, sep = '\t'))
	df %>% mutate_at(vars('full_stack_state'), funs(as.character)) # change the full_stack_state from factor to character
	return(df)
}

calculate_gap_rows_among_state_groups <- function(state_annot_df){ 
	state_group_ordered_by_appearance <- unique(state_annot_df$state_type) # list of different state groups, ordered by how they appear in the heatmap from top to bottom
	count_df <- state_annot_df %>% count(state_type)
	count_df <- count_df[match(state_group_ordered_by_appearance, count_df$state_type),] # order the rows such that the state_type are ordered based on state_group_ordered_by_appearance
	results <- cumsum(count_df$n) # cumulative sum of the count for groups of states, which will be used to generate the gaps between rows of the heatmaps
	return(results)
}

get_column_annot_df_by_ct <- function(colnames_count_df){
	ct_column_df <- data.frame(ct_segment = colnames_count_df)
	ct_column_df['ct_segment'] <- as.character(ct_column_df$ct_segment) # conver to character
	ct_column_df['ct'] <- sapply(ct_column_df$ct_segment, function(x){unlist(strsplit(x, '_'))[1]})
	ct_column_df <- ct_column_df %>% dplyr::left_join(metadata, by = c("ct" = "CT_NAME")) # combine data and group into cell groups 
	rownames(ct_column_df) <- as.character(ct_column_df$ct_segment)
	ct_column_df <- ct_column_df %>% select(c('GROUP'))
	return(ct_column_df)	
}

get_column_annot_df_by_group <- function(colnames_count_df){
	ct_column_df <- data.frame(ct_segment = colnames_count_df)
	ct_column_df['ct_segment'] <- as.character(ct_column_df$ct_segment) # conver to character
	ct_column_df['GROUP'] <- sapply(ct_column_df$ct_segment, function(x){unlist(strsplit(x, '_'))[1]})
	ct_column_df$GROUP <- recode(ct_column_df$GROUP, 'Blood...T.cell' = 'Blood & T-cell', 'HSC...B.cell' = 'HSC & B-cell', 'Sm..Muscle' = 'Sm. Muscle', 'ES.deriv' = 'ES-deriv')
	ct_column_df['GROUP'] <- as.character(ct_column_df$GROUP)
	rownames(ct_column_df) <- as.character(ct_column_df$ct_segment)
	ct_column_df <- ct_column_df %>% select(c('GROUP'))
	return(ct_column_df)
}

pheatmap_sample_ct_state_to_represent_full_stack_state <- function(count_sample_region_fn, save_fn, fS_state_annot_fn, states_to_plot, by_ct_or_by_group){
	# by_ct_or_by_group: "ct" or 'group'
	count_df <- read_count_df(count_sample_region_fn) # rows: full_stack states, columns: cell types and the number of sample segments that we got for them
	state_annot_df <- read_full_stack_state_annot_df(fS_state_annot_fn, states_to_plot) # here state_annot_df refers to the annotation of full stack states
	states_rows_to_select <- paste0('E', state_annot_df$state) # list of states that have been rearranged. This will be used as column names that we will select from states_to_select
	count_df <- count_df %>% slice(match(states_rows_to_select, full_stack_state)) # only choose the states that are specified by users, and also ordered based on the state group that they belong
	print(head(count_df%>%select(colnames(count_df)[c(1,2,3)])))
	count_df <- count_df %>% select(- c('full_stack_state'))
	rownames(count_df) <- state_annot_df$mneumonics # add rownames of states, so that the pheatmap function can create plots with the full_stack states on it
	# create the df for annotation of rows in our heatmap
	rownames(state_annot_df) <- state_annot_df$mneumonics
	state_annot_df <- state_annot_df %>% select(state_type)
	# create the df for annotation of columns in our heatmap
	if (by_ct_or_by_group == "ct"){
		ct_column_df <- get_column_annot_df_by_ct(colnames(count_df))
	}
	else { # group
		ct_column_df <- get_column_annot_df_by_group(colnames(count_df))
	}
	# right now columns are of the form <ct>_<index of the sample segment>
	# calculate the indices of rows where we want to create gaps among groups of states
	gap_row_indices <- calculate_gap_rows_among_state_groups(state_annot_df)
	NUM_CT_SPEC_MODEL_STATE <- 25
	state_value_ranges <- seq(1,NUM_CT_SPEC_MODEL_STATE,1)
	color_vector_25state <- get_25state_color_vector()
	hm_no_grouping <- pheatmap(count_df, color = color_vector_25state, fontsize = 5, annotation_col = ct_column_df, annotation_row = state_annot_df, annotation_colors = list(GROUP = CELL_GROUP_COLOR_CODE, state_type = STATE_COLOR_DICT), gaps_row = gap_row_indices, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, fontsize_col = 3, angle_col = 90 , cellheight = 10, cellwidth = 0.1, filename = paste0(save_fn))
	return (hm_no_grouping)
}


args = commandArgs(trailingOnly=TRUE)
NUM_MANDATORY_ARGS <- 4
if (length(args) < 3)
{
	stop("wrong command line argument format", call.=FALSE)
}
count_sample_region_fn <- args[1] # output of count_sample_region_per_full_stack_state.py, also input to this program
fS_state_annot_fn <- args[2] # where annotations of full-stack states, the state group and state index by groups are stored
save_fn <- args[3] # where the figures should be stored
states_to_plot <- as.integer(args[NUM_MANDATORY_ARGS:length(args)])
by_ct_or_by_group <- args[NUM_MANDATORY_ARGS]
if (by_ct_or_by_group != 'ct' && by_ct_or_by_group != 'group'){
	stop('argument by_ct_or_by_group can only take 2 values: ct or group', call. = FALSE)
}
if (length(args) == NUM_MANDATORY_ARGS){
	states_to_plot <- seq(1, 100)
} # if the user does not specify states to plots, then we would assume that we have to plot everything
print(states_to_plot)
pheatmap_sample_ct_state_to_represent_full_stack_state(count_sample_region_fn, save_fn, fS_state_annot_fn, states_to_plot,by_ct_or_by_group)
print("Done")