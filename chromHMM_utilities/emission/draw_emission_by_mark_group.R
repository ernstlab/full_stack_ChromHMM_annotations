# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

source("../helper_function.R") # some useful functions: get_mark_ct, get_chrom_mark_name, paste_col_name, and to get the metadata dataframe containing information about different cell types, provided by ROADMAP
library(pheatmap)
CHROM_MARK_COLOR_CODE <- c('class1_acetyl' = '#E5EAE7', 'H3K4me3'= '#F13712', 'H3K4me2'= '#FFA500', 'H3K4me1'= '#EDF732', 'H2A.Z'= '#C3D59C', 'H3K27me3'= '#A6A6A4', 'H3K79me1'= '#9DCEA8', 'H3K79me2'= '#377A45', 'H3K9ac'= '#F2A626', 'H3K36me3'= '#49AC5E', 'H4K20me1'= '#6AD1C8', 'DNase'= '#DBE680', 'H3K9me3'= '#677BF6', 'H3K27ac'= '#F7CB4D')
ANATOMY_COLOR_CODE = c('BLOOD'= '#e34a33', 'ESC_DERIVED'= '#fdbb84', 'ESC'= '#fef0d9', 'BRAIN'= '#ffffd4', 'LUNG'= '#006837', 'SKIN'= '#78c679', 'MUSCLE'= '#c2e699', 'IPSC'= '#7a0177', 'GI_STOMACH' = '#f768a1', 'HEART'= '#fbb4b9', 'BREAST' = '#d7b5d8', 'GI_INTESTINE' = '#253494', 'GI_RECTUM' = "#2c7fb8", 'GI_COLON' = "#a1dab4", 'FAT' = "#54278f", 'LIVER' = '#9e9ac8', 'VASCULAR' = '#f2f0f7', 'PANCREAS'= '#1c9099', 'STROMAL_CONNECTIVE' = '#bdc9e1', 'THYMUS' = "#f6eff7", 'PLACENTA'= '#3182bd', 'GI_DUODENUM' = "#636363", 'CERVIX'= '#cccccc', 'BONE'= "#051C34", 'ADRENAL'= '#55A4F4', 'KIDNEY'= '#F4556C', 'MUSCLE_LEG' = '#9B7A7F', 'OVARY' = '#F2FA10', 'SPLEEN'= '#11FA10', 'GI_ESOPHAGUS' = '#10F7FA', 'NaN'= "#040404")
CELL_GROUP_COLOR_CODE = c('Neurosph'= '#FFD924', 'HSC & B-cell'= '#678C69', 'Mesench'= '#B65C73', 'Brain'= '#C5912B', 'Adipose'= '#AF5B39', 'Thymus'= '#DAB92E', 'Sm. Muscle'= '#F182BC', 'IMR90'= '#E41A1C', 'Myosat'= '#E67326', 'iPSC'= '#69608A', 'Muscle'= '#C2655D', 'Digestive'= '#C58DAA', 'ESC'= '#924965', 'Epithelial'= '#FF9D0C', 'Heart'= '#D56F80', 'ENCODE2012'= '#000000', 'ES-deriv'= '#4178AE', 'Other'= '#999999', 'Blood & T-cell'= '#55A354', 'NA' = 'black')
STATE_COLOR_DICT = c('quescient' = "#ffffff", "HET" = '#b19cd9', 'acetylations' = '#fffacd', 'enhancers'= '#FFA500', 'ct_spec enhancers' = '#FFA500', 'transcription' = '#006400', 'weak transcription' = '#228B22', 'weak enhancers' = '#ffff00', 'transcribed and enhancer' = '#ADFF2F', 'exon' = '#3cb371', 'promoters' = '#ff4500', 'TSS' = '#FF0000', 'weak promoters' = '#800080', 'polycomb repressed' = '#C0C0C0', 'others' = '#fff5ee', 'znf' = '#7fffd4', 'DNase' = '#fff44f')
#############
### FUNCTIONS TO GET DATA FROM CHROMHMM EMISSION FORMAT INTO DATA FRAME
#############
# return emisison_df with attached metadata about the marks. Original data taken from ChromHMM emission matrix format
# mark_group, S<state_index>, ct, chrom_mark, GROUP, COLOR

get_emission_df_by_mark_group <- function(emission_fn){
	emission_df <- as.data.frame(read.table(emission_fn, header = FALSE, sep = "\t"))
	# first line of mark names are also part of the data, not the header
	emission_df <- as.data.frame(t(emission_df)) # transpos--> rows: marks, columns: states
	emission_colnames <- c("mark_group", sapply(emission_df[1,-1], paste_col_name, "S")) # mark_group, S1 --> S100
	emission_df <- emission_df[-1,] # get rid of the first rows which are column names or somethings like that
	colnames(emission_df) <- emission_colnames # assign column names
	emission_df <- emission_df %>% mutate_all(.funs = as.character)%>% mutate_at(vars(starts_with("S")), funs(as.double)) %>% mutate_at(vars("mark_group"), funs(as.character)) 
	# have to change everything to character first because otherwise converting from factor to numeric will convert to indices of levels of factor
	return(emission_df) # rows: marks, columns: states. There are slight variations because we also have to annotate the marks so columns will not just be states.
}

read_state_annot_df <- function(state_annot_fn, states_to_plot){
	annot_df <- as.data.frame(read.csv(state_annot_fn, header = TRUE, stringsAsFactors = FALSE))
	tryCatch({
		annot_df <- annot_df %>% rename('state_type' = 'Group')
		}, error = function (e) {message("tried to change column names in annot_df and nothing worth worrying happened")})
	annot_df <- annot_df %>% filter(state %in% states_to_plot)
	annot_df <- annot_df %>% arrange(state_order_by_group) # order rows based on the index that we get from state_order_by_group column
	return(annot_df)	
}

pheatmap_emission_no_grouping <- function(emission_df, save_fn, states_to_plot){
	uniq_chrom_mark <- c('H3K9me3', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H4K20me1', 'H3K79me2', 'H3K79me1', 'H3K36me3', 'DNase', 'H3K27ac', 'H3K27me3', 'H2A.Z', 'H3K9ac', 'H2AK9ac', 'H4K12ac', 'H2BK20ac', 'H3K56ac', 'H4K5ac', 'K3BK15ac', 'H2BK12ac', 'H3K14ac', 'H4K91ac', 'H2AK5ac', 'H2BK120ac', 'H2BK5ac', 'H3K18ac', 'H3K23ac', 'H3K4ac', 'H2BK15ac', 'H4K8ac', 'class1_acetyl', 'H3T11ph', 'H3K23me2', 'H3K9me1') # put chromatin marks in an order that puts all the acetylation marks together --> we are putting the rarely studied acetylation marks together, and some other less well studied marks together #sort(unique(emission_df$chrom_mark)) # list of unique chromatin marks
	# c("ENCODE2012", "ES-deriv", "ESC", "Blood & T-cell", "Digestive", "Other", "Brain", "HSC & B-cell", "Epithelial", "iPSC", "Muscle", "Heart", "IMR90", "Sm. Muscle", "Mesench", "Thymus", "Neurosph", "Adipose", "Myosat")
	plot_df <- data.frame(matrix(ncol = ncol(emission_df), nrow = 0))
	colnames(plot_df) <- colnames(emission_df)
	for (chrM in uniq_chrom_mark){
	  this_chrM_df <- emission_df %>% filter(mark_group == chrM)
	  plot_df <- bind_rows(plot_df, this_chrM_df)
	}
	emission_dist <- plot_df %>% select(paste0('S', states_to_plot))
	rownames(emission_dist) <- plot_df$mark_group # add rownames of marks, so that the pheatmap function can create plots with the chrom mark names on it
	chrom_mark_df <- data.frame(chrom_mark = plot_df$mark_group)
	chrom_mark_df['chrom_mark'] <- as.character(chrom_mark_df$chrom_mark) # conver to character
	rownames(chrom_mark_df) <- as.character(plot_df$mark_group)
	emission_value_ranges <- seq(0,1, 0.01)
	hm_no_grouping <- pheatmap(t(emission_dist), breaks = emission_value_ranges, fontsize = 5, annotation_col = chrom_mark_df, annotation_colors = list(chrom_mark = CHROM_MARK_COLOR_CODE), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, fontsize_col = 3, angle_col = 90 , cellheight = 10, filename = paste0(save_fn))
	return (hm_no_grouping)
}

calculate_gap_rows_among_state_groups <- function(state_annot_df){ 
	state_group_ordered_by_appearance <- unique(state_annot_df$state_type) # list of different state groups, ordered by how they appear in the heatmap from top to bottom
	count_df <- state_annot_df %>% count(state_type)
	count_df <- count_df[match(state_group_ordered_by_appearance, count_df$state_type),] # order the rows such that the state_type are ordered based on state_group_ordered_by_appearance
	results <- cumsum(count_df$n) # cumulative sum of the count for groups of states, which will be used to generate the gaps between rows of the heatmaps
	return(results)
}

pheatmap_emission_annot_state_group <- function(emission_df, save_fn, state_annot_fn, states_to_plot){
	uniq_chrom_mark <- c('H3K9me3', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H4K20me1', 'H3K79me2', 'H3K79me1', 'H3K36me3', 'DNase', 'H3K27ac', 'H3K27me3', 'H2A.Z', 'H3K9ac', 'H2AK9ac', 'H4K12ac', 'H2BK20ac', 'H3K56ac', 'H4K5ac', 'K3BK15ac', 'H2BK12ac', 'H3K14ac', 'H4K91ac', 'H2AK5ac', 'H2BK120ac', 'H2BK5ac', 'H3K18ac', 'H3K23ac', 'H3K4ac', 'H2BK15ac', 'H4K8ac', 'class1_acetyl', 'H3T11ph', 'H3K23me2', 'H3K9me1') # put chromatin marks in an order that puts all the acetylation marks together --> we are putting the rarely studied acetylation marks together, and some other less well studied marks together #sort(unique(emission_df$chrom_mark)) # list of unique chromatin marks
	plot_df <- data.frame(matrix(ncol = ncol(emission_df), nrow = 0))
	colnames(plot_df) <- colnames(emission_df) # rows marks columns states and other annotation of the marks
	for (chrM in uniq_chrom_mark){ # rearrange such that all expeirments of the same marks are put next to each other
	  this_chrM_df <- emission_df %>% filter(mark_group == chrM)
	  plot_df <- bind_rows(plot_df, this_chrM_df) # rearrange the marks
	}
	state_annot_df <- read_state_annot_df(state_annot_fn, states_to_plot)
	states_columns_to_select <- paste0('S', state_annot_df$state) # list of states that have been rearranged. This will be used as column names that we will select from states_to_select
	emission_dist <- plot_df %>% select(states_columns_to_select) # only choose the states that are specified by users, and also ordered based on the state group that they belong
	colnames(emission_dist) <- state_annot_df$mneumonics
	rownames(emission_dist) <- plot_df$mark_group # add rownames of marks, so that the pheatmap function can create plots with the chrom mark names on it
	# create the df for annotation of rows in our heatmap
	rownames(state_annot_df) <- state_annot_df$mneumonics
	state_annot_df <- state_annot_df %>% select(state_type)
	# create the df for annotation of columns in our heatmap
	chrom_mark_df <- data.frame(chrom_mark = plot_df$mark_group)
	chrom_mark_df['chrom_mark'] <- as.character(chrom_mark_df$chrom_mark) # conver to character
	rownames(chrom_mark_df) <- as.character(plot_df$mark_group)
	# calculate the indices of rows where we want to create gaps among groups of states
	gap_row_indices <- calculate_gap_rows_among_state_groups(state_annot_df)
	# cell_group_df <- data.frame(cell_group = plot_df$GROUP)
	emission_value_ranges <- seq(0,1, 0.01)
	hm_no_grouping <- pheatmap(t(emission_dist), breaks = emission_value_ranges, fontsize = 5, annotation_col = chrom_mark_df, annotation_row = state_annot_df, annotation_colors = list(chrom_mark = CHROM_MARK_COLOR_CODE, state_type = STATE_COLOR_DICT), gaps_row = gap_row_indices, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, fontsize_col = 3, angle_col = 90 , cellheight = 10, filename = paste0(save_fn))
	return (hm_no_grouping)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3)
{
	stop("wrong command line argument format", call.=FALSE)
}
emission_fn <- args[1] # output of calculate_emission_by_group.py , also input to this program, should be ./output_emission/emissions_by_mark_group.txt
state_annot_fn <- args[2] # where annotations of states, the state group and state index by groups are stored
save_fn <- args[3] # where the figures should be stored. 
states_to_plot <- as.integer(args[4:length(args)])
if (length(args) == 3){
	states_to_plot <- seq(1, 100)
} # if the user does not specify states to plots, then we would assume that we have to plot everything
print(states_to_plot)
emission_df <- get_emission_df_by_mark_group(emission_fn) 
print("Done getting emission_df")
pheatmap_emission_annot_state_group(emission_df, save_fn, state_annot_fn, states_to_plot) 