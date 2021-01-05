# this program will take as input the ovelap enrichment output from ChromHMM and draw a figure that show such data. The resulting plot will be red-based, and each column will have its own color scale and enrihcment values will also be denoted in the plot.
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3)
{
        stop("wrong command line argument format", call.=FALSE)
}
input_fn <- args[1] # output of ChromHMM OverlapEnrichment, also input to this program
output_fn <- args[2] # where the figure should be stored
state_annot_fn <- args[3] # where the data of states' characterization  are stored (state groups, ordering of states based on the group, etc.)

library(tidyverse)
library(tidyr)
library(dplyr)
library(pheatmap)
STATE_COLOR_DICT = c('quescient' = "#ffffff", "HET" = '#b19cd9', 'acetylations' = '#fffacd', 'enhancers'= '#FFA500', 'ct_spec enhancers' = '#FFA500', 'transcription' = '#006400', 'weak transcription' = '#228B22', 'weak enhancers' = '#ffff00', 'transcribed and enhancer' = '#ADFF2F', 'exon' = '#3cb371', 'promoters' = '#ff4500', 'TSS' = '#FF0000', 'weak promoters' = '#800080', 'polycomb repressed' = '#C0C0C0', 'others' = '#fff5ee', 'znf' = '#7fffd4', 'DNase' = '#fff44f')

read_state_annot_df <- function(state_annot_fn){
	annot_df <- as.data.frame(read.csv(state_annot_fn, header = TRUE, stringsAsFactors = FALSE))
	tryCatch({
		annot_df <- annot_df %>% rename('state_type' = 'Group')
		}, error = function (e) {message("tried to change column names in annot_df and nothing worth worrying happened")})
	annot_df <- annot_df %>% arrange(state_order_by_group) # order rows based on the index that we get from state_order_by_group column
	return(annot_df)	
}

calculate_gap_rows_among_state_groups <- function(state_annot_df){
	state_group_ordered_by_appearance <- unique(state_annot_df$state_type) # list of different state groups, ordered by how they appear in the heatmap from top to bottom
	count_df <- state_annot_df %>% count(state_type)
	count_df <- count_df[match(state_group_ordered_by_appearance, count_df$state_type),] # order the rows such that the state_type are ordered based on state_group_ordered_by_appearance
	results <- cumsum(count_df$n) # cumulative sum of the count for groups of states, which will be used to generate the gaps between rows of the heatmaps
	return(results)
}


draw_transition_matrix <- function(transition_fn, save_fig_fn, state_annot_fn){
	state_annot_df <- read_state_annot_df(state_annot_fn)
	df <- read.csv(transition_fn, sep = '\t', header = TRUE)
	df <- df %>% rename("state" = "state..from.to...Emission.order.")
	df <- df %>% select(-'state')	
	df <- df %>% slice(state_annot_df$state) # rearrange the states so that it follows the states from state_annot_df with states of the same group are put
	colnames(df) <- sapply(colnames(df), function(x){return(substr(x, 2, nchar(x)))}) # because currently all the colnames are like 'X1', 'X2' --> 'X100'
	df <- df %>% select(state_annot_df$state) # rearrange the columns so that the states of the same group are put next to each other.
	# change the colnames and rownames so that it will match with the annotation that we give when we draw the heatmap
	colnames(df) <- state_annot_df$mneumonics
	rownames(df) <- state_annot_df$mneumonics
	# prepare the state_annot_df for drawing
	rownames(state_annot_df) <- state_annot_df$mneumonics
	state_annot_df <- state_annot_df %>% select(c('state_type'))
	break_list <- seq (0, 1, by = 0.01)
	state_gap_indices <- calculate_gap_rows_among_state_groups(state_annot_df) # calculate the indices of lines that we will break the heatmap into such that states of the same group will be in the same block
	pheatmap(df, cluster_rows = FALSE, breaks = break_list, cluster_cols = FALSE, annotation_row = state_annot_df, annotation_col = state_annot_df, annotation_colors = list(state_type = STATE_COLOR_DICT), gaps_row = state_gap_indices,fontsize = 4.5, scale = 'none', , filename = save_fig_fn)
}

# call function for this file
draw_transition_matrix(transition_fn = input_fn, save_fig_fn = output_fn, state_annot_fn = state_annot_fn)
