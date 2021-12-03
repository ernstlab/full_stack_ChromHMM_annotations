# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

library("tidyverse")
library("tidyr")
library("dplyr")
library('ggplot2')
library(pheatmap)
# STATE_COLOR_DICT = c('quescient' = "#ffffff", "HET" = '#b19cd9', 'acetylations' = '#fffacd', 'enhancers'= '#FFA500', 'ct_spec enhancers' = '#FFA500', 'transcription' = '#006400', 'weak transcription' = '#228B22', 'weak enhancers' = '#ffff00', 'transcribed and enhancer' = '#ADFF2F', 'exon' = '#3cb371', 'promoters' = '#ff4500', 'TSS' = '#FF0000', 'weak promoters' = '#800080', 'bivalent promoters' = '#7030A0', 'polycomb repressed' = '#C0C0C0', 'others' = '#fff5ee', 'znf' = '#7fffd4', 'DNase' = '#fff44f')
get_rid_of_bedGZ_context_name <- function(enrichment_name){
	if (endsWith(enrichment_name, '.bed.gz')){
		result <- substring(enrichment_name, 1, nchar(enrichment_name) - 7)
		return(result)
	}
	return(enrichment_name)
}

get_overlap_enrichment_df <- function(enrichment_fn){
	if (endsWith(enrichment_fn, '.xlsx') | endsWith(enrichment_fn, '.xls')) {
		print("READING AN EXCEL FILE")
		enrichment_df <- as.data.frame(read_excel(enrichment_fn, sheet = 1, col_names = TRUE))
	} # if this is a excel file
	else{
		print ("READING A TEXT FILE")
		enrichment_df <- as.data.frame(read_tsv(enrichment_fn, col_names = TRUE))
	}
	tryCatch({
			enrichment_df <- enrichment_df %>% rename("state" = "state..Emission.order.", "percent_in_genome" = "Genome..", "state" = "state (Emission order)", "state (Emission order)" = "state") # rename("percent_in_genome" = "percent_in_genome")
		}, error = function(e) {message("Nothing bad happended")})
	score_colnames <- colnames(enrichment_df)[2:ncol(enrichment_df)]
	score_colnames <- sapply(score_colnames, get_rid_of_bedGZ_context_name)
	score_colnames <- c('state', score_colnames)
	colnames(enrichment_df) <- score_colnames
	num_state <- nrow(enrichment_df) 
	tryCatch({
	  enrichment_df <- enrichment_df %>% select(-"percent_in_genome")
	}, error = function(e) {message("There is not a percent_in_genome column. Nothing bad happend")})
	return(enrichment_df)
}

read_state_annot_df <- function(){
  state_annot_fn <- '~/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/figures/supp_excel/state_annotations_processed.csv'
	annot_df <- as.data.frame(read.csv(state_annot_fn, header = TRUE, stringsAsFactors = FALSE))
	tryCatch({
		annot_df <- annot_df %>% rename('state_type' = 'Group')
		}, error = function (e) {message("tried to change column names in annot_df and nothing worth worrying happened")})
	annot_df <- annot_df %>% arrange(state_order_by_group) # order rows based on the index that we get from state_order_by_group column
	STATE_COLOR_DICT <- as.character(annot_df$color)
	names(STATE_COLOR_DICT) <- as.character(annot_df$state_type)	
	STATE_COLOR_DICT <- STATE_COLOR_DICT[!duplicated(names(STATE_COLOR_DICT))]
	return_obj <- list(annot_df = annot_df, STATE_COLOR_DICT = STATE_COLOR_DICT)
	return(return_obj)
}

calculate_gap_rows_among_state_groups <- function(state_annot_df){ 
	state_group_ordered_by_appearance <- unique(state_annot_df$state_type) # list of different state groups, ordered by how they appear in the heatmap from top to bottom
	count_df <- state_annot_df %>% count(state_type)
	count_df <- count_df[match(state_group_ordered_by_appearance, count_df$state_type),] # order the rows such that the state_type are ordered based on state_group_ordered_by_appearance
	results <- cumsum(count_df$n) # cumulative sum of the count for groups of states, which will be used to generate the gaps between rows of the heatmaps
	return(results)
}

draw_enrichment_plot_no_color_scale <- function (enrichment_fn, save_fig_fn) {
	# draw enrichment plot where each column (each enrichment context) is on its own color scale
	annot_obj <- read_state_annot_df() # state: one-based state indices, arranged based on the groups of states, long annotations, short annotations, state_type, color, itemRgb, state_ordered_by_group: 0 --> 99 already, mneumonics in the format <state_type><state_type_index_one_based>
	annot_df <- annot_obj$annot_df
	STATE_COLOR_DICT <- annot_obj$STATE_COLOR_DICT
	enrichment_df <- get_overlap_enrichment_df(enrichment_fn) # columns: enrichment contex, rows: states
	enrichment_colnames <-  colnames(enrichment_df)[-1] # get rid of the first item in enrichment_colnames because we want to get rid of the 'state'
	# rearranged_colnames <- paste("consHMM_state", seq(1, 100) , sep = "")
	enrichment_df <- enrichment_df %>% slice(annot_df$state) # rearrange the states based on the ordered based on the order from annot_df
	enrichment_df <- as.data.frame(enrichment_df[, -1]) # get rid of the first column (states) so that the enrichment data is just the enrichment values. This is nescessary to that we can get input for pheatmap function later
	colnames(enrichment_df) <- enrichment_colnames
	rownames(enrichment_df) <- annot_df$mneumonics # row names and columns names are the names of x and y axis tick labels in the figure
	display_num_df <- round(enrichment_df, digits = 1)
	# calculate the indices of rows where we want to create gaps among groups of states
	gap_row_indices <- calculate_gap_rows_among_state_groups(annot_df)
	# prepare the annot_df for plotting as the row_annotation in the pheatmap function	
	rownames(annot_df) <- annot_df$mneumonics
	annot_df <- annot_df %>% select(c('state_type')) 
	#pheatmap(enrichment_df, cluster_rows = FALSE, cluster_cols = TRUE, fontsize = 4.5, scale = 'none', display_numbers = display_num_df, number_format = '%.1f', color = colorRampPalette(c('white', 'red'))(300), filename = save_fig_fn, cellwidth = 10, lengend = FALSE) 
	pheatmap(enrichment_df, cluster_rows = FALSE, cluster_cols = FALSE,annotation_row = annot_df, gaps_row = gap_row_indices, annotation_colors = list(state_type = STATE_COLOR_DICT), fontsize = 8, scale = 'none', color = colorRampPalette(c('white', 'blue'))(300), filename = save_fig_fn, cellwidth = 9, cellheight = 8, lengend = TRUE) 
	print(paste('DONE! Figure is saved at:', save_fig_fn))	
}
state_annot_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/figures/supp_excel/state_annotations_processed.csv' 
# tss_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/neighborhood_enrichment/hg19/neighborhood_enrichment_tss_hg19.txt'
# save_tss_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/neighborhood_enrichment/hg19/neighborhood_enrichment_tss_hg19.png'
# draw_enrichment_plot_no_color_scale(tss_fn, save_tss_fn)
# tes_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/neighborhood_enrichment/hg19/neighborhood_enrichment_tes_hg19.txt'
# save_tes_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/neighborhood_enrichment/hg19/neighborhood_enrichment_tes_hg19.png'
# draw_enrichment_plot_no_color_scale(tss_fn, save_tss_fn)
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2)
{
	stop("wrong command line argument format", call.=FALSE)
}
neighborhood_enrichment_fn <- args[1] # output of ChromHMM NeighborhoodEnrichment, also input to this program
save_fn <- args[2] # where the figures should be stored
draw_enrichment_plot_no_color_scale(neighborhood_enrichment_fn, save_fn)