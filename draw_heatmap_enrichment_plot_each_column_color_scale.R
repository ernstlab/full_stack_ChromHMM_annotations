# this program will take as input the ovelap enrichment output from ChromHMM and draw a figure that show such data. The resulting plot will be red-based, and each column will have its own color scale and enrihcment values will also be denoted in the plot.
#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# if (length(args) != 2)
# {
#         stop("wrong command line argument format", call.=FALSE)
# }
# input_fn <- args[1] # output of ChromHMM OverlapEnrichment, also input to this program
# output_fn <- args[2] # where the figure should be stored

library(tidyverse)
library(tidyr)
library(dplyr)
library(pheatmap)
library(readxl)
GENOME_CONTEXT_COLOR = c("CpGIsland" = '#73C6B6' ,"RefSeqExon" = '#D2B4DE', "RefSeqGene" = '#F1948A' ,"RefSeqTES" = '#F7DC6F', "RefSeqTSS" = '#2980B9', "RefSeqTSS2kb" = '#27AE60', "laminB1lads" = '#D7DBDD', 'clean_assembly_gap' = '#8E44AD', 'znf_genes' = '#7B241C')
STATE_COLOR_DICT = c('quescient' = "#ffffff", "HET" = '#b19cd9', 'acetylations' = '#fffacd', 'enhancers'= '#FFA500', 'ct_spec enhancers' = '#FFA500', 'transcription' = '#006400', 'weak transcription' = '#228B22', 'weak enhancers' = '#ffff00', 'transcribed and enhancer' = '#ADFF2F', 'exon' = '#3cb371', 'promoters' = '#ff4500', 'TSS' = '#FF0000', 'weak promoters' = '#800080', 'polycomb repressed' = '#C0C0C0', 'others' = '#fff5ee', 'znf' = '#7fffd4', 'DNase' = '#fff44f')

get_rid_of_bedGZ_context_name <- function(enrichment_name){
	if (endsWith(enrichment_name, '.bed.gz')){
		result <- substring(enrichment_name, 1, nchar(enrichment_name) - 7)
		return(result)
	}
	return(enrichment_name)
}

read_state_annot_df <- function(){
	state_annot_fn <- '~/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/figures/supp_excel/state_annotations_processed.csv'
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

get_overlap_enrichment_df <- function(enrichment_fn){
	if (endsWith(enrichment_fn, '.xlsx') | endsWith(enrichment_fn, '.xls')) {
		print("READING AN EXCEL FILE")
		enrichment_df <- as.data.frame(read_excel(enrichment_fn, sheet = 1, col_names = TRUE))
	} # if this is a excel file
	else{
		print ("READING A TEXT FILE")
		enrichment_df <- as.data.frame(read.table(enrichment_fn, sep = "\t", header = TRUE))
	}
	tryCatch({
			enrichment_df <- enrichment_df %>% rename("state" = "state..Emission.order.", "percent_in_genome" = "Genome..") # rename("percent_in_genome" = "percent_in_genome")
		}, error = function(e) {message("Nothing bad happended")})
	score_colnames <- colnames(enrichment_df)[3:ncol(enrichment_df)]
	score_colnames <- sapply(score_colnames, get_rid_of_bedGZ_context_name)
	score_colnames <- c(colnames(enrichment_df)[1:2], score_colnames)
	colnames(enrichment_df) <- score_colnames
	num_state <- nrow(enrichment_df) - 1
	enrichment_df <- enrichment_df[1:num_state, ]  # get rid of the last row because it does not  contain infrormation about enrichment but rather information about summary of different components
	enrichment_df <- enrichment_df %>% select(-"percent_in_genome")
	enrichment_df$state <- as.integer(as.character(enrichment_df$state)) # convert from factor to integer
	return(enrichment_df)
}

draw_enrichment_plot_no_color_scale <- function (enrichment_fn, save_fig_fn) {
	# draw enrichment plot where each column (each enrichment context) is on its own color scale
	enrichment_df <- get_overlap_enrichment_df(enrichment_fn) # columns: enrichment contex, rows: states
	enrichment_colnames <-  colnames(enrichment_df)[-1]
	# rearranged_colnames <- paste("consHMM_state", seq(1, 100) , sep = "")
	# enrichment_df <- enrichment_df %>% select(c(colnames(enrichment_df)[1], rearranged_colnames)) # rearranged the enrichment_df so that the columns are consHMM_state1 -> consHMM_state100
	# colnames(enrichment_df) <- c('state', paste("consHMM_", seq(1, 100), sep = "")) # change enrichment contexts names from consHMM_state1 to consHMM_1
	enrichment_df <- as.data.frame(enrichment_df[, c(12, seq(2,11))]) # get rid of the first column (states) so that the enrichment data is just the enrichment values. This is nescessary to that we can get input for pheatmap function later
	state_annot_df <- read_state_annot_df()
	enrichment_df <- enrichment_df %>% slice(state_annot_df$state) # order state based on the group, based on the state_annot_df
	rownames(enrichment_df) <- state_annot_df$mneumonics
	colnames(enrichment_df) <- enrichment_colnames
	display_num_df <- round(enrichment_df, digits = 1)
	# create the df for annotationof rows in our heatmap
	rownames(state_annot_df) <- state_annot_df$mneumonics
	state_annot_df <- state_annot_df %>% select(state_type)
	# cell_group_df <- data.frame(cell_group = plot_df$GROUP)
	# calculate the gaps between groups of states
	gap_row_indices <- calculate_gap_rows_among_state_groups(state_annot_df)
	pheatmap(enrichment_df, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 4.5, scale = 'none', display_numbers = display_num_df, number_format = '%.1f', annotation_row = state_annot_df, annotation_colors = list(state_type = STATE_COLOR_DICT), gaps_row = gap_row_indices, color = colorRampPalette(c('white', 'red'))(300), filename = save_fig_fn, cellwidth = 10, lengend = FALSE) 
	print(paste('DONE! Figure is saved at:', save_fig_fn))	
}

draw_enrichment_plot_each_column_color_scale <- function (enrichment_fn, save_fig_fn) {
	# draw enrichment plot where each column (each enrichment context) is on its own color scale
	enrichment_df <- get_overlap_enrichment_df(enrichment_fn) # columns: enrichment contex, rows: states
	enrichment_colnames <- colnames(enrichment_df)[-1]
	state_list <- enrichment_df$state
	enrichment_df <- as.data.frame(enrichment_df[, -1]) # get rid of the first column (states) so that the enrichment data is just the enrichment values. This is nescessary to that we can get input for pheatmap function later
	colnames(enrichment_df) <- enrichment_colnames
	rownames(enrichment_df) <- state_list # row names and columns names are the names of x and y axis tick labels in the figure
	maxCol <- sapply(enrichment_df, max, na.rm = TRUE) # max enrichment value for each enrichment context
	minCol <- sapply(enrichment_df, min, na.rm = TRUE) # min enrichment value for each enrichment context
	rangeCol <- maxCol - minCol # range of enrichment values fro each enrichemnt context

	# now onto building the scaled enrichment heatmap where each column is on its own color scale 
	hm_df <- data.frame(matrix(ncol = ncol(enrichment_df), nrow = nrow(enrichment_df)))
	colnames(hm_df) <- colnames(enrichment_df)
	rownames(hm_df) <- state_list
	for (cn in colnames(hm_df)) { # for each column name, i.e. for each enrichment context
		this_min = minCol[cn]
		this_range = rangeCol[cn]
		hm_df[cn] <- (enrichment_df[cn] - this_min) / this_range
	}
	display_num_df <- round(enrichment_df, digits = 1) # this df is used to display the value of enrichments on the plot.
	pheatmap(hm_df, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 4.5, scale = 'none', display_numbers = display_num_df, number_format = '%.1f', color = colorRampPalette(c('white', 'red'))(300), filename = save_fig_fn, cellwidth = 10, legend = FALSE) # draw the heatmap and save it into a file
	print(paste('DONE! Figure is saved at:', save_fig_fn))
}

draw_genome_context_enrichment_plot <- function (enrichment_fn, save_fig_fn) {
	enrichment_df <- get_overlap_enrichment_df(enrichment_fn) # columns: enrichment contex, rows: states
	enrichment_colnames <- colnames(enrichment_df)[-1]
	state_list <- enrichment_df$state
	enrichment_df <- as.data.frame(enrichment_df[, -1]) # get rid of the first column (states) so that the enrichment data is just the enrichment values. This is nescessary to that we can get input for pheatmap function later
	colnames(enrichment_df) <- enrichment_colnames
	rownames(enrichment_df) <- state_list # row names and columns names are the names of x and y axis tick labels in the figure

	maxCol <- sapply(enrichment_df, max, na.rm = TRUE) # max enrichment value for each enrichment context
	minCol <- sapply(enrichment_df, min, na.rm = TRUE) # min enrichment value for each enrichment context
	rangeCol <- maxCol - minCol # range of enrichment values fro each enrichemnt context

	# now onto building the scaled enrichment heatmap where each column is on its own color scale 
	hm_df <- data.frame(matrix(ncol = ncol(enrichment_df), nrow = nrow(enrichment_df)))
	colnames(hm_df) <- colnames(enrichment_df)
	rownames(hm_df) <- state_list
	for (cn in colnames(hm_df)) { # for each column name, i.e. for each enrichment context
		this_min = minCol[cn]
		this_range = rangeCol[cn]
		hm_df[cn] <- (enrichment_df[cn] - this_min) / this_range
	}
	display_num_df <- round(enrichment_df, digits = 1) # this df is used to display the value of enrichments on the plot.

	max_fold_enrich_context <- colnames(enrichment_df)[apply(enrichment_df, 1, which.max)]
	rowAnnot_df <- data.frame(max_context = max_fold_enrich_context)
	rownames(rowAnnot_df) <- state_list
	colAnnot_df <- data.frame(enrichment_context = enrichment_colnames)
	rownames(colAnnot_df) <- enrichment_colnames
	pheatmap(hm_df, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 4.5, scale = 'none', display_numbers = display_num_df, number_format = '%.1f', color = colorRampPalette(c('white', 'red'))(300), filename = save_fig_fn, cellwidth = 10, legend = FALSE, annotation_row = rowAnnot_df, annotation_col = colAnnot_df, annotation_legend = F, annotation_names_col = F, annotation_colors = list(enrichment_context = GENOME_CONTEXT_COLOR, max_context = GENOME_CONTEXT_COLOR)) # draw the heatmap and save it into a file
} 
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2)
{
	stop("wrong command line argument format", call.=FALSE)
}
enrichment_fn <- args[1] # output of ChromHMM OverlapEnrichment, also input to this program
save_fn <- args[2] # where the figures should be stored
# draw_enrichment_plot_no_color_scale(enrichment_fn, save_fn)
# call function for this file
# draw_enrichment_plot_no_color_scale(enrichment_fn = input_fn, save_fig_fn = output_fn)
draw_genome_context_enrichment_plot(input_fn, output_fn)