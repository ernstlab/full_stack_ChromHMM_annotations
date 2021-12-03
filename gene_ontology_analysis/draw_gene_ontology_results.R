library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggpubr)
input_folder <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/gene_ontology_analysis/GREAT/combined_all_states/'
output_folder <- file.path(input_folder, 'plots')
GO_TYPE_LIST <- c('GO_BP', 'GO_CC', 'GO_MF')#, 'HP')
state_group_list <- c('DNase', 'HET', 'TSS', 'acetylations', 'enhancers', 'exon', 'others', 'polycomb_repressed', 'promoters', 'quescient', 'transcribed_and_enhancer', 'transcription', 'weak_enhancers', 'weak_promoters', 'weak_transcription', 'znf')
STATE_COLOR_DICT = c('quescient' = "#ffffff", "HET" = '#b19cd9', 'acetylations' = '#fffacd', 'enhancers'= '#FFA500', 'ct_spec enhancers' = '#FFA500', 'transcription' = '#006400', 'weak transcription' = '#228B22', 'weak enhancers' = '#ffff00', 'transcribed and enhancer' = '#ADFF2F', 'exon' = '#3cb371', 'promoters' = '#ff4500', 'TSS' = '#FF0000', 'weak promoters' = '#800080', 'polycomb repressed' = '#C0C0C0', 'others' = '#fff5ee', 'znf' = '#7fffd4', 'DNase' = '#fff44f')
GO_TYPE_COLOR_DICT = c('GO_BP' = '#FF5733', 'GO_CC' = '#FFC300', 'GO_MF' = "#3498DB") #, 'HP' = '#A569BD')

read_state_annot_df <- function(){
	state_annot_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/figures/supp_excel/state_annotations_processed.csv'
	annot_df <- as.data.frame(read.csv(state_annot_fn, header = TRUE, stringsAsFactors = FALSE))
	tryCatch({
		annot_df <- annot_df %>% rename('state_type' = 'Group')
		}, error = function (e) {message("tried to change column names in annot_df and nothing worth worrying happened")})
	annot_df <- annot_df %>% arrange(state_order_by_group) # order rows based on the index that we get from state_order_by_group column
	return(annot_df)	
}

read_one_input_fn <- function(GO_TYPE, group_name){
	fn <- file.path(input_folder, GO_TYPE, paste0(group_name, '_topGO.txt.gz'))
	df <- as.data.frame(read.table(fn, header = TRUE, sep = '\t', stringsAsFactors = F))
	df$go_type <- GO_TYPE
	return(df)
}

calculate_gap_rows <- function(row_annot_df, colname_to_count){
	go_group_ordered_by_appearance <- unique(row_annot_df[colname_to_count]) # list of different state groups, ordered by how they appear in the heatmap from top to bottom
	count_df <- row_annot_df %>% count(!!as.name(colname_to_count))
	results <- cumsum(count_df$n) # cumulative sum of the count for groups of states, which will be used to generate the gaps between rows of the heatmaps
	return(results)
}

get_all_states_max_FoldEnrich_df <- function(output_folder, go_type){
	df <- read_one_input_fn(go_type, 'all_states')
	df <- df %>% group_by(state) %>% slice_max(FoldEnrich, n = 1, with_ties = F) %>% ungroup()
	df$GO <- paste0(df$ID, ":", df$Desc)
	df$neglogFdrQ <- -log10(df$FdrQ)
	df$neglogFdrQ[df$neglogFdrQ > 100] <- 100 # any p values that are less than e-100 will be capped at e-100, including the pvalues 0.0
	# annotation of the GO terms
	df <- df %>% select(c('GO', 'state', 'neglogFdrQ'))
	w <- as.data.frame(pivot_wider(df, names_from = state, values_from = neglogFdrQ))
	rownames(w) <- w$GO
	w <- w %>% select(-c('GO'))
	# now, w will have rows: GO, columns: state. Next, rearrange the states such that they appear as how they are supposed to TSS1--> TSS2 etc. Then, we will also create an annotation of the GO terms that will be used for pheatmap coloring of rows and columns later
	# annotation of the states
	state_annot_plot_df <- data.frame(state = colnames(w))
	state_annot_plot_df <- state_annot_plot_df %>% merge(read_state_annot_df(), by = 'state') %>% arrange(state_order_by_group) %>% select(c('mneumonics', 'state_type', 'state'))
	w <- w %>% select(state_annot_plot_df$state) # rearrange the columns in w so that the states are ordered they way they appear in the paper
	w <- w %>% select(state_annot_plot_df$state) # rearrange the columns in w so that the states are ordered they way they appear in the paper
	colnames(w) <- state_annot_plot_df$mneumonics
	rownames(state_annot_plot_df) <- state_annot_plot_df$mneumonics
	state_annot_plot_df <- state_annot_plot_df %>% select('state_type')
	# now prepare the matrix to show the magnitude of the p values
	text_df <- round(w, digits = 1)
	text_df[is.na(text_df)] <- '' # by now all the 
	text_df[text_df == 100] <- '>=100'
	gap_col_indices <- calculate_gap_rows(state_annot_plot_df, 'state_type')
	save_fn <- file.path(output_folder, paste0(go_type, '_allStates_topGO.png'))
	p <- pheatmap(w, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = text_df, fontsize = 9, scale = 'none', na_col = '#F5F1F8', color = colorRampPalette(c('white', 'red'))(300), annotation_col = state_annot_plot_df, annotation_colors = list(state_type = STATE_COLOR_DICT, go_type = GO_TYPE_COLOR_DICT), gaps_col =  gap_col_indices, filename = save_fn) # draw the heatmap and save it into a file	
	return(p)
}

create_heatmap_one_input_file <- function(group_name, output_folder){
	df_list <- lapply(GO_TYPE_LIST, read_one_input_fn, group_name = group_name) # a list of data frame each corrsponding to one GO_TYPE
	df <- bind_rows(df_list)
	df$GO <- paste0(df$ID, ":", df$Desc)
	df$neglogFdrQ <- -log10(df$FdrQ)
	df$neglogFdrQ[df$neglogFdrQ > 100] <- 100 # any p values that are less than e-100 will be capped at e-100, including the pvalues 0.0
	# annotation of the GO terms
	go_annot_df <- df %>% select(c('GO','go_type')) %>% distinct()
	df <- df %>% select(c('GO', 'state', 'neglogFdrQ'))
	w <- as.data.frame(pivot_wider(df, names_from = state, values_from = neglogFdrQ))
	rownames(w) <- w$GO
	w <- w %>% select(-c('GO'))
	go_annot_plot_df <- data.frame(uniq_GO = rownames(w))
	go_annot_plot_df <- go_annot_plot_df %>% left_join(go_annot_df, by = c('uniq_GO' = 'GO'))
	rownames(go_annot_plot_df) <- go_annot_plot_df$uniq_GO
	go_annot_plot_df <- go_annot_plot_df %>% select(-c('uniq_GO'))
	# now, w will have rows: GO, columns: state. Next, rearrange the states such that they appear as how they are supposed to TSS1--> TSS2 etc. Then, we will also create an annotation of the GO terms that will be used for pheatmap coloring of rows and columns later
	# annotation of the states
	state_annot_plot_df <- data.frame(state = colnames(w))
	state_annot_plot_df <- state_annot_plot_df %>% merge(read_state_annot_df(), by = 'state') %>% arrange(state_order_by_group) %>% select(c('mneumonics', 'state_type', 'state'))
	w <- w %>% select(state_annot_plot_df$state) # rearrange the columns in w so that the states are ordered they way they appear in the paper
	colnames(w) <- state_annot_plot_df$mneumonics
	rownames(state_annot_plot_df) <- state_annot_plot_df$mneumonics
	state_annot_plot_df <- state_annot_plot_df %>% select('state_type')
	# now prepare the matrix to show the magnitude of the p values
	text_df <- round(w, digits = 1)
	text_df[is.na(text_df)] <- '' # by now all the 
	text_df[text_df == 100] <- '>=100'
	gap_row_indices <- calculate_gap_rows(go_annot_plot_df, 'go_type')
	save_fn <- file.path(output_folder, paste0(group_name, '_topGO.png'))
	p <- pheatmap(w, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = text_df, fontsize = 9, scale = 'none', na_col = '#F5F1F8', color = colorRampPalette(c('white', 'red'))(300), annotation_row = go_annot_plot_df, annotation_colors = list(state_type = STATE_COLOR_DICT, go_type = GO_TYPE_COLOR_DICT), gaps_row =  gap_row_indices, filename = save_fn) # draw the heatmap and save it into a file	
	return(p)
}
state_group_list <- c('DNase', 'HET', 'TSS', 'acetylations', 'enhancers', 'exon', 'others', 'polycomb_repressed', 'promoters', 'quescient', 'transcribed_and_enhancer', 'transcription', 'weak_enhancers', 'weak_promoters', 'weak_transcription', 'znf')
for (state_group in c('znf')){
	create_heatmap_one_input_file(state_group, output_folder)
}