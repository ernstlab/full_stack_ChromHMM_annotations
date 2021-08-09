library("tidyverse")
library("tidyr")
library("dplyr")
library('ggplot2')
library(pheatmap)
library(reshape2)
STATE_COLOR_DICT = c('quescient' = "#ffffff", "HET" = '#b19cd9', 'acetylations' = '#fffacd', 'enhancers'= '#FFA500', 'ct_spec enhancers' = '#FFA500', 'transcription' = '#006400', 'weak transcription' = '#228B22', 'weak enhancers' = '#ffff00', 'transcribed and enhancer' = '#ADFF2F', 'exon' = '#3cb371', 'promoters' = '#ff4500', 'TSS' = '#FF0000', 'weak promoters' = '#800080', 'polycomb repressed' = '#C0C0C0', 'others' = '#fff5ee', 'znf' = '#7fffd4', 'DNase' = '#fff44f')
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
	}, error = function(e) {message("There is not a percent_in_genome column")})
	enrichment_df <- as.data.frame(t(enrichment_df))
	colnames(enrichment_df) <- enrichment_df[1,] # first row is also header
	enrichment_df <- enrichment_df[-1,] # get rid of the first row
	enrichment_df$coord <- as.integer(rownames(enrichment_df))
	colnames(enrichment_df)[-length(colnames(enrichment_df))] <- sapply(colnames(enrichment_df)[-length(colnames(enrichment_df))], function (x){paste0('S',x)}) # convert columns from state_index to "S<state_index>". exclude that last column because that's coord
	enrichment_df <- melt(enrichment_df, id = c('coord')) 
	colnames(enrichment_df) <- c('coord', 'state', 'enrichment')
	enrichment_df$state <- as.character(enrichment_df$state)
	return(enrichment_df)
}

get_state_type <- function(){
  state_annot_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/state_annotations.csv'
  state_annot_df <- as.data.frame(read_csv(state_annot_fn, col_names = T))
  state_annot_df$state <- sapply(as.character(state_annot_df$state), function (x){paste0('S', x)})
  return(state_annot_df)
}

draw_2DLine_neighborhood_enrichment <- function (enrichment_fn, save_fig_fn) {
	# draw enrichment plot where each column (each enrichment context) is on its own color scale
	enrichment_df <- get_overlap_enrichment_df(enrichment_fn) # columns: enrichment contex, rows: states
	annot_df <- get_state_type()
	enrichment_df <- enrichment_df %>% left_join(annot_df, by = 'state')
	color_dict <- enrichment_df$color
	names(color_dict) <- enrichment_df$state 
	#color_dict has keys: states, values: the coloro corresponding to that state
	ggplot(data = enrichment_df, aes(x = coord, y = enrichment, group = state)) +
	geom_line(aes(color = state)) +
	scale_colour_manual(values = color_dict) +
	theme_bw() +
	theme(legend.position = 'None')
	ggsave(save_fig_fn)
	print(paste('DONE! Figure is saved at:', save_fig_fn))	
}

tss_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/neighborhood_enrichment/hg19/neighborhood_enrichment_tss_hg19.txt'
save_tss_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/neighborhood_enrichment/hg19/2DLine_neighborhood_enrichment_tss_hg19.png'
draw_2DLine_neighborhood_enrichment(tss_fn, save_tss_fn)
tes_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/neighborhood_enrichment/hg19/neighborhood_enrichment_tes_hg19.txt'
save_tes_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/neighborhood_enrichment/hg19/2DLine_neighborhood_enrichment_tes_hg19.png'
draw_2DLine_neighborhood_enrichment(tes_fn, save_tes_fn)