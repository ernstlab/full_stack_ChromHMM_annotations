#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2)
{
        stop("wrong command line argument format", call.=FALSE)
}
input_fn <- args[1] # output of ChromHMM OverlapEnrichment, also input to this program
output_folder <- args[2] # where the figures should be stored
dir.create(output_folder, recursive = TRUE) # create the output folder
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)
# source('/Users/vuthaiha/Desktop/window_hoff/source/summary_analysis/draw_overlap_enrichment_heatmap.R')
annot_fn <- '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/state_annotations.csv'
print(paste("The file of states' annotation is: ", annot_fn))
annot_df <- as.data.frame(read.csv(annot_fn, sep = ',', header = TRUE))
STATE_COLOR_DICT <- c('quescient' = "white", "HET" = '#DDA0DD', 
        'acetylations' = 'lemonchiffon', 
        'enhancers'= 'orange', 
        'transcription' = '#306754', 
        'weak enhancers' = 'yellow', 
        'transcribed and enhancer' = 'greenyellow', 
        'exon' = 'mediumseagreen', 
        'promoters' = 'orangered', 
        'TSS' = 'red', 
        'weak promoters' = 'purple', 
        'polycomb repressed' = '#C0C0C0',
        'others' = '#E7D0CA', 
        'znf' = 'aquamarine', 
        'DNase' = '#fff44f'
        )

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
	state_col_index <- which(colnames(enrichment_df) %in% c('state')) # get the index of the column that contains 'state'. We will only take into account the columsn from state leftward of the enrichment_fn
	enrichment_df <- enrichment_df[1:num_state, state_col_index[1]:ncol(enrichment_df)]  # get rid of the last row because it does not  contain infrormation about enrichment but rather information about summary of different components, Also get rid of any columns preceeding the 'state' column
	enrichment_df <- enrichment_df %>% select(-"percent_in_genome")
	return(enrichment_df)
}

get_enrichment_df_combined_state_annot <- function(input_fn){
	df <- get_overlap_enrichment_df(input_fn)
	df$state <- as.integer(as.character(df$state)) # change from factor to integer type, for better plot later
	combined_df <- annot_df %>% left_join(df, by = c('state' = 'state')) #state from the enrichment df is the same as states from the state_annot df
	return (combined_df)
}

draw_barplot_one_enrichment_context  <- function(combined_df, save_fn, enr_cont_name) {
	num_states <- nrow(combined_df) # number of chromHMM states
	ggplot(combined_df, aes_string(x = 'state', y = enr_cont_name, fill = 'state_type')) +
	geom_bar(stat = 'identity') + 
	theme_bw() +
	scale_x_discrete(name = 'state', limits = seq(1, num_states)) + 
	# coord_cartesian(ylim = c(0, 2))  +
	#scale_y_continuous(name = 'enrichment', breaks = c(0, 0.5, 1, 1.5, 2)) +
	geom_hline(yintercept = 1.0) +
	theme(axis.text.x = element_text(size = 3.6, angle = 90), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 6), legend.position = 'bottom', legend.title = element_blank()) +
	scale_fill_manual('state_type', values = STATE_COLOR_DICT)
	ggsave(save_fn, width = 6, height = 3, dpi = 300) # save the figure that we just created
}

create_barPlots_enrichment <- function(input_fn, output_folder){
	# get the data frame with state, percent_in_genome, state_type, other enrichment context columns
	combined_df <- get_enrichment_df_combined_state_annot(input_fn)
	# get the list of enrichment contexts that we have data for, and so we will create figures corresponding to each of these enrichment contexts
	enr_cont_list <- setdiff(colnames(combined_df), c('state', 'state_annotation', 'state_type', 'color')) 
	for (enrichment_context in enr_cont_list){
		save_fn <- file.path(output_folder, paste0(enrichment_context, '.png'))
		draw_barplot_one_enrichment_context(combined_df, save_fn, enrichment_context)
		print (paste("Done drawing figure for ", enrichment_context))
	}
}

# call function for this file
create_barPlots_enrichment(input_fn, output_folder)
