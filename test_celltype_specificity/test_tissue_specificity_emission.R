library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

emission_fn <- "../ChromHMM_utilities/emission/emissions_100.txt" # TO BE UPDATED: this should be changed to the path to the file where you store the .txt emission matrix file, outputted from ChromHMM LearnModel and we also provided in the github for the paper  
output_folder <- "./output/"  # TO BE UPDATED: This should be changed to the path to the folder where you want the output of testing tissue specificity to be at.
dir.create(output_folder)

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

test_one_mark <- function(this_chrom_mark){
	this_mark_df <- emission_df %>% filter(chrom_mark == this_chrom_mark) # get only emission probabilities of experiments associated with the chromatin mark that we are trying to test
	distinct_ana <- sort(unique(this_mark_df$GROUP)) # unique GROUP of the cell types being questioned
	result_df <- data.frame(matrix(ncol = length(distinct_ana), nrow = num_states))
	colnames(result_df) <- distinct_ana
	rownames(result_df) <- paste0('S', seq(1, num_states))
	# result_df: rows: states
	# columns: each of the different cell groups/ GROUP that we are testing 
	for (ana in distinct_ana){
		this_ana_df <- this_mark_df %>% filter(GROUP == ana)
		other_ana_df <- this_mark_df %>% filter(GROUP != ana)
		for (state_index in seq(1, num_states)){
			x <- as.numeric(this_ana_df[paste0('S', state_index)][,1])
			y <- as.numeric(other_ana_df[paste0('S', state_index)][,1])
			t <- wilcox.test(x, y, alternative = 'greater')
			result_df[state_index, ana] <- t$p.value
			total_number_test <<- total_number_test + 1
		}
	}
	save_fn <- paste0(output_folder, this_chrom_mark, "_group_mann_whitney_test.txt")
	write.table(result_df, save_fn, sep = '\t', quote = FALSE)
}

emission_df <- get_emission_df_from_chromHMM_emission(emission_fn)
CHROM_MARK_TO_TEST = c('H3K9me3', 'H3K4me1', 'H3K4me3', 'H3K36me3', 'H3K27me3', 'H3K27ac', 'H3K9ac', 'DNase')
num_states = 100
total_number_test <<- 0
for (this_chrom_mark in CHROM_MARK_TO_TEST){
	test_one_mark(this_chrom_mark)
	cat("Finish testing for chrom_mark: ", this_chrom_mark, "\n")
}
cat("total_number_test: ", total_number_test, '\n')