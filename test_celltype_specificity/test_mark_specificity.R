library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
source("/Users/vuthaiha/Desktop/window_hoff/source/summary_analysis/draw_emission_matrix_functions.R")

emission_fn <- "/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/emissions_100.txt"
output_folder <- "/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/chrom_mark_spec_test/"
dir.create(output_folder)
emission_df <- get_emission_df_from_chromHMM_emission(emission_fn)
num_states = 100
total_number_test <<- 0
test_chrom_mark_specificity <- function(){
	distinct_marks <- sort(unique(emission_df$chrom_mark)) # unique GROUP of the cell types being questioned
	result_df <- data.frame(matrix(ncol = length(distinct_marks), nrow = num_states))
	colnames(result_df) <- distinct_marks
	rownames(result_df) <- paste0('S', seq(1, num_states))
	# result_df: rows: states
	# columns: each of the different cell groups/ GROUP that we are testing 
	for (chromM in distinct_marks){
		this_chromM_df <- emission_df %>% filter(chrom_mark == chromM)
		other_chromM_df <- emission_df %>% filter(chrom_mark != chromM)
		for (state_index in seq(1, num_states)){
			x <- as.numeric(this_chromM_df[paste0('S', state_index)][,1])
			y <- as.numeric(other_chromM_df[paste0('S', state_index)][,1])
			t <- wilcox.test(x, y, alternative = 'greater')
			result_df[state_index, chromM] <- t$p.value
			total_number_test <<- total_number_test + 1
		}
	}
	save_fn <- paste0(output_folder, "chrom_mark_group_mann_whitney_test.txt")
	write.table(result_df, save_fn, sep = '\t', quote = FALSE)
	cat("Done! Results are saved in ", save_fn, "\n")
}

test_chrom_mark_specificity()
cat("total_number_test: ", total_number_test, '\n')
