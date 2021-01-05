library("tidyverse")
library("tidyr")
library("dplyr")
library('ggplot2')


draw_models_auc <- function(auc_folder, save_fn){
	auc_fn_list <- list.files(auc_folder, pattern = 'auc_*', full.names = TRUE)
	total_df <- data.frame()
	for (fn in auc_fn_list) {
		df <- read.table(fn, sep = '\t', header = FALSE)
		total_df <- rbind(total_df, df)
	}
	colnames(total_df) <- c('model', 'auc')
	p <- ggplot(total_df, aes(x = model, y = auc)) + 
	geom_boxplot() + 
	theme_bw()+ 
	theme(axis.text.x = element_text(size = 5, angle = 90)) 
	ggsave(save_fn)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2)
{
        stop("wrong command line argument format", call.=FALSE)
}
auc_folder <- args[1] # folder that contains all the auc files that we care about
save_fn <- args[2] # where the figures should be stored
draw_models_auc(auc_folder, save_fn)