library(ggplot2)
library(dplyr)
library(tidyr)

get_genome_context <- function(fn){
  result <- unlist(strsplit(fn, 'auc_'))[2]
  result <- unlist(strsplit(result, '.txt'))[1]
  return(result)
}

get_one_auc_df <- function(fn, context_name){
  df <- as.data.frame(read.csv(fn, header = FALSE, sep = '\t', stringsAsFactors = FALSE))
  colnames(df) <- c('model', 'auc')
  df$context <- rep(context_name, nrow(df))
  return(df)
}


get_genome_context <- function(fn){
  result <- unlist(strsplit(fn, 'auc_'))[2]
  result <- unlist(strsplit(result, '.txt'))[1]
  return(result)
}

draw_compare_auc_plots <- function(auc_folder, save_fn){
	auc_fn_list = list.files(path = auc_folder, pattern = '^auc_*') # list of files starting with auc_ in this folder
	genome_context = sapply(auc_fn_list, get_genome_context) # get the names of enrichment contexts that we are trying to compare among models based on the names of auc files
	auc_fn_list <- paste0(auc_folder, auc_fn_list)
	all_df <- data.frame(model = character(0), auc = numeric(0), context = character(0)) # df to store the auc values for all model that is not the full-stack model
	full_df <- data.frame(model = character(0), auc = numeric(0), context = character(0)) # df to store the auc values for full-stack models in all the genomic contexts
	for (i in 1:length(auc_fn_list)){
		df <- get_one_auc_df(auc_fn_list[i], genome_context[i])
		this_full <- df[1,]
		all_df <- rbind(all_df, df[-1,])
		full_df <- rbind(full_df, this_full)
	}
	p <- ggplot() +
	geom_boxplot(data = all_df, aes(x = context, y = auc)) +
	geom_point(data = full_df, aes(x = context, y = auc), color = 'red', size = 7) + 
	theme_bw() + 
	theme(legend.position = 'None', axis.text.x = element_text(angle = 90)) 
	ggsave(save_fn)
}

################ GETTING THE COMMAND LINE ARGUMENTS #####################
print ('This code create plots of the AUC of cell type speicific model vs full-stack model in terms of recovering genomic contexts of interests inside the auc_folder arugments. This code assumes that the first line of each of these files always shows the AUC of the full-stack model')
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2)
{
        stop("wrong command line argument format", call.=FALSE)
}
auc_folder <- args[1] # folder where all the auc data is stored
save_fn <- args[2] # where the figure will be saved. Should have tail .png
draw_compare_auc_plots(auc_folder, save_fn)
