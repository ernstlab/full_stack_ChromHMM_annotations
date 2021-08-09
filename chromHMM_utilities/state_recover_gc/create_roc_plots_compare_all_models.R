library(ggplot2)
library(dplyr)
library(tidyr)
source('~/Desktop/window_hoff/source/summary_analysis/draw_emission_matrix_functions.R')
ANATOMY_COLOR_CODE = c('full' = '#000000', 'BLOOD'= '#e34a33', 'ESC_DERIVED'= '#fdbb84', 'ESC'= '#fef0d9', 'BRAIN'= '#ffffd4', 'LUNG'= '#006837', 'SKIN'= '#78c679', 'MUSCLE'= '#c2e699', 'IPSC'= '#7a0177', 'GI_STOMACH' = '#f768a1', 'HEART'= '#fbb4b9', 'BREAST' = '#d7b5d8', 'GI_INTESTINE' = '#253494', 'GI_RECTUM' = "#2c7fb8", 'GI_COLON' = "#a1dab4", 'FAT' = "#54278f", 'LIVER' = '#9e9ac8', 'VASCULAR' = '#f2f0f7', 'PANCREAS'= '#1c9099', 'STROMAL_CONNECTIVE' = '#bdc9e1', 'THYMUS' = "#f6eff7", 'PLACENTA'= '#3182bd', 'GI_DUODENUM' = "#636363", 'CERVIX'= '#cccccc', 'BONE'= "#5A2B2B", 'ADRENAL'= '#55A4F4', 'KIDNEY'= '#F4556C', 'MUSCLE_LEG' = '#9B7A7F', 'OVARY' = '#F2FA10', 'SPLEEN'= '#11FA10', 'GI_ESOPHAGUS' = '#10F7FA', 'NaN'= "#040404")
CELL_GROUP_COLOR_CODE = c('Adipose' = "#54278f", 'Blood & T-cell' = '#e34a33', 'Brain' = '#ffffd4', 'Digestive' = '#253494', 'ENCODE2012' = '#040404', 'Epithelial' = '#11FA10', 'ES-deriv' = '#fdbb84', 'ESC' = '#fef0d9', 'Heart' = '#fbb4b9', 'HSC & B-cell' = '#F2FA10', 'IMR90' = '#F4556C', 'iPSC' ='#7a0177', 'Mesench' = "#051C34", 'Muscle' = '#c2e699', 'Myosat' = '#9B7A7F', 'Neurosph' = '#10F7FA', 'Other' = '#006837', 'Sm. Muscle' = '#78c679', 'Thymus' = "#f6eff7")
##########################################################################
############ ALL THE FUNCTIONS ###########################################

get_enrichment_context_name <- function(roc_fn) {
	roc_fn <- unlist(strsplit(roc_fn, "/")) %>% last() # last file name in a string of full fill path, each component separated by '/'
	return(substr(roc_fn, 1, nchar(roc_fn) - 15)) # get rid of the '.bed.gz_roc.csv' tail
}

get_roc_df <- function(roc_fn){
	df <- read.table(roc_fn, sep = ',', header = FALSE, fill = TRUE)
	rownames(df) <- df$V1
	df <- df %>% select(c(-'V1')) # first column is already assigned the rownames of the df
	df <- as.data.frame(t(df))
	return(df)
}

get_max_auc_df <- function(auc_fn){
	auc_df <- read.table(auc_fn, sep = "\t", header = FALSE) 
	colnames(auc_df) <- c('ct', 'auc')
	full_auc <- auc_df[1,'auc']
	auc_df <- auc_df[-1,] # get rid of the first row, which is the full stack model auc. We can skip it because we will draw it anyway.
	auc_df <- auc_df %>% dplyr::left_join(metadata, by = c("ct" = "CT_NAME"))
	max_auc_df <- auc_df %>% group_by(ANATOMY) %>% filter(auc == max(auc)) #For each group of tissue type, get the cell type whose model yields the highest auc
	return(list('max_auc_df' = max_auc_df, 'full_auc' = full_auc))
}

get_bound_auc_df <- function(auc_fn){
  auc_df <- read.table(auc_fn, sep = "\t", header = FALSE) 
  colnames(auc_df) <- c('ct', 'auc')
  full_auc <- auc_df[1,'auc']
  auc_df <- auc_df[-1,] # get rid of the first row, which is the full stack model auc. We can skip it because we will draw it anyway.
  auc_df <- auc_df %>% dplyr::left_join(metadata, by = c("ct" = "CT_NAME"))
  result_df <- as.data.frame(matrix(nrow = 0, ncol = 4))
  colnames(result_df) <- c('ct', 'auc', 'auc_type', 'ANATOMY')
  min_df <- auc_df %>% filter(auc == min(auc)) %>% select(c('ct', 'auc', 'ANATOMY'))
  min_df['auc_type'] <- c('min_auc')
  max_df <- auc_df %>% filter(auc == max(auc)) %>% select(c('ct', 'auc', 'ANATOMY'))
  max_df['auc_type'] <- c('max_auc')
  med_df <- auc_df %>% arrange(abs(auc - median(auc))) %>% slice(1) %>% select(c('ct', 'auc', 'ANATOMY'))
  med_df['auc_type'] <- c('med_auc')
  result_df <- bind_rows(min_df, max_df, med_df)
  return(list('sum_df' = result_df, 'full_auc' = full_auc)) # ct, auc, ANATOMY, auc_type (min_auc, med_auc, max_auc)
}

draw_auc_plot_rep_tissue_model <- function(roc_fn, auc_fn, save_fn){
	roc_df <- get_roc_df(roc_fn) # full_false_pos, full_true_pos, <ct>_false_pos, <ct>_true_pos
	auc_data <- get_max_auc_df(auc_fn) # ct, auc, GROUP, COLOR, TYPE, Epig_name, ANATOMY
	max_auc_df <- auc_data$max_auc_df
	full_auc <- auc_data$full_auc
	uniq_tissue_list <- sort(unique(max_auc_df$ANATOMY)) # unique tissue types
	# first draw the full stack data fp and tp
	full_df <- data.frame('full_false_pos' = roc_df$full_false_pos, 'full_true_pos' = roc_df$full_true_pos, 'ana' = rep('full', length(roc_df$full_false_pos)))
	full_auc_label <- paste("full_auc ==" , full_auc)
	p <- ggplot() + 
	geom_line(data = full_df, aes(x=full_false_pos, y = full_true_pos, color  = ana)) + 
	scale_color_manual(values = ANATOMY_COLOR_CODE)+
	theme_bw() +
	annotate('text', x = 0.1, y = 0.9, label = full_auc_label, parse = TRUE)

	for (ana in uniq_tissue_list){
		this_ana_max_auc_ct <- (max_auc_df %>% filter(ANATOMY == ana) %>% select('ct'))
		this_ana_max_auc_ct <- (this_ana_max_auc_ct$ct)[1]
		fp_colname <- paste0(this_ana_max_auc_ct, '_false_pos')
		tp_colname <- paste0(this_ana_max_auc_ct, '_true_pos')
		num_state_this_ct_model <- length(roc_df[fp_colname])
		this_ct_df <- data.frame('fp' = roc_df[fp_colname], 'tp' = roc_df[tp_colname], 'ana' = rep(ana, num_state_this_ct_model)) %>% drop_na()
		colnames(this_ct_df) <- c('fp', 'tp', 'ana')
		p <- p + 
		geom_line(data = this_ct_df, aes(x = fp, y = tp, color = ana), alpha = 0.3) +
		scale_colour_manual(values = ANATOMY_COLOR_CODE)
	}
	enrichment_context_name <- get_enrichment_context_name(roc_fn) # get the name of the enrichment context thatwe are looking at, so that we can write the title of the graph
	print(enrichment_context_name)
	p <- p +
	ggtitle(enrichment_context_name) +
	scale_x_continuous(name = 'false_pos') +
	scale_y_continuous(name = 'true_pos') + 
	theme(legend.position="bottom")
	ggsave(save_fn)
}

get_ct_names_from_header <- function(header){
	# E003_false_pos --> E003
	# return(unlist(strsplit(header, '_'))[1])
  return (unlist(strsplit(header, '\\_false_|\\_true_| '))[1])
}

draw_auc_plot_bounds_tissue_model <- function(roc_fn, auc_fn, save_fn){
  roc_df <- get_roc_df(roc_fn) # full_false_pos, full_true_pos, <ct>_false_pos, <ct>_true_pos
  auc_data <- get_bound_auc_df(auc_fn) # ct, auc, GROUP, COLOR, TYPE, Epig_name, ANATOMY
  sum_auc_df <- auc_data$sum_df
  print(sum_auc_df)
  full_auc <- auc_data$full_auc
  # first draw the full stack data fp and tp
  full_df <- data.frame('full_false_pos' = roc_df$full_false_pos, 'full_true_pos' = roc_df$full_true_pos, 'ana' = rep('full', length(roc_df$full_false_pos)))
  full_auc_label <- paste("full_auc ==" , full_auc)
  p <- ggplot() + 
    geom_line(data = full_df, aes(x=full_false_pos, y = full_true_pos, color  = 'full')) + 
    scale_color_manual(values = c('full' = 'grey0')) +
    theme_bw() +
    annotate('text', x = 0.2, y = 0.9, label = full_auc_label, parse = TRUE) 
  min_auc_ct <- sum_auc_df %>% filter(auc_type == 'min_auc') %>% select('ct')
  max_auc_ct <- sum_auc_df %>% filter(auc_type == 'max_auc') %>% select('ct')
  med_auc_ct <- sum_auc_df %>% filter(auc_type == 'med_auc') %>% select('ct')
  plot_title <- get_enrichment_context_name(save_fn)
  print(min_auc_ct)
  p <- p +
    geom_line(data = roc_df, aes_string(x = paste0(min_auc_ct, '_false_pos'), y = paste0(min_auc_ct, '_true_pos'), na.rm = TRUE, color = shQuote("min_ct_model_auc")), alpha = 0.7) + 
    geom_line(data = roc_df, aes_string(x = paste0(max_auc_ct, '_false_pos'), y = paste0(max_auc_ct, '_true_pos'), na.rm = TRUE, color = shQuote('max_ct_model_auc')), alpha = 0.7) + 
    geom_line(data = roc_df, aes_string(x = paste0(med_auc_ct, '_false_pos'), y = paste0(med_auc_ct, '_true_pos'), na.rm = TRUE, color = shQuote('med_ct_model_auc')), alpha = 0.8) + 
    scale_color_manual(values = c("full" = 'grey0', "min_ct_model_auc" = 'red', "max_ct_model_auc" = 'blue', "med_ct_model_auc" = 'green')) +
    theme(legend.position = 'bottom') +
    ggtitle(plot_title)
  ggsave(save_fn)
  return(p)
}

draw_roc_plot_cell_type_blurred <- function(roc_fn, auc_fn, save_fn){
	roc_df <- get_roc_df(roc_fn) # full_false_pos, full_true_pos, <ct>_false_pos, <ct>_true_pos
	plot_title <- get_enrichment_context_name(save_fn)
	auc_data <- get_bound_auc_df(auc_fn) # ct, auc, GROUP, COLOR, TYPE, Epig_name, ANATOMY
	full_auc <- auc_data$full_auc

	# first draw the full stack data fp and tp
	full_df <- data.frame('full_false_pos' = roc_df$full_false_pos, 'full_true_pos' = roc_df$full_true_pos, 'ana' = rep('full', length(roc_df$full_false_pos)))
	full_auc_label <- paste("full_auc ==" , full_auc)
	p <- ggplot() + 
	geom_line(data = full_df, aes(x=full_false_pos, y = full_true_pos, color  = 'full')) + 
	scale_color_manual(values = c('full' = 'grey0')) +
	theme_bw() +
	annotate('text', x = 0.2, y = 0.9, label = full_auc_label, parse = TRUE) 
	ct_list <- sapply(colnames(roc_df), FUN = get_ct_names_from_header) 
	ct_list <- unique(ct_list)[-1] # get the list of cell types that we have in the roc data file, get rid of the first element because it's 'full'
	for (ct in ct_list){
		p <- p + 
		geom_line(data = roc_df, aes_string(x = paste0(ct, '_false_pos'), y = paste0(ct, '_true_pos'), na.rm = TRUE), color = "blue", alpha = 0.1)
	}
	p <- p + 
	theme(legend.position = 'None')+
    ggtitle(plot_title)
    ggsave(save_fn)
    return(p)
}

draw_roc_plot_each_model_one_color <- function(roc_fn, auc_fn, save_fn){
	roc_df <- get_roc_df(roc_fn) # full_false_pos, full_true_pos, <ct>_false_pos, <ct>_true_pos
	plot_title <- get_enrichment_context_name(save_fn)
	auc_data <- get_bound_auc_df(auc_fn) # ct, auc, GROUP, COLOR, TYPE, Epig_name, ANATOMY
	ct_list <- sapply(colnames(roc_df), FUN = get_ct_names_from_header) 
	color_palette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")
	plot_df <- data.frame(false_pos = numeric(), true_pos = numeric(), model = character())
	for (ct in ct_list){
		this_ct_df <- roc_df %>% select(c(paste0(ct, '_false_pos'), paste0(ct, '_true_pos')))
		colnames(this_ct_df) <- c('false_pos', 'true_pos')
		this_ct_df['model'] <- rep(ct, nrow(this_ct_df))
		plot_df <- rbind(plot_df, this_ct_df)
	}
	p <- ggplot()  + theme_bw() +
		geom_line(data = plot_df, aes(x = false_pos, y = true_pos, na.rm = TRUE, color = model)) +
		theme(legend.position = 'bottom') +
		ggtitle(plot_title)
	ggsave(save_fn)
}
# this program will takes as input the file that has auc data of all cell type specific model's auc and the full stack model's auc. It also take in roc_fn: the file that contains information about the true and false positive rates of different models. 
# It will draw a plots of roc lines: full stack model will be black line. And all other tissue types' models will be represented by the cell type of that tissue that has the most auc. The plot will be stored in save_fn
#!/usr/bin/env Rscript
################ GETTING THE COMMAND LINE ARGUMENTS #####################
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3)
{
        stop("wrong command line argument format", call.=FALSE)
}
roc_fn <- args[1] # roc_fn is the output of /u/home/h/havu73/project-ernst/source/model_evaluation/create_roc_curve_overlap_enrichment.py
auc_fn <- args[2] # the output of /u/home/h/havu73/project-ernst/source/model_evaluation/create_roc_curve_overlap_enrichment.py
save_fn <- args[3] # where the figure will be saved


################ CALLING FUNCTIONS FOR THIS PROGRAM ######################
# draw_roc_plot_each_model_one_color(roc_fn, auc_fn, save_fn)
draw_roc_plot_cell_type_blurred(roc_fn, auc_fn, save_fn)
##########################################################################

