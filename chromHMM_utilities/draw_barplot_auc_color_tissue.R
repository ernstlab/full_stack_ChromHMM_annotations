library("tidyverse")
library("tidyr")
library("dplyr")
library('ggplot2')
library('plotly')
metadata_fn <- "/Users/vuthaiha/Desktop/window_hoff/data/roadmap_epigenome/roadmap_celltype_metadata.csv"
metadata_df <- as.data.frame(read.csv(metadata_fn, sep = ","))
metadata_df <- metadata_df %>% select("Epigenome.ID..EID.", "GROUP", "COLOR", "TYPE", "Epigenome.name..from.EDACC.Release.9.directory.", "ANATOMY", "TYPE") %>% rename("CT_NAME" = "Epigenome.ID..EID.") %>% rename("Epig_name" = "Epigenome.name..from.EDACC.Release.9.directory.")
ANATOMY_COLOR_CODE = c('BLOOD'= '#e34a33', 'ESC_DERIVED'= '#fdbb84', 'ESC'= '#fef0d9', 'BRAIN'= '#ffffd4', 'LUNG'= '#006837', 'SKIN'= '#78c679', 'MUSCLE'= '#c2e699', 'IPSC'= '#7a0177', 'GI_STOMACH' = '#f768a1', 'HEART'= '#fbb4b9', 'BREAST' = '#d7b5d8', 'GI_INTESTINE' = '#253494', 'GI_RECTUM' = "#2c7fb8", 'GI_COLON' = "#a1dab4", 'FAT' = "#54278f", 'LIVER' = '#9e9ac8', 'VASCULAR' = '#f2f0f7', 'PANCREAS'= '#1c9099', 'STROMAL_CONNECTIVE' = '#bdc9e1', 'THYMUS' = "#f6eff7", 'PLACENTA'= '#3182bd', 'GI_DUODENUM' = "#636363", 'CERVIX'= '#cccccc', 'BONE'= "#051C34", 'ADRENAL'= '#55A4F4', 'KIDNEY'= '#F4556C', 'MUSCLE_LEG' = '#9B7A7F', 'OVARY' = '#F2FA10', 'SPLEEN'= '#11FA10', 'GI_ESOPHAGUS' = '#10F7FA', 'NaN'= "#040404")

draw_auc_compare_models_one_context <- function(auc_fn, save_fn, enrichment_name){
  auc_df <- as.data.frame(read.table(auc_fn, sep = "\t"))
  auc_df <- auc_df %>% rename("model_name" = "V1", "auc" = "V2")
  full_model_auc <- auc_df[1,2]
  auc_df <- auc_df[2:nrow(auc_df),]
  auc_df <- inner_join(auc_df, metadata_df, by = c("model_name"= "CT_NAME"))
  uniq_anatomy <- sort(unique(auc_df$ANATOMY))
  plot_df <- data.frame(matrix(ncol = ncol(auc_df), nrow = 0))
  colnames(plot_df) <-  colnames(auc_df)
  for (anatomy in uniq_anatomy){ # get all the cell types model of the same anatomy to be placed next to each other
    this_ana_df <- auc_df %>% filter(ANATOMY == anatomy)
    plot_df <- bind_rows(plot_df, this_ana_df)
  }
  # the following line is to tell R that it needs to repsect the order of x axies as it is in the plot_df
  plot_df$model_name <- factor(plot_df$model_name, levels = plot_df$model_name[order(plot_df$ANATOMY)])
  p <- ggplot(aes(x= model_name, y = auc, fill = ANATOMY), data = plot_df) +
    theme_bw() + 
    geom_bar(stat="identity") + 
    scale_fill_manual(name = "tissue_type", values = ANATOMY_COLOR_CODE) +
    geom_abline(slope=0, intercept=full_model_auc,  col = "black",lty=2) +
    coord_cartesian(ylim=c(0.5, 1))  + 
    labs(fill = 'tissue_type') + 
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.text = element_text(size = 5), legend.title = element_blank(), legend.position = 'bottom') +
    ggtitle(enrichment_name)
  gg <- ggplotly(p)
  ggsave(filename = save_fn)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3)
{
        stop("wrong command line argument format", call.=FALSE)
}
auc_fn <- args[1] # output of ChromHMM OverlapEnrichment, also input to this program
save_fn <- args[2] # where the figures should be stored
enrichment_name <- args[3] # name of the enrichment context, so that it can be title of the figure
draw_auc_compare_models_one_context(auc_fn, save_fn, enrichment_name)