auc_folder=$1
echo $auc_folder
draw_code=/Users/vuthaiha/Desktop/window_hoff/source/chromHMM_utilities/create_roc_plots_compare_all_models.R

for auc_fn in $auc_folder/auc_*
do 
	context=$(echo $auc_fn | awk -F'/' '{print $NF}' | awk -F'auc_' '{print $2}' | awk -F'.' '{print $1}')
	save_fn=${auc_folder}/${context}_bounds_roc.png
	roc_fn=${auc_folder}/${context}_roc.csv
	Rscript ${draw_code} ${roc_fn} ${auc_fn} ${save_fn}
done