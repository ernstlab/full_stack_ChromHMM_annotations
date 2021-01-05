import pandas as pd 
import numpy as np 
import os
import sys
import glob
import chromHMM_utilities_common_functions_helper as helper

def get_summary_auc_statistics_full_vs_others(auc_folder, output_fn):
	auc_fn_list = glob.glob(auc_folder + "/auc_*.txt") # list of auc files for different contexts
	genome_context_list = list(map(lambda x: x.split('/')[-1].split('auc_')[-1].split('.txt')[0], auc_fn_list)) # get teh name of the genomic contexts corresponding to the 
	auc_df_list = list(map(lambda x: pd.read_csv(x, header = None, sep = '\t'), auc_fn_list)) # list of auc_df corresponding to the auc_fn_list. Two columns: 0, 1
	auc_df_list = list(map(lambda x: x.set_index(0), auc_df_list)) # set the first column, which is model names as the index. Now 1 column: 1
	full_auc_list = list(map(lambda x: x.loc['full'][1], auc_df_list)) # get the auc of the full-stack model for each of the genomic context
	auc_df_list = list(map(lambda x: x.drop(index = 'full'), auc_df_list)) # drop the line corresponding to the full-stack model auc. Now we are trying to get data corresponding to other models
	result_df = pd.DataFrame(columns = ['entity'] + genome_context_list)
	result_df.loc[0] = ['full_auc'] + full_auc_list
	result_df.loc[1] = ['mean_others'] + list(map(lambda x: np.mean(x[1]), auc_df_list))
	result_df.loc[2] = ['max_others'] + list(map(lambda x: np.max(x[1]), auc_df_list))
	result_df.loc[3] = ['min_others'] + list(map(lambda x: np.min(x[1]), auc_df_list))
	result_df.loc[4] = ['std_others'] + list(map(lambda x: np.std(x[1]), auc_df_list))
	result_df.to_csv(output_fn, header = True, index = False, sep = ',')

def main():
	if len(sys.argv) != 3:
		usage()
	auc_folder = sys.argv[1]
	helper.check_dir_exist(auc_folder)
	output_fn = sys.argv[2]
	helper.create_folder_for_file(output_fn)
	print ("Done checking command line arguments")
	get_summary_auc_statistics_full_vs_others(auc_folder, output_fn)
	print ("Done!")

def usage():
	print ("python create_summary_auc_full_others.py")
	print ("auc_folder: where all the files 'auc_<context>.txt' are stored")
	print ('output_fn: where summary data of full-stack vs. cell-type-specific models\'s AUC are stored')
	exit(1)

main()
