import sys
import os
import numpy as np 
import pandas as pd 
import replicate_helper as helper
import glob

def get_data_all_ct(average_data_folder):
	ct_fn_list = glob.glob(average_data_folder + "/E*_avg_exp_tss_relative.gz") # list of cell types that we have data of gExp around the TSS
	ct_list = list(map(lambda x: (x.split('/')[-1]).split('_')[0], ct_fn_list))
	ct_df_list = list(map(lambda x: pd.read_csv(x, header = 0, sep = '\t'), ct_fn_list))
	return ct_df_list

def main():
	if len(sys.argv) != 3:
		usage()
	average_data_folder = sys.argv[1]
	helper.check_dir_exist(average_data_folder)
	output_fn = sys.argv[2]
	helper.create_folder_for_file(output_fn)
	print ("Done getting command line argument")
	ct_df_list = get_data_all_ct(average_data_folder)
	# now average over all the dataframes
	avg_df = pd.concat(ct_df_list).groupby(level=0).mean() # columns: states, rows: genomic bins relative to the tss
	avg_df.to_csv(output_fn, header = True, index = False, sep = '\t', compression = 'gzip')
	print("Done!")

def usage():
	print ("python average_over_ct_gExp_tss_relative.py")
	print ("average_data_folder: where the data are stored for each cell type in the file name <ct>_avg_exp_tss_relative.gz")
	print ("output_fn: where data that we average over all the cell types are stored")
	exit(1)
main()