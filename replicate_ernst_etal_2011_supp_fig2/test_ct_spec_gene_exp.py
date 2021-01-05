import replicate_helper as helper 
from scipy.stats import mannwhitneyu
import pandas as pd 
import numpy as np 
import sys
import os
NUM_FULL_STACK_STATES = 100

def get_metadata():
	meta_fn = '../../ROADMAP_aligned_reads/ROADMAP_metadata_july2013.csv'
	df = pd.read_csv(meta_fn, header = 0, index_col = None, sep = ',')
	df = df.rename(columns = {'Epigenome ID (EID)' : 'ct'})
	df = df[['ct', 'GROUP', 'ANATOMY']]
	return df

def get_state_annot_df():
	state_annot_fn = '../../ROADMAP_aligned_reads/chromHMM_model/model_100_state/figures/supp_excel/state_annotations_processed.csv'
	state_annot_df = pd.read_csv(state_annot_fn, sep = ',', header = 0)
	state_annot_df = state_annot_df[['state', 'color', 'mneumonics', 'state_order_by_group']]
	return state_annot_df	

def get_state_avg_ge_df(state_avg_gene_exp_fn):
	state_avg_ge_df = pd.read_csv(state_avg_gene_exp_fn, header = 0, index_col = None, sep = '\t')
	state_avg_ge_df['state'] = state_avg_ge_df['state'].apply(lambda x: 'S' + x[1:])
	state_avg_ge_df = state_avg_ge_df.set_index('state')
	state_avg_ge_df = state_avg_ge_df.transpose()
	state_avg_ge_df = state_avg_ge_df.reset_index()
	state_avg_ge_df = state_avg_ge_df.rename(columns = {'index' : 'ct'})
	meta_df = get_metadata()
	state_avg_ge_df = state_avg_ge_df.merge(meta_df, how = 'left', on = 'ct')
	return state_avg_ge_df

def test_ct_spect(state_avg_ge_df, group):
	group_df = state_avg_ge_df.groupby(group)
	results_df = pd.DataFrame(columns = group_df.groups.keys()) # columns are the groups of cell types that we detect in the data
	for g_i, g_name in enumerate(group_df.groups.keys()): # loop through each group of cell types
		not_groups = group_df.groups.keys()[:g_i] + group_df.groups.keys()[(g_i+1):] # names of cell groups that are not this group
		not_group_df = pd.concat(group_df.get_group(g) for g in not_groups) # get the df that has data for all the groups that are not the particular group that we are testing
		this_group_df = group_df.get_group(g_name) # get the data frame corresponding to this group
		this_group_result = pd.Series()
		for state in range(NUM_FULL_STACK_STATES):
			state_name = 'S' + str(state + 1) # S1 --> S100
			t = mannwhitneyu(this_group_df[state_name], not_group_df[state_name], use_continuity = False, alternative = 'greater')
			this_group_result[state_name] = t.pvalue 
		results_df[g_name] = this_group_result # append the p values correspond to this group to the results_df
	results_df = results_df.reset_index()
	results_df = results_df.rename(columns = {'index' : 'state'})
	results_df['state'] = results_df['state'].apply(lambda x: int(x[1:])) # convert S1 to 1 as an int
	total_test_num = NUM_FULL_STACK_STATES * len(group_df.groups.keys())
	return results_df, total_test_num

def color_significant_pval(pval_row, threshold):
	# pval_row: a rows of p-values
	white_color = '#ffffff' # white
	blue_color = '#85BCE5' # light blue
	red_color = '#FF7F7F' # light red
	results = pd.Series(['background-color: %s' % white_color for x in pval_row])
	results.index = pval_row.index
	# change colors to blue if below the thresholds
	below_threshold_indices = (pval_row <= threshold)
	results[below_threshold_indices] = 'background-color: %s' % blue_color
	results[pval_row ==  pval_row.min()] = 'background-color: %s' % red_color
	return results

def color_state_annotation(row_data, index_to_color):
	results = [""] * len(row_data) # current paint all the cells in the rows with nothing, no format yet
	state_annot_color = row_data['color']
	results[index_to_color] = 'background-color: %s' % state_annot_color # the third cell from the left is the state annotation cells
	return results

def paint_results_excel(test_df, test_num, output_fn, sheet_name):
	# test_df adn test_num are the output of test_ct_spect
	ALPHA = 0.01
	threshold = ALPHA / float(test_num)
	columns_to_paint_pval = test_df.columns[1:] # names of columns that we will paint for the p values calculations
	print(columns_to_paint_pval)
	state_annot_df = get_state_annot_df()
	test_df = test_df.merge(state_annot_df, how = 'left', left_on = 'state', right_on = 'state', suffixes = ('_x', '_y'))
	test_df = test_df.sort_values(by = 'state_order_by_group')
	colored_df = test_df.style.apply(lambda x: color_significant_pval(x, threshold), axis = 1, subset = pd.IndexSlice[:, columns_to_paint_pval]) #exclude coloring the first column which is state annotation
	mneumonics_index_in_row = test_df.columns.get_loc('mneumonics') # column index, zero-based of the menumonics entries, which we will use to paint the right column later
	colored_df = colored_df.apply(lambda x: color_state_annotation(x, mneumonics_index_in_row), axis = 1)
	return colored_df

def main():
	if len(sys.argv) != 3: 
		usage()
	state_avg_gene_exp_fn = sys.argv[1]
	helper.check_file_exist(state_avg_gene_exp_fn)
	output_fn = sys.argv[2]
	helper.create_folder_for_file(output_fn)
	print ("Done getting command line arguments")
	state_avg_ge_df = get_state_avg_ge_df(state_avg_gene_exp_fn)
	print ("Done getting input")
	group_test_df, group_test_num = test_ct_spect(state_avg_ge_df, 'GROUP')
	anatomy_test_df, ana_test_num = test_ct_spect(state_avg_ge_df, 'ANATOMY')
	print ("Done calculating p-values")
	colored_group_df = paint_results_excel(group_test_df, group_test_num, output_fn, 'GROUP')
	colored_ana_df = paint_results_excel(anatomy_test_df, ana_test_num, output_fn, 'ANATOMY')
	writer = pd.ExcelWriter(output_fn, engine = 'xlsxwriter')
	colored_group_df.to_excel(writer, sheet_name = 'GROUP')
	colored_ana_df.to_excel(writer, sheet_name = 'ANATOMY')
	writer.save() 
	print("Done paiting excel results")
	
def usage():
	print ("python test_ct_spec_gene_exp.py")
	print ("state_avg_gene_exp_fn: one average per state per cell type")
	print ("output_fn: results of excel with p-values and highlighted significant cells")
	exit(1)
main()