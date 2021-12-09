# after we already have an input file that shows the full stack states and the corresponding 25-state states in 127 cell types, we will count, for each full stack state, the number (out of num_segment_per_state from sample_region_for_state_representation.py) of occurrences of cell type specific states that occur in the cell type specific models that overlap with regions associated with the full-stack states.
import os 
import sys
import pandas as pd 
import helper
import numpy as np 
NUM_STATE = 100
NUM_CELL_TYPE_SPEC_STATE = 25 # because this fiel is inside the folder associated with 25-state cell type specific system. 
NUM_NON_CT_COLUMNS = 4 # the first NUM_NON_CT_COLUMNS in bed_map_df are chrom, start, end, full_stack
state_annot_fn = '/u/home/h/havu73/project-ernst/data/roadmap_epigenome/25_imputed12marks_model_downloaded/state_annot.txt'
def read_ct_annot_df():
	roadmap_ct_annot_fn = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/ROADMAP_metadata_july2013.csv'
	ct_annot_df = pd.read_csv(roadmap_ct_annot_fn, header = 0, index_col = None, sep = ',')
	ct_annot_df = ct_annot_df.rename(columns = {'Epigenome ID (EID)': 'eid'})
	ct_annot_df = ct_annot_df[['eid', 'GROUP', 'COLOR', 'ANATOMY']] 
	ct_annot_df = ct_annot_df.sort_values('GROUP') # sort by gruups so that we can later arrange the data in specific order
	return ct_annot_df

def output_values_sorted_by_ct(output_fn, bed_map_df, num_columns_results):
	# bed_map_df is the grouped bed_map_df from count_sample_region_each_full_stack_state
	result_df = pd.DataFrame(columns = range(num_columns_results))
	result_colnames = ['full_stack_state']
	for state, state_df in bed_map_df: # loop through each group of full-stack state
		row_result = [state] # first element is the full_stack state
		state_df = state_df[state_df.columns[NUM_NON_CT_COLUMNS:]] # rid of the first NUM_NON_CT_COLUMNS columns (chrom, start, end, full_stack)
		state_df = state_df.applymap(lambda x: int(x[1:])) # convert state annotation E1 --> 1
		for ct in state_df.columns: # each ct that we already arranged in the right order
			this_ct_state_list = np.sort(state_df[ct]) # ascending sort of the state numbers in the 100 places that we sampled for this full_stack state
			row_result += list(this_ct_state_list) # append these cells into the row, each group fo cells appended corresponds to a sorted list of state numbers in the ct-spec models
			if state == 'E1': # only need to add this once
				result_colnames += list(map(lambda x: ct + '_' + str(x), range(state_df.shape[0]))) # ct_0 --> ct_99
		result_df.loc[result_df.shape[0]] = row_result
	result_df.columns = result_colnames
	result_df.to_csv(output_fn, header = True, index = False, sep = '\t', compression = 'gzip')

def get_dictionary_cellgroup_eid(ct_annot_df):
	ct_grouped_df = ct_annot_df.groupby('GROUP')
	result = {}
	for group, cg_df in ct_grouped_df:
		result[group] = list(cg_df.eid)
	return result # key: cell group, values, list of eid that are in this group

def output_values_sorted_by_group(output_fn, bed_map_df, ct_annot_df, num_columns_results): 
	# bed_map_df is the grouped bed_map_df from count_sample_region_each_full_stack_state
	result_df = pd.DataFrame(columns = range(num_columns_results))
	result_colnames = ['full_stack_state']
	cellgroup_eid_dict = get_dictionary_cellgroup_eid(ct_annot_df) # key: cell group, values, list of eid that are in this group
	print (cellgroup_eid_dict)
	cell_group_list = np.unique(ct_annot_df.GROUP)
	for state, state_df in bed_map_df: # loop through each group of full-stack state
		row_result = [state] # first element is the full_stack state
		state_df = state_df[state_df.columns[NUM_NON_CT_COLUMNS:]] # rid of the first NUM_NON_CT_COLUMNS columns (chrom, start, end, full_stack)
		state_df = state_df.applymap(lambda x: int(x[1:])) # convert state annotation E1 --> 1
		for group in cell_group_list:
			this_group_state_list = []
			for ct in cellgroup_eid_dict[group]:
				this_group_state_list += list(state_df[ct])
			this_group_state_list = np.sort(this_group_state_list) # sort the state values that are in the cell types in each group. This list has length: num_ct_in_this_group * 100
			row_result += list(this_group_state_list)
			if state == 'E1': # we only need to add this once
				result_colnames += list(map(lambda x: group + "_"+ str(x), range(len(this_group_state_list)))) # <group>_1 --> <group>_<final index>
		result_df.loc[result_df.shape[0]] = row_result
	result_df.columns = result_colnames
	result_df.to_csv(output_fn, header = True, index = False, sep = '\t', compression = 'gzip')

def count_sample_region_each_full_stack_state(bed_map_fn, output_folder):
	ct_annot_df = read_ct_annot_df() # eid, GROUP, COLOR, ANATOMY
	bed_map_df = pd.read_csv(bed_map_fn, header = 0, index_col = None, sep = '\t')
	# rearrange the columns such that cell types of the same group are together
	num_segment_per_state = bed_map_df.shape[0] / NUM_STATE
	# chrom, start, end, full_stack, E001, E002, etc. 
	input_ct_list = bed_map_df.columns[NUM_NON_CT_COLUMNS:] # list of 127 ct 
	ct_annot_df = ct_annot_df[ct_annot_df['eid'].isin(input_ct_list)] # filter out rows that are not in the list of input ct
	bed_map_df = bed_map_df[list(bed_map_df.columns[:NUM_NON_CT_COLUMNS]) + list(ct_annot_df['eid'])] # rearrange columns
	bed_map_df = bed_map_df.groupby('full_stack')
	num_ct_total = ct_annot_df.shape[0] # 127
	num_columns_results = int(1 + num_ct_total * num_segment_per_state) # number of fields in the result dataframe
	ct_grouped_fn = os.path.join(output_folder, 'sample_state_sorted_by_ct.txt.gz')
	output_values_sorted_by_ct(ct_grouped_fn, bed_map_df, num_columns_results)
	group_grouped_fn = os.path.join(output_folder, 'sample_state_sorted_by_group.txt.gz')
	output_values_sorted_by_group(group_grouped_fn, bed_map_df, ct_annot_df, num_columns_results)
	return 

def main():
	if len(sys.argv) != 3:
		usage()
	bed_map_fn = sys.argv[1]
	helper.check_file_exist(bed_map_fn)
	output_folder = sys.argv[2]
	helper.check_dir_exist(output_folder)
	print('Done getting command line arguments')
	count_sample_region_each_full_stack_state(bed_map_fn, output_folder)
	print('Done!')
	return 

def usage():
	print ('python count_sample_region_each_full_stack_state.py')
	print ('bed_map_fn: output of the step where we do a bunch of bedtools map, each map correspoidng to one cell type')
	print ('output_folder: where we save the data of the counts. Two files: sample_state_sorted_by_ct.txt.gz and sample_state_sorted_by_group.txt.gz')
	exit(1)
main()