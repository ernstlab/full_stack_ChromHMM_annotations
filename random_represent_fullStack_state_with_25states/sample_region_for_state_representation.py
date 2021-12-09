# this file will look into the genome wide segmentation file of full stack state, and then it will sample, for each state, an equal number of 200bp regions that belong to that state. In the end, we get a bed files showing the coordinate of all those regions. 
import os 
import sys
import pandas as pd 
import helper
import numpy as np 
NUM_STATE = 100

def get_one_segment_from_merged_segment(selected_state_df): 
	"""
	selected_state_df : comtain rows that have num_bins > 1. We will pick one bin out of the multiple bins that are contained in each row. R
	"""
	selected_state_df['picked_segment_index'] = selected_state_df.apply(lambda x: (np.random.choice(range(x['num_bins']), 1))[0], axis = 1) # for each row, pick a number from 0 --> num_bins -1 where we present as the picked bin for this segment
	selected_state_df['start'] = selected_state_df['start'] + selected_state_df['picked_segment_index'] * helper.NUM_BP_PER_BIN
	selected_state_df['start'] = pd.to_numeric(selected_state_df['start'], downcast = 'integer')
	selected_state_df['end'] = selected_state_df['start'] + helper.NUM_BP_PER_BIN
	return selected_state_df

def sample_regions_equal_state_segmentation(semgent_fn, num_segment_per_state, output_fn, seed):
	segment_df = pd.read_csv(semgent_fn, header = None, sep = '\t', index_col = None)
	segment_df.columns = ['chrom', 'start', 'end', 'state']
	segment_df = segment_df[~segment_df.chrom.isin(['chrY', 'chrM'])] # avoid these two chromosomes because they are not annotated in cell-type-spec 25-state models
	result_df = pd.DataFrame(columns = segment_df.columns)
	segment_df = segment_df.groupby('state')
	for state, state_df in segment_df: # loop through each group
		state_df = state_df.reset_index(drop = True)
		np.random.seed(seed)
		row_indices_to_sample = np.random.choice(state_df.index, num_segment_per_state, replace = False)
		state_df = state_df.loc[row_indices_to_sample] # pick the places that we
		state_df = state_df.reset_index(drop = True)
		state_df['num_bins'] = (state_df['end'] - state_df['start']) / helper.NUM_BP_PER_BIN
		state_df['num_bins'] = pd.to_numeric(state_df['num_bins'], downcast = 'integer')
		state_df = get_one_segment_from_merged_segment(state_df) # this is now fixed such that the each row only contains one bin representing the segment 
		state_df = state_df[['chrom', 'start', 'end', 'state']]
		result_df = result_df.append(state_df)
	result_df.to_csv(output_fn, header = False, index = False, sep = '\t', compression = 'gzip')
	return result_df

def main():
	if len(sys.argv) != 5:
		usage()
	semgent_fn = sys.argv[1]
	helper.check_file_exist(semgent_fn)
	num_segment_per_state = helper.get_command_line_integer(sys.argv[2])
	output_fn = sys.argv[3]
	helper.create_folder_for_file(output_fn)
	seed = helper.get_command_line_integer(sys.argv[4])
	print("Done getting command line arguments")
	sample_regions_equal_state_segmentation(semgent_fn, num_segment_per_state, output_fn, seed)
	print ('Done!')

def usage():
	print("python sample_regions_for_state_representation.py")
	print('semgent_fn: the genome-wide segmentation file as input')
	print('num_segment_per_state: number of segments per state that we want to sample')
	print('output_fn: where we store the data of full-stack segmentation that we would like to sample from. 4 columns: 3 bed file columns and full-stack state')
	print('seed: the seed for random number generator')
	exit(1)
main()