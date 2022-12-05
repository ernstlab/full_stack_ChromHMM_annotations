import pandas as pd 
import helper 
import numpy as np 
import os 
import argparse
parser = argparse.ArgumentParser(description = 'This code takes in a segmentation file and reaarranges the rows such that if consecutive genomic bins are annotated as the same state, they will be combined. At the same time, the output file will be sorted by chrom, start_bp such that it will be readily available to use for downstream analysis with bedtools or pybedtools')
parser.add_argument('--segment_bed_fn', type = str, required = True,
	help = 'input segment_bed_fn')
parser.add_argument('--output_fn', type = str, required = True, help = 'output_fn')
args = parser.parse_args()
print(args)
helper.check_file_exist(args.segment_bed_fn)
helper.create_folder_for_file(args.output_fn)

def get_compressed_state_segment_one_chrom_manual(chrom_segment_df, chrom):
	# combine multiple rows with the same state into a row so that we can save space on the computer
	# chrom_segment_df has columns 'chrom', 'start', 'end', 'state'
	chrom_segment_df = chrom_segment_df.reset_index(drop = True)
	chrom_segment_df = chrom_segment_df.sort_values(['start']) # sort by ascending start
	result_df = pd.DataFrame(columns = ['chrom', 'start', 'end', 'state'])
	current_start_bin_index = 0
	current_start = chrom_segment_df.loc[current_start_bin_index, 'start']
	current_end = chrom_segment_df.loc[current_start_bin_index, 'end']
	current_state = chrom_segment_df.loc[current_start_bin_index, 'state']
	for index in range(1, chrom_segment_df.shape[0]): # skip the first one because it's already reported into the current data
		if (current_state == chrom_segment_df.loc[index, 'state']) and current_end == chrom_segment_df.loc[index, 'start']: 
			# if in the same state and the previous row and the current row are continous segments on the genome (no gap)
			current_end = chrom_segment_df.loc[index, 'end']
		else: # change of state or there is a gap in the segment --> report into the result_df
			add_row = [chrom, current_start, current_end, current_state]
			result_df.loc[result_df.shape[0]] = add_row
			current_state = chrom_segment_df.loc[index, 'state']
			current_start = chrom_segment_df.loc[index, 'start']
			current_end = chrom_segment_df.loc[index, 'end']
	add_row = [chrom, current_start, current_end, current_state]
	result_df.loc[result_df.shape[0]] = add_row
	# now we are done producing
	return result_df

def get_compressed_state_segment_one_chrom_auto(chrom_df, chrom):
	chrom_df = chrom_df.sort_values('start')
	chrom_df['state'] = chrom_df['state'].apply(lambda x: int(x[1:])) # convert E1 --> 1
	# now we will calculate the distance between each line's end and the following line's start
	chrom_df.loc[:, 'next_start'] = [np.NaN] + list(chrom_df['end'][1:]) 
	chrom_df.loc[:, 'diff_with_next_row'] = chrom_df['next_start'] - chrom_df['end']
	# now we will merge the rows of df such that the consecutive rows with the same value for state, and that they are consecutive in terms of genomic coordinates will be grouped into one row. References: https://stackoverflow.com/questions/26911851/how-to-use-pandas-to-find-consecutive-same-data-in-time-series
	chrom_df.loc[:, 'start_grp'] = (chrom_df['diff_with_next_row'] != 0).cumsum()
	chrom_df.loc[:, 'state_grp'] = (chrom_df['state'].diff(periods=1) != 0).astype('int').cumsum() 
	group_grp = chrom_df.groupby(['start_grp', 'state_grp'])
	result_df = pd.DataFrame({'start': group_grp.start.first(), 'end': group_grp.end.last(), 'state': group_grp.state.first()})
	result_df.loc[:, 'chrom'] = chrom
	result_df['state'] = result_df['state'].apply(lambda x: 'E{}'.format(x))
	result_df = result_df [['chrom', 'start', 'end', 'state']]
	return result_df 

def compress_segmentation(segment_bed_fn, output_fn):
	# segment_df should have 4 columns: chrom, start, end, state
	segment_df = pd.read_csv(segment_bed_fn, header = None, sep = '\t', index_col = None)
	segment_df.columns = ['chrom', 'start', 'end', 'state']
	print('Done reading in data from segment_bed_fn')
	group_df = segment_df.groupby('chrom')
	result_df_list = []
	for chrom, chrom_segment_df in segment_df.groupby('chrom'):
		compressed_chrom_df = get_compressed_state_segment_one_chrom_auto(chrom_segment_df, chrom)
		result_df_list.append(compressed_chrom_df)
		print('Done for chromosome: {}'.format(chrom))
	result_df = pd.concat(result_df_list, ignore_index = True)
	result_df.to_csv(output_fn, header = False, index = False, sep = '\t', compression = 'gzip')
	return 


compress_segmentation(args.segment_bed_fn, args.output_fn)