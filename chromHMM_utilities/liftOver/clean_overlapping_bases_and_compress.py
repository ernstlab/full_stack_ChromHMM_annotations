import helper 
import numpy as np 
import pandas as pd 
import os 
import argparse
import pybedtools 
import compress_segmentation_data as compress
def convert_index_to_consecutive_segment(idx_series, chrom):
	'''
	If input is pd.Series([1, 2, 3, 5, 6, 7, 10, 11, 15]), this function shoudl return 
	start  end
0      1    3
1      5    7
2     10   11
3     15   15
	Function was written by ChatGPT and verified by Ha
	'''
	idx_series = pd.Series(idx_series)
	diffs = idx_series.diff() > 1
	# Step 2: Label groups
	group_labels = diffs.cumsum()
	# Step 3: Aggregate groups to find start and end
	grouped = idx_series.groupby(group_labels).agg(['min', 'max']).reset_index(drop=True)
	grouped.columns = ['start', 'end']		
	grouped['end'] = grouped['end'] + 1
	# because bedtools index is 0-based, half-open interval [start,end) while pandas series indexing is closed interval [start,end], this line is needed to convert the pandas series idx system to bedtools's idx system
	# this line is basically trying to convert from the pandas series indexing system to the bedtools indexing system
	grouped['chrom'] = chrom
	grouped = grouped[['chrom', 'start', 'end']]
	return grouped 


def clean_and_compress_segment(segment_bed_fn, output_fn):
	'''
	Assumption is that the segment_bed_fn is from only one chromosome
	This file will try to find the uniquely mapped bases in the end assembly, and compress the resulting segmentation
	'''

	df = pd.read_csv(segment_bed_fn, header = None, sep = '\t', index_col = None)
	df.columns = ['chrom', 'start', 'end', 'state']
	chrom = np.unique(df['chrom'])
	assert len(chrom) == 1, 'Function clean_segment_file assumes segment_bed_fn is from only one chromosome, which is not the case'
	chrom = chrom[0]
	df = df.sort_values('start').reset_index(drop=True)  # may be redudant, but we want to sort the start coordinates 
	chromStart = df.loc[0, 'start']
	chromEnd = df.loc[df.index[-1], 'end']
	counter = pd.Series(0, index = range(chromStart, chromEnd), dtype=int)
	for idx, row in df.iterrows():
		try:
			counter.loc[row['start']:(row['end']-1)] += 1
			# because bedtools index is 0-based, half-open interval [start,end) while pandas series indexing is closed interval [start,end], this line is needed to make sure that the endpoints do not get counted twice
		except: 
			print(row['start'], row['end'])
	# now, filter for regions where the count is 1, meaning they are uniquely mapped region
	uniq_counter = counter[counter==1]
	# now, put consecutive regions that are uniquely mapped into one line in a bed frame
	uniq_df = convert_index_to_consecutive_segment(uniq_counter.index, chrom) # columns chrom, start, end
	# number of bp spanning this uniq_df is exactly the number of bp in uniq_counter
	# now we will intersect the uniq-mapped region with the chromatin state assignments of the original 200-bp mappings
	rawbed = pybedtools.BedTool.from_dataframe(df)
	uniqbed = pybedtools.BedTool.from_dataframe(uniq_df)
	intersection = rawbed.intersect(uniqbed)
	intersection = intersection.to_dataframe()
	intersection.rename(columns = {'name': 'state'}, inplace=True)
	# this has passed the test that the number of bps in intersection is exactly the number of bp in uniq_counter
	# now we will try to convert the intersection into more condensed format such that consecutive region with the same state annotation will be condensed into one line
	intersection = compress.get_compressed_state_segment_one_chrom_auto(intersection, chrom) 
	# this function has passed the test that the number of bp in the compressed intersection is exactly the same as the number of bp in uniqbed
	intersection.to_csv(output_fn, header = False, index = False, sep = '\t', compression ='gzip')
	return





	


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'This code will take one single chromosome\'s segmentation file that has just been lifted from one assembly to another, and the code will figure out if there are any bases that are mapped multiple times regardless of whether those correspond to the same state. This code will get rid of all those overlapping region, and return a cleaned-up bed file')
	parser.add_argument('--segment_bed_fn', type = str, required = True,
		help = 'input segment_bed_fn')
	parser.add_argument('--output_fn', type = str, required = True, help = 'output_fn')
	args = parser.parse_args()
	helper.check_file_exist(args.segment_bed_fn)
	helper.create_folder_for_file(args.output_fn)
	print('clean_overlapping_bases_and_compress.py: Done getting command line argument')
	clean_and_compress_segment(args.segment_bed_fn, args.output_fn)
	print('Done!')
