import pandas as pd 
import numpy as np 
import sys
import os
import chromHMM_utilities_common_functions_helper as helper
def get_total_num_segments_200bp (segment_fn):
	segment_df = pd.read_csv(segment_fn, header = None, sep = '\t')
	segment_df.columns = ['chrom', 'start', 'end', 'state']
	segment_df['num_bins'] = (segment_df.end - segment_df.start) / 200
	total_bins = np.sum(segment_df.num_bins)
	return total_bins, segment_df

def get_total_num_segments_from_liftOver (segment_fn):
	segment_df = pd.read_csv(segment_fn, header = None, sep = '\t')
	segment_df.columns = ['chrom', 'start', 'end', 'state']
	total_bins = int(segment_df.shape[0]) # num_rows is the same as the number of semgnets. This is based on the assumption that the segment_fn file is converted from a segment file that's 200bp
	return total_bins, segment_df

def convert_coordinates_to_int (input_fn, output_fn):
	df = pd.read_csv(input_fn, header = None, sep = '\t')
	df.columns = ['chr', 'start', 'end', 'state']
	df['start'] = df.start.astype('int64')
	df['end'] = df.end.astype('int64')
	df.to_csv(output_fn, sep = '\t', header = False, index = False, compression = 'gzip')
	print ("Done converting data from " + input_fn + ' to ' + output_fn)

def get_subset_rows_and_save(big_df, indices, save_fn):
	df = (big_df.loc[indices])[['chrom', 'start', 'end', 'state']] 
	df['start'] = df.start.astype('int64')
	df['end'] = df.end.astype('int64')
	df.to_csv(save_fn, compression = 'gzip', index = False, header = False, sep = '\t')
	print ("Done saving files " + save_fn)

def main():
	if len(sys.argv) != 4:
		usage()
	segment_fn = sys.argv[1]
	helper.check_file_exist(segment_fn)
	fract_gene_to_sample = float(sys.argv[2])
	output_folder = sys.argv[3]
	helper.make_dir(output_folder)
	print ("Done getting command_line arguments")
	total_bins, segment_df = get_total_num_segments_from_liftOver(segment_fn) #15478457
	print ("Done getting the total number of bins: " + str(total_bins))
	num_bins_to_pick = int(float(total_bins) * fract_gene_to_sample)
	train_indices = np.random.choice(range(int(total_bins)), num_bins_to_pick, replace = False) # pick without replacement
	train_indices = np.sort(train_indices)
	train_fn = os.path.join(output_folder, 'train_segments.bed.gz')
	get_subset_rows_and_save(segment_df, train_indices, train_fn)
	test_indices = np.setdiff1d(np.array(range(int(total_bins))), train_indices) # indices of lines that we will take as the test segments
	test_indices =  np.sort(test_indices)
	test_fn = os.path.join(output_folder, 'test_segments.bed.gz')
	get_subset_rows_and_save(segment_df, test_indices, test_fn)
	print ("Done getting test and train segmentation data")

def usage():
	print ("python sample_random_segments.py")
	print ("segment_fn: segments should be converted to the 200bp format")
	print ("fraction of genome to sample")
	print ("output_folder: where bed file output to get the sampled segments")
	exit(1)

main()