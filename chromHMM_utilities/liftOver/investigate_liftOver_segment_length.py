import pandas as pd 
import numpy as np 
import argparse
import seaborn as sns
import helper

def record_count_of_non200bp_segment(bed_fn, output_fn= None):
	print ('record_count_of_non200bp_segment:', bed_fn)
	try:
		bed_df = pd.read_csv(bed_fn, header = None, index_col = None, sep = '\t')
	except:
		print('Could not open file {} because it is most likely empty')
		return pd.DataFrame() # return an empty dataframe
	bed_df.columns = ['chrom', 'start', 'end', 'state']
	bed_df['diff'] = bed_df['end'] - bed_df['start']
	outlier_df = (bed_df[bed_df['diff'] != 200]).copy()
	length, counts = np.unique(outlier_df['diff'], return_counts = True)
	result_df = pd.DataFrame({'length': pd.Categorical(length), 'num_segments': counts })
	if output_fn != None:
		result_df.to_csv(output_fn, header = True, index=False, sep=',')
	return result_df


def record_count_of_non200bp_segment_from_list(bed_fn_list, output_fn):
	result_df_list = list(map(lambda x: record_count_of_non200bp_segment(x), bed_fn_list))
	result_df = pd.concat(result_df_list, ignore_index=True)
	result_df = result_df.groupby('length', as_index=False)['num_segments'].sum()
	if output_fn != None:
		result_df.to_csv(output_fn, header = True, index=False, sep=',')
	return result_df

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "This code look into the liftOver results before and after we have gotten rid of segments in dest_assembly (hg19 --> hg38 then dest_assembly is hg38) that are mapped from multiple org_assembly (hg19--> hg38 then org_assembly is hg19). Normally, we would like 200-bp from org_assembly to be mapped exactly to another 200-bp segment in dest_assembly.")
	parser.add_argument('--end_segment_raw', required = False, type = str, help='Output of liftOver when we input the org_segment_200bp_fn to ucsc genome browser\'s liftOver')
	parser.add_argument('--end_segment_no_overlap_from_org_list', nargs='+', required = True, type = str, help= 'After liftOver result in end_segment_raw, and we go an extra step to filter out all the segments that were mapped from multiple segments in org_assembly .')
	parser.add_argument('--output_prefix', required = True, type = str, help = 'Prefix to the files where we will store the csv files of the count of segments that are of different length than 200bp in the end_assembly')
	args = parser.parse_args()
	helper.check_file_exist(args.end_segment_raw)
	map(helper.check_file_exist, args.end_segment_no_overlap_from_org_list)
	helper.create_folder_for_file(args.output_prefix)
	print ('investigate_liftOver_segment_length.py: Done getting command line argument')
	raw_segment_count_fn = args.output_prefix + '_non200bp_raw_segments.csv'
	record_count_of_non200bp_segment(args.end_segment_raw, raw_segment_count_fn)
	nonOverlap_segment_count_fn = args.output_prefix+ '_non200bp_noOverlap_segments.csv'
	record_count_of_non200bp_segment_from_list(args.end_segment_no_overlap_from_org_list, nonOverlap_segment_count_fn)
	print('Done!')
