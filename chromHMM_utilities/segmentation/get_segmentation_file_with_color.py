#import seaborn as sns
import pandas as pd 
import numpy as np 
import os
import sys
import time
import chromHMM_utilities_common_functions_helper as helper
	
def get_rgb_format_right(rgb):
	# convert from (255, 245, 238) to 255,245,238
	# numbers = (rgb[1:-1]).split(',') # get rid of () 
	numbers = (rgb[:]).split(',')
	numbers = list(map(lambda x: str(int(x)), numbers))
	return ",".join(numbers)

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    numbers = tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
    numbers = list(map(lambda x: str(x), numbers))
    return ",".join(numbers)

def get_standard_ucsc_format_bed_file(segment_fn, state_annot_fn, output_fn, track_name, chrom, start_bp, end_bp):
	segment_df = pd.read_csv(segment_fn, sep = '\t', header = None)
	segment_df.columns = ['chrom', 'start_bp', 'end_bp', 'state']
	# filter the segmentation data for the region that we are interested in 
	if chrom != None and start_bp != None and end_bp != None:
		segment_df = segment_df[segment_df['chrom'] == chrom]
		segment_df = segment_df[segment_df['start_bp'] >= start_bp]
		segment_df = segment_df[segment_df['end_bp'] <= end_bp]
	# get state annotation data
	state_df = pd.read_csv(state_annot_fn, sep = ',', header = 0)
	state_df['state'] = (state_df['state']).apply(lambda x: 'E' + str(x))
	segment_df = pd.merge(segment_df, state_df, how = 'left', left_on = 'state', right_on = 'state', left_index = False, right_index = False)
	print(segment_df.columns)
	segment_df = segment_df [['chrom', 'start_bp', 'end_bp', 'mneumonics', 'color']]
	# print(segment_df.head())
	(nrow, ncol) = segment_df.shape
	print(segment_df.head())
	segment_df['color'] = (segment_df['color']).apply(hex_to_rgb)
	segment_df['score'] = ['1'] * nrow
	segment_df ['strand'] = ['.'] * nrow
	segment_df['thickStart'] = segment_df['start_bp']
	segment_df['thickEnd'] = segment_df['start_bp']
	# segment_df['blockCount'] = ['.'] * nrow
	# segment_df['blockSizes'] = ['.'] * nrow
	# segment_df['blockStarts'] = ['.'] * nrow
	segment_df = segment_df.rename(columns = {'start_bp' : 'chromStart', 'end_bp' : 'chromEnd', 'mneumonics' : 'name'})
	segment_df = segment_df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'color']]
	outF = open(output_fn, 'a')
	header_comment = "track name=\"" + track_name + "\" description=\"\" visibility=1 itemRgb=\"On\"\n"
	outF.write(header_comment)
	segment_df.to_csv(outF, sep = '\t', header = False, index = False)
	outF.close()


def main():
	print(sys.argv)
	print (len(sys.argv))
	if len(sys.argv) != 8: 
		usage()
	segment_fn = sys.argv[1]
	helper.check_file_exist(segment_fn)
	state_annot_fn = sys.argv[2]
	helper.check_file_exist(state_annot_fn)
	output_fn = sys.argv[3]
	try:
		os.remove(output_fn) # remove to replace later with a new one
	except: # if the file does not exist, that's okay, move on!
		pass
	helper.create_folder_for_file(output_fn)
	track_name = sys.argv[4]
	try:
		chrom = sys.argv[5]
	except:
		chrom = None
	try:
		start_bp = int(sys.argv[6])
	except:
		start_bp = None
	try:
		end_bp = int(sys.argv[7])
	except:
		end_bp = None
	if chrom==None or start_bp==None or end_bp==None:
		chrom = None
		start_bp = None
		end_bp = None
	print ("Done reading command line arguments")
	get_standard_ucsc_format_bed_file(segment_fn, state_annot_fn, output_fn, track_name, chrom, start_bp, end_bp)
	print ("Done!")
	
def usage():
	print ("python get_segmentation_file_with_color.py")
	print ("segment_fn: where the segmentation data is stored")
	print ("state_annot_fn: where the data of annotations of states (colors, annotations, etc.) are stored")
	print ("output_fn: where the full bed file used for ucsc track, readable from igv should be saved. Note: not gzip file, have to gzip separately")
	print ("track_name: the name of the track that will be on display when we upload onto uscs genome browser")
	print ("chrom: chromosome where we want to get the segment data for")
	print ("start_bp: start_bp where we want to get the segment data for")
	print ("end_bp: end_bp where we wwant to get the segment data for")
	exit(1)

main()
