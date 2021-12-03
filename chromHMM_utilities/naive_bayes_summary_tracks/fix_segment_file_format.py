# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
When our segment is of the from: 
chr11	72233800	72234000	E41
chr11	72234000	72234200	E42
chr11	72234200	72234400	E59
chr11	72234400	72234600	E58
chr11	72234600	72234800	E57
chr11	72234800	72235200	E58
chr11	72235200	72235600	E56
we want it to become 
chr11.72	233800	234000	E41
chr11.72	234000	234200	E42
chr11.72	234200	234400	E59
chr11.72	234400	234600	E58
chr11.72	234600	234800	E57
chr11.72	234800	235200	E58
chr11.72	235200	235600	E56

Tested correct on 08/31/2018
zcat genome_100_segments.bed.gz |awk '{print $4}' > raw_states.txt
awk '{print $4}' genome_100_segments_divided.bed >fixed_states.txt
diff raw_states.txt fixed_states.txt --> empty
"""
SEGMENT_LENGTH = 200
NUM_BASE_PER_BIN_FILE = 1000000 # number of base pari covered in a binary data file
import numpy as np 
import pandas as pd
import naive_bayes_helpers as nbh
import os
import sys
import string
import time
start_time = time.clock()
def append_subscript(chrom_name):
	return ".".join([chrom_name[0], str(chrom_name[1])])

def create_data_frame_from_file(fn, this_header):
	"""
	Open a zipped or unzipped file. Return the pandas file object
	"""
	try:
		# zip file
		df = pd.read_table(fn, compression = "gzip", header = this_header)
	except: 
		# non zip file
		df = pd.read_table(fn, header = this_header)
	return df

def num_section_separated_by_dot(word):
	return len(word.split("."))

def find_chrom_subscript(start_index):
	return int(start_index / NUM_BASE_PER_BIN_FILE)

def get_segment_dataframe(segment_fn):
	"""
	Put the data from segmentation (output of ChromHMM) into a pandas data frame with more information
	"""
	# get the data into a data frame
	segment_df = create_data_frame_from_file(segment_fn, this_header = None)
	# check what kind of format of chromosome regions we are in 
	num_sec_name = segment_df[0].apply(num_section_separated_by_dot)
	if (num_sec_name == 1).all(): # if all the chromosome region names are of the form "chrsomethingsomething"	
		# fix the name of chromosome so that now we have correct chromosome name. This is necessary when I processed data from Jason with stacked marks
		chrom_subscript = segment_df[1].apply(find_chrom_subscript)
		chrom_name = zip(segment_df[0], chrom_subscript)
		chrom_name = map(append_subscript, chrom_name)
		segment_df[0] = chrom_name
	elif (num_sec_name == 1).all(): # if all the chromosome region names are of the form "chrsomethingsomething.somethingsomething"
		pass
	else: 
		print "The names of chromsome regions in segmentation files are not consistent and maybe not of the known format"
		exit(1)

	# transform the start and end index to be mod 1000000
	segment_df[segment_df.columns[1:3]] %= NUM_BASE_PER_BIN_FILE
	# calculate the number of occurences of each state
	# segment_df["num_occ"] = (segment_df[2] - segment_df[1]) / SEGMENT_LENGTH
	return segment_df

def main():
	if len(sys.argv) != 3:
		usage()
	segment_fn = sys.argv[1]
	output_segment_fn = sys.argv[2]
	if not os.path.isfile(segment_fn):
		print "segment file not found"
		usage()
	nbh.create_folder_for_file(output_segment_fn)
	print "Done processing input argment after: " + str(time.clock() - start_time)
	segment_df = get_segment_dataframe(segment_fn)
	print "Done processing the segment data after: " + str(time.clock() - start_time)
	segment_df.to_csv(output_segment_fn, sep = '\t', header = False, index = False)

def usage():
	print "python fix_segment_file_format.py "
	print "segment_fn"
	print "output_segment_fn"
	exit(1)
main()