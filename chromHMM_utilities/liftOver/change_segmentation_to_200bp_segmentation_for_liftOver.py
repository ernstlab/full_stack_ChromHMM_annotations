# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import os
import sys
import time
import gzip
start_time = time.clock()
NUM_BP_PER_BIN = 200
CHR_COL_IND = 0 # column index that contains information about the chromosome
START_COL_IND = 1
END_COL_IND = 2
STATE_COL_IND = 3
def open_file(fn): 
	"""
	Open a zipped or unzipped file. Return the open file object
	"""
	if fn[-3:] == ".gz":
	        # zip file
	        F = gzip.open(fn, 'rb')
	else: 
	        # non zip file
	        F = open(fn, 'r')
	return F

def split_segmentation_into_200bp(org_seg_fn, output_fn):
	orgF = open_file(org_seg_fn)
	outF = gzip.open(output_fn, 'wb')
	for line in orgF:
		line_data = line.strip().split() # chr , start_pos, end_pos, chromatin state
		org_end_bp = int(line_data[END_COL_IND])
		org_start_bp = int(line_data[START_COL_IND])
		num_bins_this_line = (org_end_bp - org_start_bp) / NUM_BP_PER_BIN
		for bin_ind in range(num_bins_this_line): 
			this_bin_start = org_start_bp + bin_ind * NUM_BP_PER_BIN
			this_bin_end = this_bin_start + NUM_BP_PER_BIN
			outF.write(line_data[CHR_COL_IND] + "\t" + str(this_bin_start) + "\t" + str(this_bin_end) + "\t" + line_data[STATE_COL_IND] + "\n")
	outF.close()
	orgF.close()

def main():
	if len(sys.argv) != 3:
		usage()
	org_seg_fn = sys.argv[1]
	if not os.path.isfile(org_seg_fn):
		print "org_seg_fn: "+ org_seg_fn + " does not exit"
		usage()
	output_fn = sys.argv[2]
	print "Done getting command line argument after: " + str(time.clock() - start_time)
	split_segmentation_into_200bp(org_seg_fn, output_fn)
	print "Done after: " + str(time.clock() - start_time)

def usage():
	print "python change_semgentation_to_200bp_segmentation_for_liftOver.py"
	print "org_seg_fn"
	print "output_fn"
	exit(1)

main()