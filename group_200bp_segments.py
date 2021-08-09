# this file is tested to do the correct work on 01172020
# it does the inverse of change_segmentation_to_200bp_segmentation_for_liftOver.py
import os
import sys
import gzip
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

def group_200bp_segments(org_seg_fn, output_fn):
	orgF = open_file(org_seg_fn)
	outF = gzip.open(output_fn, 'wb')	
	orgF.readline() # skip the first line
	current_chrom = ""
	current_state = ""
	current_start = 0
	current_end = 0
	for line in orgF:
		line_data = line.strip().split() # chr , start_pos, end_pos, chromatin state
		org_chrom = line_data[CHR_COL_IND]
		org_end_bp = int(line_data[END_COL_IND])
		org_start_bp = int(line_data[START_COL_IND])
		org_state = line_data[STATE_COL_IND]
		if org_chrom != current_chrom:
			if current_chrom != "": # no need to write if it's ""
				write_data = current_chrom + "\t" + str(current_start) + "\t" + str(current_end) + "\t" + current_state
				outF.write(write_data + "\n")
			current_chrom = org_chrom
			current_state = org_state
			current_start = org_start_bp
			current_end = org_end_bp 
			continue # change of chromosome
		# same chromosome and different state
		if org_state != current_state:
			write_data = current_chrom + "\t" + str(current_start) + "\t" + str(current_end) + "\t" + current_state
			outF.write(write_data + "\n")
			current_start = org_start_bp
			current_state = org_state
			current_end = org_end_bp
			continue # change of state
		# same chromosome and same state
		if org_start_bp == current_end: # it's a segment that is of the same state and also is a continuation of the previous segment 
			current_end = org_end_bp
		else: # same chromsome, same state, but it's a segment that has a gap with the previous segment. Therefore, we will write the previous segment down and update the current segment
			write_data = current_chrom + "\t" + str(current_start) + "\t" + str(current_end) + "\t" + current_state
			outF.write(write_data + "\n")
			current_state = org_state
			current_start = org_start_bp
			current_end = org_end_bp
	write_data = current_chrom + "\t" + str(current_start) + "\t" + str(current_end) + "\t" + current_state
	outF.write(write_data + "\n")
	outF.close()
	orgF.close()

def main():
	if len(sys.argv) != 3:
		usage()
	org_seg_fn = sys.argv[1]
	if not os.path.isfile(org_seg_fn):
		print ("org_seg_fn: "+ org_seg_fn + " does not exit")
		usage()
	output_fn = sys.argv[2]
	print ("Done getting command line argument")
	group_200bp_segments(org_seg_fn, output_fn)
	print ("Done!")

def usage():
	print ("python group_200bp_segments.py")
	print ("org_seg_fn")
	print ("output_fn")
	exit(1)

main()