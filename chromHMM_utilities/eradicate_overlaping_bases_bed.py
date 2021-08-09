print "There is a potential that there are overlapping regions in a CDS or noncoding bed file"
print "It means that the end of one line in a bed file is within the range of the next line"
print "We want to identify that and then fix that"

import sys
import os
import time
import gzip
start_time = time.clock()

def create_folder_for_file(fn):
	last_slash_index = string.rfind(fn, '/')
	if last_slash_index != -1: # path contains folder
		make_dir(fn[:last_slash_index])
	return 

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

def get_num_uniq_bases (input_rearranged_fn, output_fn):
	inputF = open_file(input_rearranged_fn)
	outF = gzip.open(output_fn, 'wb')
	# initialize these values so that the first line will definitely be reported into the output file
	current_chromosome = ""
	current_end_bp = -1
	total_bp = 0 # total number of uniq base pairs in this bed file
	for line in inputF: 
		line = line.strip()
		line_data = line.split("\t")
		this_chrom = line_data[0]
		this_start_bp = int(line_data[1])
		this_end_bp = int(line_data[2])
		if this_chrom != current_chromosome:
			outF.write(line + "\n")
			total_bp += (this_end_bp - this_start_bp)
		else: # same chromosome
			if this_start_bp > current_end_bp:
				outF.write(line + "\n")
				total_bp += (this_end_bp - this_start_bp)
			else: # start of this line is within the range of the previous line's start_bp - end_bp
				if current_end_bp == this_end_bp: # there is no point in reporting this region where the start and end base pairs are the same
					continue 
				if current_end_bp > this_end_bp: # this region is actually a part of the previous region, so we skip it
					continue
				line_data[1] = str(current_end_bp) # replace the starting bp of this line and report it later
				outF.write("\t".join(line_data) + "\n")
				total_bp += ((this_end_bp - this_start_bp) - (current_end_bp - this_start_bp))
		current_chromosome = this_chrom
		current_end_bp = this_end_bp
	inputF.close()
	outF.close()
	print "In total: There are: " + str(total_bp) + " base pairs in this bed file"

def main():
	if len(sys.argv) != 3:
		usage()
	input_rearranged_fn = sys.argv[1]
	if not os.path.isfile(input_rearranged_fn):
		print "input_rearranged_fn: " + input_rearranged_fn + " does not exist"
		usage()
	output_fn = sys.argv[2] 
	print "Done getting command line arguments after: " + str(time.clock() - start_time)
	get_num_uniq_bases (input_rearranged_fn, output_fn)
	print "Done getting uniq bases after: " + str(time.clock() - start_time)

def usage():
	print "pypy fix_cds_non_coding_region_overlap.py"
	print "input_rearranged_fn: this file should be the output of reorganize_non_coding_GENCODE.sh"
	print "output_fn: the file with no overlapping non coding / cds regions' data is stored"
	exit(1)
main()
