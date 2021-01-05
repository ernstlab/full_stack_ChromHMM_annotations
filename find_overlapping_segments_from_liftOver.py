# find segments that map to the same place on the genome when we perform some liftOver job
# it only works one chromosome at a time
# it only report the invalid lines and not fix the problem
# Therefore, this code is a little useless
import os
import sys
import string
import gzip
CHROM_INDEX = 0
START_INDEX = 1
END_INDEX = 2
STATE_INDEX = 3
###### some helper functions, can also be found at chromHMM_utilities_common_functions_helper.py #######
def make_dir(directory):
	try:
		os.makedirs(directory)
	except:
		print 'Folder' + directory + ' is already created'



def check_file_exist(fn):
	if not os.path.isfile(fn):
		print "File: " + fn + " DOES NOT EXISTS"
		exit(1)
	return 

def check_dir_exist(fn):
	if not os.path.isdir(fn):
		print "Directory: " + fn + " DOES NOT EXISTS"
		exit(1)
	return 
	
def create_folder_for_file(fn):
	last_slash_index = string.rfind(fn, '/')
	if last_slash_index != -1: # path contains folder
		make_dir(fn[:last_slash_index])
	return 
###### end of some helper functions, can also be found at chromHMM_utilities_common_functions_helper.py #######
def report_overlapping_bases(input_fn, output_fn):
	inF = gzip.open(input_fn, mode = 'rb')
	outF = gzip.open(output_fn, 'wb')
	line = inF.readline() # read the first line 
	invalid_lines = []
	line_data = line.strip().split()
	current_start = int(line_data[START_INDEX])
	current_end = int(line_data[END_INDEX])
	current_line_index = 0
	for line in inF:
		line_data = line.strip().split()
		this_start =  int(line_data[START_INDEX])
		this_end = int(line_data[END_INDEX])
		if this_start < current_end:
			invalid_lines.append(current_line_index)
			invalid_lines.append(current_line_index + 1)
		current_line_index += 1
		current_start = this_start
		current_end = this_end
	invalid_lines = list(set(invalid_lines)) # get only the unique invalid lines
	invalid_lines = map(lambda x: str(x), invalid_lines) # convert to string so that we can write to file			
	outF.write('\n'.join(invalid_lines))
	inF.close()
	outF.close()
	

def main():
	if len(sys.argv) != 3: 
		usage()
	input_fn = sys.argv[1]
	check_file_exist(input_fn)
	output_fn = sys.argv[2]
	create_folder_for_file(output_fn)
	print "Done getting command line argument"
	report_overlapping_bases(input_fn, output_fn)
	print "Done!"

def usage():
	print "python find_overlapping_segments_from_liftOver.py"
	print "input_fn: the results of liftOver. Assumption: this file only contains infomation for one chromsome. In other words, the chromosome column in this file is exactly similar across different rows. This file should be sorted by increasing starting bp"
	print "output_fn: the corrected segmentation file. Any segments that map to the same places with any other segments will be eradicated"
	exit(1)
main()