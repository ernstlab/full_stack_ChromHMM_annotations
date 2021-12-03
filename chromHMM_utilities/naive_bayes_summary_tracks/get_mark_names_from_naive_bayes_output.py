# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import sys
import os
import naive_bayes_helpers as nbh
import time
start_time = time.clock()
def strip(word):
	return word.strip()

def get_chosen_mark_names_from_output(mark_rank_fname, mark_names_list, output_fName):
	rankF = open(mark_rank_fname, 'r')
	outF = open(output_fName, 'w')
	for line in rankF:
		line = line.strip()
		mark_index = int(line.split()[0])
		outF.write(mark_names_list[mark_index] + "\n")
	rankF.close()
	outF.close()

def get_all_mark_names_from_binary_data(binary_data_folder):
	files_in_this_folder = os.listdir(binary_data_folder) # this is based on the assumption that all files in this folder are just binary data files
	binF = nbh.open_file(os.path.join(binary_data_folder, files_in_this_folder[0]))
	binF.readline() # first line is just the genomic position, we do not care
	headers = binF.readline() # second line is the the names of the marks
	mark_names_list = headers.strip().split()
	binF.close()
	return mark_names_list

def main():
	if len(sys.argv) != 4: usage()
	else:
		mark_rank_fname = sys.argv[1]
		binary_data_folder = sys.argv[2]
		if not os.path.isdir(binary_data_folder): usage()
		output_fName = sys.argv[3]
		mark_names_list = get_all_mark_names_from_binary_data(binary_data_folder)
		print "Got all mark names after: " + str(time.clock() - start_time)
		get_chosen_mark_names_from_output(mark_rank_fname, mark_names_list, output_fName)
		print "Got chosen mark names after: " + str(time.clock() - start_time)

def usage():
	print "Wrong input man!"
	print "python get_mark_names_from_output.py \
	<mark_rank_fname> \
	<binary_data_folder> \
	<output_fName>"
	exit(1)
main()