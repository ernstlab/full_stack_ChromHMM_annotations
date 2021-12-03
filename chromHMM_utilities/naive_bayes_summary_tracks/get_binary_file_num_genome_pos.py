# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
This file contains function that will take as argument path to the binary files. Then, it will report to a file the number of genomic positions in each file
"""

import os
import gzip
import time
import sys
import multiprocessing as mp
import naive_bayes_helpers as nbh
import numpy as np
start_time = time.clock()

def get_chrom_name_from_bin_fn(bin_fn): 
	"""
	Assuming that the bin_fn is like this: genome_chr2.211_binary.txt.gz
	Chrom name that is returned in this function will be like thisL chr2.211
	"""
	return bin_fn.split("_")[1]

def get_length_one_binary_data(bin_folder, bin_fn, output):
	full_bin_fn = os.path.join(bin_folder, bin_fn)
	chrom_name = get_chrom_name_from_bin_fn(bin_fn) # genome_chr2.211_binary.txt.gz --> chr2.211
	bin_data = np.loadtxt(full_bin_fn, skiprows = 2) # first two lines are meta data
	(num_genome_pos, num_marks) = bin_data.shape
	output.put((chrom_name, num_genome_pos, num_marks))
	print "Done: " + chrom_name
 
def get_binary_file_length_all_binary_files(binary_folder, output_fn, num_processes):
	binary_file_list = os.listdir(binary_folder)
	outF = open(output_fn, 'w')
	# initialze parallel processes to find the length of each binary file
	num_bin_file = len(binary_file_list)
	num_processing_round = int(num_bin_file / num_processes) + 1
	for i in range (num_processing_round):
		output = mp.Queue()
		start_index = i * num_processes
		end_index = min(num_bin_file, (i + 1) * num_processes)
		processes = [mp.Process(target = get_length_one_binary_data, args = (binary_folder, binary_file_list[file_index], output)) for file_index in range(start_index, end_index)]
		for p in processes:
			p.start()
		for p in processes:
			p.join()
			(chrom_name, num_genome_pos, num_marks) = output.get()
			# write out the data
			outF.write(chrom_name  + "\t" + str(num_genome_pos) + "\t" + str(num_marks) + "\n")
	outF.close()

def main():
	if len(sys.argv) != 4:
		usage()
	else:
		try:
			# get input data
			binary_folder = sys.argv[1]
			output_fn = sys.argv[2]
			num_processes = int(sys.argv[3])
			# Checking input 
			if (not os.path.isdir(binary_folder)):
				print "Something wrong with the binary folder"
				usage()
			try: 
				nbh.create_folder_for_file(output_fn)
			except: 
				pass
		except: 
			print "Wrong input format"
			usage()

	print "Done getting command line input: " + str(time.clock() - start_time)
	get_binary_file_length_all_binary_files(binary_folder, output_fn, num_processes)
	print "Done after: " + str(time.clock() - start_time)

def usage():
	print "python (dont use pypy fro this file) get_binary_file_num_genome_pos.py  <binary folder>  <output fn where each binary file's number of genomic positions are reported> <number of parallel processes used to run this job>"	
	exit(1)

main()