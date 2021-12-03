# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from collections import Counter
import os
import gzip
import time
import sys
# import pandas as pd
# import numpy as np
# import functools
import naive_bayes_helpers as nbh
SEGMENT_LENGTH = nbh.SEGMENT_LENGTH
NUM_BASE_PER_BIN_FILE = nbh.NUM_BASE_PER_BIN_FILE # number of base pari covered in a binary data file
start_time = time.clock()


def get_one_prior_probability(state, segment_df, total_positions):
	"""
	Returns the number of genomics positions that are of the state
	"""
	state_df = segment_df[segment_df[3] == state]
	num_pos_this_state = float(sum(state_df['num_occ']))
	return num_pos_this_state / float(total_positions)

"""
def get_prior_probabilities(segment_fn, output_fn): # this is the more organized version of finding prior probabilities. It cannot beat the speed of prior probabilities calculated using pypy
	# Get the correct form of segment data frame 
	
	# chromosome region (chr10_10) --> should be consistent with the file name in binary data, start index (mod NUM_BASE_PER_BIN_FILE), end index (mod NUM_BASE_PER_BIN_FILE), state, num_occ: number of 200bp windows that covers that spans this states ((end_index - start_index) / SEGMENT_LENGTH) 
	
	# get segmentation data
	segment_df = nbh.get_segment_dataframe(segment_fn)
	# np array of unique states
	states = list(segment_df[3].unique())
	total_positions = sum(segment_df["num_occ"])
	prior_probabilities = map(functools.partial(get_one_prior_probability, segment_df = segment_df, total_positions =  total_positions), states)
	pp_dict = dict(zip(states, prior_probabilities))
	outputF = open(output_fn, 'w')
	for state in pp_dict:
		outputF.write(state + "\t" + str(pp_dict[state]) + "\n")
	outputF.close()
"""



def get_prior_probabilities(segment_fn, output_fn):
	segmentF = nbh.open_file(segment_fn)
	outputF = open(output_fn, 'w')
	state_counter = Counter([])
	start_time = time.clock()
	for line in segmentF:
		[chrom, start_index, end_index, state] = line.strip().split()
		num_appearances = (int(end_index) - int(start_index)) / SEGMENT_LENGTH
		state_counter.update({state:num_appearances})
	total_occurences = sum(state_counter.values())
	print "Time passed: " + str(time.clock() - start_time)
	for state in state_counter:
		outputF.write(state + "\t" + str(float(state_counter[state]) / float(total_occurences)) + "\n")
	outputF.close()
	segmentF.close()




def get_conditional_probabilities(segment_fn, binary_folder, output_folder, num_marks):
	binary_files = os.listdir(binary_folder)
	start_time = time.clock()
	state_counter, segment_dictionary, chr_name_list  = nbh.get_segment_data_into_dictionary(segment_fn, binary_files)
	for chrom_name in chr_name_list:
		if len(segment_dictionary[chrom_name]) != 5000:
			print chrom_name + "\t" + str(len(segment_dictionary[chrom_name]))
	print segment_dictionary.keys()[:10]
	num_states = len(state_counter)
	print "Done with processing the segment file: num_states: " + str(num_states)
	print "Time passed: " + str(time.clock() - start_time)
	num_marks_in_states = [([0] * num_marks) for i in range(num_states)]
	for i, bin_fn in enumerate(binary_files):
		chrom_name = chr_name_list[i]
		print "Binary file: " + bin_fn
		state_annotation = segment_dictionary[chrom_name]

		bin_fn = os.path.join(binary_folder, bin_fn)
		binF = gzip.open(bin_fn, 'rb')
		binF.readline()
		binF.readline()
		line_index = 0
		# loop through each line in the file
		for line in binF:
			this_200bp_marks = map(lambda x: int(x), line.strip().split())
			this_group_state = state_annotation[line_index]
			# update the number of marks that are present in this state
			for j in range(num_marks):
				num_marks_in_states[this_group_state][j] += this_200bp_marks[j]
			line_index += 1
		print "processed file: " + bin_fn
		print "time passed: " + (str(time.clock() - start_time))
		binF.close()

	state_total_fn = os.path.join(output_folder, "state_total.txt")
	state_totalF = open(state_total_fn, 'w')
	for state in state_counter:
		state_totalF.write(str(state) + "\t" + str(state_counter[state]) + "\n")
	state_totalF.close()
	print "Done finding the total of states"
	print "time passed: " + str(time.clock() - start_time)
	cond_prob_file = open(output_folder + "/conditional_probabilities.txt", 'w')
	for i in range(num_states):
		for j in range (num_marks):
			if num_marks_in_states[i][j] == 0:
				num_marks_in_states[i][j] = 1
			cond_prob = float(num_marks_in_states[i][j]) / float(state_counter[i] + 1)
			cond_prob_file.write(str(cond_prob) + "\t")
		cond_prob_file.write("\n")
	cond_prob_file.close()
	print "Done calculating conditional probabilities and writing out things"
	print "time passed: " + str(time.clock() - start_time)


def get_num_marks(binary_folder, binary_file_list_fn):
	"""
	This functiononly looks at one binary data file and it will see how many marks are investigated.
	Return the number of marks
	"""
	bin_fileF = open(binary_file_list_fn, 'r')
	bin_fn = os.path.join(binary_folder, bin_fileF.readline().strip())  # just look at the first file in the list
	binF = gzip.open(bin_fn, 'rb')
	# The first two lines are header lines so we will skip those
	binF.readline()
	binF.readline()
	num_marks = len(binF.readline().strip().split())
	return num_marks

def main():
	if len(sys.argv) != 5:
		usage()
	else:
		try:
			# get input data
			segment_fn = sys.argv[1]
			binary_folder = sys.argv[2]
			binary_file_list_fn = sys.argv[3]
			output_folder = sys.argv[4]
			# Checking input 
			if (not os.path.isfile(segment_fn) or not os.path.isdir(binary_folder) or not os.path.isfile(binary_file_list_fn)):
				print "Something wrong here"
				usage()
			try: 
				os.makedirs(output_folder)
			except: 
				pass
		except: 
			print "Wrong input format"
			usage()

	print "Done getting command line input: " + str(time.clock() - start_time)
	get_prior_probabilities(segment_fn = segment_fn, output_fn = os.path.join(output_folder , "prior_probabilities.txt"))
	print "Done getting prior probabilities: " + str(time.clock() - start_time)
	num_marks = get_num_marks(binary_folder, binary_file_list_fn)
	print "Done getting num_marks: " + str(num_marks) + "  "+ str(time.clock() - start_time)
	get_conditional_probabilities(segment_fn, binary_folder, output_folder, num_marks)
	print "Done getting conditional_probabilities: " + str(time.clock() - start_time)
	
def usage():
	print "python get_naive_bayes_data.py <segmentation file> <binary folder> <binary_file list file> <output folder>"	
	exit(1)

main()