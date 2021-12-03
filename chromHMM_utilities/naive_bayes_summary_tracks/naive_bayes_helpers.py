"""
This file contains functions that are used in get_naive_bayes_data.py and naive_bayes_version2.py
This is specifically designed for running Jason's 
"""
import os
import gzip
from collections import Counter
import string
# import pandas as pd
# import numpy as np
SEGMENT_LENGTH = 200
NUM_BASE_PER_BIN_FILE = 1000000 # number of base pari covered in a binary data file


def make_dir(directory):
	try:
		os.makedirs(directory)
	except:
		print 'Folder' + directory + ' is already created'

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


def get_segment_data_into_dictionary(segment_fn, binary_files):
	chr_name_list = []
	state_counter = Counter([])

	for i, bin_fn in enumerate(binary_files):
		name_list = bin_fn.split("_")
		chrom_name =  name_list[1] # this is specifically designed for the binary data file name that follows the format:genome_chr10.0_binary.txt.gz
		#".".join([name_list[0], name_list[1]]), 
		chr_name_list.append(chrom_name)
	chr_name_counter = Counter(chr_name_list)
	segmentF = open_file(segment_fn)
	segment_dictionary = {}
	for line in segmentF:
		[chrom, start_index, end_index, state] = line.strip().split()
		# chrom_subscript = int(int(start_index) / NUM_BASE_PER_BIN_FILE)  # subscipt starts from 0
		# chrom = ".".join([chrom, str(chrom_subscript)])
		if chr_name_counter[chrom] > 0: # if this state is is the list of binary data file that we want to process
			state_index = int(state[1:]) - 1 # because raw state would be like 'E100'
			num_appearances = (int(end_index) - int(start_index)) / SEGMENT_LENGTH
			state_counter.update({state_index:num_appearances})
			if chrom in segment_dictionary:
				(segment_dictionary[chrom]).extend([state_index] * num_appearances)
			else: 
				(segment_dictionary[chrom]) = [state_index] * num_appearances
	return state_counter, segment_dictionary, chr_name_list

def get_input_file_list(input_file_list_fn):
	inF = open(input_file_list_fn, 'r')
	input_file_list = []
	for line in inF:
		input_file_list.append(line.strip())
	return input_file_list

def get_marks_index(mark_include_fn, num_marks_to_include): # example file:  /u/home/h/havu73/project-ernst/naive_bayes/output_cluster_02252018/naive_bayes2_small80_cmd_425_02252018.txt
# Mark index is from 0
	markF = open(mark_include_fn, 'r')
	marks_to_include = []
	for line in markF:
		marks_to_include.append(int(line.strip().split()[0]))
		if len(marks_to_include) == num_marks_to_include: # We only take the first num_marks_to_include marks
			break
	markF.close()
	return marks_to_include