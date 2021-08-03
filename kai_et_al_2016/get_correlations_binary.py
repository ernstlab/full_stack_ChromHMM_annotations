import sys
import os
import numpy as np
import gzip
import math
from scipy.stats.stats import pearsonr
def calculate_correlation_matrix(num_marks, num_genome_postions, sum_x, sum_x_squared, sum_xy, outDir):
	new_num_marks = num_marks
	for i in range(num_marks): # we will not calculate the marks that are all zeros
		if sum_x[i] == 0:
			print "Mark " + str(i) + " contains all 0s! "
			new_num_marks -= 1
	print "The total number of marks that we can calculate correlations now is " + str(new_num_marks)
	correlations = np.ones((new_num_marks, new_num_marks))
	real_i = 0
	for i in range(num_marks):
		real_j = 1
		if sum_x[i] != 0:
			n_avg_x_squared = float(np.power(sum_x[i], 2.0)) / float(num_genome_postions)
			for j in range(i+1, num_marks):
				if sum_x[j] != 0:
					n_avg_x_avg_y = float(sum_x[i] * sum_x[j]) / float(num_genome_postions)
					n_avg_y_squared = float(np.power(sum_x[j], 2.0)) / float(num_genome_postions)
					nominator = sum_xy[i,j] - n_avg_x_avg_y
					denominator = np.power(sum_x_squared[i]- n_avg_x_squared, 0.5) * \
					np.power(sum_x_squared[j] - n_avg_y_squared, 0.5)
					corr = nominator / denominator
					correlations[real_i, real_i + real_j] =  corr
					correlations[real_j + real_i,real_i] = corr
					real_j += 1
			real_i += 1

	np.savetxt(os.path.join(outDir, "correlations_binary.txt"), correlations)
	correlations = np.absolute(correlations)
	np.savetxt(os.path.join(outDir, "abs_correlations_binary.txt"), correlations)

"""
num_snp = 3
num_samples = 4
sum_x = np.array([10, 25, 20])
sum_x_squared = np.array([30, 171, 120])
sum_xy = np.array([[0, 67, 60], [0, 0, 134], [0, 0, 0]])
corrF = open("trial_to_be_deleted.txt", 'w')
correlations = calculate_correlation_matrix(num_snp, num_samples, sum_x, sum_x_squared, sum_xy, corrF)
corrF.close()
"""
def read_report_mark_correlation(inputFolder, input_file_list, outputFName):
	data_files = [os.path.join(inputFolder, fn) for fn in input_file_list]
	certain_num_marks = 0 
	sum_x = None
	sum_x_squared = None
	sum_xy = None
	num_genome_postions = 0 
	for i, file in enumerate(data_files): # process each file in the folder
		print "Processing file: " + file
		try:
			binary_data = np.loadtxt(file, skiprows=2)
		except:
			print "could not open file " + file
			continue
		if i == 0: # the first time we process the files, gotta set up the number of marks, and the data structures 
			print binary_data.shape
			genome_postions, certain_num_marks = binary_data.shape
			sum_x = np.zeros((certain_num_marks,))
			sum_x_squared = np.zeros((certain_num_marks,))
			sum_xy = np.zeros((certain_num_marks, certain_num_marks))
		(genome_postions, num_marks) = binary_data.shape
		num_genome_postions += genome_postions
		assert num_marks == certain_num_marks, "The number of marks in file " + file + " is weird: " + str(num_marks)
		for j in range(num_marks):
			sum_xy[j,j] = 1
			sum_x[j] += np.sum(binary_data[:,j])
			sum_x_squared[j] += np.sum(np.square(binary_data[:,j]))
			for k in range(j+1, num_marks):
				this_xy = np.sum(np.multiply(binary_data[:,j], binary_data[:,k]))
				sum_xy[j,k] += this_xy
				sum_xy[k,j] += this_xy
	calculate_correlation_matrix(num_marks, num_genome_postions, sum_x, sum_x_squared, sum_xy, outputFName)
	return num_marks

def write_mark_names(num_marks, mark_name_fn):
	outF = open(mark_name_fn, 'w')
	for i in range(num_marks):
		outF.write(str(i) + "\n")
	outF.close()

def get_input_file_list(input_file_list_fn):
	inF = open(input_file_list_fn, 'r')
	input_file_list = []
	for line in inF:
		input_file_list.append(line.strip())
	return input_file_list

def main():
	if len(sys.argv) != 4:
		print "Wrong input man!"
		print "python get_correlation_between_marks.py \
		<input folder> \
		<output file>"
		exit(1)
	else:
		inputFolder = sys.argv[1]
		input_file_list_fn = sys.argv[2]
		outDir = sys.argv[3]
		try: 
			os.makedirs(outDir) # create output directory recursively
			print "Created output directory recursively"
		except: 
			print "Output directory already exists"

		mark_name_fn = os.path.join(outDir, "mark_names.txt")
		input_file_list = get_input_file_list(input_file_list_fn)
		print "Done getting input file list"
		num_marks = read_report_mark_correlation(inputFolder, input_file_list, outDir)
		print "Done calculating the correlations"
		write_mark_names(num_marks, mark_name_fn)
		print "Done writing down mark names"

main()
