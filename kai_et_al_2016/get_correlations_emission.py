import sys
import os
import numpy as np
from scipy.stats.stats import pearsonr
def read_report_correlations_emission(inputFName, outDir):
	"""
	Fixed to add calculate the absolute values of correlations this time
	Given the emission matrix output from ChromHMM training we will calculate the correlations between marks given its emission vector, outputed to outputFName
	"""
	# get the number of columns in this file, which is related to the number of marks
	with open(inputFName) as f:
		f.readline() # we skip the first line which contains mark names
		ncols = len(f.readline().strip().split())
		f.close()
	print "NUmber of columsn is " + str(ncols)
    ### Now only calculating emission
	emission = np.loadtxt(inputFName, skiprows=1, usecols = range(1,ncols)) # skip the first column and load data until the end of line
	(num_states, num_marks) = emission.shape 
	correlations = np.zeros((num_marks, num_marks))
	for i in range(num_marks):
		correlations[i,i] = 1
		for j in range(i+1, num_marks):
			(corr, p_value) = pearsonr(emission[:,i], emission[:,j])
			correlations[i,j] = corr
			correlations[j,i] = corr
	np.savetxt(os.path.join(outDir, "correlations_emission.txt") ,correlations)
	correlations = np.absolute(correlations)
	np.savetxt(os.path.join(outDir, "abs_correlations_emission.txt"), correlations)
	return num_marks

def read_report_absolute_correlations(inputFName, outputFName):
	"""
	This function was created due to some stupid mistakes during the research process. I fogot to calcuate the absolute values of correlations and printed out the correlations only. This fucntion takes the correlations matrix and return its absolute into output file
	"""
	correlations = np.loadtxt(inputFName)
	correlations = np.absolute(correlations)
	np.savetxt(outputFName, correlations)

def write_mark_names(num_marks, mark_name_fn):
	outF = open(mark_name_fn, 'w')
	for i in range(num_marks):
		outF.write(str(i) + "\n")
	outF.close()

def main():
	if len(sys.argv) != 3:
		print "Wrong input man!"
		print "python get_correlation_between_marks.py \
		<input folder> \
		<output file>"
		exit(1)
	else:
		inputFName = sys.argv[1]
		outDir = sys.argv[2]
		try: 
			os.makedirs(outDir) # create output directory recursively
			print "Created output directory recursively"
		except: 
			print "Output directory already exists"

		corr_out_fn = os.path.join(outDir, "correlations_emission.csv")
		mark_name_fn = os.path.join(outDir, "mark_names.txt")
		print "Calculating correlations between marks based on emission matrix"
		num_marks =  read_report_correlations_emission(inputFName, outDir)
		print "Writing mark names"
		write_mark_names(num_marks, mark_name_fn)
		print "Done"

main()

def read_report_state_diagonal_correlation_output_kai(inputFName, outputFName, num_states):
	"""
	Given the input fiel name as the file that contains confusion matrix from ChromHMM EvalSubset, we extract info of diagonal entries
	"""
	state_correlation_matrix = np.loadtxt(inputFName, skiprows=2, usecols = range(1, num_states + 1))
	outF = open(outputFName, 'w')
	for i in range(num_states):
		outF.write("E"+ str(i) + "\t" + str(state_correlation_matrix[i,i]) + "\n")
	outF.close()

