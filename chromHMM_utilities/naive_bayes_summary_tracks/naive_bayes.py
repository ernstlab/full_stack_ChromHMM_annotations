import math
import sys
import time
import os
from collections import Counter
import gzip
import naive_bayes_helpers as nbh
start_time = time.clock()

SEGMENT_LENGTH = nbh.SEGMENT_LENGTH
NUM_BASE_PER_BIN_FILE = nbh.NUM_BASE_PER_BIN_FILE # number of base pari covered in a binary data file
"""
state_counter: a counter dictionary [state] --> number of times that this state got appeared in the segmentation
state_recovery_counter : [state][mark]: how many times that this staete got recovered when this mark is included in the chosen marks
chosen_marks: counter dictionary [mark index] --> 1 or 0
binary_data_dict: dictionary of 2D tables: key: chomosome name, values: talbe of binary d atae [position][mark index]
prior_dict: dictionary {key: state index, values: prior probabilities of that state}
cond_prob and inverse_cond_prob: 2D table [state_index][mark index]
segment_dictionary: diction {key: genomeregion name} {value: list of state at each position in the region}
chrom_name_list has the same order as binary_files.
	Ex: binary_files: chr1_1_binary.txt.gz, chr1_2_binary.txt.gz
	then chrom_name_list will be chr1.1 and chr1.2
	we need chrom_name_list so that the names of files as keys to binary_data_dict will be the same as the keys in segment_data_dict
current_likelihood_dict: key: region on the genome, values: lists of lists: [position in the genomic region] [state index]
"""
def get_state_mark_cond_prob(cond_prob_fn):
	cond_probF = open(cond_prob_fn, 'r')
	cond_prob = []
	inverse_cond_prob = []
	for line in cond_probF:
		probs = map(lambda x: float(x), line.strip().split())
		log_probs = []
		inverse_log_probs = []
		for x in probs: 
			if x == 0.0: 
				log_probs.append(float("-inf"))
				inverse_log_probs.append(0.0)
			else:
				log_probs.append(math.log(x, 10))
				inverse_log_probs.append(math.log(1.0 - x, 10))
		cond_prob.append(log_probs)
		inverse_cond_prob.append(inverse_log_probs)
	return cond_prob, inverse_cond_prob

def get_prior_probabilities(prior_prob_fn):
	prior_dict = {}
	priorF = open(prior_prob_fn, 'r')
	for line in priorF:
		[state, prob] = line.strip().split()
		state = int(state[1:]) - 1 # sample state: 'E23'
		prob = math.log(float(prob), 10.0)
		prior_dict[state] = prob
	return prior_dict

def most_likely_state(current_likelihood, marks, cond_prob, inverse_cond_prob, mark_index, NUM_STATES, NUM_MARKS):
	
	max_log_likelihood = float("-inf")
	max_state = None # The index of the state that is most likely from this
	for state_index in range(NUM_STATES):
		if marks[mark_index] == 1: ## this mark is present at the position
			this_obs_given_state = current_likelihood[state_index] + cond_prob[state_index][mark_index]
		else: # this mark is not present at the position
			this_obs_given_state = current_likelihood[state_index] + inverse_cond_prob[state_index][mark_index]
		if this_obs_given_state > max_log_likelihood:
			max_log_likelihood = this_obs_given_state
			max_state = state_index
	return max_state, max_log_likelihood


def find_best_mark(chosen_marks, state_recovery_counter, state_counter, NUM_STATES, NUM_MARKS):
	max_recovery = 0 
	max_mark = None
	for mark_index in range(NUM_MARKS):
		if chosen_marks[mark_index] == 0: # is this mark has not been chosen before
			this_mark_recovery = 0
			for state_index in range(NUM_STATES):
				try:
					this_mark_recovery += float(state_recovery_counter[state_index][mark_index]) / \
					float(state_counter[state_index])
				except: 
					pass
			if this_mark_recovery > max_recovery:
				max_recovery = this_mark_recovery
				max_mark = mark_index
	return max_mark, max_recovery

def combined_most_likely_state(current_likelihood, marks, cond_prob, inverse_cond_prob, mark_index, marks_to_consider, prior):
	max_log_likelihood = float("-inf")
	max_state = None # The index of the state that is most likely from this
	max_log_v1 = float("-inf")
	max_state_v1 = None
	for state_index in range(NUM_STATES):
		log_likelihood_v1 = prior[state_index]
		for m, i in enumerate(marks_to_consider):
			if marks[m] == 1: # this mark is present at the position
				log_likelihood_v1 += cond_prob[state_index][m]
			else: # this mark is not present at the position
				log_likelihood_v1 += inverse_cond_prob[state_index][m]

		if marks[mark_index] == 1: ## this mark is present at the position
			this_obs_given_state = current_likelihood[state_index] + cond_prob[state_index][mark_index]
		else: # this mark is not present at the position
			this_obs_given_state = current_likelihood[state_index] + inverse_cond_prob[state_index][mark_index]
		if this_obs_given_state != log_likelihood_v1:
			print "when we consider mark: " + str(mark_index)
			print "marks_to_consider: " + str(marks_to_consider)
			print "presence or not (v2): " + str(marks[mark_index])
			print "presence or not (v1): " + str([marks[m] for m in marks_to_consider])
			print "state: " + str(state_index)
			print "current_likelihood: " + str(current_likelihood[state_index])
			print "cond_prob(v2): "+ str(cond_prob[state_index][mark_index])
			print "inverse_cond_prob(v2): "+ str(inverse_cond_prob[state_index][mark_index])
			print "cond_prob(v1): "+ str([cond_prob[state_index][m] for m in marks_to_consider])
			print "inverse_prob(v1): "+ str([inverse_cond_prob[state_index][m] for m in marks_to_consider])
			print "prior: "+ str(prior[state_index])
			print 
		if this_obs_given_state > max_log_likelihood:
			max_log_likelihood = this_obs_given_state
			max_state = state_index
	return max_state, max_log_likelihood, max_state_v1, max_log_v1

def one_round_naive_base_to_find_best_mark(current_likelihood_dict, chosen_marks, segment_dictionary, cond_prob, inverse_cond_prob, prior_dict,\
 state_counter, binary_data_dict, chrom_name_list, NUM_STATES, NUM_MARKS):
	chosen_marks_list = list(chosen_marks.elements())
	state_recovery_counter = [([0] * NUM_MARKS) for i in range(NUM_STATES)] # [state index] [mark index]
	for chrom_name in chrom_name_list:
		binary_data = binary_data_dict[chrom_name]
		this_chrom_current_likelihood_table = current_likelihood_dict[chrom_name]
		assert len(binary_data) == len(this_chrom_current_likelihood_table), "This chromosome: " + chrom_name + \
		" does not have consistent lenth: " + str(len(binary_data)) + "   " + str(len(this_chrom_current_likelihood_table))
		this_chrom_state_list = segment_dictionary[chrom_name]
		for i, position in enumerate(binary_data):
			this_pos_state = this_chrom_state_list[i]
			this_pos_current_likelihood = this_chrom_current_likelihood_table[i]
			for mark_index in range(NUM_MARKS):
				if chosen_marks[mark_index] == 0:# this mark has not been included before
					marks_to_consider = chosen_marks_list + [mark_index]
					max_state, max_log_likelihood = most_likely_state(this_pos_current_likelihood, position, cond_prob, inverse_cond_prob, mark_index, NUM_STATES, NUM_MARKS)
					if max_state == this_pos_state: # the naive bayes gives the same results of state assignment as chromHMM
						state_recovery_counter[max_state][mark_index] += 1

	max_mark, max_recovery = find_best_mark(chosen_marks, state_recovery_counter, state_counter, NUM_STATES, NUM_MARKS)
	return max_mark, max_recovery

def initialize_likelihood(chrom_name_list, prior_dict, cond_prob, inverse_cond_prob, binary_data_dict, NUM_STATES, NUM_MARKS, included_marks):
	current_likelihood_dict = {} # key: binary file, value: [position][state index]
	state_prior_list = []
	num_included_marks = len(included_marks)
	for state_index in range(NUM_STATES):
		state_prior_list.append(prior_dict[state_index])
	for i, chrom_name in enumerate(chrom_name_list): # for each binary file
		# print "initialize_likelihood: chrom_name:" + chrom_name
		this_chrom_prior_table = [] # [position][state index]
		num_pos = len(binary_data_dict[chrom_name])
		for j in range (num_pos): # for each genomic position in this binary file
			# print "num_pos: " + str(j)
			# print "state_prior_list"
			# print state_prior_list
			this_pos_current_likelihood = list(state_prior_list) # I need the list thing because otherwise what we do after this will also change the state_prior_list
			for k, mark_index in enumerate(included_marks):
				# print "mark index included:" + str(mark_index)
				# print " binary_data_dict[chrom_name][j][mark_index]: " + str(binary_data_dict[chrom_name][j][mark_index])
				if binary_data_dict[chrom_name][j][mark_index] == 1: # if this mark is present in this position
					for state_index in range(NUM_STATES):
						# print "		state_index " + str(state_index)
						# print "		this_pos_current_likelihood[state_index] : " + str(this_pos_current_likelihood[state_index])
						# print "		cond_prob[state_index][mark_index]: " + str(cond_prob[state_index][mark_index])
						this_pos_current_likelihood[state_index] += cond_prob[state_index][mark_index]
						# print "		after calculation: this_pos_current_likelihood[state_index] : " + str(this_pos_current_likelihood[state_index])
				else: # fi this mark is not present in this position
					for state_index in range(NUM_STATES): 
						# print "		state_index " + str(state_index)
						# print "		this_pos_current_likelihood[state_index] : " + str(this_pos_current_likelihood[state_index])
						# print "		iverse_cond_prob[state_index][mark_index]: " + str(inverse_cond_prob[state_index][mark_index])
						this_pos_current_likelihood[state_index] += inverse_cond_prob[state_index][mark_index]
						# print "		after calculation: this_pos_current_likelihood[state_index] : " + str(this_pos_current_likelihood[state_index])

			# print "this_pos_current_likelihood:"
			# print this_pos_current_likelihood
			this_chrom_prior_table.append(this_pos_current_likelihood)
		# print "this_chrom_prior_table"
		# print this_chrom_prior_table
		current_likelihood_dict[chrom_name] = this_chrom_prior_table
	return current_likelihood_dict

def write_included_marks_into_output_file(output_file, included_marks):
	for mark in included_marks:
		output_file.write(str(mark) + "\n")
		output_file.flush()

def naive_bayes(binary_data_dict, chrom_name_list, num_marks_to_pick, cond_prob, inverse_cond_prob, prior_dict, \
	state_counter, segment_dictionary, output_file, NUM_STATES, NUM_MARKS, included_marks):

	chosen_marks = Counter(included_marks)
	print "Starting to run naive bayse ..."
	current_likelihood_dict = initialize_likelihood(chrom_name_list, prior_dict, cond_prob, inverse_cond_prob, binary_data_dict, NUM_STATES, NUM_MARKS, included_marks)
	# HA: debug
	# print "current_likelihood_dict"
	# print current_likelihood_dict
	print "Done with getting the current_likelihood_dict table: time passed: " + str(time.clock() - start_time)
	for i in range(num_marks_to_pick):
		#print "finding mark " + str(i)
		max_mark, max_recovery = one_round_naive_base_to_find_best_mark(current_likelihood_dict, chosen_marks, segment_dictionary, cond_prob, \
			inverse_cond_prob, prior_dict, state_counter, binary_data_dict, chrom_name_list, NUM_STATES, NUM_MARKS)
		# update the current likelihood dictionary
		for chrom_name,this_current_likelihood in current_likelihood_dict.iteritems():
			new_chrom_current_likelihood = []
			for position_index in range(len(this_current_likelihood)):
				new_position_likelihood = []
				if binary_data_dict[chrom_name][position_index][max_mark] == 1: 
					for state_index in range(NUM_STATES):
						new_position_likelihood.append(this_current_likelihood[position_index][state_index] + cond_prob[state_index][max_mark])

				else:
					for state_index in range(NUM_STATES):
						new_position_likelihood.append(this_current_likelihood[position_index][state_index] + inverse_cond_prob[state_index][max_mark])
				new_chrom_current_likelihood.append(new_position_likelihood)
			current_likelihood_dict[chrom_name] = new_chrom_current_likelihood


		output_file.write(str(max_mark) + "\t" + str(max_recovery) + "\n")
		output_file.flush()
		chosen_marks.update([max_mark])
		print "Found mark "+ str(max_mark) + "   " + str(max_recovery)
		print "time passed: " + str(time.clock() - start_time)
	output_file.close()

def get_binary_data(binary_input_folder, binary_files, chrom_name_list):
	"""
	chrom_name_list has the same order as binary_files.
	Ex: binary_files: chr1_1_binary.txt.gz, chr1_2_binary.txt.gz
	then chrom_name_list will be chr1.1 and chr1.2
	we need chrom_name_list so that the names of files as keys to binary_data_dict will be the same as the keys in segment_data_dict
	"""
	full_binary_fn = [os.path.join(binary_input_folder, x) for x in binary_files]
	full_binary_files = [gzip.open(bin_fn, 'rb') for bin_fn in full_binary_fn]
	binary_data_dict = {}
	for i, binF in enumerate(full_binary_files):
		file_data = []
		binF.readline()
		binF.readline()
		line_index = 0
		# loop through each line in the file
		for line in binF:
			this_200bp_marks = map(lambda x: int(x), line.strip().split())
			file_data.append(this_200bp_marks)
		binary_data_dict[chrom_name_list[i]] = file_data
		binF.close()		
	return binary_data_dict

def main():
	if len(sys.argv) < 8:
		usage()
	cond_prob_fn = sys.argv[1]
	prior_fn = sys.argv[2]
	binary_input_folder = sys.argv[3]
	binary_file_list_fn = sys.argv[4]
	segment_fn = sys.argv[5]
	output_file = open(sys.argv[6], 'w', 0) # no buffering of files
	output_file.write("Hello!")
	num_marks_to_pick = int(sys.argv[7])
	if len(sys.argv) == 10:
		existing_mark_fn = sys.argv[8]
		num_existing_mark = int(sys.argv[9])
		included_marks = nbh.get_marks_index(existing_mark_fn, num_existing_mark)
		print "Get existing marks after " + str(time.clock() - start_time)
		print "Marks that are included are " 
		print included_marks
		write_included_marks_into_output_file(output_file, included_marks)
	else: 
		print "There are no included marks"
		included_marks = []

	binary_files = nbh.get_input_file_list(binary_file_list_fn)
	
	cond_prob, inverse_cond_prob = get_state_mark_cond_prob(cond_prob_fn)
	# print "cond_prob"
	# print cond_prob
	# print "inverse_cond_prob"
	# print inverse_cond_prob
	NUM_STATES = len(cond_prob)
	NUM_MARKS = len(cond_prob[0])
	print "get all conditional probabilites after " + str(time.clock() - start_time)

	prior_dict = get_prior_probabilities(prior_fn)
	# print "prior_dict"
	# print prior_dict
	assert len(prior_dict) == NUM_STATES, "prior probabilites dictionary has weird length: " + str(len(prior_dict))
	print "Get all probabilites after " + str(time.clock() - start_time)
	
	state_counter, segment_dictionary, chrom_name_list = nbh.get_segment_data_into_dictionary(segment_fn, binary_files)
	# print "binary_files"
	# print binary_files
	# print "chrom_name_list"
	# print chrom_name_list
	# print
	# print "segment_dictionary:"
	# print segment_dictionary
	# print "Get segmentation data after " + str(time.clock() - start_time)
	# print "Chrom nam list: ", chrom_name_list
	
	binary_data_dict = get_binary_data(binary_input_folder, binary_files, chrom_name_list)
	# print "binary_data_dict"
	# print binary_data_dict
	print "There are "+ str(len(chrom_name_list)) + " chromosomes files" 
	print "Get binary data after " + str(time.clock() - start_time)
	

	#### HA:starting from here: allow existing marks to be incldue. Remember to write existing marks into output file
	naive_bayes(binary_data_dict, chrom_name_list, num_marks_to_pick, cond_prob, inverse_cond_prob, prior_dict, \
	state_counter, segment_dictionary, output_file, NUM_STATES, NUM_MARKS, included_marks)

def usage():
	print "python naive_bayes_version2.py <cond_prob_fn> <prior_fn> <binary_input_folder> <binary_file_list_fn> <segment_fn> <output_file> <num_marks_to_pick>"
	print "This program will run naive bayes and find the best marks based on naive bayes procedure. Mark indices will be written into output file name"
	print "Example:"
	print "python naive_bayes_version2.py ../../chrom_mark_selection_project/naive_bayes/trial_naive_bayes/conditional_probabilities.txt ../../chrom_mark_selection_project/naive_bayes/trial_naive_bayes/prior_probabilities.txt ../../chrom_mark_selection_project/naive_bayes/trial_naive_bayes/bin_data/ ../../chrom_mark_selection_project/naive_bayes/trial_naive_bayes/bin_file_name ../../chrom_mark_selection_project/naive_bayes/trial_naive_bayes/segments.bed ../../chrom_mark_selection_project/naive_bayes/trial_naive_bayes/nb_output 1 ../../chrom_mark_selection_project/naive_bayes/trial_naive_bayes/chosen_marks 1"
	print "Another example:"
	print "pypy naive_bayes_version2.py ../../chrom_mark_selection_project/naive_bayes/conditional_probabilities.txt  ../../chrom_mark_selection_project/naive_bayes/prior_probabilities.txt /u/project/ernst/ernst/STACKEDCHROMHMM/BINARY_ROADMAP_NARROWPEAKCONSOLIDATED ../../chrom_mark_selection_project/binary_files_list ../../chrom_mark_selection_project/jason_original_model/genome_100_segments.bed ../../chrom_mark_selection_project/naive_bayes/top_160Marks_included_marks  80 ../../chrom_mark_selection_project/naive_bayes/decrease_40fold/top_80Marks 80"
	exit(1)
main()


