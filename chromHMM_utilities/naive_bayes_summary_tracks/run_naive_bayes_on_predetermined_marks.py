import math
import sys
import time
import os
from collections import Counter
import gzip
NUM_MARKS = 4649
NUM_STATES = 100
SEGMENT_LENGTH = 200
"""
state_counter: a counter dictionary [state] --> number of times that this state got appeared in the segmentation
state_recovery_counter : [state][mark]: how many times that this staete got recovered when this mark is included in the chosen marks
chosen_marks: counter dictionary [mark index] --> 1 or 0
binary_data_dict: dictionary of 2D tables: key: chomosome name, values: talbe of binary d atae [position][mark index]
prior_dict: dictionary {key: state index, values: prior probabilities of that state}
cond_prob and inverse_cond_prob: 2D table [state_index][mark index]
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

def most_likely_state(current_likelihood, marks, cond_prob, inverse_cond_prob, mark_index):
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

def get_segment_data_into_dictionary(segment_fn, binary_files):
	chrom_name_list = []
	state_counter = Counter([])

	for i, bin_fn in enumerate(binary_files):
		name_list = bin_fn.split("_")
		chrom_name = ".".join([name_list[0], name_list[1]])
		chrom_name_list.append(chrom_name)
	chr_name_counter = Counter(chrom_name_list)
	segmentF = open(segment_fn, 'r')
	segment_dictionary = {}
	for line in segmentF:	
		[chrom, start_index, end_index, state] = line.strip().split()
		if chr_name_counter[chrom] > 0: # if this state is is the list of binary data file that we want to process
			state_index = int(state[1:]) - 1 # because raw state would be like 'E100'
			num_appearances = (int(end_index) - int(start_index)) / SEGMENT_LENGTH
			state_counter.update({state_index:num_appearances})
			if chrom in segment_dictionary:
				(segment_dictionary[chrom]).extend([state_index] * num_appearances)
			else: 
				(segment_dictionary[chrom]) = [state_index] * num_appearances
	return state_counter, segment_dictionary, chrom_name_list

def find_best_mark(chosen_marks, state_recovery_counter, state_counter):
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

def calculate_recovery(state_recovery_counter, state_counter):
	recovery = 0
	for state_index in range(NUM_STATES):
		try:
			recovery += float(state_recovery_counter[state_index]) / float(state_counter[state_index])
		except: 
			pass
	return recovery

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

def one_round_naive_base_to_find_best_mark(current_likelihood_dict, mark_index, segment_dictionary, cond_prob, inverse_cond_prob, prior_dict,\
 state_counter, binary_data_dict, chrom_name_list):
	state_recovery_counter = [0] * NUM_STATES
	for chrom_name in chrom_name_list:
		binary_data = binary_data_dict[chrom_name]
		this_chrom_current_likelihood_table = current_likelihood_dict[chrom_name]
		assert len(binary_data) == len(this_chrom_current_likelihood_table), "This chromosome: " + chrom_name + \
		" does not have consistent lenth: " + str(len(binary_data)) + "   " + str(len(this_chrom_current_likelihood_table))
		this_chrom_state_list = segment_dictionary[chrom_name]
		for i, position in enumerate(binary_data):
			this_pos_state = this_chrom_state_list[i]
			this_pos_current_likelihood = this_chrom_current_likelihood_table[i]
			max_state, max_log_likelihood = most_likely_state(this_pos_current_likelihood, position, cond_prob, inverse_cond_prob, mark_index)
			if max_state == this_pos_state:
				state_recovery_counter[max_state] += 1

	this_mark_recovery = calculate_recovery(state_recovery_counter, state_counter)
	return mark_index, this_mark_recovery

def initialize_likelihood(chrom_name_list, prior_dict, binary_data_dict):
	current_likelihood_dict = {}
	state_prior_list = []
	for state_index in range(NUM_STATES):
		state_prior_list.append(prior_dict[state_index])
	for i, chrom_name in enumerate(chrom_name_list):
		this_chrom_prior_table = [] # [position][state index]
		num_pos = len(binary_data_dict[chrom_name])
		for j in range (num_pos):
			this_chrom_prior_table.append(state_prior_list)
		current_likelihood_dict[chrom_name] = this_chrom_prior_table
	return current_likelihood_dict


def naive_bayes(binary_data_dict, chrom_name_list, mark_indices, cond_prob, inverse_cond_prob, prior_dict, \
	state_counter, segment_dictionary, output_file, start_time):
	print "Starting to run naive bayse ..."
	current_likelihood_dict = initialize_likelihood(chrom_name_list, prior_dict, binary_data_dict)
	print "Done with getting the current_likelihood_dict table: time passed: " + str(time.clock() - start_time)
	for mark_index in mark_indices:
		#print "finding mark " + str(i)
		index, recovery = one_round_naive_base_to_find_best_mark(current_likelihood_dict, mark_index, segment_dictionary, cond_prob, \
			inverse_cond_prob, prior_dict, state_counter, binary_data_dict, chrom_name_list)
		# update the current likelihood dictionary
		for chrom_name,this_current_likelihood in current_likelihood_dict.iteritems():
			new_chrom_current_likelihood = []
			for position_index in range(len(this_current_likelihood)):
				new_position_likelihood = []
				if binary_data_dict[chrom_name][position_index][mark_index] == 1: 
					for state_index in range(NUM_STATES):
						new_position_likelihood.append(this_current_likelihood[position_index][state_index] + cond_prob[state_index][mark_index])

				else:
					for state_index in range(NUM_STATES):
						new_position_likelihood.append(this_current_likelihood[position_index][state_index] + inverse_cond_prob[state_index][mark_index])
				new_chrom_current_likelihood.append(new_position_likelihood)
			current_likelihood_dict[chrom_name] = new_chrom_current_likelihood

		output_file.write(str(mark_index) + "\t" + str(recovery) + "\n")
		output_file.flush()
		print "Found mark "+ str(mark_index) + "   " + str(recovery)
		print "time passed: " + str(time.clock() - start_time)
	output_file.close()

def get_binary_data(binary_input_folder, binary_files, chrom_name_list):
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

def get_mark_indices(marks_fn):
	results = []
	marksF = open(marks_fn, 'r')
	for line in marksF:
		index = int(line.strip()[5:])
		results.append(index)
	marksF.close()
	return results

def main():
	cond_prob_fn = sys.argv[1]
	prior_fn = sys.argv[2]
	binary_input_folder = sys.argv[3]
	segment_fn = sys.argv[4]
	output_file = open(sys.argv[5], 'w', 0) # no buffering of files
	output_file.write("Hello!")
	marks_fn = sys.argv[6]
	start_time = time.clock()
	binary_files = os.listdir(binary_input_folder)
	
	cond_prob, inverse_cond_prob = get_state_mark_cond_prob(cond_prob_fn)
	assert len(cond_prob) == NUM_STATES and len(inverse_cond_prob) == NUM_STATES, "Length of cond_prob and inverse_cond_prob are wierd: " + str(len(cond_prob)) + " " + str(len(inverse_cond_prob))
	for i in range(NUM_STATES):
		assert len(cond_prob[i]) == NUM_MARKS and len(inverse_cond_prob[i]) == NUM_MARKS, "Length of cond_prob and inverse_cond_prob at state " + str(i) + " are weird: " + str(len(cond_prob[i])) + " " + str(len(inverse_cond_prob[i]))
	print "get all conditional probabilites after " + str(time.clock() - start_time)

	prior_dict = get_prior_probabilities(prior_fn)
	assert len(prior_dict) == NUM_STATES, "prior probabilites dictionary has weird length: " + str(len(prior_dict))
	print "Get all probabilites after " + str(time.clock() - start_time)
	
	state_counter, segment_dictionary, chrom_name_list = get_segment_data_into_dictionary(segment_fn, binary_files)
	print "Get segmentation data after " + str(time.clock() - start_time)
	print "Chrom nam list: ", chrom_name_list
	
	binary_data_dict = get_binary_data(binary_input_folder, binary_files, chrom_name_list)
	print "There are "+ str(len(chrom_name_list)) + " chromosomes files" 
	print "Get segmentation data after " + str(time.clock() - start_time)

	mark_indices = get_mark_indices(marks_fn)
	print "There are "+ str(len(mark_indices)) + " mark indices to go through" 
	print "Get mark indices data after " + str(time.clock() - start_time)
	naive_bayes(binary_data_dict, chrom_name_list, mark_indices, cond_prob, inverse_cond_prob, prior_dict, \
	state_counter, segment_dictionary, output_file, start_time)
main()

