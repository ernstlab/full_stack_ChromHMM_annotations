print ("Given the file that has information about semgentation, gene expression and TSS regions' boundaries, we want to be able to aggregate data about the average gene expression in each state, along each bin defined through its distance from the TSS")
print ("")

import sys
import os
import gzip
import replicate_helper as rhelp
import pandas as pd 
import numpy as np 
import multiprocessing as mp
TSS_BIN_INDEX = rhelp.HALF_WINDOW_SIZE / rhelp.BIN_SIZE # the index of the bin that contains the tss in the window region that we consider
def get_cell_type_list(celltype_list_fn):
	'''
	get list of cell types that are present in the gene_exp_TSS_fn. The order of cell types in the two files should be the same. 
	return: list --> cell type names
	'''
	NUM_LINE_TO_SKIP_CELL_TYPE_FILE = 1 # number of the first few lines that need to be skipped when we look at the cell type file and get the list of cell types that correspond to the cell types in gene_exp_TSS_fn
	ctF = rhelp.open_file(celltype_list_fn)
	ct_list = ctF.readlines()[NUM_LINE_TO_SKIP_CELL_TYPE_FILE:] # we exclude the first NUM_LINE_TO_SKIP_CELL_TYPE_FILE because those are lines that are not related to the cell types
	ct_list = list(map(lambda x: x.strip(), ct_list)) # for each line, strip the trailing spaces
	ctF.close()
	return ct_list

def get_start_end_index_no_strand_exp(strand, tss_bp, start_bp_this_segment, end_bp_this_segment):
	left_bp_tss = rhelp.round_down_to_the_hundred(tss_bp - rhelp.HALF_WINDOW_SIZE) # the smaller index of the tss region, does notmean it's the start of the tss region because the start is determined also by the strand that the tss is on.
	right_bp_tss = rhelp.round_up_to_the_hundred(tss_bp + rhelp.HALF_WINDOW_SIZE) # the larger index 
	# calculate start and end indices of regardless of the strand
	start_index_no_strand = (start_bp_this_segment - left_bp_tss) / rhelp.BIN_SIZE
	end_index_no_strand = (end_bp_this_segment - left_bp_tss) / rhelp.BIN_SIZE
	# start and end indices when we take into account the strand
	if strand == 1:
		start_index_to_report = int(start_index_no_strand)
		end_index_to_report = int(end_index_no_strand)
	if strand == -1: # when strand is -1, the place to add gene expression will be flipped. For the plot if the TSS is on the negative strand, then those positions with greater coordinates would go to the left of the TSS in the plot and those with less than the TSS would go the right of the TSS, which is the opposite if it is on the positive strand.
		start_index_to_report = int(rhelp.TOTAL_NUM_BINS_INCLUDING_TSS - end_index_no_strand) 
		end_index_to_report = int(rhelp.TOTAL_NUM_BINS_INCLUDING_TSS - start_index_no_strand)
	return start_index_to_report, end_index_to_report # these are the indices of places in the array where we have to add the gene expression of this gene to all the cell type's total gene expression that overlap with this state.

def get_gene_exp_one_state_by_distance_from_TSS(gene_exp_TSS_folder, ct_list, output_folder, one_based_state_index, header_list): 
	# 1. Get the input data of the positions that overlap with the states, and the gene expression of genes whose TSS-relative distance overlap with the state
	this_state_gene_exp_fn = os.path.join(gene_exp_TSS_folder, "state_" + str(one_based_state_index) + "_gene_exp.gz")		
	this_state_all_data_df = pd.read_csv(this_state_gene_exp_fn, sep = '\t', header = None)
	this_state_all_data_df.columns = header_list
	this_state_all_data_df[ct_list] = np.log(this_state_all_data_df[ct_list] + 1) # TRANSFORM EXPRESSION DATA BY ADDING 1 AND TAKING THE LOG
	# 2. Declare the output data frame
	this_state_exp_df = pd.DataFrame(0, index = np.arange(rhelp.TOTAL_NUM_BINS_INCLUDING_TSS), columns = ct_list) # rows: genomic positions relative to the TSS (TSS is always in the middle), columns: different cell types, entries: avg_gene expression that this state fall into, given that the state at each position relative to the TSS of the gene
	# 3. Declare an important data array
	num_gene_overlap_state_each_bin = np.zeros(rhelp.TOTAL_NUM_BINS_INCLUDING_TSS) # number of genes whose regions from the TSS this state overlap with. We need to get the average gene expression at each position later

	for row_index, row_data in this_state_all_data_df.iterrows(): # loop through each row of this state , aka each place that this state overlap with
		tss_bp = row_data['tss']
		this_strand = row_data['strand'] 
		start_bp_this_segment = max(row_data['start_bp_segment'], row_data['start_tss_region']) 
		end_bp_this_segment = min(row_data['end_bp_segment'], row_data['end_tss_region'])
		this_gene_exp = row_data[ct_list] # --> numpy array of gene expression of this gene in these cell types
		start_index_to_report, end_index_to_report = get_start_end_index_no_strand_exp(this_strand, tss_bp, start_bp_this_segment, end_bp_this_segment) # --> indices of rows of the this_state_exp_df where we would need to add the gene expression of this gene
		# now add the gene expression to each position that we need to add
		for add_row_index in range(start_index_to_report, end_index_to_report): 
			this_state_exp_df.loc[add_row_index] = (this_state_exp_df.loc[add_row_index]).copy() + this_gene_exp 
		# now add 1  to each position that we just added the gene expression, meaning that we just found a gene whose TSS_relative positions overlap with the state that we are looking at
		num_gene_overlap_state_each_bin[start_index_to_report:end_index_to_report] += 1

	# now we are done processing all the genes. We have to get the average gene expression for all the positions
	num_gene_overlap_state_each_bin[num_gene_overlap_state_each_bin == 0] = 1 # any position that does not overlap with any gene will have denominator changed to 1, so as to avoid division by 0
	this_state_exp_df = this_state_exp_df.divide(num_gene_overlap_state_each_bin, axis = 'rows') # get the average gene expression at each position that this state overlap with 50kb window around annotated TSSs.
	# Finally, save the file of results related to this state
	this_state_save_fn = os.path.join(output_folder, "state_" + str(one_based_state_index) + "_avg_exp_tss_relative.gz")
	this_state_exp_df.to_csv(this_state_save_fn, index = False, header = True, sep = '\t') # rows: genomic positions relative to the TSS (TSS is always in the middle), columns: different cell types, entries: avg_gene expression that this state fall into, given that the state at each position relative to the TSS of the gene
	return

def get_gene_exp_by_state_by_distance_from_TSS_all_states(gene_exp_TSS_folder, output_folder, ct_list, state_list_this_process, num_state):
	FIRST_FEW_HEADERS = ['chr', 'start_bp_segment', 'end_bp_segment', 'state', 'chr_again', 'start_tss_region', 'end_tss_region', 'gene_id', 'strand', 'tss'] # the first set of start and end coordinate correspond to segmentation data, the second set corresponds to the TSS regions, this is the output of bedtools intersect
	CHROM_LIST=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
	# put the gene expression data into data frame in pandas
	# gene_exp_df = pd.read_csv(gene_exp_TSS_fn, sep = '\t', header = None)
	header_list = FIRST_FEW_HEADERS + ct_list
	# gene_exp_df.columns = header_list
	# print "Done getting in the data from gene_exp_TSS_folder after: " + str(time.clock() - start_time)
	# now declare all the dataframes for all cell types. These data frames should be stored inside 
	state_header = list(map(lambda x: "state" + str(x+1), range(num_state))) # --> ['state1', 'state2, 'state3',....]
	for state_index in state_list_this_process:
		one_based_state_index = state_index + 1	
		get_gene_exp_one_state_by_distance_from_TSS(gene_exp_TSS_folder, ct_list, output_folder, one_based_state_index, header_list) # Will process and print out a data frame with rows: genomic positions relative to the TSS (TSS is always in the middle), columns: different cell types, entries: avg_gene expression that this state fall into, given that the state at each position relative to the TSS of the gene
		# the data frame will later be opened and copied into better forms that are cell type specific.
		print ("Done processing data for state: " + str(one_based_state_index))
	return

def distribute_states_for_processes(num_state, num_processes, output_folder):
	'''
	Given the number of states and number of processes, try to divide even the state indices to each process. That way, each process has a list of states that it will process. 
	return: list of list : process_index --> list of states that this process should work on
	'''
	remaining_states = range(num_state) # find states that have not generated all the output files yet
	# for state_i in range(num_state):
	# 	fn = os.path.join(output_folder, "state_" + str(state_i + 1) + "_avg_exp_tss_relative.gz")
	# 	if not os.path.isfile(fn):
	# 		remaining_states.append(state_i)
	num_state_per_process = int(len(remaining_states) / num_processes)
	process_state_list = []
	for process_index in range(num_processes - 1):
		start_state_index_this_process = process_index * num_state_per_process
		end_state_index_this_process = (process_index + 1) * num_state_per_process
		process_state_list.append(remaining_states[start_state_index_this_process:end_state_index_this_process])
	last_start_state_index = (num_processes - 1) * num_state_per_process
	last_end_state_index = num_state
	process_state_list.append(remaining_states[last_start_state_index:])
	return process_state_list

def combine_data_from_state_to_cell_type(output_folder, num_state, ct_list):
	# we got the average gene expression at TSS-relative positions for all cell types. Now, we need to add this data to the right data frame of each cell type. For each cell-type data frame, we will need to add the average gene expression data that we just got from this state	
	header_list_for_ct_exp_df = list(map(lambda x: 'state_' + str(x+1), range(num_state)))
	celltype_df_list = []
	for ct_index, ct in enumerate(ct_list):
		celltype_df_list.append(pd.DataFrame(0, index = np.arange(rhelp.TOTAL_NUM_BINS_INCLUDING_TSS), columns = header_list_for_ct_exp_df)) # for each cell type, create a data frame the will report the final results for each cell type. rows: genomic positions relative to the TSS (TSS is always in the middle), columns: different chromHMM states
	for state_index in range(num_state):
		one_based_state_index = state_index + 1
		this_state_gene_exp_fn = os.path.join(output_folder, "state_" + str(one_based_state_index) + "_avg_exp_tss_relative.gz")
		this_state_gene_exp_df = pd.read_csv(this_state_gene_exp_fn, header = 0, sep = '\t')
		this_state_header = 'state_' + str(one_based_state_index)
		for ct_index, ct in enumerate(ct_list):
			(celltype_df_list[ct_index])[this_state_header] = this_state_gene_exp_df[ct] # copy data from the state df for this cell type, to this ct df in this state
		print ("Done copying data for state: " + str(one_based_state_index))
	# now save each of the ct df
	for ct_index, ct in enumerate(ct_list):
		this_ct_exp_fn = os.path.join(output_folder, ct + "_avg_exp_tss_relative.gz")
		(celltype_df_list[ct_index]).to_csv(this_ct_exp_fn, index = False, header = True, sep = '\t')
	print ("Done saving all cell type df ")

def process_gene_exp_data_multiple_processes(gene_exp_TSS_folder, output_folder, ct_list, num_state, num_processes):
	process_state_list = distribute_states_for_processes(num_state, num_processes, output_folder) # list of list : process_index --> list of states that this process should work on
	# get_gene_exp_by_state_by_distance_from_TSS_all_states(gene_exp_TSS_folder, output_folder, ct_list, process_state_list[0], num_state)
	processes = []
	# call each process
	for process_index in range(num_processes):
		state_list_this_process = process_state_list[process_index]
		p = mp.Process(target = get_gene_exp_by_state_by_distance_from_TSS_all_states, args = (gene_exp_TSS_folder, output_folder, ct_list, state_list_this_process, num_state))
		processes.append(p)
		p.start()
		print ("Started process: " + str(process_index))
	# now wait for processes to finish:
	for p in processes:
		p.join()
	print ("Done getting all states' expression data")

	# # now we just have to combine the data from all states into data from all cell types
	combine_data_from_state_to_cell_type(output_folder, num_state, ct_list)

def main():
	if len(sys.argv) != 6:
		usage()
	gene_exp_TSS_folder = sys.argv[1]
	if not os.path.isdir(gene_exp_TSS_folder):
		print ("gene_exp_TSS_folder: " + gene_exp_TSS_folder + " DOES NOT EXIST")
		usage()
	output_folder = sys.argv[2] # rows: states, columns: bins based on distance from the TSS, entries: average gene expression of the genes that this state overlaps with, as the bin that we are looking at, relative to the distance from the TSS (note: we have specific rules to apply whether the gene is on negative or positive strands)
	rhelp.make_dir(output_folder)	
	celltype_list_fn = sys.argv[3]
	if not os.path.isfile(celltype_list_fn):
		print ("celltype_list_fn: " + celltype_list_fn + " DOES NOT EXIST")
		usage()
	num_state = int(sys.argv[4])
	num_processes = int(sys.argv[5])
	print ("Done getting command lines")
	ct_list = get_cell_type_list(celltype_list_fn)
	print ("Done get_cell_type_list ")
	process_gene_exp_data_multiple_processes(gene_exp_TSS_folder, output_folder, ct_list, num_state, num_processes)

def usage():
	print ('python get_gene_exp_by_state_by_distance_from_TSS.py')
	print ("gene_exp_TSS_folder: result of step 3 of the pipeline written in pipeline_replicate_ernst_etal_2011_supp_fig2.sh")
	print ('output_folder: folder where all output_fn will be stored. Each output fn is corresponding to one cell type\'s data frame. rows: states, columns: bins based on distance from the TSS, entries: average gene expression of the genes that this state overlaps with, as the bin that we are looking at, relative to the distance from the TSS (note: we have specific rules to apply whether the gene is on negative or positive strands)')
	print ('celltype_list_fn: fn that contains cell types that correspond to the order presetned in the input gene_exp_TSS_folder. It is assumed that (1) the file is new-line-separated, (2) the first line is gene_id, so we will skip the first line')
	print ('num_state')
	print ('num_processes')
	print ("")
	exit(1)

main()