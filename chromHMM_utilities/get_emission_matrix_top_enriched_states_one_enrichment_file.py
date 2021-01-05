import pandas as pd 
import numpy as np 
import sys
import os
import chromHMM_utilities_common_functions_helper as chromHelp
NUM_TOP_STATE = 5 # for each enrichment context, how many most enriched states do we pick 
def get_enrichment_df(enrichment_fn):
	enrichment_df = pd.read_csv(enrichment_fn, sep = "\t", header = 0)
	enrichment_df = enrichment_df.rename(columns = {u'state (Emission order)': 'state', u'Genome %' : 'percent_in_genome'})
	(num_state, num_enr_cont) = (enrichment_df.shape[0] - 1, enrichment_df.shape[1] - 1) 
	percent_genome_of_cont = enrichment_df.iloc[num_state, 2:]
	enrichment_df = enrichment_df.loc[:(num_state - 1)]
	# num_enr_cont : number of enrichment context (TSS, exon, etc.)
	return enrichment_df, num_state, num_enr_cont, percent_genome_of_cont

def get_emission_df(emission_fn): # this function only works for emission matrix that is the output of ChromHMM
	emission_df = pd.read_csv(emission_fn, sep = '\t', header = 0)
	emission_df = emission_df.rename(columns = {'state (Emission order)': 'state'}) # change the column name so that it gets easier to understand 
	return emission_df

def get_one_genome_context_top_states(df, col_name, top_num):
	sorted_gc, sorted_index = (pd.Index(df[col_name])).sort_values\
	(ascending = False, return_indexer = True)
	return sorted_index[:top_num] # state_indices here are all zero-based

def get_emission_matrix_top5_state(enrichment_fn, emission_fn, output_excel_fn):
	enrichment_df, num_state, num_enr_cont, percent_genome_of_cont = get_enrichment_df(enrichment_fn)
	enrichment_cont_list = enrichment_df.columns[2:] # list of enrichment context that we have enrichment data for
	joint_top_states = []
	for enrichment_cont in enrichment_cont_list:
		top_state_this_cont = get_one_genome_context_top_states(enrichment_df, enrichment_cont, NUM_TOP_STATE) # get the list of top states most enriched with this enrichment context
		joint_top_states += list(top_state_this_cont)
	joint_top_states = np.unique(joint_top_states) # get the set of top enriched states to include
	# get emission data
	emission_df = get_emission_df(emission_fn)
	# now get the emission data that include only the top most enriched states
	top_states_df = emission_df.loc[joint_top_states, :] # joint_top_states contain zero-based state indices
	top_states_df['state'] = (top_states_df['state']).astype(str).astype(int) # convert the state column from object to int64
	# now change the colname of state back to its original form, so that we have an emission matrix that looks exactly like what ChromHMM would give us
	top_states_df = top_states_df.rename(columns = {'state' : 'state (Emission order)'})
	# now save the subsetted emission matrix into a text file, so that we can do further analysis later
	top_states_df.to_csv(output_excel_fn, sep = '\t', header = True, index = False)
	return 

def main() :
	if len(sys.argv) != 4:
		usage()
	enrichment_fn = sys.argv[1]
	chromHelp.check_file_exist(enrichment_fn)
	emission_fn = sys.argv[2]
	chromHelp.check_file_exist(emission_fn)
	output_excel_fn = sys.argv[3]
	chromHelp.create_folder_for_file(output_excel_fn)
	print "Done getting command line arguments and getting the state annotation data frame"
	# get the excel file that we want 
	top_states_df = get_emission_matrix_top5_state(enrichment_fn, emission_fn, output_excel_fn)
	print "Done getting the subsetted emission matrix with top enriched states"


def usage():
	print "python get_emission_matrix_top_enriched_states_one_enrichment_file.py"
	print "enrichment_fn"
	print "emission_fn"
	print "output_excel_fn"
	exit(1)

main()