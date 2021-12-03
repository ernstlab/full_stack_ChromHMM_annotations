# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import pandas as pd 
import numpy as np 
import seaborn as sns
import sys
import os
import helper as chromHelp
MAX_NUM_TOP_STATE = 10 # for each enrichment context, how many most enriched states do we pick 
clustered_repeat_elements = ['acro_family.hg19', 'rRNA_class.hg19', 'rRNA_family.hg19', 'tRNA_class.hg19', 'telo_family.hg19', 'Merlin_family.hg19', 'Unknown?_class.hg19', 'Unknown?_family.hg19', 'Deu_family.hg19', 'PiggyBac?_family.hg19', 'L1?_family.hg19', 'scRNA_class.hg19', 'scRNA_family.hg19', 'snRNA_class.hg19', 'snRNA_family.hg19', 'TcMar_family.hg19', 'hAT?_family.hg19', 'Helitron?_family.hg19', 'DNA_family.hg19', 'TcMar?_family.hg19', 'CR1_family.hg19', 'Unknown_class.hg19', 'Unknown_family.hg19', 'ERVL?_family.hg19', 'L2_family.hg19', 'MIR_family.hg19', 'ERVL_family.hg19', 'Gypsy?_family.hg19', 'Gypsy_family.hg19', 'LTR_family.hg19', 'hAT-Tip100_family.hg19', 'hAT-Blackjack_family.hg19', 'hAT_family.hg19', 'Alu_family.hg19', 'SINE_class.hg19', 'repeats.gz', 'L1_family.hg19', 'LINE_class.hg19', 'PiggyBac_family.hg19', 'TcMar-Mariner_family.hg19', 'TcMar-Tigger_family.hg19', 'MuDR_family.hg19', 'RTE_family.hg19', 'TcMar-Tc2_family.hg19', 'Helitron_family.hg19', 'RC_class.hg19', 'DNA_class.hg19', 'hAT-Charlie_family.hg19', 'LINE?_class.hg19', 'Penelope?_family.hg19', 'LTR?_class.hg19', 'LTR?_family.hg19', 'SINE?_class.hg19', 'SINE?_family.hg19', 'Dong-R4_family.hg19', 'DNA?_class.hg19', 'DNA?_family.hg19', 'RTE-BovB_family.hg19', 'SINE_family.hg19', 'ERV_family.hg19', 'RNA_class.hg19', 'RNA_family.hg19', 'srpRNA_class.hg19', 'srpRNA_family.hg19', 'Other_class.hg19', 'Other_family.hg19', 'ERVK_family.hg19', 'ERV1_family.hg19', 'ERVL-MaLR_family.hg19', 'LTR_class.hg19', 'Simple_repeat_class.hg19', 'Simple_repeat_family.hg19', 'tRNA_family.hg19', 'Low_complexity_class.hg19', 'Low_complexity_family.hg19', 'Satellite_family.hg19', 'Satellite_class.hg19', 'centr_family.hg19']

clustered_repeat_elements = map(lambda x: x + '.bed.gz' if (x != 'repeats.gz') else x, clustered_repeat_elements)

def fix_one_enr_context_name(enr_cont):
	fixed_enr_cont = enr_cont # by default, no need to change anything
	if 'non_coding' in enr_cont or 'whole_genome' in enr_cont: # that means enr_cont can be FATHMM_non_coding_top_0.01, so we want FATHMM instead
		enr_cont_list = enr_cont.split('_')
		enr_cont_list = enr_cont_list[:-4] # get rid of the last 4 items because 
		fixed_enr_cont = '_'.join(enr_cont_list)
	if '_class.hg19.bed.gz' in enr_cont: # repeat classes
		fixed_enr_cont = enr_cont.split('_class.hg19.bed.gz')[0]
	return fixed_enr_cont

def get_enrichment_df(enrichment_fn):
	enrichment_df = pd.read_csv(enrichment_fn, sep = "\t", header = 0)
	enrichment_df = enrichment_df.rename(columns = {u'state (Emission order)': 'state', u'Genome %' : 'percent_in_genome'})
	(num_state, num_enr_cont) = (enrichment_df.shape[0] - 1, enrichment_df.shape[1] - 2) # the last row is the base, the first two columns are states and percent_in_genome 
	enr_cont_list = list(map(fix_one_enr_context_name, enrichment_df.columns[2:])) # fix the names of the enrichment context 
	enrichment_df.columns = list(enrichment_df.columns[:2]) + enr_cont_list
	percent_genome_of_cont = enrichment_df.iloc[num_state, 2:]
	enrichment_df = enrichment_df.loc[:(num_state - 1)] # get rid of the last row which is a row that indicates the percentage of genome each enrichment context occupies
	# num_enr_cont : number of enrichment context (TSS, exon, etc.)
	return enrichment_df, num_state, num_enr_cont, percent_genome_of_cont


def get_one_genome_context_top_states(df, col_name, top_num):
	sorted_gc, sorted_index = (pd.Index(df[col_name])).sort_values\
	(ascending = False, return_indexer = True)
	return sorted_index[:top_num] # state_indices here are all zero-based

def color_top_five(column, num_top_state):
	rank = column.rank(ascending = False)
	color_format = []
	for row_index, row_rank in enumerate(rank):
		if row_rank > MAX_NUM_TOP_STATE or row_rank > num_top_state: # greater than 5
			color_format += ['background-color: white']
			continue
		if row_rank == 1:
			color_format += ['background-color: #E74C3C'] 
		elif row_rank == 2: 
			color_format += ['background-color: #EC7063']
		elif row_rank == 3: 
			color_format += ['background-color: #F1948A']
		elif row_rank == 4: 
			color_format += ['background-color: #F5B7B1']
		elif row_rank == 5: 
			color_format += ['background-color: #FADBD8']
		elif row_rank > 5 and row_rank <=10:
			color_format += ['background-color: #faeae8']
	return color_format

def color_state_annotation(row_data, index_to_color):
	results = [""] * len(row_data) # current paint all the cells in the rows with nothing, no format yet
	state_annot_color = row_data['color']
	results[index_to_color] = 'background-color: %s' % state_annot_color # the third cell from the left is the state annotation cells
	return results

def get_colored_excel_top5_state(enrichment_fn, output_excel_fn, state_annot_df, num_top_state):
	enrichment_df, num_state, num_enr_cont, percent_genome_of_cont = get_enrichment_df(enrichment_fn)
	enrichment_cont_list = enrichment_df.columns[2:] # list of enrichment context that we have enrichment data for
	
	joint_top_states = []
	for enrichment_cont in enrichment_cont_list:
		top_state_this_cont = get_one_genome_context_top_states(enrichment_df, enrichment_cont, num_top_state) # get the list of top states most enriched with this enrichment context
		joint_top_states += list(top_state_this_cont)
	joint_top_states = np.unique(joint_top_states) # get the set of top enriched states to include
	# now get the enrichment data that include only the top most enriched states
	top_states_df = enrichment_df.loc[joint_top_states, :] # joint_top_states contain zero-based state indices
	top_states_df['state'] = (top_states_df['state']).astype(str).astype(int) # convert the state column from object to int64
	# join with state information 
	top_states_df = top_states_df.merge(state_annot_df, how = 'left', left_on = 'state', right_on = 'state', suffixes = ('_x', '_y'))
	top_states_df = top_states_df.sort_values(by = 'state_order_by_group')
	print(top_states_df.head())
	print ("Done getting the data of the most enriched states")
	# create a painted excel file from this
	colored_df = top_states_df.style.apply(lambda x: color_top_five(x, num_top_state), subset = pd.IndexSlice[:, top_states_df.columns[2 : 2 + num_enr_cont]]) # color each of the columsn corresponding to the enrichment contexts, color the rank of the states
	mneumonics_index_in_row = top_states_df.columns.get_loc('mneumonics') # column index, zero-based of the menumonics entries, which we will use to paint the right column later
	colored_df = colored_df.apply(lambda x: color_state_annotation(x, mneumonics_index_in_row), axis = 1)
	# annot_df = pd.DataFrame(columns = top_states_df.columns)
	# annot_df.loc[annot_df.shape[0]] = ['',''] + list(map(lambda x: x.split('_')[0], annot_df.columns[2:-3])) + ['','','']
	# annot_df.loc[annot_df.shape[0]] = ['',''] + list(map(lambda x: '_'.join(x.split('_')[1:-1]), annot_df.columns[2:-3])) + ['','','']
	# save to excel 
	writer = pd.ExcelWriter(output_excel_fn, engine='xlsxwriter')
	colored_df.to_excel(writer, sheet_name='top_states')
	# annot_df.to_excel(writer, sheet_name = 'annot_columns')
	writer.save()
	print ("Done styling the excel file")
	return top_states_df
	
def main() :
	if len(sys.argv) != 4:
		usage()
	enrichment_fn = sys.argv[1]
	chromHelp.check_file_exist(enrichment_fn)
	output_excel_fn = sys.argv[2]
	chromHelp.create_folder_for_file(output_excel_fn)
	num_top_state = chromHelp.get_command_line_integer(sys.argv[3])
	# if (num_top_state > MAX_NUM_TOP_STATE):
	# 	print "num_top_state has be be at most: " + str(MAX_NUM_TOP_STATE)
	# 	usage()
	# get data about state annotation
	state_annot_fn = sys.argv[4]
	chromHelp.check_file_exist(state_annot_fn)
	print ('Done getting command line arguments')
	state_annot_df = pd.read_csv(state_annot_fn, sep = ',', header = 0)
	state_annot_df = state_annot_df[['state', 'color', 'mneumonics', 'state_order_by_group']]
	print ("Done getting command line arguments and getting the state annotation data frame")
	# get the excel file that we want 
	top_states_df = get_colored_excel_top5_state(enrichment_fn, output_excel_fn, state_annot_df, num_top_state)


def usage():
	print ("python get_top_enriched_states_one_enrichment_file.py")
	print ("enrichment_fn")
	print ("output_excel_fn")
	print ("num_top_state: number of states that we want to collect for each enrichment context. The states that get selected have to be in the top states in at least one annotation in the file")
	print ('state_annot_fn: annotations of the states')
	exit(1)

main()