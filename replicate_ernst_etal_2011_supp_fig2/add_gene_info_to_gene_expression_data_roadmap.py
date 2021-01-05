import sys
import os
import numpy as np 
import pandas as pd 
import replicate_helper as helper

#### README ####
'''
This file will take as input the raw gene expression data from road map with rows as genes and columns as cell types. 
It also take the gene information: locations of genes and stuffs like that
It will combine the two data into a format that looks like a bed file with the first three columns being chromosom, start coord, end coord 
The following columns will be information that we already have in the raw gene expression file. The output of this file should be put into a bedtools intersect command later to get the states that intersect with each of the genes
'''
SEGMENT_LENGTH = 200
def get_combined_state_gene_exp_df(raw_exp_fn, gene_info_fn):
	gene_info_df = pd.read_csv(gene_info_fn,  sep = "\t", header = None) 
	gene_info_df.columns = ['chr_segmentation', 'start_bp_segment', 'end_bp_segment', 'segment', 'chr_gene', 'start_bp_gene', 'end_bp_gene', 'strand', 'gene_id'] # column names of tables produced in the step 3 of pipeline_get_avg_gene_exp
	gene_info_df['gene_length'] = gene_info_df['end_bp_gene'] - gene_info_df['start_bp_gene']
	gene_info_df['gene_length_segments'] = gene_info_df['gene_length'] / SEGMENT_LENGTH
	gene_info_df['com_start_bp'] = gene_info_df[['start_bp_segment', 'start_bp_gene']].max(axis = 1) #start of the region we are interested is the start of either the gene or the segment of chromatin state. Whichever is greater is the start of the intersection. 
	gene_info_df['com_end_bp'] = gene_info_df[['end_bp_segment', 'end_bp_gene']].min(axis = 1) # same argument as above with the end of the region 
	gene_info_df['num_segments'] = (gene_info_df['com_end_bp'] - gene_info_df['com_start_bp']) / SEGMENT_LENGTH # number of segments that are shared between the genes and the chromatin state
	gene_info_df['num_segments'] = gene_info_df['num_segments'].apply(np.ceil) # round up the number of segments
	# Number of segments is the number of basepair in the region we are interested in / number of basepairs per segment bin
	gene_info_df = gene_info_df[[u'segment', u'strand', u'gene_id', u'gene_length', u'com_start_bp', u'com_end_bp', u'num_segments', 'gene_length_segments']] # get rid of redundant data
	raw_exp_df = pd.read_csv(raw_exp_fn, sep = '\t', header = 0, index_col = False)
	com_df = raw_exp_df.merge(gene_info_df, left_on = 'gene_id', right_on = 'gene_id', suffixes = (False, False)) # combine gene information df with raw_exp_df so that now we have a big tables with gene_exp in all cell types, and the information of segments of genes that are associated with different chromatin states. Example list of columns: u'gene_id', u'E000', u'E003', u'E004', u'E005', u'E006', u'E007',u'E011', u'E012', u'E013', u'E016', u'E024', u'E027', u'E028', u'E037', u'E038', u'E047', u'E050', u'E053', u'E054', u'E055', u'E056', u'E057', u'E058', u'E059', u'E061', u'E062', u'E065', u'E066', u'E070', u'E071', u'E079', u'E082', u'E084', u'E085', u'E087', u'E094', u'E095', u'E096', u'E097', u'E098', u'E100', u'E104', u'E105', u'E106', u'E109', u'E112', u'E113', u'E114', u'E116', u'E117', u'E118', u'E119', u'E120', u'E122', u'E123', u'E127', u'E128', u'segment', u'strand', u'gene_length', u'com_start_bp', u'com_end_bp', u'num_segments', u'gene_length_segments'
	return(com_df)

def get_correct_output_format(raw_exp_fn, gene_info_fn, output_folder, num_states):
	df = get_combined_state_gene_exp_df(raw_exp_fn, gene_info_fn) # result is a df that contains both gene expresions and the coordinate of the shared segment between the chromtatin state and the gene
	# segment: chromatin state 
	# gene_length: exact gene_length
	# com_start_bp: the max of the segment start and gene exact start bp
	# com_end_bp: the min of the segment end and gene exact end bp
	# num_segments: number of segments of intersection between genes and the chromatin states, calculated using com_start_bp and com_end_bp
	# gene_length_segments: number of segments that spand this gene, equal exact gene_length / SEGMENT_LENGTH, floored.
	print ("Done getting and pre-processing data" )
	# get the list of column names that starts with "E", which signifies that they are columns of gene expression in different cell types
	is_gene_exp_columns = list(map(lambda x: x.startswith("E") and x != 'E000', df.columns))
	gene_exp_columns =  df.columns[is_gene_exp_columns]
	# declare the data frame where we will store the average gene expression for each state in each cell types --> rows: states, columns: cell types
	res_index = list(map(lambda x: "E" + str(x+1), range(num_states))) # result index: states: E1 --> E100
	bin_uniform_result_df = pd.DataFrame(data = None, columns = gene_exp_columns, index = res_index) # Rows; states, columns: different cell types of the form Exxx. Sum over regions (num_bin of this region * gene_exp of this region) / sum over regions(num_bin of this state)
	bp_uniform_result_df = pd.DataFrame(data = None, columns = gene_exp_columns, index = res_index) # Sum over regions (num_bp of this region * gene_exp of this region) / sum over regions(num_bp of this state)
	for state_index in range(num_states):
		one_based_state_index = state_index + 1
		this_state_org_df = (df[df['segment'] == 'E' + str(one_based_state_index)]) # get the rows from the df that combines data of chromHMM states and gene expression that correspond to the specific state that we are looking at
		bin_unif_state_df = pd.DataFrame() # rows: positions on the genome that corresponds to this state. columns: cell types
		bp_unif_state_df = pd.DataFrame()
		for ct_col in gene_exp_columns: # loop through each column of gene expression in each cell type
			bin_unif_state_df[ct_col] = np.log(this_state_org_df[ct_col] + 1) * this_state_org_df['num_segments']
			bp_unif_state_df[ct_col] = np.log(this_state_org_df[ct_col] + 1) * this_state_org_df['num_segments'] / this_state_org_df['gene_length_segments']
		# now, bin_unif_state_df has columns: ct gene expressions, rows: different regions in the genome that correspond to the state we are looking at, cells: expressions * number of segments in each regions. Next, we have to calculate the weighted average gene expression, bin-uniformly for each state in each cell type
		total_segments = np.sum(this_state_org_df['num_segments']) # total number segments overlapping any genes AND this specific state
		total_weighted_segments = np.sum(this_state_org_df['num_segments'] / this_state_org_df['gene_length_segments'])
		avg_exp_all_ct_bin_unif = bin_unif_state_df.sum(axis = 0) # column sum
		avg_exp_all_ct_bin_unif = avg_exp_all_ct_bin_unif / total_segments # Sum over regions (num_bin of this region * gene_exp of this region) / sum over regions(num_bin of this region)
		avg_exp_all_ct_bp_unif = bp_unif_state_df.sum(axis = 0) # column sum
		avg_exp_all_ct_bp_unif = avg_exp_all_ct_bp_unif / total_weighted_segments # Sum over regions (num_bin of this region * gene_exp of this region) / sum over regions(num_bin of this region / num_bins associated with the gene in this region)

		# put data into the right place in result_df
		bin_uniform_result_df.loc['E' + str(one_based_state_index)] = avg_exp_all_ct_bin_unif # index of each row is actually 'E' + state_index
		bp_uniform_result_df.loc['E' + str(one_based_state_index)] = avg_exp_all_ct_bp_unif
		print ("Done processing data for state: " + str(one_based_state_index))
	# save data
	bin_uniform_result_df['state'] = bin_uniform_result_df.index
	bp_uniform_result_df['state'] = bp_uniform_result_df.index
	bin_uniform_result_df.to_csv(output_folder + "/segment_uniform_avg_exp.txt", sep ='\t', header = True, index = False)
	bp_uniform_result_df.to_csv(output_folder + "/bp_uniform_avg_exp.txt", sep = '\t', header = True, index = False)



def main():
	if len(sys.argv) != 5:
		usage()
	raw_exp_fn = sys.argv[1]
	if not os.path.isfile(raw_exp_fn):
		print ("raw_exp_fn: " + raw_exp_fn + " is NOT VALID")
		usage()
	gene_info_fn = sys.argv[2] # this file contains data about the genes coordinates, strand and the chromHMM state associated with each segment in the gene.
	if not os.path.isfile(gene_info_fn):
		print ("gene_info_fn: " + gene_info_fn + " is NOT VALID")
		usage()
	output_folder = sys.argv[3]
	helper.make_dir(output_folder)
	try:
		num_states = int(sys.argv[4])
	except:
		print ("num_states is NOT VALID")
		usage()
	#### DEBUGGING CODE ####
	# df = get_combined_state_gene_exp_df(raw_exp_fn, gene_info_fn)
	# print(df.head())
	# df = df.groupby('segment')
	# for state in range(100):
	# 	state_name = 'E' + str(state + 1)
	# 	this_state_org_df = df.get_group(state_name)

	########################
	get_correct_output_format(raw_exp_fn, gene_info_fn, output_folder, num_states)

def usage():
	print ("raw_exp_fn: raw gene expression data from roadmap")
	print ("gene_info_fn: information about genes' coordinates and strands and etc. from roadmap")
	print ('output_fn: chromosome, start coord, end coord, gene names, expressions in different cell types')
	print ("number of states")
	exit(1)

main()