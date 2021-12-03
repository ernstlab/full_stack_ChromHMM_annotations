import pandas as pd 
import numpy as np 
import os
import sys
import helper
import glob
GO_TYPE_LIST = ['GO Biological Process', 'GO Cellular Component', 'GO Molecular Function', 'Human Phenotype']
GO_TYPE_MNE_LIST = ['GO_BP', 'GO_CC', 'GO_MF', 'HP']
def read_one_GREAT_output_file(fn):
	df = pd.read_csv(fn, comment = "#", sep = '\t', index_col = None, header = None, skiprows = 1)
	header_line = "# Ontology,ID,Desc,BinomRank,BinomP,BinomBonfP,BinomFdrQ,RegionFoldEnrich,ExpRegions,ObsRegions,GenomeFrac,SetCov,HyperRank,HyperP,HyperBonfP,HyperFdrQ,GeneFoldEnrich,ExpGenes,ObsGenes,TotalGenes,GeneSetCov,TermCov,Regions,Genes"
	header_list = header_line.split(',')
	header_list[0] = header_list[0][2:]
	df.columns = header_list
	return df

def get_one_state_one_GO_type_GO_list(df, GO_type, rank_filter_list):
	go_df = df[df['Ontology'] == GO_type]
	go_df = go_df[go_df['BinomRank'].isin(rank_filter_list)]
	id_list_str = '|'.join(go_df['ID'])
	desc_list_str = '|'.join(go_df['Desc'])
	fdr_p_list_str = '|'.join(list(map(lambda x: str(x), go_df['BinomFdrQ'])))
	enrich_list_str = '|'.join(list(map(lambda x: str(x), go_df['RegionFoldEnrich'])))
	return [id_list_str, desc_list_str, fdr_p_list_str, enrich_list_str]

def combine_result_with_state_annotation(result_df):
	state_annot_fn = '../../ROADMAP_aligned_reads/chromHMM_model/model_100_state/figures/supp_excel/state_annotations_processed.csv'
	state_annot_df = pd.read_csv(state_annot_fn, header = 0, index_col = None, sep = ',')
	state_annot_df = state_annot_df[['state', 'state_order_by_group', 'mneumonics', 'Group']]
	result_df['state'] = result_df['state'].apply(lambda x: int(x)) # so that we can match with the state_annot_df soon
	result_df = result_df.merge(state_annot_df, how = 'left', left_on = 'state', right_on = 'state')
	result_df = result_df.sort_values('state_order_by_group')
	return result_df

def get_plot_data_one_state_group_one_go_type(state_group_df, go_type, state_group_name, output_folder):
	# state_group_df corresponds to the df with each row corresponding to one state in this group 
	result_df = pd.DataFrame(columns = ['ID', 'Desc', 'state', 'FdrQ', 'FoldEnrich'])
	go_id_colname = go_type + "_" + 'ID'
	go_desc_colname = go_type + "_" + "Desc"
	go_fdrq_colname = go_type + "_" + "FdrQ"
	go_fold_colname = go_type + "_" + "FoldEnrich"
	for index, row in state_group_df.iterrows():
		this_row_df = pd.DataFrame() # an empty data frame that corresponds to this row
		go_id = row[go_id_colname]
		this_row_df['ID'] = go_id.split('|')
		go_desc = row[go_desc_colname]
		this_row_df['Desc'] = go_desc.split('|')
		go_fdrq = row[go_fdrq_colname]
		this_row_df['FdrQ'] = str(go_fdrq).split('|')
		go_fold = row[go_fold_colname]
		this_row_df['FoldEnrich'] = str(go_fold).split('|')
		this_row_df['state'] = [row['state']] * this_row_df.shape[0] # a list of values equal to this state, and the length 
		result_df = result_df.append(this_row_df, ignore_index = True)
	state_group_name = '_'.join(state_group_name.split())
	save_fn = os.path.join(output_folder, state_group_name + "_topGO.txt.gz")
	result_df.to_csv(save_fn, header = True, index = False, sep = '\t',  compression = 'gzip')
	return

def prepare_data_for_plotting(all_state_df, output_folder):
	# all_state_df corresponds to the result_df in function combine_GO_output_all_states
	grouped_state_df = all_state_df.groupby('Group')
	group_list = grouped_state_df.groups.keys()
	for GO_type in GO_TYPE_MNE_LIST:
		this_GO_output_folder = os.path.join(output_folder, GO_type)
		helper.make_dir(this_GO_output_folder)
		get_plot_data_one_state_group_one_go_type(all_state_df, GO_type, 'all_states', this_GO_output_folder)
	for state_group_name in group_list:
		df = grouped_state_df.get_group(state_group_name)
		for GO_type in GO_TYPE_MNE_LIST:
			this_GO_output_folder = os.path.join(output_folder, GO_type)
			this_group_df = get_plot_data_one_state_group_one_go_type(df, GO_type, state_group_name, this_GO_output_folder)
	return 

def combine_GO_output_all_states(all_state_folder, num_max_GO_terms, output_folder):
	rank_filter_list = list(map(lambda x: x+1, range(num_max_GO_terms)))
	all_state_fn_list = glob.glob(all_state_folder + "/state_E*_GO.tsv")
	all_state_list = list(map(lambda x: int(x.split('/')[-1].split('_')[1][1:]), all_state_fn_list))
	all_state_df_list = list(map(read_one_GREAT_output_file, all_state_fn_list))
	columns_list = ['state']
	for go_type in GO_TYPE_MNE_LIST:
		columns_list += list(map(lambda x: go_type + "_" + x, ['ID', 'Desc', 'FdrQ', 'FoldEnrich']))
	result_df = pd.DataFrame(columns = np.arange(len(columns_list))) # keep the column names as numbers to ensure correct appending of rows
	for index, state in enumerate(all_state_list):
		state_df = all_state_df_list[index]
		go_bp_results = get_one_state_one_GO_type_GO_list(state_df, 'GO Biological Process', rank_filter_list)
		go_cc_results = get_one_state_one_GO_type_GO_list(state_df, 'GO Cellular Component', rank_filter_list)
		go_mf_results = get_one_state_one_GO_type_GO_list(state_df, 'GO Molecular Function', rank_filter_list)
		hp_results = get_one_state_one_GO_type_GO_list(state_df, 'Human Phenotype', rank_filter_list)
		row_result = pd.Series([str(state)] + go_bp_results + go_cc_results + go_mf_results + hp_results)
		result_df = result_df.append(row_result, ignore_index = True)
	result_df.columns = columns_list
	result_df = combine_result_with_state_annotation(result_df)
	# result_df = pd.read_csv(os.path.join(output_folder, 'top' + str(num_max_GO_terms) + '_GO_all_states.txt'), header = 0, index_col = None, sep = '\t')
	print('Done combining data for all states')
	# this is to create a bunch of files for plotting the data
	prepare_data_for_plotting(result_df, output_folder)
	print('Done getting data to make different plots')
	# this is to write down all the data that we have aggregated
	output_fn = os.path.join(output_folder, 'top' + str(num_max_GO_terms) + '_GO_all_states.txt.gz')
	result_df.to_csv(output_fn, header = True, index = False, sep = '\t' , compression = 'gzip')
	return 

def main():
	if len(sys.argv) != 4:
		usage()
	all_state_folder = sys.argv[1]
	helper.check_dir_exist(all_state_folder)
	num_max_GO_terms = helper.get_command_line_integer(sys.argv[2])
	output_folder = sys.argv[3]
	helper.make_dir(output_folder)
	print("Done getting command line arguments")
	combine_GO_output_all_states(all_state_folder, num_max_GO_terms, output_folder)
	print("Done!")
	return 

def usage():
	print("python process_GREAT_output_all_states.py")
	print("all_state_folder: where output of GREAT for all states are stored")
	print("num_max_GO_terms: number of maximum GO terms that we would like to collect")
	print("output_folder")
	exit(1)
main()