import pandas as pd 
import numpy as np 
import os
import sys
import helper
import glob
import matplotlib.pyplot as plt
import seaborn as sns
def read_one_GREAT_output_file(fn):
	df = pd.read_csv(fn, comment = "#", sep = '\t', index_col = None, header = None, skiprows = 1)
	header_line = "# Ontology,ID,Desc,BinomRank,BinomP,BinomBonfP,BinomFdrQ,RegionFoldEnrich,ExpRegions,ObsRegions,GenomeFrac,SetCov,HyperRank,HyperP,HyperBonfP,HyperFdrQ,GeneFoldEnrich,ExpGenes,ObsGenes,TotalGenes,GeneSetCov,TermCov,Regions,Genes"
	header_list = header_line.split(',')
	header_list[0] = header_list[0][2:]
	df.columns = header_list
	df = df[['Ontology', 'ID', 'Desc', 'BinomRank', 'BinomP', 'RegionFoldEnrich', 'ExpRegions', 'ObsRegions','GenomeFrac', 'SetCov','GeneFoldEnrich', 'ExpGenes', 'ObsGenes', 'TotalGenes','GeneSetCov', 'TermCov','Genes']]
	return df

def get_state_name_dict(state_annot_fn):
	state_annot_df = pd.read_csv(state_annot_fn, header = 0, index_col = None, sep = ',')
	state_annot_df = state_annot_df[['state', 'state_order_by_group', 'mneumonics', 'Group']]
	return dict(zip(state_annot_df.state, state_annot_df.mneumonics))

def draw_barplot_one_state(df, GO_type_list, rank_filter_list, save_fn, num_terms_to_plot, plot_title):
	'''
	Right now the step of separating rank_fileter_list and num_terms_to_plot feel redundant, but we keep it for now in case we want to query the data differently in the future. 
	'''
	go_df = df[df['Ontology'].isin(GO_type_list)]
	min_nonzero_p = np.min(go_df.BinomP[go_df.BinomP != 0]) # find this number first before we go into filter for only the top enriched terms
	go_df = go_df[go_df['BinomRank'].isin(rank_filter_list)]
	go_df['full_stack_state'] = plot_title
	plot_df = pd.DataFrame(columns = ['ID', 'Desc', 'BinomP'])
	num_terms_to_plot_per_go = int(num_terms_to_plot / len(GO_type_list))
	for go_type in GO_type_list:
		this_go_df = go_df[go_df['Ontology'] == go_type]
		this_go_df = this_go_df.nsmallest(num_terms_to_plot_per_go, 'BinomRank')[['ID', 'Desc', 'BinomP']]
		plot_df = plot_df.append(this_go_df)
	plot_df['BinomP'] = plot_df['BinomP'].replace(0, min_nonzero_p) # replace the p values of 0 with the smallest non zero values
	plot_df['negLog10_p'] = - np.log10(plot_df['BinomP'])
	plot_df = plot_df.sort_values(by = ['negLog10_p'], ascending = False)
	plot_df['GO_term'] = plot_df['ID'] + ' ' + plot_df['Desc']
	fig, ax = plt.subplots()
	sns.barplot(x='negLog10_p', y = 'GO_term', data = plot_df, ax = ax, color = '#338DFF').tick_params(axis = 'x', rotation = 45)
	ax.set_title(plot_title) # blue color
	ax.set_position([0.6, 0.15,0.36, 0.75])
	fig.set_size_inches(12.3, 5)
	fig.savefig(save_fn) 
	return go_df



def draw_barplot_all_states(all_state_folder, num_max_GO_terms, output_folder, num_terms_to_plot, state_annot_fn):
	state_name_dict = get_state_name_dict(state_annot_fn)
	rank_filter_list = list(map(lambda x: x+1, range(num_max_GO_terms)))	
	all_state_fn_list = glob.glob(all_state_folder + "/state_E*_GO.tsv.gz")
	# state_list = list(map(lambda x: x+1, range(100))) #[43, 44, 45, 46, 49, 98, 99]#[47, 48]#[50, 51, 52, 97]
	# all_state_fn_list = list(map(lambda x: os.path.join(all_state_folder, 'state_E' + str(x) + '_GO.tsv.gz'), state_list))
	all_state_list = list(map(lambda x: int(x.split('/')[-1].split('_')[1][1:]), all_state_fn_list))
	all_state_df_list = list(map(read_one_GREAT_output_file, all_state_fn_list))
	GO_type_list = ['GO Biological Process', 'GO Molecular Function'] #['GO Biological Process', 'GO Cellular Component', 'GO Molecular Function']
	processed_go_df_list = []
	for index, state in enumerate(all_state_list):
		state_df = all_state_df_list[index]
		plot_title = state_name_dict[state]
		save_fn = os.path.join(output_folder, 'S_'+ plot_title + '_barplot.pdf')
		top_go_df = draw_barplot_one_state(state_df, GO_type_list, rank_filter_list, save_fn, num_terms_to_plot, plot_title)
		processed_go_df_list.append(top_go_df)
		print('Done for state {}'.format(state))
	all_state_go_df = pd.concat(processed_go_df_list)
	save_fn = os.path.join(output_folder, 'top10_GO_terms_all_states.txt.gz')
	all_state_go_df.to_csv(save_fn, header = True, index = False, sep = '\t', compression = 'gzip')
	return 
	
def main():
	if len(sys.argv) != 6:
		usage()
	all_state_folder = sys.argv[1]
	helper.check_dir_exist(all_state_folder)
	num_max_GO_terms = helper.get_command_line_integer(sys.argv[2])
	output_folder = sys.argv[3]
	helper.make_dir(output_folder)
	num_terms_to_plot = helper.get_command_line_integer(sys.argv[4])
	state_annot_fn = sys.argv[5]
	helper.check_file_exist(state_annot_fn)
	print("Done getting command line arguments")
	draw_barplot_all_states(all_state_folder, num_max_GO_terms, output_folder, num_terms_to_plot, state_annot_fn)
	print("Done!")
	return 

def usage():
	print("python draw_barplot_per_state.py: for each state, draw a barplot of the top N GO terms")
	print("all_state_folder: where output of GREAT for all states are stored")
	print("num_max_GO_terms: number of maximum GO terms that we would like to collect per state")
	print("output_folder")
	print('num_terms_to_plot: number of GO terms to include into the plots')
	print('state_annot_fn')
	exit(1)
main()