import pandas as pd 
import numpy as np 
import os
import sys
import helper
import glob
from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator
import matplotlib.pyplot as plt
def read_one_GREAT_output_file(fn):
	df = pd.read_csv(fn, comment = "#", sep = '\t', index_col = None, header = None, skiprows = 1)
	header_line = "# Ontology,ID,Desc,BinomRank,BinomP,BinomBonfP,BinomFdrQ,RegionFoldEnrich,ExpRegions,ObsRegions,GenomeFrac,SetCov,HyperRank,HyperP,HyperBonfP,HyperFdrQ,GeneFoldEnrich,ExpGenes,ObsGenes,TotalGenes,GeneSetCov,TermCov,Regions,Genes"
	header_list = header_line.split(',')
	header_list[0] = header_list[0][2:]
	df.columns = header_list
	return df

def draw_wordCloud_one_state(df, GO_type_list, rank_filter_list, save_fn):
	go_df = df[df['Ontology'].isin(GO_type_list)]
	go_df = go_df[go_df['BinomRank'].isin(rank_filter_list)]
	id_list_str = '|'.join(go_df['ID'])
	desc_list_str = ' '.join(go_df['Desc'])
	fdr_p_list_str = '|'.join(list(map(lambda x: str(x), go_df['BinomFdrQ'])))
	enrich_list_str = '|'.join(list(map(lambda x: str(x), go_df['RegionFoldEnrich'])))
	wc = WordCloud().generate(desc_list_str)
	wc.to_file(save_fn) # save the wordcloud to file
	return

def draw_wordCloud_all_states(all_state_folder, num_max_GO_terms, output_folder):
	rank_filter_list = list(map(lambda x: x+1, range(num_max_GO_terms)))	
	all_state_fn_list = glob.glob(all_state_folder + "/state_E*_GO.tsv.gz")
	all_state_list = list(map(lambda x: int(x.split('/')[-1].split('_')[1][1:]), all_state_fn_list))
	all_state_df_list = list(map(read_one_GREAT_output_file, all_state_fn_list))
	GO_type_list = ['GO Molecular Function'] #['GO Biological Process', 'GO Cellular Component', 'GO Molecular Function']
	for index, state in enumerate(all_state_list):
		state_df = all_state_df_list[index]
		save_fn = os.path.join(output_folder, 'S'+str(state) + '_wordcloud.png')
		draw_wordCloud_one_state(state_df, GO_type_list, rank_filter_list, save_fn)
		print('Done for state {}'.format(state))
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
	draw_wordCloud_all_states(all_state_folder, num_max_GO_terms, output_folder)
	print("Done!")
	return 

def usage():
	print("python draw_wordCloud_per_state.py: for each state, draw a wordcloud of the top N GO terms")
	print("all_state_folder: where output of GREAT for all states are stored")
	print("num_max_GO_terms: number of maximum GO terms that we would like to collect per state")
	print("output_folder")
	exit(1)
main()