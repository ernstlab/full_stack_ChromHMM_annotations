import pandas as pd 
import numpy as np 
import os
import sys
import chromHMM_utilities_common_functions_helper as helper
def get_mark_group_dictionary(emission_df):
	result = {} # keys: mark, value: group that the mark belongs to
	CLASS1_ACETYLATIONS = ['H2AK5ac', 'H2AK9ac', 'H2BK120ac', 'H2BK12ac', 'H2BK15ac', 'H2BK20ac', 'H2BK5ac', 'H3K14ac', 'H3K18ac', 'H3K23ac', 'H3K4ac', 'H3K56ac', 'H4K12ac', 'H4K5ac', 'H4K8ac', 'H4K91ac']
	OTHERS_MARKS = ['H3K23me2', 'H3K9me1', 'H3T11ph']
	for mark in np.unique(emission_df['chrom_mark']):
		if mark in CLASS1_ACETYLATIONS:
			result[mark] = 'class1_acetyl'
		elif mark in OTHERS_MARKS:
			result[mark] = 'other'
		else: 
			result[mark] = mark
	return result

def get_full_emission_df():
	emission_fn = '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/emissions_100.txt'
	emission_df = pd.read_csv(emission_fn, header = 0, sep = '\t') # columns: marks, rows: states
	emission_df = emission_df.rename(columns = {'state (Emission order)' : 'state'})
	emission_df = emission_df.transpose() # columns: states, rows: marks
	emission_df = emission_df.drop(index = 'state') # drop the first row cuz we will make column names soon
	num_state = emission_df.shape[1] # currently the number of columns is the number of state
	emission_df = emission_df.reset_index() # current index is the mark names, reset to make it a separate column
	emission_df.columns = ['mark'] + map(lambda x: str(x+1), range(num_state))
	emission_df['ct'] = (emission_df['mark']).apply(lambda x: x.split('-')[0])
	emission_df['chrom_mark'] = (emission_df['mark']).apply(lambda x: x.split('-')[1])
	mark_group_dict = get_mark_group_dictionary(emission_df)
	# get the annotation of marks
	mark_annot_fn = '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/ROADMAP_metadata_july2013.csv'
	annot_df = pd.read_csv(mark_annot_fn, header = 0, sep = ',')
	annot_df = annot_df.rename(columns = {'Epigenome ID (EID)' : 'ct'})
	annot_df = annot_df[['ct', 'GROUP', 'ANATOMY']] # only get the ones needed to join with emission_df
	emission_df = emission_df.merge(annot_df, left_on = 'ct', right_on = 'ct')
	emission_df['mark_group'] = emission_df['chrom_mark'].apply(lambda x: mark_group_dict[x])
	return emission_df, num_state

def get_emission_by_mark_df (emission_df, num_state):
	mark_emission_df = emission_df[['mark_group'] + map(lambda x: str(x+1), range(num_state))] # get only the states adn the chrom mark to make this emission_df
	mark_emission_df = mark_emission_df.groupby('mark_group').mean() # --> columns: state, rows: groups of marks, values: average emission for experiments associated with each of the group of marks
	mark_emission_df = mark_emission_df.reset_index()
	mark_emission_df = mark_emission_df.transpose() # rows: states, columns: mark_group
	mark_emission_df = mark_emission_df.reset_index() # reset index after the transpose job
	mark_emission_df.columns = mark_emission_df.loc[0] # column names is now ['makr_group', 'DNase', 'H2A.Z', etc.]
	mark_emission_df = mark_emission_df.drop(index = [0]) 
	mark_emission_df.columns = ['state'] + list(mark_emission_df.columns[1:]) 
	save_fn = '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/emission_results/good_plots/emissions_by_mark_group.txt'
	mark_emission_df.to_csv(save_fn, header = True, index = False ,sep = '\t')
	return mark_emission_df

def get_emission_by_group_df (emission_df, num_state):
	group_emission_df = emission_df[['GROUP'] + map(lambda x: str(x+1), range(num_state))]
	group_emission_df = group_emission_df.groupby('GROUP').mean() # --> columns: state, rows: groups of cell types, values: average emission for experiments associated with each of the group of cell types
	group_emission_df = group_emission_df.reset_index()
	group_emission_df = group_emission_df.transpose()
# rows: states, columns: cell type grou
	group_emission_df = group_emission_df.reset_index() # reset index after the transpose job
	group_emission_df.columns = group_emission_df.loc[0] # column names is now ['makr_group', 'DNase', 'H2A.Z', etc.]
	group_emission_df = group_emission_df.drop(index = [0]) 
	group_emission_df.columns = ['state'] + list(group_emission_df.columns[1:]) 
	save_fn = '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/emission_results/good_plots/emissions_by_cell_type_group.txt'
	group_emission_df.to_csv(save_fn, header = True, index = False ,sep = '\t')
	return group_emission_df

def get_emission_by_anatomy_df (emission_df, num_state):
	ana_emission_df = emission_df[['ANATOMY'] + map(lambda x: str(x+1), range(num_state))]
	ana_emission_df = ana_emission_df.groupby('ANATOMY').mean() # --> columns: state, rows: groups of cell types, values: average emission for experiments associated with each of the group of cell types
	ana_emission_df = ana_emission_df.reset_index()
	ana_emission_df = ana_emission_df.transpose()
# rows: states, columns: cell type grou
	ana_emission_df = ana_emission_df.reset_index() # reset index after the transpose job
	ana_emission_df.columns = ana_emission_df.loc[0] # column names is now ['makr_group', 'DNase', 'H2A.Z', etc.]
	ana_emission_df = ana_emission_df.drop(index = [0]) 
	ana_emission_df.columns = ['state'] + list(ana_emission_df.columns[1:]) 
	save_fn = '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/emission_results/good_plots/emissions_by_ana_group.txt'
	ana_emission_df.to_csv(save_fn, header = True, index = False, sep = '\t')
	return ana_emission_df

emission_df, num_state = get_full_emission_df()
get_emission_by_mark_df (emission_df, num_state)
get_emission_by_group_df(emission_df, num_state)
get_emission_by_anatomy_df(emission_df, num_state)