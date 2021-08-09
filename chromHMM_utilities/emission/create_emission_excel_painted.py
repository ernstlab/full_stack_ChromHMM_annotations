import pandas as pd
import seaborn as sns
import numpy as np 
import os
import sys
sys.path.append('/u/home/h/havu73/project-ernst/source/characterize_full_stack_model_states')
import characterize_model_helper as cmh
CHROM_MARK_COLOR_CODE = {'H2BK12ac': '#E5EAE7', 'H2AK5ac': '#E5EAE7', 'H2BK20ac': '#E5EAE7', 'H3K27ac': '#F7CB4D', 'H2BK5ac': '#E5EAE7', 'H4K12ac': '#E5EAE7', 'H4K5ac': '#E5EAE7', 'H3K4me3': '#F13712', 'H3K4me2': '#FFA500', 'H3K4me1': '#EDF732', 'H2A.Z': '#C3D59C', 'H3K9me1': '#F3AEAE', 'H3K23ac': '#E5EAE7', 'H4K91ac': '#E5EAE7', 'H3K4ac': '#E5EAE7', 'H3T11ph': '#F3AEAE', 'H3K23me2': '#F3AEAE', 'H2BK120ac': '#E5EAE7', 'H3K27me3': '#A6A6A4', 'H3K79me1': '#9DCEA8', 'H3K79me2': '#377A45', 'H3K9ac': '#F2A626', 'H3K18ac': '#E5EAE7', 'K3BK15ac': '#E5EAE7', 'H3K36me3': '#49AC5E', 'H4K20me1': '#6AD1C8', 'H3K56ac': '#E5EAE7', 'DNase': '#DBE680', 'H3K9me3': '#677BF6', 'H2AK9ac': '#E5EAE7', 'H3K14ac': '#E5EAE7', 'H4K8ac': '#E5EAE7', 'NA' : 'black', 'H2BK15ac' : '#E5EAE7'}
ANATOMY_COLOR_CODE = {'BLOOD': '#e34a33', 'ESC_DERIVED': '#fdbb84', 'ESC': '#fef0d9', 'BRAIN': '#ffffd4', 'LUNG': '#006837', 'SKIN': '#78c679', 'MUSCLE': '#c2e699', 'IPSC': '#7a0177', 'GI_STOMACH' : '#f768a1', 'HEART': '#fbb4b9', 'BREAST' : '#d7b5d8', 'GI_INTESTINE' : '#253494', 'GI_RECTUM' : "#2c7fb8", 'GI_COLON' : "#a1dab4", 'FAT' : "#54278f", 'LIVER' : '#9e9ac8', 'VASCULAR' : '#f2f0f7', 'PANCREAS': '#1c9099', 'STROMAL_CONNECTIVE' : '#bdc9e1', 'THYMUS' : "#f6eff7", 'PLACENTA': '#3182bd', 'GI_DUODENUM' : "#636363", 'CERVIX': '#cccccc', 'BONE': "#051C34", 'ADRENAL': '#55A4F4', 'KIDNEY': '#F4556C', 'MUSCLE_LEG' : '#9B7A7F', 'OVARY' : '#F2FA10', 'SPLEEN': '#11FA10', 'GI_ESOPHAGUS' : '#10F7FA', 'NaN': "#040404"}
CELL_GROUP_COLOR_CODE = {'Neurosph': '#FFD924', 'HSC & B-cell': '#678C69', 'Mesench': '#B65C73', 'Brain': '#C5912B', 'Adipose': '#AF5B39', 'Thymus': '#DAB92E', 'Sm. Muscle': '#F182BC', 'IMR90': '#E41A1C', 'Myosat': '#E67326', 'iPSC': '#69608A', 'Muscle': '#C2655D', 'Digestive': '#C58DAA', 'ESC': '#924965', 'Epithelial': '#FF9D0C', 'Heart': '#D56F80', 'ENCODE2012': '#000000', 'ES-deriv': '#4178AE', 'Other': '#999999', 'Blood & T-cell': '#55A354', 'NA' : 'black'}
MARK_NAME_BASED = 0
GROUP_BASED = 1
NUM_STATE = 100
def color_mark_names(val):
	if val == "":
	    color = CHROM_MARK_COLOR_CODE['NA']
	else:
	    color = CHROM_MARK_COLOR_CODE[val]
	return 'background-color: %s' % color

def color_anatomy_names(val):
	if val == "":
	    color = ANATOMY_COLOR_CODE['NaN']
	else:
	    color = ANATOMY_COLOR_CODE[val]
	return 'background-color: %s' % color

def color_cell_group_names(val):
	if val == "":
		color = CELL_GROUP_COLOR_CODE['NA']
	else:
		color = CELL_GROUP_COLOR_CODE[val]
	return 'background-color: %s' % color

def color_one_emission_df(emission_df, df_type_index): # df_type_index is either 0 or 1, corresponding to mark_name_based or group_based
	COLOR_FUNCTION_LIST = [color_mark_names, color_cell_group_names]
	COLNAME_TO_COLOR_LIST = ['chrom_mark', 'GROUP']
	colored_df = emission_df.style.applymap(COLOR_FUNCTION_LIST[df_type_index], subset = pd.IndexSlice[:, [COLNAME_TO_COLOR_LIST[df_type_index]]])
	cm = sns.light_palette("red", as_cmap=True)
	state_colnames_list = map(lambda x: 'S' + str(x + 1), range(NUM_STATE))
	colored_df = colored_df.background_gradient(subset = pd.IndexSlice[:,state_colnames_list], cmap = cm)
	return colored_df


def get_emission_df_marks_ordered_for_plotting(emission_fn): # order the marks so that experiments associated with the same chromatin marks are put next to each other
	emission_df, NUM_STATE, num_marks = cmh.get_emission_df_right(emission_fn) # rows: marks, columns: state and marks' parameter
	# merge with metadata so that we know more information about the experiments (cell group, tissue types, etc.)
	meta_df = cmh.get_metada(cmh.DEFAULT_METADATA_FN)
	emission_df = emission_df.merge(meta_df, left_on = "ct", right_on = "CT_NAME")
	uniq_chrom_mark = ['H3K9me3', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H4K20me1', 'H3K79me2', 'H3K79me1', 'H3K36me3', 'DNase', 'H3K27ac', 'H3K27me3', 'H2A.Z', 'H3K9ac', 'H2AK9ac', 'H4K12ac', 'H2BK20ac', 'H3K56ac', 'H4K5ac', 'K3BK15ac', 'H2BK12ac', 'H3K14ac', 'H4K91ac', 'H2AK5ac', 'H2BK120ac', 'H2BK5ac', 'H3K18ac', 'H3K23ac', 'H3K4ac', 'H2BK15ac', 'H4K8ac', 'H3T11ph', 'H3K23me2', 'H3K9me1'] # put chromatin marks in an order that puts all the acetylation marks together --> we are putting the rarely studied acetylation marks together, and some other less well studied marks together #sort(unique(emission_df$chrom_mark)) # list of unique chromatin marks
	plot_df = pd.DataFrame(columns = emission_df.columns) # create an empty df that have the same order of columns as emission df
	for mark in uniq_chrom_mark: # we will include marks that are associated with the marks listed in order in uniq_chrom_mark one at a time
		this_chrom_mark_df = emission_df[emission_df['chrom_mark'] == mark]
		plot_df = plot_df.append(this_chrom_mark_df)
	columns_to_collect = ['mark_name', 'chrom_mark'] + map(lambda x: x + 1, range(NUM_STATE))
	plot_df = plot_df[columns_to_collect]  # get only columns of interests
	plot_df.columns = ['mark_name', 'chrom_mark'] + map(lambda x: 'S' + str(x + 1), range(NUM_STATE))
	#print(plot_df.head())
	print "Inside get_emission_df_marks_ordered_for_plotting"
	plot_df = color_one_emission_df(plot_df, MARK_NAME_BASED)
	return plot_df

def get_emission_df_cell_groups_ordered_for_plotting(emission_fn):
	emission_df, NUM_STATE, num_marks = cmh.get_emission_df_right(emission_fn) # rows: marks, columns: state and marks' parameter 
	meta_df = cmh.get_metada(cmh.DEFAULT_METADATA_FN) # merge with metadata so that we know more information about the experiments (cell group, tissue types, etc.)
	emission_df = emission_df.merge(meta_df, left_on = "ct", right_on = "CT_NAME")
	cell_group_count = emission_df['GROUP'].value_counts() # get uniq cell groups, placed in descending order #  a Pandas series with index being the cell group names and values being the count of the indices. All stored in descending order
	plot_df = pd.DataFrame(columns = emission_df.columns) # create an empty df that have the same order of columns as emission df
	for group in cell_group_count.index: 
		this_group_df = emission_df[emission_df['GROUP'] == group]
		plot_df = plot_df.append(this_group_df)
	columns_to_collect = ['mark_name', 'GROUP'] + map(lambda x: x + 1, range(NUM_STATE))
	plot_df = plot_df[columns_to_collect]  # get only columns of interests
	plot_df.columns = ['mark_name', 'GROUP'] + map(lambda x: 'S' + str(x + 1), range(NUM_STATE))
	plot_df = color_one_emission_df(plot_df, df_type_index = GROUP_BASED)
	return plot_df

def get_emission_df_excel(emission_fn, save_fn):
	mark_ordered_df = get_emission_df_marks_ordered_for_plotting(emission_fn)
	group_ordered_df = get_emission_df_cell_groups_ordered_for_plotting(emission_fn)
	writer = pd.ExcelWriter(save_fn, engine='xlsxwriter')
	mark_ordered_df.to_excel(writer, sheet_name = 'mark_ordered')
	group_ordered_df.to_excel(writer, sheet_name = 'group_ordered')
	writer.save()


emission_fn = "/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/emissions_100.txt"
save_fn = "/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/emissions_100.xlsx"

get_emission_df_excel(emission_fn, save_fn)