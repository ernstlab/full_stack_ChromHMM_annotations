import pandas as pd
import numpy as np
import sys
import os
import time 
import characterize_model_helper as cmh

start_time = time.clock()
CHROM_MARK_COLOR_CODE = {'H2BK12ac': '#E5EAE7', 'H2AK5ac': '#E5EAE7', 'H2BK20ac': '#E5EAE7', 'H3K27ac': '#F7CB4D', 'H2BK5ac': '#E5EAE7', 'H4K12ac': '#E5EAE7', 'H4K5ac': '#E5EAE7', 'H3K4me3': '#F13712', 'H3K4me2': '#FFA500', 'H3K4me1': '#EDF732', 'H2A.Z': '#C3D59C', 'H3K9me1': '#F3AEAE', 'H3K23ac': '#E5EAE7', 'H4K91ac': '#E5EAE7', 'H3K4ac': '#E5EAE7', 'H3T11ph': '#F3AEAE', 'H3K23me2': '#F3AEAE', 'H2BK120ac': '#E5EAE7', 'H3K27me3': '#A6A6A4', 'H3K79me1': '#9DCEA8', 'H3K79me2': '#377A45', 'H3K9ac': '#F2A626', 'H3K18ac': '#E5EAE7', 'K3BK15ac': '#E5EAE7', 'H3K36me3': '#49AC5E', 'H4K20me1': '#6AD1C8', 'H3K56ac': '#E5EAE7', 'DNase': '#DBE680', 'H3K9me3': '#677BF6', 'H2AK9ac': '#E5EAE7', 'H3K14ac': '#E5EAE7', 'H4K8ac': '#E5EAE7', 'NA' : 'black', 'H2BK15ac' : '#E5EAE7'}
ANATOMY_COLOR_CODE = {'BLOOD': '#e34a33',\
'ESC_DERIVED': '#fdbb84', \
'ESC': '#fef0d9', \
'BRAIN': '#ffffd4', \
'LUNG': '#006837', \
'SKIN': '#78c679', \
'MUSCLE': '#c2e699', \
'IPSC': '#7a0177', \
'GI_STOMACH' : '#f768a1', \
'HEART': '#fbb4b9', \
'BREAST' : '#d7b5d8', \
'GI_INTESTINE' : '#253494', \
'GI_RECTUM' : "#2c7fb8", \
'GI_COLON' : "#a1dab4", \
'FAT' : "#54278f", \
'LIVER' : '#9e9ac8', \
'VASCULAR' : '#f2f0f7', \
'PANCREAS': '#1c9099', \
'STROMAL_CONNECTIVE' : '#bdc9e1', \
'THYMUS' : "#f6eff7", \
'PLACENTA': '#3182bd', \
'GI_DUODENUM' : "#636363", \
'CERVIX': '#cccccc', \
'BONE': "#051C34", \
'ADRENAL': '#55A4F4', \
'KIDNEY': '#F4556C', \
'MUSCLE_LEG' : '#9B7A7F', \
'OVARY' : '#F2FA10', \
'SPLEEN': '#11FA10', \
'GI_ESOPHAGUS' : '#10F7FA', \
'NaN': "#040404"}
CELL_GROUP_COLOR_CODE = {'Neurosph': '#FFD924', 'HSC & B-cell': '#678C69', 'Mesench': '#B65C73', 'Brain': '#C5912B', 'Adipose': '#AF5B39', 'Thymus': '#DAB92E', 'Sm. Muscle': '#F182BC', 'IMR90': '#E41A1C', 'Myosat': '#E67326', 'iPSC': '#69608A', 'Muscle': '#C2655D', 'Digestive': '#C58DAA', 'ESC': '#924965', 'Epithelial': '#FF9D0C', 'Heart': '#D56F80', 'ENCODE2012': '#000000', 'ES-deriv': '#4178AE', 'Other': '#999999', 'Blood & T-cell': '#55A354', 'NA' : 'black'}

def get_ranked_mark_name (row_data):
	sorted_rank = row_data.sort_values(ascending = False) # 1 --> n
	return pd.Series(sorted_rank.index) # get the index whcih is the list of experiments, ordered from least emitted to most emitted in each state. Convert to pd.Series so that we can concatenate into data frame later.

def get_ranked_exp_df (fn) : # fn should be emission_fn
	df = pd.read_csv(fn, header = 0, sep = '\t') # get emission df
	df = df.rename(columns = {'state (Emission order)' : 'state'}) # one column is state, others  are all experiment names
	df.index = df['state'] # state column becomes the index of the dataframe
	df = df[df.columns[1:]] # get rid of state column so that we can rank the emission probabilities of each state
	rank_df = df.rank(axis = 1) # index: state, columns: experiments ordered just like in df. Cells: rank of each experiment within each state, 1: least emission --> n: highest emission
	rank_df = rank_df.apply(get_ranked_mark_name, axis = 1) # index; state, columns: rank 1 --> n most emitted experiment to least emitted experiments
	rank_df.columns = map(lambda x: 'r_' + str(x + 1), range(len(rank_df.columns)))
	return rank_df

def color_anatomy_names(val):
	if val == "":
		color = ANATOMY_COLOR_CODE['NaN']
	else:
		color = ANATOMY_COLOR_CODE[val]
	return 'background-color: %s' % color

def color_cell_group_names(val):
	if val == "":
		color = CELL_GROUP_COLOR_CODE["NA"]
	else:
		color = CELL_GROUP_COLOR_CODE[val]
	return 'background-color: %s' % color

def color_mark_names(val):
	if val == "":
		color = CHROM_MARK_COLOR_CODE['NaN']
	else:
		color = CHROM_MARK_COLOR_CODE[val]
	return 'background-color: %s' % color

def get_painted_excel_ranked_exp(rank_df, output_fn, num_top_marks, celltype_ana_dict, celltype_cell_group_dict):
	ct_df = rank_df.applymap(lambda x: x.split('-')[0]) # E118
	ct_df = ct_df.applymap(lambda x: celltype_cell_group_dict[x]) # convert the cell types to the tissue types
	ct_df.reset_index(inplace = True)
	chrom_mark_df = rank_df.applymap(lambda x: x.split('-')[1]) # H3K9me3
	chrom_mark_df.reset_index(inplace = True)
	columns_to_paint = map(lambda x: 'r_' + str(x + 1), range(num_top_marks))
	ct_df = ct_df[['state'] + columns_to_paint]
	chrom_mark_df = chrom_mark_df[['state'] + columns_to_paint]
	colored_ct_df = ct_df.style.applymap(color_cell_group_names, subset = columns_to_paint)
	colored_chrom_mark_df = chrom_mark_df.style.applymap(color_mark_names, subset = columns_to_paint)
	# save file
	writer = pd.ExcelWriter(output_fn, engine='xlsxwriter')
	colored_ct_df.to_excel(writer, sheet_name = 'cell_group')
	colored_chrom_mark_df.to_excel(writer, sheet_name = 'chrom_mark')
	writer.save()
	print "Done saving data into " + output_fn + " after: " + str(time.clock() - start_time)

def main():
	if len(sys.argv) != 5:
		usage()
	emission_fn = sys.argv[1]
	output_fn = sys.argv[2]
	try: 
		num_top_marks = int(sys.argv[3])
	except:
		print 'num_top_marks is NOT VALID'
		usage()
	metadata_fn = sys.argv[4]
	if not os.path.isfile(metadata_fn):
		print "****************** ATTENTION ****************"
		print "metadata_fn: " + metadata_fn + " does not exist. We will use the default metadata_fn: " + cmh.DEFAULT_METADATA_FN
		print "****************** END OF ATTENTION ****************"
		metadata_fn = cmh.DEFAULT_METADATA_FN
	meta_df = cmh.get_metada(metadata_fn) 
	print "Got meta_df after: " + str(time.clock() - start_time)
	celltype_ana_dict = dict(zip(meta_df.CT_NAME, meta_df.ANATOMY)) # keys: cell types, values: anatomy of the cell types
	celltype_cell_group_dict = dict(zip(meta_df.CT_NAME, meta_df.GROUP))
	print "Done getting command line argument after: " + str(time.clock() - start_time)
	rank_df = get_ranked_exp_df(emission_fn) # index: states, columns: rank of experiments r_1: highest emitted, r_n: lowest emitted --> cells: names of experiments. Example: E118-H3K9me3
	print "Done getting ranked data after: " + str(time.clock() - start_time)
	get_painted_excel_ranked_exp(rank_df, output_fn, num_top_marks, celltype_ana_dict, celltype_cell_group_dict)



def usage():
	print "python rank_chrom_mark_celltype_emission.py"
	print "emission_fn"
	print "output_fn: execl file where the output data of ranked cell type and chrom marks are stored"
	print "num_top_marks: number of top marks that we want to report"
	print "metadata_fn: file name of meta data of different reference epigenomes from ROADMAP, optional"
	exit(1)

main()
# fn = '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/emissions_100.txt'
# df = pd.read_csv(fn, header = 0, sep = '\t')
# df = df.rename(columns = {'state (Emission order)' : 'state'})
# df.index = df['state']
# df = df[df.columns[1:]]
# rank_df.max(axis = 1)
# mark_names = df.columns
# chrom_mark = map(lambda x: x.split('-')[1], mark_names)
# ct_list = map(lambda x: x.split('-')[0], mark_names)

# def get_ranked_mark_name (row_data):
# 	sorted_rank = row_data.sort_values(ascending = False)
# 	return pd.Series(sorted_rank.index)

# res_df = rank_df.apply(get_ranked_mark_name, axis = 1)
# res_df.columns = map(lambda x: 'r_' + str(x + 1), range(len(res_df.columns)))
# ct_df = res_df.applymap(lambda x: x.split('-')[0])

# chrom_mark_df = res_df.applymap(lambda x: x.split('-')[1])
# chrom_mark_df.reset_index(inplace = True)