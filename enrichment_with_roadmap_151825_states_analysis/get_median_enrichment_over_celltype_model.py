import pandas as pd
import seaborn as sns
import numpy as np
import os
import sys
sys.path.append('/Users/vuthaiha/Desktop/window_hoff/source/characterize_full_stack_model_states')
import characterize_model_helper as cmh
import time
from glob import glob
import string
start_time = time.clock()
TWENTY_FIVE_STATE_ANNOTATION = ['1_TssA','2_PromU','3_PromD1','4_PromD2', \
'5_Tx5p','6_Tx','7_Tx3p','8_TxWk',\
'9_TxReg','10_TxEnh5p','11_TxEnh3p','12_TxEnhW',\
'13_EnhA1','14_EnhA2','15_EnhAF',\
'16_EnhW1','17_EnhW2','18_EnhAc','19_DNase',\
'20_ZNF/Rpts',\
'21_Het',\
'22_PromP',\
'23_PromBiv',\
'24_ReprPC',\
'25_Quies']
TWENTY_FIVE_STATE_COLOR = {1: 'Red', 2: 'Orange Red', 3: 'Orange Red', 4: 'Orange Red', 5: 'Green', 6: 'Green', 7: 'Green', 8: 'Lighter Green', 9: 'Electric Lime', 10: 'Electric Lime', 11: 'Electric Lime', 12: 'Electric Lime', 13: 'Orange', 14: 'Orange', 15: 'Orange', 16: 'Yellow', 17: 'Yellow', 18: 'Yellow', 19: 'Lemon', 20: 'Aquamarine', 21: 'Light Purple', 22: 'Pink', 23: 'Dark Purple', 24: 'Gray', 25: 'White'}

# ['red','red','red','red',\
# 'green','green','green','#228B22',\
# '#adff2f','#adff2f','#adff2f','#adff2f',\
# 'orange','orange','orange',\
# 'yellow','yellow','yellow','yellow',\
# '#40E0D0',\
# '#9370db',\
# '#FFC0CB',\
# 'purple',\
# 'silver',\
# 'white']

ANATOMY_COLOR_CODE = {'BLOOD': '#e34a33', 'ESC_DERIVED': '#fdbb84', 'ESC': '#fef0d9', 'BRAIN': '#ffffd4', 'LUNG': '#006837', 'SKIN': '#78c679', 'MUSCLE': '#c2e699', 'IPSC': '#7a0177', 'GI_STOMACH' : '#f768a1', 'HEART': '#fbb4b9', 'BREAST' : '#d7b5d8', 'GI_INTESTINE' : '#253494', 'GI_RECTUM' : "#2c7fb8", 'GI_COLON' : "#a1dab4", 'FAT' : "#54278f", 'LIVER' : '#9e9ac8', 'VASCULAR' : '#f2f0f7', 'PANCREAS': '#1c9099', 'STROMAL_CONNECTIVE' : '#bdc9e1', 'THYMUS' : "#f6eff7", 'PLACENTA': '#3182bd', 'GI_DUODENUM' : "#636363", 'CERVIX': '#cccccc', 'BONE': "#051C34", 'ADRENAL': '#55A4F4', 'KIDNEY': '#F4556C', 'MUSCLE_LEG' : '#9B7A7F', 'OVARY' : '#F2FA10', 'SPLEEN': '#11FA10', 'GI_ESOPHAGUS' : '#10F7FA', 'NaN': "#E5B8E8", 'NA' : '#E5B8E8'}
CELL_GROUP_COLOR_CODE = {'Neurosph': '#FFD924', 'HSC & B-cell': '#678C69', 'Mesench': '#B65C73', 'Brain': '#C5912B', 'Adipose': '#AF5B39', 'Thymus': '#DAB92E', 'Sm. Muscle': '#F182BC', 'IMR90': '#E41A1C', 'Myosat': '#E67326', 'iPSC': '#69608A', 'Muscle': '#C2655D', 'Digestive': '#C58DAA', 'ESC': '#924965', 'Epithelial': '#FF9D0C', 'Heart': '#D56F80', 'ENCODE2012': '#000000', 'ES-deriv': '#4178AE', 'Other': '#999999', 'Blood & T-cell': '#55A354', 'NA' : 'black'}

meta_df = cmh.get_metada(cmh.DEFAULT_METADATA_FN)  # get the default metada with the tissue type for each of the cell type
# column names: ["CT_NAME", "Epig_name", "GROUP", "TYPE", "ANATOMY"]
meta_df['color'] = (meta_df['GROUP']).map(CELL_GROUP_COLOR_CODE)
def get_celltype_color_map():
	CELLTYPE_COLOR_MAP = pd.Series(meta_df.color.values, index = meta_df.CT_NAME).to_dict() # convert two columns in to a dictionary: keys: cell type, values: color corresponding to that tissue
	CELLTYPE_COLOR_MAP['NA'] = '#E5B8E8'
	return CELLTYPE_COLOR_MAP

CELLTYPE_COLOR_MAP = get_celltype_color_map()

def get_rid_of_stupid_file_tail(context_name):
	if context_name.endswith('.bed.gz'):
		return(context_name[:-7])
	else:
		return(context_name)

def color_cell_type_names(val):
	if val == "":
		color = CELLTYPE_COLOR_MAP['NA']
	else:
		color = CELLTYPE_COLOR_MAP[val]
	return 'background-color: %s' % color

def color_tissue_names(val):
	if val == "":
		color = CELL_GROUP_COLOR_CODE['NA']
	else:
		color = CELL_GROUP_COLOR_CODE[val]
	return 'background-color: %s' % color

def get_one_enrichment_ct_model_df(fn, num_state_ct_model):
	df = pd.read_csv(fn, sep = '\t', header = 0)
	df = df.rename(columns = {'state (Emission order)' : 'state', 'Genome %': 'percent_in_genome'}) # rename some columns so that it is easier to write column names later
	state_colName_list = map(lambda x: 'state_' + str(x + 1) + ".bed.gz", range(num_state_ct_model))
	df = df[['state', 'percent_in_genome'] + state_colName_list]
	return df

def get_median_enrichment_one_genome_context(all_ct_df_list, num_state_ct_model, state_context_colName): 
	# state_context_colName: name of the column that we are looking at so that we know what column to look at for each data frame
	this_context_df = pd.DataFrame({'state' : (all_ct_df_list[0])['state']})
	# collect data from each of the cell type specific model --> columsn: cell type specific model, rows: 100 full stack states
	for df_index, df in enumerate(all_ct_df_list):
		this_context_df['ct' + str(df_index)] = df[state_context_colName]
	this_context_df['median_enrichment'] = (this_context_df.drop('state', axis = 1)).median(axis = 1)
	return np.array(this_context_df['median_enrichment'])

def get_max_enrichment_one_genome_context(all_ct_df_list, all_ct_name_list, num_state_ct_model, state_context_colName): 
	# state_context_colName: name of the column that we are looking at so that we know what column to look at for each data frame
	this_context_df = pd.DataFrame({'state' : (all_ct_df_list[0])['state']})
	# collect data from each of the cell type specific model --> columsn: cell type specific model, rows: 100 full stack states
	for df_index, df in enumerate(all_ct_df_list):
		ct_name = all_ct_name_list[df_index] # all_ct_fn_list and all_ct_df_list are responsive to each other, i.e. each element in each list corresponds to the same cell type. We tried zipping teh two list, but that gave a bug message. so this code here is not the most graceful.
		this_context_df[ct_name] = df[state_context_colName]
	this_context_df['max_enrichment'] = (this_context_df.drop('state', axis = 1)).max(axis = 1)
	this_context_df['context_max_enrichment'] = (this_context_df.drop(['state', 'max_enrichment'], axis = 1)).idxmax(axis = 1)
	return this_context_df[['max_enrichment', 'context_max_enrichment']]

def get_25_state_annot(state):
	# 'state_1' --> ''1_TssA''
	state_index = int(state.split('_')[1]) - 1 # zero-based
	return TWENTY_FIVE_STATE_ANNOTATION[state_index]

def get_TWENTY_FIVE_STATE_COLOR(state):
	# '1_TssA' --> 'red'
	state_index = int(state.split('_')[0]) # one-based
	color = TWENTY_FIVE_STATE_COLOR[state_index]
	return 'background-color: %s' % color

def paint_excel_max_enrichment(enrichment_df, max_epig_df, output_folder):
	# here enrichment_df is actually max_enrichment_df from get_all_ct_model_max_enrichment
	cm = sns.light_palette("red", as_cmap=True)
	(num_state, num_enr_cont) = (enrichment_df.shape[0] - 1, enrichment_df.shape[1] - 1)

	enr_cont_list = enrichment_df.columns[2:]
	enr_cont_list = map(get_rid_of_stupid_file_tail, enr_cont_list) # fix the name of the enrichment context. If it contains the tail '.bed.gz' then we get rid of it
	enrichment_df.columns = list(enrichment_df.columns[:2]) + list(enr_cont_list) # now change the column names of the enrichment contexts. If the contexts' names contain '.bed.gz' then get rid of it.
	max_epig_df.columns = map(get_rid_of_stupid_file_tail, max_epig_df.columns)
	no_state_df = enrichment_df[enrichment_df.columns[2:]] # only get the data of enrichment contexts for now, don't consider state and percent_in_genome

	enrichment_df['max_fold_context'] = no_state_df.apply(lambda x: x.idxmax(), axis = 1) # name of the context that are most enriched in this state
	enrichment_df['max_enrich'] = no_state_df.apply(lambda x: x.max(), axis = 1) # the value of the max fold enrichment in this state
	max_epig_df['max_fold_context'] = enrichment_df['max_fold_context']  #copy the data of the max_fold_context to the max_epig_df so that we can later get the information of the maximum enrichment epig for each state
	max_epig_df['max_enrich_ct'] = max_epig_df.apply(lambda x: x[x['max_fold_context']], axis = 1) # get the cell type that has the state that are most enriched in this state
	enrichment_df['max_enrich_ct'] = max_epig_df['max_enrich_ct']
	enrichment_df = pd.merge(enrichment_df, meta_df, how = 'left', left_on = 'max_enrich_ct', right_on = 'CT_NAME') # now we will have these additional columns: Epig_name, GROUP, TYPE, ANATOMY
	enrichment_df = enrichment_df[['state', 'percent_in_genome'] + list(enr_cont_list) + ['max_fold_context', 'max_enrich', 'max_enrich_ct', 'GROUP']] # select only the columns that we care about
	enrichment_df['max_fold_context'] = map(get_25_state_annot ,enrichment_df['max_fold_context']) # convert the state that is max enriched here with the right state annotation of the 25 state system. This is specific for the 25 state enrichment analysis
	colored_df = enrichment_df.style.background_gradient(subset = pd.IndexSlice[:(num_state - 1), [enr_cont_list[0]]], cmap = cm) # color the enrichment values into a red-white gradient
	for enr_cont in enr_cont_list[1:]:
	    colored_df = colored_df.background_gradient(subset = pd.IndexSlice[:(num_state - 1), [enr_cont]], cmap = cm)
	colored_df = colored_df.applymap(get_TWENTY_FIVE_STATE_COLOR, subset = pd.IndexSlice[:, ['max_fold_context']]) # color the 25-state that is most enriched in each state
	colored_df = colored_df.applymap(color_cell_type_names, subset = pd.IndexSlice[:, ['max_enrich_ct']]) # color the tissue type of the cell types whose 25-state is most enriched in each state
	colored_df = colored_df.applymap(color_tissue_names, subset = pd.IndexSlice[:, ['GROUP']]) # color the tissue type of the cell types whose 25-state is most enriched in each state

	# now color the max_epig_df
	colored_max_epig_df = max_epig_df.style.applymap(color_cell_type_names, subset = pd.IndexSlice[0:len(max_epig_df.index) - 2, enr_cont_list]) # paint each of the cell so that it corresponds to the tissue color of the cell type that has the max enrichment value. We skip coloring the last row, which is the Base row. In fact, we can ignore this row.
	print "Done coloring"
	# now save into two sheets
	output_fn = os.path.join(output_folder + "max_enrichment_all_ct.xlsx")
	writer = pd.ExcelWriter(output_fn, engine='xlsxwriter')
	colored_df.to_excel(writer, sheet_name = 'max_enrichment')
	colored_max_epig_df.to_excel(writer, sheet_name = 'max_epig_enrichment')
	writer.save()


def get_all_ct_model_max_enrichment(input_folder, output_folder, num_state_ct_model):
	all_ct_folders = glob(input_folder + "/E*/")
	all_ct_fn_list = map(lambda x: x + "/overlap_enrichment.txt", all_ct_folders)
	all_ct_name_list = map(lambda x: (x.split("/"))[-2], all_ct_folders)
	all_ct_df_list = map(lambda x: get_one_enrichment_ct_model_df(x, num_state_ct_model), all_ct_fn_list)
	# up until here, we have finished getting all the data from all the overlap enrichment files of all cell types into data frame --> all_ct_df_list [index] --> data frame of the cell type
	# create a data frame that will report the median enrichment over all the cell types
	max_enrichment_df = pd.DataFrame({'state': (all_ct_df_list[0])['state'], 'percent_in_genome' : (all_ct_df_list[0])['percent_in_genome']}) # this will report for each state in the full stack model, and for each state in the cell type specific model, what are the highest enrichment in each state-state combination. 
	max_enrichment_df = max_enrichment_df[['state', 'percent_in_genome']] # some how when we create the data frame it moves state to be the second column instead of the first one, so this line of code is there to fix that problem, not a big deal
	max_epig_df = pd.DataFrame({'state': (all_ct_df_list[0])['state'], 'percent_in_genome' : (all_ct_df_list[0])['percent_in_genome']})
	max_epig_df = max_epig_df[['state', 'percent_in_genome']]

	# now loop through each of the context that we care about
	enrichment_context_list = (all_ct_df_list[0]).columns[2:]
	for enr_context in enrichment_context_list:
		this_context_max_enrichment_df = get_max_enrichment_one_genome_context(all_ct_df_list, all_ct_name_list, num_state_ct_model, enr_context)
		max_enrichment_df[enr_context] = this_context_max_enrichment_df['max_enrichment']
		max_epig_df[enr_context] = this_context_max_enrichment_df['context_max_enrichment']
	print "Done getting all the necessary data after: " + str(time.clock() - start_time)
	# max_enrichment_save_fn = os.path.join(output_folder, 'max_enrichment_all_ct.txt')
	# max_enrichment_df.to_csv(max_enrichment_save_fn, sep = '\t', index = False, header = True)
	# max_epig_save_fn = os.path.join(output_folder, 'max_epig_name')
	# max_epig_df.to_csv(max_epig_save_fn + ".txt", sep = '\t', index = False, header = True)
	paint_excel_max_enrichment(max_enrichment_df, max_epig_df, output_folder)
	print "Done getting the max_enrichment_df after: " + str(time.clock() - start_time)

def get_all_ct_model_median_enrichment(input_folder, output_folder, num_state_ct_model):
	output_fn = os.path.join(output_folder, "median_enrichment_all_ct.txt")
	all_ct_folders = glob(input_folder + "/E*/")
	all_ct_fn_list = map(lambda x: x + "/overlap_enrichment.txt", all_ct_folders)
	all_ct_df_list = map(lambda x: get_one_enrichment_ct_model_df(x, num_state_ct_model), all_ct_fn_list)
	# up until here, we have finished getting all the data from all the overlap enrichment files of all cell types into data frame --> all_ct_df_list [index] --> data frame of the cell type
	# create a data frame that will report the median enrichment over all the cell types
	median_enrichment_df = pd.DataFrame({'state': (all_ct_df_list[0])['state'], 'percent_in_genome' : (all_ct_df_list[0])['percent_in_genome']})
	# now loop through each of the context that we care about
	enrichment_context_list = (all_ct_df_list[0]).columns[2:]
	for enr_context in enrichment_context_list:
		median_enrichment_df[enr_context] = get_median_enrichment_one_genome_context(all_ct_df_list, num_state_ct_model, enr_context)
	median_enrichment_df.to_csv(output_fn, sep = '\t', index = False, header = True)
	print "Done getting the median_enrichment_df after: " + str(time.clock() - start_time)

def main():
	if len(sys.argv) != 4:
		usage()
	input_folder = sys.argv[1]
	if not os.path.isdir(input_folder):
		print "input_folder is not a valid directory"
		usage()
	output_folder = sys.argv[2]
	cmh.make_dir(output_folder)
	try: 
		num_state_ct_model = int(sys.argv[3]) # number of state in the cell type specific model, 15 25 or 18
	except:
		print "num_state_ct_model is not valid"
		usage() 
	print "Done getting command line arguments after: " + str(time.clock() - start_time)
	#get_all_ct_model_median_enrichment(input_folder, output_folder, num_state_ct_model)
	get_all_ct_model_max_enrichment(input_folder, output_folder, num_state_ct_model)

def usage():
	print "python get_median_enrichment_over_celltype_model.py"
	print "input_folder: where thare are multiple subfolders, each containing enrichment data for different cell type specific model"
	print "output_fn: Where the median fold enrichment across all cell types for each state is reported"
	print "num_state_ct_model: number of states in the cell type specific models"
	exit(1)

main()