import pandas as pd
import seaborn as sns
import numpy as np
import os
import sys
sys.path.append('/Users/vuthaiha/Desktop/window_hoff/source/characterize_full_stack_model_states')
import characterize_model_helper as cmh
import time
import glob
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
# orange red : #ff4500
# lighter orange red: #ff7d4d
# medium purple: #9370DB
# lemon: #fff44f
# darker lemon: #ffc84f
# electric lime: #CEFA05
# light purple: #b19cd9
# lighter green: #90ee90

TWENTY_FIVE_STATE_COLOR = {1: 'Red', 2: '#ff7d4d', 3: '#ff7d4d', 4: '#ff7d4d', 5: 'Green', 6: 'Green', 7: 'Green', 8: '#90ee90', 9: '#CEFA05', 10: '#CEFA05', 11: '#CEFA05', 12: '#CEFA05', 13: 'Orange', 14: 'Orange', 15: 'Orange', 16: 'Yellow', 17: 'Yellow', 18: 'Yellow', 19: '#ffc84f', 20: '#7fffd4', 21: '#b19cd9', 22: 'Pink', 23: '#9370DB', 24: 'Gray', 25: 'White'}

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
CELLTYPE_CELL_GROUP_MAP = pd.Series(meta_df.GROUP.values, index = meta_df.CT_NAME).to_dict() # keys: cell type name, values: cell group
print type(CELLTYPE_CELL_GROUP_MAP)
def get_rid_of_stupid_file_tail(context_name):
	if context_name.endswith('.bed.gz'):
		return(context_name[:-7])
	else:
		return(context_name)

def color_cell_group_names(val):
	if val == "":
		color = CELL_GROUP_COLOR_CODE['NA']
	else:
		color = CELL_GROUP_COLOR_CODE[val]
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
	df = df[['state', 'percent_in_genome'] + state_colName_list] # get the data frame to display columns in the expected order
	df.columns = ['state', 'percent_in_genome'] + TWENTY_FIVE_STATE_ANNOTATION # rename the columns so that instead of 'state1.bed.gz' we have '1_TssA'
	df['max_enrichment'] = (df.drop(['state', 'percent_in_genome'], axis = 1)).max(axis = 1) # find the maximum enrichment values in this row 
	df['max_enrichment_context'] = (df.drop(['state', 'percent_in_genome'], axis = 1)).idxmax(axis = 1)
	(nrow, ncol) = df.shape
	df = df.drop(nrow - 1) # drop the last row, which is the 'Base' row with percentage that each enrichment context occupies the genome
	return df

def get_25_state_annot(state):
	# 'state_1' --> ''1_TssA''
	state_index = int(state.split('_')[1]) - 1 # zero-based
	return TWENTY_FIVE_STATE_ANNOTATION[state_index]

def get_TWENTY_FIVE_STATE_COLOR(state):
	# '1_TssA' --> 'red'
	state_index = int(state.split('_')[0]) # one-based
	color = TWENTY_FIVE_STATE_COLOR[state_index]
	return 'background-color: %s' % color

def color_rank_ct_df(rank_ct_df):
	color_colname_list = rank_ct_df.columns[1:]
	rank_cell_group_df = (rank_ct_df[color_colname_list]).applymap(lambda x: CELLTYPE_CELL_GROUP_MAP[x]) # change the cell type data into the respective cell 
	rank_cell_group_df['state'] = rank_ct_df['state'] 
	rank_cell_group_df = rank_cell_group_df[['state'] + list(color_colname_list)] # without the list casting, the colnames would be staterank1, staterank2, etc.
 	colored_df = rank_cell_group_df.style.applymap(color_cell_group_names, subset = pd.IndexSlice[:, color_colname_list]) # only color the ranked cell types
	return colored_df

def color_rank_25_state_df(rank_25_state_df):
	color_colname_list = rank_25_state_df.columns[1:]
	colored_df = rank_25_state_df.style.applymap(get_TWENTY_FIVE_STATE_COLOR, subset = pd.IndexSlice[:, color_colname_list]) # only color the ranked 25-state
	return colored_df


def get_max_enrichment_25_state_df_ordered_cell_group(max_enrich_25_state_all_ct_df, max_enrich_value_all_ct_df):
	"""
	max_enrich_25_state_all_ct_df: rows: full stack states, columns: state of the 25-state system that is most enriched with the full-stack states WITHIN a cell type. Cell type are the column names  
	max_enrich_value_all_ct_df: Same data order as max_enrich_25_state_all_ct_df, except this dataframe stores the enrihcment values themselves
	--> a plot where cell types are juxtaposed based on the cell groups they belong to. The output dataframe will be colored properly
	"""
	cell_group_df = pd.DataFrame(columns = ['cell_type']) # rows: the cell types that are used in this analysis (colnames of max_enrich_25_state_all_ct_df), columns : 
	cell_group_df['cell_type'] = max_enrich_25_state_all_ct_df.columns[1:] # we skip the first column because that's 'state'
	cell_group_df = pd.merge(cell_group_df, meta_df, how = 'left', left_on = 'cell_type', right_on = 'CT_NAME') # merge the cell types so that we can get the information that we want
	# Now let's get the count of the number of cell types that are of each cell groups, and then get the unique cell groups, ordered by descending counts
	cell_group_count = cell_group_df['GROUP'].value_counts() # pandas Series: index: cellgroups, values: counts
	unique_cell_groups = cell_group_count.index
	# list of cell types, arranged such that those of the same cell_groups are juxtaposed
	rearranged_cell_types = []
	rearranged_cell_groups = []
	for cg in unique_cell_groups:
		this_cg_df = cell_group_df[cell_group_df['GROUP'] == cg] # filter out rows with this cellgroup
		rearranged_cell_types += list(this_cg_df['cell_type'])
		rearranged_cell_groups += list(this_cg_df['GROUP'])
	# now we will rearrange the columns of max_enrich_25_state_all_ct_df and max_enrich_value_all_ct_df based on the order that we just got from rearranged_cell_types
	max_enrich_25_state_all_ct_df = max_enrich_25_state_all_ct_df[['state'] + rearranged_cell_types]
	max_enrich_value_all_ct_df = max_enrich_value_all_ct_df[['state'] + rearranged_cell_types]
	cell_group_row_index = max_enrich_25_state_all_ct_df.shape[0]
	max_enrich_25_state_all_ct_df.loc[cell_group_row_index] = [''] + rearranged_cell_groups # add one more row that specific the cell_group of each of the cell types
	max_enrich_value_all_ct_df.loc[cell_group_row_index] = [''] + rearranged_cell_groups 
	colored_25_state_df = max_enrich_25_state_all_ct_df.style.applymap(color_cell_group_names, subset = pd.IndexSlice[cell_group_row_index, :]) # color the row that contain the cell_group of cell types
	colored_25_state_df = colored_25_state_df.applymap(get_TWENTY_FIVE_STATE_COLOR, subset = pd.IndexSlice[:(cell_group_row_index - 1) , colored_25_state_df.columns[1:]]) # color the states of the 25-system that is most enriched in each of the full-stack states
	return colored_25_state_df


def get_rank_25_state_df(rank_df, max_enrich_25_state_all_ct_df, max_enrich_value_all_ct_df, output_fn):
	(num_state, num_ct) = rank_df.shape # nrow, ncol
	rank_25_state_df = pd.DataFrame(columns = ['state'] + map(lambda x: 'rank' + str(x + 1),  range(num_ct))) # create an empty df: state, rank_1, rank_2, ... --> rank 25-state types
	rank_ct_df = pd.DataFrame(columns = ['state'] + map(lambda x: 'rank' + str(x + 1),  range(num_ct))) # create an empty df: state, rank_1, rank_2, ... --> rank cell types by the corresponding state enrichments
	rank_enrichment_value_df = pd.DataFrame(columns = ['state'] + map(lambda x: 'rank' + str(x + 1),  range(num_ct)))
	ct_list = rank_df.columns # list of ct as ordered by the rank_df, which is useful for the creation of rank_ct_df
	rank_25_state_df['state'] = map(lambda x: x+1, range(num_state)) # full-stack states
	rank_ct_df['state'] = map(lambda x: x + 1, range(num_state))
	rank_enrichment_value_df['state'] = map(lambda x: x+1, range(num_state))
	for row_index, row in rank_df.iterrows(): # loop through each full-stack state
		full_stack_state = row_index + 1
		ranked_ct_this_row = (row.sort_values()).index
		rank_ct_df.loc[row_index] = [str(full_stack_state)] + list(ranked_ct_this_row)
		
		# get this full-stack 25-state data 
		this_row_25_state_data = max_enrich_25_state_all_ct_df.iloc[row_index, 1:] # skip the first column because that's the state column, get the data for this full-stack state
		ranked_25_state_this_row = this_row_25_state_data[ranked_ct_this_row] # reorder the max_enriched 25_state based on the level of enrichments
		rank_25_state_df.loc[row_index] = [str(full_stack_state)] + list(ranked_25_state_this_row)
		
		# get this full-stack enrichment values across cell types
		this_row_enrichment_val_data = max_enrich_value_all_ct_df.iloc[row_index, 1:] # skip the first column because that's the state column, get the data for this full-stack state
		ranked_enrichment_value_this_row = this_row_enrichment_val_data[ranked_ct_this_row] # reorder the enrichment values, descending
		rank_enrichment_value_df.loc[row_index] = [full_stack_state] + list(ranked_enrichment_value_this_row)

	# color cell type
	ct_colored_df = color_rank_ct_df(rank_ct_df)
	# color 25-state
	colored_25_state_df = color_rank_25_state_df(rank_25_state_df)
	# color 25_state_cell_group_df: df where cell types of the same cell group are juxtaposed
	colored_25_state_cell_group_df = get_max_enrichment_25_state_df_ordered_cell_group(max_enrich_25_state_all_ct_df, max_enrich_value_all_ct_df)

	# save the excel
	writer = pd.ExcelWriter(output_fn, engine='xlsxwriter')
	ct_colored_df.to_excel(writer, sheet_name = 'ranked_cell_type')
	colored_25_state_df.to_excel(writer, sheet_name = 'ranked_25_state')
	colored_25_state_cell_group_df.to_excel(writer, sheet_name = 'cell_group_25_state')
	rank_enrichment_value_df.to_excel(writer, sheet_name = 'ranked_enrichment_value')
	writer.save()



def get_ranked_max_enriched_25_state(input_folder, output_folder, num_state_ct_model):
	output_fn = os.path.join(output_folder, "max_enriched_states_all_cell_types.xlsx")
	all_ct_folders = glob.glob(input_folder + "/E*/")
	all_ct_list = map(lambda x: x.split('/')[-2], all_ct_folders) # from '/path/to/E129/' to 'E129' 
	all_ct_fn_list = map(lambda x: x + "/overlap_enrichment.txt", all_ct_folders)
	all_ct_df_list = map(lambda x: get_one_enrichment_ct_model_df(x, num_state_ct_model), all_ct_fn_list)

	max_enrich_value_all_ct_df = pd.DataFrame({'state' : (all_ct_df_list[0])['state']}) # create a data frame with only a column called state that are the states in the full-stack model. This data frame store the enrichment values
	max_enrich_25_state_all_ct_df = pd.DataFrame({'state' : (all_ct_df_list[0])['state']}) # this data frame stores the names of states that are most enriched with each of the full-stack state in each of the cell type that we look at
	for ct_index, ct in enumerate(all_ct_list):
		this_ct_fn = all_ct_fn_list[ct_index] # get the overlap_enrichment file associated with this cell type
		this_ct_enrichment_df = get_one_enrichment_ct_model_df(this_ct_fn, num_state_ct_model) # get the enrichment data frame associated with this cell type
		max_enrich_value_all_ct_df[ct] = this_ct_enrichment_df['max_enrichment'] # store the values of maximum enrichment in this cell type
		max_enrich_25_state_all_ct_df[ct] = this_ct_enrichment_df['max_enrichment_context'] # store the names of the state that is most enriched with each of the full-stack state for this one cell type

	# by now, we have obtained all the data from all the cell types. Now we have to order each row in the descending enrichment value
	rank_df = (max_enrich_value_all_ct_df.drop('state', axis = 1)).rank(axis = 1, ascending = False) # rank all the enrichment value within each row (full-stack state) in descending order --> descending max_enrichment across different cell type, with respect to the full-stack state that we are looking at. Each row corresponds to one full-stack state
	print "Done getting data from all cell types after: " + str(time.clock() - start_time)
	get_rank_25_state_df(rank_df, max_enrich_25_state_all_ct_df, max_enrich_value_all_ct_df, output_fn)
	print "Done ranking enrichment data across cell types after: " + str(time.clock() - start_time)


def main():
	if len(sys.argv) != 3:
		usage()
	input_folder = sys.argv[1]
	if not os.path.isdir(input_folder):
		print "input_folder is not a valid directory"
		usage()
	output_folder = sys.argv[2]
	cmh.make_dir(output_folder)
	print "Done getting command line arguments after: " + str(time.clock() - start_time)
	#get_all_ct_model_median_enrichment(input_folder, output_folder, num_state_ct_model)
	num_state_ct_model = 25 # I hard code it here, but it can be changed to be a parameter in the future
	get_ranked_max_enriched_25_state(input_folder, output_folder, num_state_ct_model)
	print "Done after: " + str(time.clock() - start_time)

def usage():
	print "python get_ranked_max_enriched_25_state.py"
	print "input_folder: where there are multiple subfolders, each containing enrichment data for different cell type specific model"
	print "output_fn: Where the 25-state with maximum fold enrichment across all 25-state states are reported for all cell types for each state is reported"
	print "num_state_ct_model must be 25. I write this code to be hard-coded work for the small models to have 25 states"
	exit(1)

main()