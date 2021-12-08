import pandas as pd
import seaborn as sns
import numpy as np
import os
import sys
import helper
from glob import glob
import string
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
TWENTY_FIVE_STATE_COLOR = {1: 'Red', \
2: '#ff4500', 3: '#ff4500', 4: '#ff4500', \
# orange red \
5: 'Green', 6: 'Green', 7: 'Green', \
8: '#90EE90', \
#Lighter green\
9: '#CCFF00', 10: '#CCFF00', 11: '#CCFF00', 12: '#CCFF00', \
# electric lime
13: 'Orange', 14: 'Orange', 15: 'Orange', \
16: 'Yellow', 17: 'Yellow', 18: 'Yellow', \
19: '#fff44f', \
# Lemon
20: '#7fffd4', \
# Aquamarine
21: '#b19cd9', \
# Light Purple
22: 'Pink', \
23: '#901068', \
# Dark purple
24: 'Gray', \
25: 'White'}

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



def get_rid_of_stupid_file_tail(context_name):
	if context_name.endswith('.bed.gz'):
		return(context_name[:-7])
	else:
		return(context_name)


def get_one_enrichment_ct_model_df(fn):
	df = pd.read_csv(fn, sep = '\t', header = 0)
	df = df.rename(columns = {'state (Emission order)' : 'state', 'Genome %': 'percent_in_genome'}) # rename some columns so that it is easier to write column names later
	context_name_list = list(map(get_rid_of_stupid_file_tail, df.columns[2:])) # get the names of enrichment contexts. If the enrichment contexts is <context>.bed.gz --> <context>
	df.columns = ['state', 'percent_in_genome'] + context_name_list
	df = df.drop([df.shape[0] - 1], axis = 0) # drop that last row, which is about the percentage of genome that each enrichment context occupies
	return df


def get_mmm_enrichment_one_genome_context(all_ct_df_list, all_ct_name_list, state_context_colName): 
	# state_context_colName: name of the column that we are looking at so that we know what column to look at for each data frame
	this_context_df = pd.DataFrame({'state' : (all_ct_df_list[0])['state']})
	# collect data from each of the cell type specific model --> columsn: cell type specific model, rows: 100 full stack states
	for df_index, df in enumerate(all_ct_df_list):
		ct_name = all_ct_name_list[df_index] # all_ct_fn_list and all_ct_df_list are responsive to each other, i.e. each element in each list corresponds to the same cell type. We tried zipping teh two list, but that gave a bug message. so this code here is not the most graceful.
		this_context_df[ct_name] = df[state_context_colName]
	this_context_df['max_enrichment'] = (this_context_df.drop('state', axis = 1)).max(axis = 1)
	this_context_df['min_enrichment'] = (this_context_df.drop('state', axis = 1)).min(axis = 1)
	this_context_df['median_enrichment'] = (this_context_df.drop('state', axis = 1)).median(axis = 1)
	return this_context_df

def get_25_state_annot(state):
	# '1' --> ''1_TssA''
	state_index = int(state) - 1 # zero-based
	return TWENTY_FIVE_STATE_ANNOTATION[state_index]

def get_TWENTY_FIVE_STATE_COLOR(state):
	# '1_TssA' --> 'red'
	state_index = int(state.split('_')[0]) # one-based
	color = TWENTY_FIVE_STATE_COLOR[state_index]
	return 'background-color: %s' % color

def paint_excel_mmm_enrichment(enrichment_df):
	# here enrichment_df is actually max_enrichment_df from get_all_ct_model_mmm_enrichment
	cm = sns.light_palette("red", as_cmap=True)
	(num_state, num_enr_cont) = (enrichment_df.shape[0], enrichment_df.shape[1] - 1)

	enr_cont_list = enrichment_df.columns[2:]
	enrichment_df['state'] = enrichment_df['state'].apply(get_25_state_annot) # 1 --> 1_TssA
	colored_df = enrichment_df.style.background_gradient(subset = pd.IndexSlice[:, enr_cont_list], cmap = cm) # color the enrichment values into a red-white gradient in the first enrichment contnext
	colored_df = colored_df.applymap(get_TWENTY_FIVE_STATE_COLOR, subset = pd.IndexSlice[:, ['state']])
	return colored_df


def get_all_ct_model_mmm_enrichment(input_folder, output_folder, num_state_ct_model):
	all_ct_folders = glob(input_folder + "/E*/")
	all_ct_fn_list = list(map(lambda x: x + "/overlap_whole.txt", all_ct_folders))
	all_ct_name_list = list(map(lambda x: (x.split("/"))[-2], all_ct_folders))
	all_ct_df_list = map(lambda x: get_one_enrichment_ct_model_df(x), all_ct_fn_list)
	all_ct_df_list = list(all_ct_df_list)
	# up until here, we have finished getting all the data from all the overlap enrichment files of all cell types into data frame --> all_ct_df_list [index] --> data frame of the cell type
	print ("Done getting all enrichment data from all cell types")
	# create a data frame that will report the median enrichment over all the cell types
	max_enrichment_df = pd.DataFrame({'state': (all_ct_df_list[0])['state']}) # this will report for each state in the 25-state model, and for each enrichment context, what are the highest enrichment in each state-context combination. 
	# now loop through each of the context that we care about
	median_enrichment_df = pd.DataFrame({'state': (all_ct_df_list[0])['state']}) # this will report for each state in the 25-state model, and for each enrichment context, what are the highest enrichment in each state-context combination. 
	# now loop through each of the context that we care about
	min_enrichment_df = pd.DataFrame({'state': (all_ct_df_list[0])['state']}) # this will report for each state in the 25-state model, and for each enrichment context, what are the highest enrichment in each state-context combination. 
	# now loop through each of the context that we care about
	enrichment_context_list = ['percent_in_genome'] + list((all_ct_df_list[0]).columns[2:]) # this is because we also want to report the max, median, min of percent_in_genome_of data from all ct-spec enrichments
	print (enrichment_context_list)
	for enr_context in enrichment_context_list:
		this_context_mmm_enrichment_df = get_mmm_enrichment_one_genome_context(all_ct_df_list, all_ct_name_list, enr_context)
		max_enrichment_df[enr_context] = this_context_mmm_enrichment_df['max_enrichment']
		min_enrichment_df[enr_context] = this_context_mmm_enrichment_df['min_enrichment']
		median_enrichment_df[enr_context] = this_context_mmm_enrichment_df['median_enrichment']		
	print ("Done getting all the necessary data!")
	max_colored_df = paint_excel_mmm_enrichment(max_enrichment_df)
	min_colored_df = paint_excel_mmm_enrichment(min_enrichment_df)
	median_colored_df = paint_excel_mmm_enrichment(median_enrichment_df)
	# now save into 3 sheets
	output_fn = os.path.join(output_folder,  "mmm_enrichment_all_ct.xlsx")
	writer = pd.ExcelWriter(output_fn, engine='xlsxwriter')
	max_colored_df.to_excel(writer, sheet_name='max')
	min_colored_df.to_excel(writer, sheet_name='min')
	median_colored_df.to_excel(writer, sheet_name='median')
	writer.save()
	print ("Done getting the mmm_enrichment_df!")

def main():
	if len(sys.argv) != 4:
		usage()
	input_folder = sys.argv[1]
	if not os.path.isdir(input_folder):
		print ("input_folder is not a valid directory")
		usage()
	output_folder = sys.argv[2]
	helper.make_dir(output_folder)
	try: 
		num_state_ct_model = int(sys.argv[3]) # number of state in the cell type specific model, 15 25 or 18
	except:
		print ("num_state_ct_model is not valid")
		usage() 
	print ("Done getting command line arguments!")
	get_all_ct_model_mmm_enrichment(input_folder, output_folder, num_state_ct_model)
	print ("Done!")

def usage():
	print ("python get_mmm_enrichments_across_ct.py") # mmm stands for min median max
	print ("input_folder: where thare are multiple subfolders, each containing enrichment data for different cell type specific model")
	print ("output_folder: Where the median fold enrichment across all cell types for each state is reported")
	print ("num_state_ct_model: number of states in the cell type specific models")
	exit(1)

main()