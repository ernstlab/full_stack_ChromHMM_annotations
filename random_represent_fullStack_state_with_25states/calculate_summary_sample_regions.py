# this file will take output of code: count_sample_region_per_full_stack_state.py and produce some summary tables of the sample regions with the different seeds. 
# one summary table will have full_stack state as rows, each column will correspond to <cell_group>_<state25>, where state25 is each of the 25 state. The table will be painted in excel and the color scale will correspond to the the color of the 25 states. It will be a beautiful excel. 
import os 
import sys
import pandas as pd 
import helper
import numpy as np 
import seaborn as sns
import glob
from scipy.stats import mannwhitneyu
NUM_FULL_STACK_state = 100
NUM_CELL_TYPE_SPEC_state = 25 # because this fiel is inside the folder associated with 25-state cell type specific system. 
NUM_NON_CT_COLUMNS = 4 # the first NUM_NON_CT_COLUMNS in bed_map_df are chrom, start, end, full_stack
def rgb2hex(rgb_str): # I took this function from online: https://stackoverflow.com/questions/3380726/converting-a-rgb-color-tuple-to-a-six-digit-code
	r,g,b = tuple(map(lambda x: int(x), rgb_str.split(',')))
	return "#{:02x}{:02x}{:02x}".format(r,g,b)

def get_state25_annot_df():
	state25_annot_fn = '../../..//data/roadmap_epigenome/25_imputed12marks_model_downloaded/state_annot.txt'
	state25_df = pd.read_csv(state25_annot_fn, header = 0, index_col = None, sep = '\t')
	state25_df['color_hex'] = state25_df['rgb'].apply(rgb2hex)
	return state25_df

def get_state25_to_mnemonic_dict(state25_df): # key: state of the form 1 --> 25, values: state mnemonic
	return dict(zip(state25_df.state, state25_df.mnemonic))

def get_full_stack_annot_df():
	fS_annot_fn = '../../..//ROADMAP_aligned_reads/chromHMM_model/model_100_state/figures/supp_excel/state_annotations_processed.csv'
	fS_annot_df = pd.read_csv(fS_annot_fn, header = 0, index_col = None, sep = ',') # state, mneumonics, Long annotations, Short Annotations, Group, color, state_order_by_group. state: 1--> 100, color: hex
	fS_annot_df = fS_annot_df[['state', 'mneumonics', 'color', 'state_order_by_group']]
	return fS_annot_df

def read_ct_annot_df():
	roadmap_ct_annot_fn = '../../..//ROADMAP_aligned_reads/ROADMAP_metadata_july2013.csv'
	ct_annot_df = pd.read_csv(roadmap_ct_annot_fn, header = 0, index_col = None, sep = ',')
	ct_annot_df = ct_annot_df.rename(columns = {'Epigenome ID (EID)': 'eid'})
	cell_group_color_dict = dict(zip(ct_annot_df.GROUP, ct_annot_df.COLOR))
	ct_annot_df = ct_annot_df[['eid', 'GROUP']] #, 'COLOR', 'ANATOMY']] 
	# ct_annot_df = ct_annot_df.sort_values('GROUP') # sort by gruups so that we can later arrange the data in specific order
	return ct_annot_df, cell_group_color_dict

def process_one_seedFile(bed_map_fn): # we process the result of one file that contains the sampled reference 
	bed_map_df = pd.read_csv(bed_map_fn, header = 0, index_col = None, sep = '\t')
	# rearrange the columns such that cell types of the same group are together
	num_segment_per_state = bed_map_df.shape[0] / NUM_FULL_STACK_state
	# chrom, start, end, full_stack, E001, E002, etc. 
	ct_colnames = pd.Series(bed_map_df.columns).str.startswith('E') # a list of True or False, True if the corrsponding column name starts with E, which is to say it signifies a cell type in ROADMAP
	input_ct_list = bed_map_df.columns[ct_colnames]
	all_state25_list = list(map(lambda x: 'E' + str(x+1), range(NUM_CELL_TYPE_SPEC_state))) # E1 --> E25
	bed_map_df = bed_map_df.groupby('full_stack') 
	result_df = pd.DataFrame(columns = ['state25', 'ct', 'prop_in_ct', 'full_stack'])
	for full_stack_state, state_df in bed_map_df:# loop through each group of full-stack state
		state_df = state_df[input_ct_list] 
		count_df = state_df[input_ct_list].apply(pd.Series.value_counts, normalize = True) # only get columns corresponding to the cell types in ROADMAP, then count the number of occurrences of each of the 25 states in each ROADMPA cell type. Then, calculate the fraction of occurrence for each of the 25 state in each of the cell type (normalize = True). There may not be all 25 states, and there will be NAN values (later replace by 0). Rows: each of the 25 states. columns: ROADMAP cell types
		count_df = count_df.fillna(0) # Rows: each of the 25 states (may not be complete all 25 states). columns: ROADMAP cell types
		missing_state_list = np.setdiff1d(all_state25_list, count_df.index)
		count_df = count_df.append(pd.DataFrame(0, columns = count_df.columns, index = missing_state_list)) # add the missing states from the 25 states, and say that the proportion of these states in the sampled regions in each ct is 0
		count_df = count_df.reset_index().rename(columns = {'index' : 'state25'}) # one more column showing the state E1 --> E25
		count_df = count_df.melt(id_vars = ['state25']).rename(columns = {'variable' : 'ct', 'value' : 'prop_in_ct'}) # 3 columns: state25, ct, prop_in_ct showing the proportion that each state25 is sampled in each ct from ROADMAP
		count_df['full_stack'] = full_stack_state # 4 columns: state25, ct, prop_in_ct , full_stack states
		result_df = result_df.append(count_df)
	return result_df

def calculate_avg_proportion_across_ct_per_group(oneSeed_prop_df, ct_annot_df):
	# ct_annot_df has 2 columns: eid and GROUP
	# oneSeed_prop_df is a df with 4 columns state25, ct, prop_in_ct, full_stack
	oneSeed_prop_df = oneSeed_prop_df.merge(ct_annot_df, how = 'left', left_on = 'ct', right_on = 'eid') # 6 columns state25, ct, prop_in_ct, full_stack, eid, GROUP
	piv_df = pd.pivot_table(oneSeed_prop_df, values = 'prop_in_ct', index = ['full_stack', 'state25'], columns = ['GROUP'], aggfunc = {'prop_in_ct' : np.mean}) # a pivot table index has two layers: full_stack and state25. Columns correspond to GROUP
	piv_df = piv_df.unstack() # df has columns with two layers: GROUP and state25, index: full_stack
	return piv_df

def color_full_stack_state(column_data, full_stack_color_dict):
	results = [''] * len(column_data)
	for row_index, value in enumerate(column_data):
		if row_index < NUM_FULL_STACK_state:
			try:	
				color = full_stack_color_dict[value]
			except:
				color = '#ffffff' # white
		else:
			color = '#ffffff' # white
		results[row_index] = 'background-color: %s' % color
	return results

def color_annot_cellGroup_state25(column_data, cell_group_color_dict, state25_color_dict): 
	results = [''] * len(column_data)
	for index, value in enumerate(column_data):
		if index == 0:
			try:
				color = cell_group_color_dict[value]
			except:
				color = '#ffffff' # white
		elif index == 1:
			try:
				color = state25_color_dict[value]
			except:
				color = '#ffffff' # white
		else:
			color = '#ffffff' # white
		results[index] = 'background-color: %s' % color
	return results

def color_avg_allSeed_group_df(avg_allSeed_group_df, state25_df):
	# test_df has exactly the same columns as avg_allSeed_group_df: mneumonics, <cell_group>_E<state25>, state (1 ,..., 100), color
	full_stack_color_dict = dict(zip(avg_allSeed_group_df.mneumonics, avg_allSeed_group_df.color))
	state25_color_dict = dict(zip(state25_df.mnemonic, state25_df.color_hex))
	state25_color_dict['Quies'] = '#3498DB' # change the color of quiescent state for better visibility
	colored_df = avg_allSeed_group_df.style.apply(lambda x: color_full_stack_state(x, full_stack_color_dict), axis = 0) # color full_stack states column
	# now onto coloring individual states
	NUM_ROWS_TO_COLOR_STATE25_WISE = NUM_FULL_STACK_state + 1  # full stack states and min and max, state25-block wise
	NUM_ROWS_TO_COLOR_THEOR_WISE = NUM_FULL_STACK_state + 3 # full stack states and theoretical min (0) and max (1)
	for state_index in range(NUM_CELL_TYPE_SPEC_state):
		state_str = 'E' + str(state_index + 1) # 0 --> E1
		state_color = state25_color_dict[state25_df.loc[state_index, 'mnemonic']]
		cm = sns.light_palette(state_color, as_cmap=True) # the color map for this state25
		colnames_to_colors = pd.Series(avg_allSeed_group_df.columns).str.endswith('_' + state_str) # a list of true or false. True if the column name ends with the state that we want to color
		colnames_to_colors = avg_allSeed_group_df.columns[colnames_to_colors] #now just a list of column names to pass into pd.IndexSlice
		for colname in colnames_to_colors:
			colored_df = colored_df.background_gradient(subset = pd.IndexSlice[:NUM_ROWS_TO_COLOR_STATE25_WISE, [colname]], cmap = cm)
	return colored_df

def helper_get_state25_from_colname(colname, state25_to_mnemonic_dict):
	s_list = colname.split('_')
	if len(s_list) < 2:
		return ''
	else: 
		return state25_to_mnemonic_dict[int(s_list[1][1:])]

def color_column_annotations(avg_df_columns, cell_group_color_dict, state25_df): #avg_df_columns should be the columns of avg_allSeed_group_df
	state25_to_mnemonic_dict = get_state25_to_mnemonic_dict(state25_df)
	state25_color_dict = dict(zip(state25_df.mnemonic, state25_df.color_hex))
	state25_color_dict['Quies'] = '#3498DB'	 # change the color of the quiescent state for better visibility
	df = pd.DataFrame(columns = avg_df_columns)
	df.loc[0] = list(map(lambda x: x.split('_')[0], avg_df_columns))
	df.loc[1] = list(map(lambda x: helper_get_state25_from_colname(x, state25_to_mnemonic_dict), avg_df_columns))
	colored_df = df.style.apply(lambda x: color_annot_cellGroup_state25(x, cell_group_color_dict, state25_color_dict), axis = 0)
	return colored_df

def calculate_state25_max_min_row(avg_allSeed_group_df):
	"""
	For each state25, there are <#cell_group> columns, each column has 100 rows corrsponding to 100 full-stack state. In this function, we will report the maximum values of proportion in across a subset of avg_allSeed_group_df with nrow = # full_stack_state and ncols  = # cell_group, corresponding to each state. The result  is a list of length <#cell_group> * 25, where every <#cell_group> entry has the same value of maximum/minimum proportion
	"""
	results_max = []
	results_min = []
	for state_index in range(NUM_CELL_TYPE_SPEC_state):
		state_str = 'E' + str(state_index + 1) # 0 --> E1
		colnames_to_this_state = pd.Series(avg_allSeed_group_df.columns).str.endswith('_' + state_str) # a list of true or false. True if the column name ends with the state that we want to color
		colnames_to_this_state = avg_allSeed_group_df.columns[colnames_to_this_state] #now just a list of column names to pass into pd.IndexSlice
		this_group_max = (avg_allSeed_group_df[colnames_to_this_state]).values.max(axis = None) # max value of the entire dataframe (this subset of the big dataframe)
		results_max += [this_group_max] * len(colnames_to_this_state)
		this_group_min = (avg_allSeed_group_df[colnames_to_this_state]).values.min(axis = None) # min value of the entire dataframe (this subset of the big dataframe)
		results_min += [this_group_min] * len(colnames_to_this_state)
	return results_max, results_min

def test_wilcoxon_group_spec(oneSeed_proportion_df_list, ct_annot_df, fS_annot_df, output_folder):
	"""
	oneSeed_proportion_df_list from count_sample_region_each_full_stack_state: each item in the list is a df with 4 columns state25, ct, prop_in_ct, full_stack. Each df corresponds to the results from one random generator seed
	ct_annot_df has columns: eid, GROUP
	Note: both full_stack and state25 are now in the form E<state_index_1_based>
	"""
	# allSeed_prop_df = pd.concat(oneSeed_proportion_df_list) # 1 big df with 4 columns state25, ct, prop_in_ct, full_stack
	allSeed_prop_df = pd.read_csv('./trial.gz', header = 0, index_col = None, sep = '\t') # this is for debug reasons
	allSeed_prop_df = allSeed_prop_df.merge(ct_annot_df, how = 'left', left_on = 'ct', right_on = 'eid')
	allSeed_prop_df = allSeed_prop_df.drop(columns = ['ct', 'eid']) # drop the unnecessary columns
	uniq_group_list = np.unique(allSeed_prop_df.GROUP) # I confirm that this list will be element-wise the same as the one we declared in function
	output_columns = []
	for state in range(NUM_CELL_TYPE_SPEC_state):
		output_columns += list(map(lambda x: x + "_E" + str(state+1), uniq_group_list)) 
	result_df = pd.DataFrame(columns = ['full_stack'] + output_columns) # we add full_stack but not include in output_columns because we want to use output_columns later when merge with fS_state_annot_df
	allSeed_prop_df = allSeed_prop_df.groupby('full_stack')
	for fS_state in range(NUM_FULL_STACK_state):
		state_str = 'E' + str(fS_state + 1)
		this_fS_df = allSeed_prop_df.get_group(state_str).groupby('state25') # this df contains all data associated with this full_stack state
		result_row = [fS_state+1] # full_stack_state, 1-based, int
		for state25 in range(NUM_CELL_TYPE_SPEC_state): # 25
			state25_str = 'E' + str(state25+1)
			this_s25_df = this_fS_df.get_group(state25_str) # state25: current state25, ct, prop_in_ct, full_stack: current fS state, GROUP
			for cell_group in uniq_group_list:
				row_bool_in_cGroup = pd.Series(this_s25_df['GROUP']).str.match(cell_group)
				x = np.array(this_s25_df['prop_in_ct'][row_bool_in_cGroup]) # in this group
				y = np.array(this_s25_df['prop_in_ct'][~row_bool_in_cGroup]) # NOT in the group
				try:  
					t = mannwhitneyu(x, y , use_continuity = False, alternative = 'greater')
					result_row.append(t.pvalue)
				except: # if the test crash, that's because all values in x and y are identical --> p value should be 1
					result_row.append(1)
		result_df.loc[result_df.shape[0]] = result_row
	result_df = result_df.merge(fS_annot_df, how = 'left', left_on = 'full_stack', right_on = 'state')
	result_df = result_df.sort_values('state_order_by_group') 
	result_df.reset_index(drop = True, inplace = True) # reset index after reordering the rows 
	result_df = result_df[['mneumonics'] + output_columns + ['state', 'color']] # rearrange the columns to be similar to avg_allSeed_group_df
	save_fn = os.path.join(output_folder, 'pval_mannU_test_across_ct.csv.gz')
	result_df.to_csv(save_fn, header = True, index = False, sep = '\t', compression = 'gzip') # debug
	return result_df

def color_significant_pvalue(column_data, pval_thres):
	results = [''] * len(column_data)
	for row_index, value in enumerate(column_data):
		if isinstance(value, float) and value <= pval_thres:
			results[row_index] = 'background-color: #A7CCE5' # blue sky
	return results

def color_pvalue_df(test_df, pval_thres):
	print("pval_thres: " + str(pval_thres))
	full_stack_color_dict = dict(zip(test_df.mneumonics, test_df.color))
	colored_df = test_df.style.apply(lambda x: color_full_stack_state(x, full_stack_color_dict), axis = 0) # color full_stack states column
	colored_df = colored_df.apply(lambda x: color_significant_pvalue(x, pval_thres), axis = 0) # color each column
	return colored_df

def report_only_significant_cells(test_df, pval_stringent_thres = 1.0E-100):
	# test_df: mneumonics, <cell_group>_<state25:Esomething>, state, color
	# pval_stringent_thres is actually much samller than the bonferroni-corrected p-value threshold
	pval_df = test_df.iloc[:,1:-2] # skip the first and the last two columns beacuse those correspond to mneumonics, state and color
	rows_pick = pval_df.le(pval_stringent_thres).any(axis = 1) # T/F for each rows that has >=1 pvalues <= pval_stringent_thres
	columns_pick = pval_df.columns[pval_df.le(pval_stringent_thres).any(axis = 0)] # list of column names for columns that has >= 1 pvalue < pval_stringent_thres	
	uniq_group_list = np.unique(list(map(lambda x: x.split('_')[0], columns_pick)))
	rearranged_colnames = [] # rearrange such that columns of the same cell groups are put next to each other
	columns_pick = pd.Series(columns_pick) # transform to pandas Series to use some useful functions
	for group in uniq_group_list:
		colnames_this_group = columns_pick[columns_pick.str.startswith(group)]
		rearranged_colnames += list(colnames_this_group)
	significant_df = test_df.loc[rows_pick, ['mneumonics'] + list(rearranged_colnames) + ['state', 'color']] # select only rows and columns with the significant p values, add the columns corresponding to mneumonics, state, color
	return significant_df

def create_excel_from_precalculated_data(output_folder, all_seed_folder):
	# This function was designed out of necessity for debuging purposes. We need this function after we have calculated the data for files pval_mannU_test_across_ct.csv.gz and avg_prop_state25_per_group.csv.gz through the function count_sample_region_each_full_stack_state, but the formatting kept being modified and we don't want to keep recalculating these files through the function count_sample_region_each_full_stack_state. So, to save time, I made this function. 
	ct_annot_df, cell_group_color_dict = read_ct_annot_df() # eid, GROUP. We need the cell_group_color_dict for later use of coloring the excel file
	state25_df = get_state25_annot_df() # state, mnemonic, description, color_string, rgb, color_hex
	fS_annot_df = get_full_stack_annot_df() # 'state', 'mneumonics', 'color', 'state_order_by_group'
	test_df = pd.read_csv(os.path.join(output_folder, 'pval_mannU_test_across_ct.csv.gz'), header = 0, index_col = None, sep = '\t') 
	avg_allSeed_group_df = pd.read_csv(os.path.join(output_folder, "avg_prop_state25_per_group.csv.gz"), header = 0, index_col = None, sep = '\t') # this is for debug
	uniq_group_list = np.unique(list(map(lambda x: x.split('_')[0],avg_allSeed_group_df.columns[1:-2]))) # this is for debug
	rearranged_colnames = []
	for state25 in range(NUM_CELL_TYPE_SPEC_state):
		rearranged_colnames += list(map(lambda x: x + '_E' + str(state25+1), uniq_group_list))
	# record the raw data, before we add some rows about the max and min values of each column
	num_columns_with_prop_values = len(rearranged_colnames) # number of columns that are associated with the numbers showing proportion of the sampled regions that are in the state25 for each cell group 
	num_columns_with_empty_values = avg_allSeed_group_df.shape[1] - (1 + num_columns_with_prop_values) # number of columns that we will fill with empty at the end when we try to append some few rows before making the colored excel file. the 1 in the code corresponds to column 'full_stack'
	state25_wise_max, state25_wise_min = calculate_state25_max_min_row(avg_allSeed_group_df)
	avg_allSeed_group_df.loc[avg_allSeed_group_df.shape[0]] = ['state25_wise_min'] + state25_wise_min + [''] * num_columns_with_empty_values # add one more row showing the min value in each block of state25 in the table
	avg_allSeed_group_df.loc[avg_allSeed_group_df.shape[0]] = ['state25_wise_max'] + state25_wise_max + [''] * num_columns_with_empty_values # add one more row showing the max value in each block of state25 in the table
	avg_allSeed_group_df.loc[avg_allSeed_group_df.shape[0]] = ['theoretical_min'] + [0] * num_columns_with_prop_values + [''] * num_columns_with_empty_values # add one more row showing the min value/ This is needed because we want to control the color of the cells in excel to be uniformly ranged (0,1) later
	avg_allSeed_group_df.loc[avg_allSeed_group_df.shape[0]] = ['theoretical_max'] + [1] * num_columns_with_prop_values + [''] * num_columns_with_empty_values # add one more row showing the max value/ This is needed because we want to control the color of the cells in excel to be uniformly ranged (0,1) later
	print ("Done getting all the necessary data for the excel")
	main_colored_df = color_avg_allSeed_group_df(avg_allSeed_group_df, state25_df)
	all_column_annot_colored_df = color_column_annotations(avg_allSeed_group_df.columns, cell_group_color_dict, state25_df)
	ALPHA = 0.01
	bonreffroni_pval_thres = ALPHA / (NUM_CELL_TYPE_SPEC_state * NUM_FULL_STACK_state * len(uniq_group_list)) # bonreffroni corrected p-value threshold
	pval_colored_df = color_pvalue_df(test_df, bonreffroni_pval_thres)
	pval_stringent_thres = 1.0E-100
	stringent_pval_df = report_only_significant_cells(test_df, pval_stringent_thres)
	significant_colored_df = color_pvalue_df(stringent_pval_df, pval_stringent_thres)
	significant_column_annot_df = color_column_annotations(stringent_pval_df.columns, cell_group_color_dict, state25_df)
	output_fn = os.path.join(output_folder,  "avg_prop_state25_per_group.xlsx")
	writer = pd.ExcelWriter(output_fn, engine='xlsxwriter')
	main_colored_df.to_excel(writer, sheet_name='avg_prop_state25_per_group')
	all_column_annot_colored_df.to_excel(writer, sheet_name = 'all_column_annot')
	pval_colored_df.to_excel(writer, sheet_name = 'pval_mwu_1s')
	significant_colored_df.to_excel(writer, sheet_name = 'pval_thres_1.0E-100')
	significant_column_annot_df.to_excel(writer, sheet_name = 'significant_column_annot')
	writer.save()
	print ("Done painting the excel")
	return

def count_sample_region_each_full_stack_state(all_seed_folder, output_folder):
	ct_annot_df, cell_group_color_dict = read_ct_annot_df() # eid, GROUP. We need the cell_group_color_dict for later use of coloring the excel file
	state25_df = get_state25_annot_df() # state, mnemonic, description, color_string, rgb, color_hex
	fS_annot_df = get_full_stack_annot_df() # 'state', 'mneumonics', 'color', 'state_order_by_group'
	bedmap_fn_list = glob.glob(all_seed_folder + 'seed_*/sample_segment_fullStack_ct25State.bed.gz')
	oneSeed_proportion_df_list = list(map(process_one_seedFile, bedmap_fn_list)) # each item in the list is a df with 4 columns state25, ct, prop_in_ct, full_stack
	allSeed_prop_df = pd.concat(oneSeed_proportion_df_list) #for debug
	test_df = test_wilcoxon_group_spec(oneSeed_proportion_df_list, ct_annot_df, fS_annot_df, output_folder)
	oneSeed_proportion_df_list = list(map(lambda x: calculate_avg_proportion_across_ct_per_group(x, ct_annot_df), oneSeed_proportion_df_list)) # list of df with columns with two layers: GROUP and state25, index: full_stack
	print("Done processing all input data")
	avg_allSeed_group_df = pd.concat(oneSeed_proportion_df_list).groupby(level = 0).mean()
	uniq_group_list = np.unique(avg_allSeed_group_df.columns.get_level_values('GROUP'))
	avg_allSeed_group_df.columns = ['_'.join(col) for col in avg_allSeed_group_df.columns.values] # now columns are just one layer of the form 
	# after the following line of code, there is one more column to this dataframe
	avg_allSeed_group_df.reset_index(drop = False, inplace = True) # full_stack, Adipose_E1  Blood & T-cell_E1  Brain_E1 --> Sm. Muscle_E25  Thymus_E25  iPSC_E25
	avg_allSeed_group_df['full_stack'] = avg_allSeed_group_df['full_stack'].apply(lambda x: int(x[1:])) # E1 --> 1
	avg_allSeed_group_df = avg_allSeed_group_df.merge(fS_annot_df, how =  'left', left_on = 'full_stack', right_on = 'state').fillna('')
	avg_allSeed_group_df = avg_allSeed_group_df.sort_values('state_order_by_group')
	avg_allSeed_group_df.reset_index(drop = True, inplace = True) # after the sorting step, the indices are not in order so we just reset it.
	avg_allSeed_group_df = avg_allSeed_group_df.drop(columns = ['full_stack', 'state_order_by_group'])
	rearranged_colnames = []
	for state25 in range(NUM_CELL_TYPE_SPEC_state):
		rearranged_colnames += list(map(lambda x: x + '_E' + str(state25+1), uniq_group_list))
	avg_allSeed_group_df = avg_allSeed_group_df[['mneumonics'] + rearranged_colnames + ['state', 'color']] # rearrange columns
	# record the raw data, before we add some rows about the max and min values of each column
	output_csv_fn = os.path.join(output_folder, "avg_prop_state25_per_group.csv.gz")
	avg_allSeed_group_df.to_csv(output_csv_fn, header = True, index = False, sep = '\t', compression = 'gzip')
	print ("Done painting the excel")
	return

def main():
	if len(sys.argv) != 3:
		usage()
	all_seed_folder = sys.argv[1]
	helper.check_dir_exist(all_seed_folder)
	output_folder = sys.argv[2]
	helper.make_dir(output_folder)
	print("Done getting command line arguments")
	count_sample_region_each_full_stack_state(all_seed_folder, output_folder)	
	create_excel_from_precalculated_data(output_folder, all_seed_folder)
	print("Done!")
	return 

def usage():
	print("python calculate_summary_sample_regions.py")
	print("all_seed_folder: where there are subfolders seed_<seed>. Inside each seed folder there should be files sample_state_sorted_by_ct.txt.gz and sample_state_sorted_by_group.txt.gz")
	print("output_folder: where we store the excel files")	
	exit(1)
main()