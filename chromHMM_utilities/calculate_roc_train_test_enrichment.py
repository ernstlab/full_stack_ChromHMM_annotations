import sys
import os
import numpy as np
from collections import Counter
import pandas as pd
from sklearn.metrics import auc
import chromHMM_utilities_common_functions_helper as helper
TRAIN_INDEX = 0
TEST_INDEX = 1
PRECISION_INDEX = 0
RECALL_INDEX = 1
PRECAL_HEADER_LIST = ["_precision", "_recall"]
PRECAL_AXIS_NAME_LIST = ["Precision", "Recall"]

TRUE_POS_INDEX = 0
FALSE_POS_INDEX = 1
TRUE_FALSE_HEADER_LIST = ["_true_pos", "_false_pos"]
TRUE_FALSE_AXIS_NAME_LIST = ["True Positive rates", "False positive rates"]


X_AXIS_INDEX = 1
Y_AXIS_INDEX = 0
ROC_INDEX = 0 
PRECISION_RECALL_INDEX = 1
DATA_FILE_SUFFIX_LIST = ["_roc.csv", "_prec_recall.csv"]
PLOT_FILE_SUFFIX_LIST = ["_roc.png", "_prec_recall.png"]
HEADER_LIST_LIST = [TRUE_FALSE_HEADER_LIST, PRECAL_HEADER_LIST]
AXIS_NAME_LIST_LIST = [TRUE_FALSE_AXIS_NAME_LIST, PRECAL_AXIS_NAME_LIST]
def get_enrichment_analysis_name_from_file_name(fn):
	"""
	Output from ChromHMM Overlap Enrichment usually denotes the enrichment analysis name as the file names existing in the coordinate directory: For example: mutation_occ1.bed.gz 
	We want the enrichment_analysis_name to be mutation_occ1
	"""
	return (fn.split("."))[0]

def open_one_enrichment_df(fn):
	df = pd.read_csv(fn, sep = '\t', header = 0)
	try:
		df = df.rename(columns={"state (Emission order)": "state", "Genome %": "percent_in_genome", "state (User order)" : "state"})
	except:
		print ("Could not convert data columns headers for data frame: " + fn)
		exit(1)
	return df

def quality_control_all_enrichment_df (data_frame_list):
	headers_list = list(map(lambda x: x.columns, data_frame_list))
	print (list(map(lambda x: len(x) , headers_list)))
	# make sure that all lenth of each of the header_list is the same
	if len(set(list(map(lambda x: len(x) , headers_list)))) != 1: # check that all enrichment files have the same headers
		print ("The length of the headers_list in all enrichment files is not the same. Exiting...")
		exit(1)
	for i in range(len(headers_list[0])): # check that all headers of all the enrichment files are the same
		if len(set(list(map(lambda x: x[i], headers_list)))) != 1:
			print ("Header indexed: " + str(i) + "  of the headers_list is not consistent")
			print (list(map(lambda x: x[i], headers_list)))
			print ("exiting ....")
			exit(1)

	if headers_list[0][0] != "state": # check that the first header is state, in all enrichment files
		print ("The first headers is not state: " + headers_list[0][0])
		print ("exiting ...")
		exit(1)
	if headers_list[0][1] != "percent_in_genome": # check that the second header is percent_in_genome, in all enrichment files
		print ("The second headers is not percent_in_genome: " + headers_list[0][0])
		print ("exiting ...")
		exit(1)
	return

def get_enrichment_data(train_test_enrich_folder_list):
	train_fn_list = list(map(lambda x: os.path.join(x, "train_overlap.txt"), train_test_enrich_folder_list))
	test_fn_list = list(map(lambda x: os.path.join(x, 'test_overlap.txt'), train_test_enrich_folder_list))
	train_df_list = list(map(open_one_enrichment_df, train_fn_list))
	test_df_list = list(map(open_one_enrichment_df, test_fn_list))
	quality_control_all_enrichment_df(train_df_list)
	quality_control_all_enrichment_df(test_df_list)
	return list(zip(train_df_list, test_df_list)) # return a list of tuples. Each item in the list: a pair of df for a model, train (0) and test(1)


def change_fract_state_be_cont_to_one (value):
	# this only happen once during calculation of CDS enrichment for 100-state E129. It's because of decimal point error
	if value > 1.0:
		return 1.0
	return value

def do_roc_analysis(this_cont_df, train_cont_percent, test_cont_percent, enr_cont_name, enrichment_model_name):
	"""
	consider 'gContext' as the analysis that we are talking about
	train_test_df_tuple: has 4 columns train_percent_in_genome, test_percent_in_genome, train_<cont_name>, test_<cont_name>
	"""
	# filter out state that is non existent on the genome
	this_cont_df = this_cont_df.loc[this_cont_df['test_percent_in_genome'] > 0].copy()
	# sort states based on decreasing enrichment of the enrichment_analysis_name, looking at the enrichment values in train data
	this_cont_df.sort_values(by='train_' + enr_cont_name, ascending = False, inplace = True)
	num_rows, num_cols = this_cont_df.shape # ncols should be 4: percent_in_genome, genomic context of interest for train and test data
	# calculate true and false positive rates. From now on all the analysis focuses on the test data, the training fold enrichment statistics were only used for ordering the states
	test_fract_genome_in_cont = test_cont_percent / 100.0 # fraction of the genome that the cont occupies
	test_fract_genome_not_cont = 1 - test_fract_genome_in_cont # fraction of the genome that the cont does not occupy
	this_cont_df['fract_in_genome'] = this_cont_df['test_percent_in_genome'] / 100.0
	this_cont_df['culm_frac_gene_in_state'] = (this_cont_df.fract_in_genome).cumsum()
	this_cont_df['fract_cont_in_state'] = this_cont_df['fract_in_genome'] * this_cont_df['test_' + enr_cont_name]
	this_cont_df['fract_state_be_cont'] = this_cont_df['test_' + enr_cont_name] * test_fract_genome_in_cont
	this_cont_df['fract_state_be_cont'] = this_cont_df['fract_state_be_cont'].apply(change_fract_state_be_cont_to_one)# this only happen once during calculation of CDS enrichment for 100-state E129. It's because of decimal point error
	this_cont_df['fract_genome_in_state_and_cont'] = this_cont_df['fract_state_be_cont'] * this_cont_df['fract_in_genome'] 
	this_cont_df['culm_fract_gene_in_state_and_cont'] = (this_cont_df.fract_genome_in_state_and_cont).cumsum()
	this_cont_df['fract_genome_in_state_not_cont'] = this_cont_df['fract_in_genome'] * (1 - this_cont_df['fract_state_be_cont'])
	this_cont_df['fract_not_cont_in_state'] = this_cont_df['fract_genome_in_state_not_cont'] / test_fract_genome_not_cont
	this_cont_df['true_pos'] = this_cont_df['culm_fract_gene_in_state_and_cont'] / test_fract_genome_in_cont
	this_cont_df['false_pos'] = (this_cont_df.fract_not_cont_in_state).cumsum()
	this_cont_df['precision'] = this_cont_df['culm_fract_gene_in_state_and_cont'] / this_cont_df['culm_frac_gene_in_state']
	this_cont_df['recall'] = this_cont_df['true_pos']
	result_df = this_cont_df[['true_pos', 'false_pos', 'precision', 'recall']]
	first_row = pd.DataFrame({'true_pos' : [0], 'false_pos' : [0], 'precision' : [0], 'recall' : [0]})
	result_df = pd.concat([first_row, result_df]).reset_index(drop = True)
	return result_df

def write_results_data(outputFolder, enrichment_model_name, result_df_list, enr_cont_name, roc_or_precall):
	ROC_INDEX = 0 
	PRECISION_RECALL_INDEX = 1
	RESULT_DF_HEADER = [['true_pos', 'false_pos'], ['precision', 'recall']]
	enr_cont_name = get_enrichment_analysis_name_from_file_name(enr_cont_name)
	num_row_list = list(map(lambda x: x.shape[0], result_df_list))
	max_num_row = max(num_row_list)
	outF = open(os.path.join(outputFolder, enr_cont_name + DATA_FILE_SUFFIX_LIST[roc_or_precall]), 'w') 
	# each model that we investigate the ROC curves has 32 lists: one list is the true positive and one list is the false positive
	for i, model_name in enumerate(enrichment_model_name):
		this_model_result_df = result_df_list[i]
		this_model_first = this_model_result_df[RESULT_DF_HEADER[roc_or_precall][0]]
		this_model_second = this_model_result_df[RESULT_DF_HEADER[roc_or_precall][1]]
		outF.write(model_name + HEADER_LIST_LIST[roc_or_precall][0])
		for j, rate in enumerate(this_model_first):
			outF.write("," + str(rate))
		num_extra_comma = max_num_row - len(this_model_first)
		outF.write("," * num_extra_comma)
		outF.write("\n")
		outF.write(model_name + HEADER_LIST_LIST[roc_or_precall][1])
		for j, rate in enumerate(this_model_second):
			outF.write("," + str(rate))
		outF.write("," * num_extra_comma)
		outF.write("\n")
	outF.close()
	print ("Done with writing out data for curves for " + \
	enr_cont_name + " after: ")

def calculate_AUC(result_df_list, enr_cont_name, outputFolder, enrichment_model_name):
	# result_df_list is a list of lists: Evey two list  contains information of true positive and false positive for a particular enrichment analysis name. This will be used to plot the ROC curves
	# to calculate AUC using the function auc, we have to provdie false positive data first, and then true positive data
	"""
	enrichment_model_name: list of names of models that we are trying to do enrichment ROC analysis on 
	"""	
	num_models = len(result_df_list)
	assert num_models == len(enrichment_model_name), "Number of models provided by user is not the same as the number of models observed from the result_df_list"
	enr_cont_name = get_enrichment_analysis_name_from_file_name(enr_cont_name)
	auc_enrichment_analysis_fn = os.path.join(outputFolder, "auc_" + enr_cont_name + ".txt")
	aucF = open(auc_enrichment_analysis_fn, 'w')
	auc_list = []
	for model_index, model_name in enumerate(enrichment_model_name): # example: CpG island, exome, etc.
		this_true_pos = (result_df_list[model_index])['true_pos']
		this_false_pos = (result_df_list[model_index])['false_pos']
		try:
			this_analysis_auc = auc(this_false_pos, this_true_pos)
		except:
			print ("FAILED HERE")
			print(model_name)
			print ((np.diff(this_false_pos)))
			print ((np.diff(this_true_pos)))
		auc_list.append(this_analysis_auc)
		aucF.write(model_name + "\t" + str(this_analysis_auc) + "\n")
	aucF.close()
	auc_list = pd.Series(auc_list)
	auc_list.index = enrichment_model_name
	print ("Max AUC: " + enr_cont_name)
	print (auc_list.idxmax())

def get_train_test_one_enr_cont_df (train_test_tuple, enr_cont_name): 
	# input: a tuple of df, for training and testing for one model. We want to extract information related to the enrichment context that we are trying to do calculations on
	# skip the last row because those are lines clarifying percent in genome that the enrichment context occupies
	train_df = ((train_test_tuple[TRAIN_INDEX])[['state', 'percent_in_genome', enr_cont_name]][:-1]).copy()
	train_df = train_df.rename(columns = {enr_cont_name: 'train_' + enr_cont_name, 'percent_in_genome' : "train_percent_in_genome", 'state' : 'train_state'})
	num_train_state = train_df.shape[0] # number of states recorded in training data set
	test_df = ((train_test_tuple[TEST_INDEX])[['state', enr_cont_name, 'percent_in_genome']][:-1]).copy()
	test_df = test_df.rename(columns = {enr_cont_name: 'test_' + enr_cont_name, 'percent_in_genome' : "test_percent_in_genome", 'state' : 'test_state'})
	num_test_state = test_df.shape[0] 
	df = train_df.merge(test_df, how = 'outer', left_on = 'train_state', right_on = 'test_state', suffixes = ('',''))
	values = {'train_percent_in_genome': 0, 'test_percent_in_genome': 0, 'train_' + enr_cont_name: 0, 'test_' + enr_cont_name: 0}
	df = df.fillna(value = values)
	train_cont_percent = (train_test_tuple[TRAIN_INDEX])[enr_cont_name][num_train_state] # percent in genoem that the enr_cont occupies, in the train enrichment data
	test_cont_percent = (train_test_tuple[TEST_INDEX])[enr_cont_name][num_test_state]# percent in genoem that the enr_cont occupies, in the test enrichment data
	return df, (train_cont_percent, test_cont_percent)

def one_enr_cont_analysis(train_test_df_list, outputFolder,  enrichment_model_name, enr_cont_name):
	'''
	train_test_df_list: list of tuple: inside each tuple contains train and test df of enrichment data for one single model
	'''
	result_df_list = [] # each item corresponds to a model that we are measuring true/ false positive rates, precision and recall
	for model_index, model_name in enumerate(enrichment_model_name):
		this_cont_df, (train_cont_percent, test_cont_percent) = get_train_test_one_enr_cont_df(train_test_df_list[model_index], enr_cont_name)
		# this_cont_df has 4 columns: train_percent_in_genome, test_percent_in_genome, train_<cont_name>, test_<cont_name>
		roc_df = do_roc_analysis(this_cont_df, train_cont_percent, test_cont_percent, enr_cont_name, enrichment_model_name) # roc_df has 4 columns: false_pos, true_pos, precision, recall
		# print roc_df.head()
		result_df_list.append(roc_df)
	write_results_data(outputFolder, enrichment_model_name, result_df_list, enr_cont_name, roc_or_precall = ROC_INDEX)
	write_results_data(outputFolder, enrichment_model_name, result_df_list, enr_cont_name, roc_or_precall = PRECISION_RECALL_INDEX)
	calculate_AUC(result_df_list, enr_cont_name, outputFolder, enrichment_model_name)



def create_roc_analysis(train_test_df_list, outputFolder, enrichment_model_name):
	"""
	enrichment_model_name: list of names of models that we are trying to do enrichment ROC analysis on 
	train_test_df_list: list of tuple: inside each tuple contains train and test df of enrichment data for one single model
	"""
	headers = list((train_test_df_list[0][0]).columns) # now that all headers of all files are the same, we take the first file's headers
	enr_cont_list = headers[2:]
	for enr_cont_name in enr_cont_list: # enrichment analysis name is the name of genomic content that we want to do analysis of enrichment. Skip the first two because it's state and percent_in_genome --> irrelevant
		one_enr_cont_analysis(train_test_df_list, outputFolder,  enrichment_model_name, enr_cont_name)
		print ("Done doing analysis for enrichment context: " + enr_cont_name)


def main():
	minimum_arvg = 3
	if len(sys.argv) < minimum_arvg:
		print ("Wrong input man!")
		print ("python create_roc_curve.py \
		<number of enrichment files to process> \
		<output folder> \
		<list of enrichment files> --> kind of optional")
		exit(1)
	else:
		num_enrich = int(sys.argv[1]) # number of models that we have enrichment data for
		outputFolder = sys.argv[2]
		helper.make_dir(outputFolder)
		if len(sys.argv) != minimum_arvg + (num_enrich) * 2:
			print ("The number of arguments are not as exact as expected according to the number of models that we are trying to compare")
			print (sys.argv)
			print ("exiting ...")
			exit(1)
		train_test_enrich_folder_list = sys.argv[minimum_arvg: (minimum_arvg + num_enrich)]
		enrichment_model_name = sys.argv[(minimum_arvg + num_enrich):]
		# Get data and do some quality control of the files
		train_test_df_list = get_enrichment_data(train_test_enrich_folder_list) # a list of tuples: each item in the list is a tuple of train and test dataframes
		print ("Get all data for enrichment after: " )
		create_roc_analysis(train_test_df_list, outputFolder, enrichment_model_name)
		print ("Done!")

def usage():
	print ("python calculate_roc_train_test_enrichment.py")
	print ("num_enrich: number of segmentations that we have enrihcment data for, i.e. number of models that we are trying to calculate ROC/ precall for")
	print ("outputFolder: where roc/precall data for each enrichment context are stored in separate files")
	print ("train_test_enrich_folder_list: each item corresponds to a folder containing files train_overlap.txt and test_overlap.txt from one model")
	print ("enrichment_model_name: list of model names. It should correspond to the order of items  in train_test_enrich_folder_list")
	exit(1)
main()
