import sys
import os
import numpy as np
import time
from collections import Counter
import pandas as pd
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import matplotlib.legend_handler as handler
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.pyplot import cm
import math
start_time = time.clock()
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

# color_list = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
def plot_results_curve(coordinate_data, outputFolder, enrichment_model_name, enrichment_analysis_name, roc_or_precall):
	assert len(coordinate_data) == len(enrichment_model_name) * 2, "The number of coordidate_data list is not consistent with the number of model we are trying to get ROC curve for: " + str(len(coordinate_data)) + "    " + str(len(enrichment_model_name))

	num_models = len(enrichment_model_name) # number of models that we have to plot
	color_list = cm.rainbow(np.linspace(0,1,num_models))
	fig = plt.figure(figsize = (6,6), dpi = 300)
	ax = fig.add_subplot(111)
	for i in range(len(enrichment_model_name)): # x: false positive, y: true positive
	# x: recall, y: precision
	# in both cases (roc curves or precision recal curves), the x axis always correspond to the list that is created later (aka: false positive or recall rates)
		ax.plot(coordinate_data[i * 2 + 1], coordinate_data[2 * i ], color = color_list[i])

	# figure legend
	legend_shapes = []
	legend_labels = []
	for i, model_name in enumerate(enrichment_model_name):
		legend_shapes.append(mpatches.Patch(color=color_list[i]))
		legend_labels.append(r''+model_name)
	num_legend_columns = int(math.ceil(float(num_models) / 4.0))
	ax.legend(legend_shapes, legend_labels, loc=9, bbox_to_anchor=(0.5 - (num_legend_columns - 1) * 0.15 ,1.2), ncol= num_legend_columns, numpoints=1, fontsize=8)

	# axis lable
	ax.set_xlabel(AXIS_NAME_LIST_LIST[roc_or_precall][X_AXIS_INDEX], fontsize=8)
	ax.set_ylabel(AXIS_NAME_LIST_LIST[roc_or_precall][Y_AXIS_INDEX], fontsize=8,)

	# figure title
	fig.suptitle(enrichment_analysis_name, fontsize =12)

	# adjust subplot
	fig.subplots_adjust(left=0.1, right =0.95, bottom = 0.1, top=0.8) 

	fig.savefig(os.path.join(outputFolder, enrichment_analysis_name + PLOT_FILE_SUFFIX_LIST[roc_or_precall]), format = "png", dpi = 300)


def get_enrichment_data(enrichment_fn_list):
	data_frame_list = [pd.read_csv(f, sep = "\t") for f in enrichment_fn_list]
	for df_index, df in enumerate(data_frame_list):
		try:
			data_frame_list[df_index] = df.rename(columns={"state (Emission order)": "state", "Genome %": "percent_in_genome", "state (User order)" : "state"})
		except:
			print "Could not convert data columns headers for data frame: " + str(df_index)
	headers_list = [list(x) for x in data_frame_list]
	# make sure that all lenth of each of the header_list is the same
	if len(set(map(lambda x: len(x) , headers_list))) != 1: # check that all enrichment files have the same headers
		print "The length of the headers_list in all enrichment files is not the same. Exiting..."
		exit(1)
	for i in range(len(headers_list[0])): # check that all headers of all the enrichment files are the same
		if len(set(map(lambda x: x[i], headers_list))) != 1:
			print "Header indexed: " + str(i) + "  of the headers_list is not consistent"
			print map(lambda x: x[i], headers_list)
			print "exiting ...."
			exit(1)

	if headers_list[0][0] != "state": # check that the first header is state, in all enrichment files
		print "The first headers is not state: " + headers_list[0][0]
		print "exiting ..."
		exit(1)
	if headers_list[0][1] != "percent_in_genome": # check that the second header is percent_in_genome, in all enrichment files
		print "The second headers is not percent_in_genome: " + headers_list[0][0]
		print "exiting ..."
		exit(1)
	return data_frame_list

def do_precision_recall_analysis(data_frame, annot_percent, enrichment_analysis_name):
	# sort states based on decreasing enrichment of the enrichment_analysis_name
	data_frame.sort_values(by=enrichment_analysis_name, ascending = False, inplace = True)
	num_rows, num_cols = data_frame.shape # ncols should be 2: percent_in_genome, genomic context of interest
	# Get the shape of the data frame

	# 1. Get the fraction of the genome that is in each state. The statistics "percent_in_genome" that we get from enrichment analysis output from ChromHMM is in the unit of percent. Therefore, we will have to divide that by 100 to get the fraction of genes that are in each state.
	data_frame['frac_gene_in_state'] = data_frame["percent_in_genome"] / 100.0

	# 2. Calculate the fraction of each state that the mark is present. frac_state_true_mark = (# SM / #State) = (#mark / # genome_pos) * ((# SM / #State) / (#mark / # genome_pos)) = frac_gene_mark_present * fold_enrichment = annot_percent / 100 * fold_enrichment
	frac_gene_mark_present = annot_percent / 100.0 
	data_frame['frac_state_true_mark'] = frac_gene_mark_present * data_frame[enrichment_analysis_name]

	# 3. Calculate the fraction of gene that is in the state and where the mark is present 
	# frac_gene_in_state_mark = (# SM / # genome_pos) = (# SM / # state) * (# state / # genome_pos) = frac_state_true_mark * frac_gene_in_state
	data_frame['frac_gene_in_state_mark'] = data_frame['frac_state_true_mark'] * data_frame['frac_gene_in_state']

	# 4. Calculate the culmulative sum of fraction of genes that are in states and also contains the mark
	data_frame['culm_frac_gene_in_state_mark'] = (data_frame['frac_gene_in_state_mark']).cumsum()

	# 5. Calculate the culmulative fractive of genes that are in states
	data_frame['culm_frac_gene_in_state'] = (data_frame['frac_gene_in_state']).cumsum()

	# 6. Precision = True positive / claimed positive 
	# = # state mark / # genome pos in state
	# = (# SM / # genome_pos) / (# genome_pos in state / # genome_pos)
	# = (frac_gene_in_state_mark) / frac_gene_in_state
	# we have to do it in a culmulative way to aggregate states' power
	data_frame['precision'] = data_frame['culm_frac_gene_in_state_mark'] / data_frame['culm_frac_gene_in_state']

	# 7. Recall = True positive / actual positives
	# = # state mark / # mark 
	# = (# state mark / # genome_pos) / (# mark / # genome_pos)
	# = frac_gene_in_state_mark / frac_gene_mark_present
	# but we have to do it in culmulative way
	data_frame['recall'] = data_frame['culm_frac_gene_in_state_mark'] / frac_gene_mark_present

	# we will return only precision and recall data vectors
	return (data_frame['precision'], data_frame['recall'])

def do_roc_analysis(this_annot_data_frame, annot_percent, enrichment_analysis_name):
	"""
	consider 'gContext' as the analysis that we are talking about
	"""
	# filter out state that is non existent on the genome
	data_frame = (this_annot_data_frame.loc[this_annot_data_frame['percent_in_genome'] > 0]).copy()
	# sort states based on decreasing enrichment of the enrichment_analysis_name
	data_frame.sort_values(by=enrichment_analysis_name, ascending = False, inplace = True)
	num_rows, num_cols = data_frame.shape # ncols should be 2: percent_in_genome, genomic context of interest
	# Get the shape of the data frame
	# calculate true positive rates
	data_frame['per_gContext_in_state'] = data_frame["percent_in_genome"] * data_frame[enrichment_analysis_name] / 100.0 # percentage of the places that are of gContext to be actually of each state
	true_pos = [0]
	current_cumulative_true_pos = 0
	for per_gContext_in_state in data_frame['per_gContext_in_state']:
		current_cumulative_true_pos += per_gContext_in_state
		true_pos.append(current_cumulative_true_pos)
	# calculate false_positive rates
	false_pos = [0]
	data_frame['percent_state_be_gContext'] = data_frame[enrichment_analysis_name] * annot_percent
	data_frame['percent_genome_in_state_not_gContext'] = data_frame["percent_in_genome"] * \
	(100 - data_frame['percent_state_be_gContext']) / 100
	percent_genome_not_gContext = data_frame['percent_genome_in_state_not_gContext'].sum()
	current_culmulative_false_pos = 0
	for percent_genome_in_state_not_gContext in data_frame['percent_genome_in_state_not_gContext']:
		current_culmulative_false_pos += (percent_genome_in_state_not_gContext / percent_genome_not_gContext)
		false_pos.append(current_culmulative_false_pos)
	return true_pos, false_pos


def write_results_data(outputFolder, enrichment_model_name, coordinate_data, enrichment_analysis_name, roc_or_precall):
	length_coord_data = map(lambda x: len(x), coordinate_data)
	max_row_length = max(length_coord_data)
	outF = open(os.path.join(outputFolder, enrichment_analysis_name + DATA_FILE_SUFFIX_LIST[roc_or_precall]), 'w') 
	# each model that we investigate the ROC curves has 32 lists: one list is the true positive and one list is the false positive
	for i, model_name in enumerate(enrichment_model_name):
		this_model_first = coordinate_data[i * 2]
		this_model_second = coordinate_data[i * 2 + 1]
		outF.write(model_name + HEADER_LIST_LIST[roc_or_precall][0])
		for j, rate in enumerate(this_model_first):
			outF.write("," + str(rate))
		num_extra_comma = max_row_length - len(this_model_first)
		outF.write("," * num_extra_comma)
		outF.write("\n")
		outF.write(model_name + HEADER_LIST_LIST[roc_or_precall][1])
		for j, rate in enumerate(this_model_second):
			outF.write("," + str(rate))
		outF.write("," * num_extra_comma)
		outF.write("\n")
	outF.close()
	print "Done with writing out data for curves for " + \
	enrichment_analysis_name + " after: " + str(time.clock() - start_time)

def calculate_AUC(roc_coordinate_data, enrichment_analysis_name, outputFolder, enrichment_model_name):
	# roc_coordinate_data is a list of lists: Evey two list  contains information of true positive and false positive for a particular enrichment analysis name. This will be used to plot the ROC curves
	# to calculate AUC using the function auc, we have to provdie false positive data first, and then true positive data
	"""
	enrichment_model_name: list of names of models that we are trying to do enrichment ROC analysis on 
	"""	
	num_models = len(roc_coordinate_data) / 2.0
	assert num_models == len(enrichment_model_name), "Number of models provided by user is not the same as the number of models observed from the roc_coordinate_data"
	enrichment_analysis_name = get_enrichment_analysis_name_from_file_name(enrichment_analysis_name)
	auc_enrichment_analysis_fn = os.path.join(outputFolder, "auc_" + enrichment_analysis_name + ".txt")
	aucF = open(auc_enrichment_analysis_fn, 'w')
	for model_index, model_name in enumerate(enrichment_model_name): # example: CpG island, exome, etc.
		this_true_pos = roc_coordinate_data[model_index * 2]
		this_false_pos = roc_coordinate_data[model_index * 2 + 1]
		this_analysis_auc = auc(this_false_pos, this_true_pos, reorder = True)
		aucF.write(model_name + "\t" + str(this_analysis_auc) + "\n")
	aucF.close()


def create_roc_analysis(data_frame_list, outputFolder, enrichment_model_name, plot_or_not):
	"""
	enrichment_model_name: list of names of models that we are trying to do enrichment ROC analysis on 
	"""
	num_enrich = len(data_frame_list) # number of enrichment files that we got
	headers = list(data_frame_list[0]) # now that all headers of all files are the same, we take the first file's headers
	num_rows_list = [data_frame_list[i].shape[0] for i in range(num_enrich)] # [model index]
	num_cols_list = [data_frame_list[i].shape[1] for i in range(num_enrich)] # [model index]
	for enrichment_analysis_name in headers[2:]: # enrichment analysis name is the name of genomic content that we want to do analysis of enrichment. Skip the first two because it's state and percent_in_genome --> irrelevant
		roc_coordinate_data = [] # lists of lists: Evey two list  contains information of true positive and false positive for a particular enrichment analysis name. This will be used to plot the ROC curves

		roc_curve_data = [(x[["percent_in_genome", enrichment_analysis_name]][:-1]) for x in data_frame_list] #skip the last row b ecause is the base row. [model_index]--> tables with 2 columns [percent_in_genome, genomic context of interest][states]
		annot_percent_list = [data_frame_list[i].loc[num_rows_list[i] - 1, enrichment_analysis_name] for i in range(num_enrich)]
		# the last row contains information about the percentage of the genome that is in the enrichment_analysis_name: ex: phastCons, TSS, etc. Get that data for each of the model that we want to plot roc curves/ precision- recall curves for.
		for model_index in range(num_enrich):
			true_pos, false_pos = do_roc_analysis((roc_curve_data[model_index]).copy(),		(annot_percent_list[model_index]), enrichment_analysis_name)
			# append the two lists to the list of list for later storage and visualization
			roc_coordinate_data.append(true_pos)
			roc_coordinate_data.append(false_pos)
		# Now, calculate and write the AUC for each of the enrichment analysis
		calculate_AUC(roc_coordinate_data, enrichment_analysis_name, outputFolder, enrichment_model_name)
		write_results_data(outputFolder, enrichment_model_name, roc_coordinate_data, enrichment_analysis_name, ROC_INDEX)
		if (plot_or_not):
			plot_results_curve(roc_coordinate_data, outputFolder, enrichment_model_name, enrichment_analysis_name, ROC_INDEX)
			print "Done with creating plot for ROC curver for " + \
			enrichment_analysis_name + " after: " + str(time.clock() - start_time)
		

def create_precision_recall_analysis(data_frame_list, outputFolder, enrichment_model_name, plot_or_not):
	"""
	enrichment_model_name: list of names of models that we are trying to do enrichment ROC analysis on 
	"""
	num_enrich = len(data_frame_list) # number of enrichment files that we got
	headers = list(data_frame_list[0]) # now that all headers of all files are the same, we take the first file's headers
	num_rows_list = [data_frame_list[i].shape[0] for i in range(num_enrich)] # [model index]
	num_cols_list = [data_frame_list[i].shape[1] for i in range(num_enrich)] # [model index]

	for enrichment_analysis_name in headers[2:]: # enrichment analysis name is the name of genomic content that we want to do analysis of enrichment. Skip the first two because it's state and percent_in_genome --> irrelevant
		prec_recall_coordinate_data = [] # lists of lists: Evey two list  contains information of precision and recall for a particular enrichment analysis name (TSS, phastCons, etc.). This will be used to plot the ROC curves
		
		enrichment_data = [(x[["percent_in_genome", enrichment_analysis_name]][:-1]) for x in data_frame_list] #skip the last row b ecause is the base row. [model_index]--> tables with 2 columns [percent_in_genome, genomic context of interest][states]
		annot_percent_list = [data_frame_list[i].loc[num_rows_list[i] - 1, enrichment_analysis_name] for i in range(num_enrich)]
		# the last row contains information about the percentage of the genome that is in the enrichment_analysis_name: ex: phastCons, TSS, etc. Get that data for each of the model that we want to plot precision- recall curves for.
		for model_index in range(num_enrich):
			# do analysis to calculate precision recall
			precision_list, recall_list = do_precision_recall_analysis((enrichment_data[model_index]).copy(), (annot_percent_list[model_index]), enrichment_analysis_name)
			# store the lists into list of lists so that we can to storage and visualization later
			prec_recall_coordinate_data.append(precision_list)
			prec_recall_coordinate_data.append(recall_list)
		# write data into a table. We cannot save data using the standard methods for saving data into csv files from pandas, because the length of lists vary depends on the number of states that each model has
		write_results_data(outputFolder, enrichment_model_name, prec_recall_coordinate_data, enrichment_analysis_name, PRECISION_RECALL_INDEX)
		if (plot_or_not):
			# plot precision recall curve
			plot_results_curve(prec_recall_coordinate_data, outputFolder, enrichment_model_name, enrichment_analysis_name, PRECISION_RECALL_INDEX)
			print "Done with creating plot for precision_recall curves for " + \
			enrichment_analysis_name + " after: " + str(time.clock() - start_time)

def main():
	minimum_arvg = 4
	if len(sys.argv) < minimum_arvg:
		print "Wrong input man!"
		print "python create_roc_curve.py \
		<number of enrichment files to process> \
		<output folder> \
		<list of enrichment files> --> kind of optional"
		exit(1)
	else:
		num_enrich = int(sys.argv[1])
		outputFolder = sys.argv[2]
		plot_or_not = int(sys.argv[3])
		if (plot_or_not not in [0,1]):
			print "plot_or_not can be NO (0) or YES (1)"
			exit(1)
		plot_or_not = [False, True][plot_or_not] # true: plot, false: don't plot
		try: 
			os.makedirs(outputFolder)
		except:
			pass
		if len(sys.argv) != minimum_arvg + (num_enrich) * 2:
			print "wrong input man!"
			print sys.argv
			print "exiting ..."
			exit(1)
		enrichment_fn_list = sys.argv[minimum_arvg: (minimum_arvg + num_enrich)]
		enrichment_model_name = sys.argv[(minimum_arvg + num_enrich):]
		# Get data and do some quality control of the files
		data_frame_list = get_enrichment_data(enrichment_fn_list)
		print "Get all data for enrichment after: " + str(time.clock() - start_time)
		create_roc_analysis(data_frame_list, outputFolder, enrichment_model_name, plot_or_not)
		print "Done with ROC curves!"
		# create_precision_recall_analysis(data_frame_list, outputFolder, enrichment_model_name, plot_or_not)
		print "Done precision recall curves!"

main()
