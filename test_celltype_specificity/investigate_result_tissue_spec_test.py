import pandas as pd
import sys
import os
import time 
import handle_emission_helper as em
start_time = time.clock()
THRESHOLD = 0.00000352112# 4.16666666667e-07
def color_significant_tests(value):
    if value < THRESHOLD:
        color = '#F1826A'
    else: 
        color  = '#F3E4E1'
    return 'background-color: %s' % color

def color_excel(result_fn, output_fn):
	df = pd.read_csv(result_fn, sep = '\t', header = 0, index_col = 0)
	df = df.style.applymap(color_significant_tests)
	df.to_excel(output_fn, engine = 'openpyxl')

def main():
	if len(sys.argv) != 4: 
		usage()
	result_fn = sys.argv[1]
	if not os.path.join(result_fn):
		print "result_fn is NOT VALID"
		usage()
	output_fn = sys.argv[2]
	em.create_folder_for_file(output_fn)
	print "Done after: " + str(time.clock() - start_time)

def summarize_one_row_one_chrom_mark_result_df (row_data):
	# row_data: pandas series corresponding to one state where indices of the series are the different cell groups or different anatomy
	significant_groups = row_data.index[row_data < THRESHOLD] # get the names of cell groups that have a p-values lower than the threshold
	return ','.join(significant_groups) # each siginificant cell group should be reported separated by ','

def summarize_one_chrom_mark_result_df(df):
	# df: one of the itme in group_test_result_df_list corresponding to results of testing cell type specificity for one chromatin mark
	# indices of df: 'S1' --> 'S100'
	# colnames of df: different cell groups or different anatomy
	return df.apply(summarize_one_row_one_chrom_mark_result_df, axis = 1)

def usage():
	print "python investigate_result_tissue_spec_test.py"
	print "result_fn: result of test_tissue_specificity_emission.R"
	print "output_fn: output with styled excel file"
	exit(1)

result_folder = "./output/" # TO BE UPDATED: change this folder to where you got the results of testing from running test_tissue_specificity_emission.R
CHROM_MARK_TO_TEST = ['H3K9me3', 'H3K4me1', 'H3K4me3', 'H3K36me3', 'H3K27me3', 'H3K27ac', 'H3K9ac', 'DNase']
num_states = 100
group_test_result_fn_list = map(lambda x: os.path.join(result_folder, x + '_group_mann_whitney_test.txt'), CHROM_MARK_TO_TEST)
print "Done getting mark names files"
group_test_result_df_list = map(lambda x: pd.read_csv(x, sep = '\t', header = 0, index_col = 0), group_test_result_fn_list)
print "Done getting mark names dataframes"
summary_df = pd.DataFrame(columns = CHROM_MARK_TO_TEST)
for mark_index, mark_name in enumerate(CHROM_MARK_TO_TEST):
	mark_df = group_test_result_df_list[mark_index]
	summary_df[mark_name] = summarize_one_chrom_mark_result_df(mark_df)
	print "Done with mark " + mark_name
summary_df.index = map(lambda x: 'S' + str(x+1), range(num_states))
summary_fn = os.path.join(result_folder, "summary_cell_group.xlsx")
summary_df.to_excel(summary_fn	)	
