import pandas as pd 
import numpy as np 
import os
import sys
import analysis_helper as helper
from scipy.stats import mannwhitneyu
emission_fn = "/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/emissions_100.txt"
meta_fn = '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/ROADMAP_metadata_july2013.csv'
output_folder = "/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/chromHMM_model/model_100_state/chrom_mark_spec_test/"
ALPHA = 0.01

def color_significant_pval(pval_row, threshold):
	# pval_row: a rows of p-values
	white_color = '#ffffff' # white
	blue_color = '#85BCE5' # light blue
	red_color = '#FF7F7F' # light red
	results = pd.Series(['background-color: %s' % white_color for x in pval_row])
	results.index = pval_row.index
	# change colors to blue if below the thresholds
	below_threshold_indices = (pval_row <= threshold)
	results[below_threshold_indices] = 'background-color: %s' % blue_color
	results[pval_row ==  pval_row.min()] = 'background-color: %s' % red_color
	return results

def get_emission_matrix_df (emission_fn, meta_fn, num_state):
	emission_df = pd.read_csv(emission_fn, header = 0, index_col = 0, sep = '\t')
	emission_df = emission_df.transpose()
	emission_df.reset_index(inplace = True)
	emission_df.columns = ['experiment'] + map(lambda x: 'S' + str(x+1), range(num_state))
	emission_df['chrom_mark'] = emission_df['experiment'].apply(lambda x: x.split('-')[1])
	emission_df['ct'] = emission_df['experiment'].apply(lambda x: x.split('-')[0])
	meta_df = pd.read_csv(meta_fn, header = 0, sep = ',')
	meta_df = meta_df.rename(columns = {'Epigenome ID (EID)' : 'ct'})
	meta_df = meta_df[['ct', 'GROUP', 'ANATOMY']]
	emission_df = pd.merge(emission_df, meta_df, how = 'left', left_on = 'ct', right_on = 'ct')
	return emission_df

def test_chrom_mark_specificity (emission_df, num_state):
	total_number_test = 0
	count_mark = emission_df.chrom_mark.value_counts()
	marks_to_test = count_mark.index[count_mark > 15]
	result_df = pd.DataFrame(columns = ['state'] + list(marks_to_test))
	# result_df: columns: marks that we want to test the significance of higher emission probabilities, rows: different states
	result_df['state'] = map(lambda x: 'S' + str(x+1), range(num_state))
	num_tests = num_state * len(marks_to_test)
	for chromM in marks_to_test:
		this_chromM_df = emission_df[emission_df['chrom_mark'] == chromM]
		other_chromM_df = emission_df[emission_df['chrom_mark'] != chromM]
		this_mark_results = []
		for state_index in range(num_state):
			x = this_chromM_df['S' + str(state_index + 1)]
			y = other_chromM_df['S' + str(state_index + 1)]
			t = mannwhitneyu(x, y, use_continuity = False, alternative = 'greater')
			this_mark_results.append(t.pvalue)
		result_df[chromM] = this_mark_results
	return result_df, num_tests

def paint_result_excel(result_df, num_tests, output_fn):
	threshold = ALPHA / float(num_tests)
	colored_df = result_df.style.apply(lambda x: color_significant_pval(x, threshold), axis = 1, subset = pd.IndexSlice[:, result_df.columns[1:]]) #exclude coloring the first column which is state annotation
	writer = pd.ExcelWriter(output_fn, engine = 'xlsxwriter')
	colored_df.to_excel(writer, sheet_name = 'mark_specificity')
	writer.save() 

def main():
	num_state = 100
	emission_df = get_emission_matrix_df(emission_fn, meta_fn, num_state)
	print "Done getting emission_df"
	result_df, num_tests = test_chrom_mark_specificity(emission_df, num_state)
	print "Done getting result_df"
	output_fn = os.path.join(output_folder, 'test_chrom_mark_specificity.xlsx')
	paint_result_excel(result_df, num_tests, output_fn)
	print "Done paint_result_excel!"
main()