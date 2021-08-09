# this file is created to produce a summary of the different cell types that may outperform full-stack state annotations in terms of recovering known genomic contexts
import sys
import os
import numpy as np
from collections import Counter
import pandas as pd
import helper
import glob

def read_one_auc_fn (fn):
	df = pd.read_csv(fn, header = None, sep = '\t', index_col = 0, squeeze = True)
	full = df['full']
	max_ct = np.max(df[1:]) # avoid the first because the first one is full-stack's auc
	min_ct = np.min(df[1:])
	median_ct = np.median(df[1:])
	ct_greater_than_full = list(df[1:].index[df[1:] > df['full']])
	ct_greater_than_full = ','.join(ct_greater_than_full)
	return full, max_ct, min_ct, median_ct, ct_greater_than_full

def process_auc_data_both_group(ct_100_folder, ct_18_folder, writer, sheet_name):
	# writer should be an excel object
	# first, check the integrity of the data
	auc18_fn_list = glob.glob(ct_18_folder + '/auc_*.txt')
	auc18_context_list = list(map(lambda x: '_'.join(x.split('/')[-1].split('.')[0].split('_')[1:]), auc18_fn_list)) # /path/to/18_state/auc_<context>.txt --> context
	auc100_fn_list = glob.glob(ct_100_folder + '/auc_*.txt')
	auc100_context_list = list(map(lambda x: '_'.join(x.split('/')[-1].split('.')[0].split('_')[1:]), auc100_fn_list)) # /path/to/18_state/auc_<context>.txt --> context
	assert set(auc18_context_list) == set(auc100_context_list), "the 100-state and 18-state contexts are not the same"
	result_df = pd.DataFrame(columns = ['full', 'context','max_18state', 'min_18state', 'median_18state', 'ct_18_better_than_full', 'max_100state', 'min_100state', 'median_100state', 'ct_100_better_than_full'])
	for context in auc18_context_list:
		s100_fn = os.path.join(ct_100_folder, 'auc_' + context + '.txt')
		s18_fn = os.path.join(ct_18_folder, 'auc_' + context + '.txt')
		full, max_ct100, min_ct100, median_ct100, ct100_greater_than_full = read_one_auc_fn(s100_fn)
		full, max_ct18, min_ct18, median_ct18, ct18_greater_than_full = read_one_auc_fn(s18_fn)
		result_df.loc[result_df.shape[0]] = [full, context, max_ct18, min_ct18, median_ct18, ct18_greater_than_full, max_ct100, min_ct100, median_ct100, ct100_greater_than_full]
	result_df.to_excel(writer, sheet_name = sheet_name)
	return result_df

output_excel = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/figures/supp_excel/summary_auc_compare_models_code_generated.xlsx'
writer = pd.ExcelWriter(output_excel, engine='xlsxwriter')
# general genome contexts
gc_100s_folder = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/overlap_enrichment/genome_context/100_state/compare_models'
gc_18s_folder = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/overlap_enrichment/genome_context/18_state/compare_models'
process_auc_data_both_group(gc_100s_folder, gc_18s_folder, writer, 'genome_context_fig3A')
# structural variants
sv_100s_folder = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/overlap_enrichment/structural_variant_hall_lab/ct_100_state/compare_models'
sv_18s_folder = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/overlap_enrichment/structural_variant_hall_lab/ct_18_state/compare_models'
process_auc_data_both_group(sv_100s_folder, sv_18s_folder, writer, 'structural_variants')
# repeat classes
repeats_100s_folder = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/overlap_enrichment/ucsc/repeats_classes_families/100_state/compare_models/all_classes/'
repeats_18s_folder = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/overlap_enrichment/ucsc/repeats_classes_families/18_state/compare_models/all_classes/'
process_auc_data_both_group(repeats_100s_folder, repeats_18s_folder, writer, 'repeat_classes')
# causaldb
causaldb_100s_folder = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/overlap_enrichment/causaldb/100_state/against_bg/compare_models'
causaldb_18s_folder = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/overlap_enrichment/causaldb/18_state/against_bg/compare_models'
process_auc_data_both_group(causaldb_100s_folder, causaldb_18s_folder, writer, 'causaldb')
# top 1% prioritized non-coding variants
vp_100s_folder = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/variant_prioritization_analysis/non_coding_0.01/100_state/against_bg/compare_models'
vp_18s_folder = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/variant_prioritization_analysis/non_coding_0.01/18_state/against_bg/compare_models'
process_auc_data_both_group(vp_100s_folder, vp_18s_folder, writer, 'prioritized_variants')
writer.save()
def main():
	if len(sys.argv) != 5:
		usage()
	ct_100_folder = sys.argv[1]
	helper.check_dir_exist(ct_100_folder)
	ct_18_folder = sys.argv[2]
	helper.check_dir_exist(ct_18_folder)
	output_excel = sys.argv[3]
	helper.create_folder_for_file(output_excel)
	sheet_name = sys.argv[4]
	print ("Done getting command line argument")

def usage():
	print ("python summarize_auc_overlap_for_publication.py")
	print ("ct_100_folder: where the compare_models folder for 100-state files are found. It should contain the auc_<genome_contex>.txt")
	print ("ct_18_folder: where the compare_models folder for 18-state files are found. It should contain the auc_<genome_contex>.txt")
	print ("output_excel: where we store the data")
	print ("sheet_name: so that we can save the excel")
	exit(1)