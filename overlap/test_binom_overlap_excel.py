import pandas as pd 
import numpy as np 
import os
import sys
import helper
from scipy.stats import binom_test
import seaborn as sns
full_segment_fn = '/u/home/h/havu73/project-ernst/ROADMAP_aligned_reads/chromHMM_model/model_100_state/hg19_segmentations/genome_100_segments.bed.gz'
ALPHA = 0.05
def color_significant_pval(pval, threshold):
	if pval <= threshold:
		color = '#85BCE5' # light blue
	else:
		color = '#ffffff' # white
	return 'background-color: %s' % color

def get_enrichment_df(enrichment_fn):
    enrichment_df = pd.read_csv(enrichment_fn, sep = "\t", header = 0)
    enrichment_df = enrichment_df.rename(columns = {u'state (Emission order)': 'state', u'Genome %' : 'percent_in_genome'})
    # now we will fix the enrichment context names to get rid of the tails of the files: AT128_cred_snps.bed.gz --> AT128
    enr_cont_list = enrichment_df.columns[2:]
    enrichment_df.columns = list(enrichment_df.columns[:2]) + list(enr_cont_list)
    (num_state, num_enr_cont) = (enrichment_df.shape[0] - 1, enrichment_df.shape[1] - 2) 
    percent_genome_of_cont = np.array(enrichment_df.iloc[num_state, 2:])
    enrichment_df = enrichment_df.loc[:(num_state - 1)]
    # num_enr_cont : number of enrichment context (TSS, exon, etc.)
    enrichment_df['frac_gene_in_state'] = enrichment_df.percent_in_genome / 100.0 # the fraction of the genome that is in the state
    return enrichment_df, num_state, num_enr_cont, percent_genome_of_cont

def calculate_p_value_one_enr_cont(df, percent_gen_this_cont, total_bins):
	# df should have three columns: state, frac_gene_in_state, <enr_cont>
	# percent_gen_this_cont : percent of the genome in this enrichment context
	this_enr_cont = df.columns[2] # the third column name is the names of the enrichment context
	num_bins_cont = total_bins * percent_gen_this_cont * (200 / 100.0) # 200: number of bases in each bin 
	df['frac_cont_in_state'] = df[this_enr_cont] * df['frac_gene_in_state'] # fraction of the genomic context that is in the state
	df['obs_num_bins_in_state_cont'] = df['frac_cont_in_state'] * num_bins_cont # number of bins that are in the state and also containing the enrichment context of interest
	df[this_enr_cont + '_binomP'] = df.apply(lambda x: binom_test(x['obs_num_bins_in_state_cont'], num_bins_cont, x['frac_gene_in_state'], alternative = 'greater'), axis = 1)
	return df[this_enr_cont + '_binomP']

def get_total_num_segments (segment_fn):
	segment_df = pd.read_csv(segment_fn, header = None, sep = '\t')
	segment_df.columns = ['chrom', 'start', 'end', 'state']
	segment_df['num_bins'] = (segment_df.end - segment_df.start) / 200
	total_bins = np.sum(segment_df.num_bins)
	return total_bins

def paint_excel_result(enrichment_df, result_df, output_fn, enr_cont_list, percent_genome_of_cont, num_state, num_enr_cont):
	num_tests = num_state * num_enr_cont # total number of tests that we conduct here
	pval_threshold = ALPHA / float(num_tests)
	print ("Bornferonni corrected p-val theshold: " + str(pval_threshold))
	pval_colnames = list(map(lambda x: x + "_binomP", enr_cont_list))
	pval_colored_df = result_df.style.applymap(lambda x: color_significant_pval(x, pval_threshold), subset = pd.IndexSlice[:, pval_colnames]) # color only the results, which excludes the state column
	enrichment_df = enrichment_df[['state', 'frac_gene_in_state'] + list(enr_cont_list)]
	enrichment_df.loc[num_state] = ['percent_genome_of_cont', 100] + list(percent_genome_of_cont)
	cm = sns.light_palette("red", as_cmap=True)
	enr_colored_df = enrichment_df.style.background_gradient(subset = pd.IndexSlice[:(num_state - 1), enr_cont_list], cmap = cm)
	writer = pd.ExcelWriter(output_fn, engine = 'xlsxwriter')
	pval_colored_df.to_excel(writer, sheet_name = 'binom_pval')
	enr_colored_df.to_excel(writer, sheet_name = 'enrichment')
	writer.save()

def calculate_pval_enrichment(total_bins, overlap_fn, output_fn):
	enrichment_df, num_state, num_enr_cont, percent_genome_of_cont = get_enrichment_df(overlap_fn)
	enr_cont_list = enrichment_df.columns[2:(2 + num_enr_cont)]
	result_df = pd.DataFrame()
	result_df['state'] = enrichment_df['state']
	for cont_index, enr_cont in enumerate(enr_cont_list):
		this_cont_df = (enrichment_df[['state', 'frac_gene_in_state', enr_cont]]).copy()
		result_df[enr_cont + "_binomP"] = calculate_p_value_one_enr_cont(this_cont_df, percent_genome_of_cont[cont_index], total_bins)
	print ("Done calculating p-values")
	paint_excel_result(enrichment_df, result_df, output_fn, enr_cont_list, percent_genome_of_cont, num_state, num_enr_cont)
	print ("Done writing to " + output_fn)

def main():
	if len(sys.argv) != 4:
		usage()
	segment_fn = sys.argv[1]
	helper.check_file_exist(segment_fn)
	overlap_fn = sys.argv[2]
	helper.check_file_exist(overlap_fn)
	output_fn = sys.argv[3]
	helper.create_folder_for_file(output_fn)
	print ("Done getting command line arguments")
	total_bins = get_total_num_segments(segment_fn) #15478457
	print ("Done getting the total number of bins: " + str(total_bins))
	calculate_pval_enrichment(total_bins, overlap_fn, output_fn)
	print ("Done!")

def usage():
	print ("python test_overlap_pval_causaldb.py")
	print ("segment_fn: the file we used to create the overlap enrichments that we will calculate p-values for this analysis")
	print ("overlap_fn: results from ChromHMM OvelapEnrichment")
	print ("output_fn: painted excel file with p p-values")
	exit(1)

main()