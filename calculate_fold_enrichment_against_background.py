import pandas as pd
import numpy as np
import sys
import os
import chromHMM_utilities_common_functions_helper as cmh

def get_rid_of_stupid_file_tail(context_name):
    if context_name.endswith('.bed.gz'):
        return(context_name[:-7])
    else:
        return(context_name)

def get_enrichment_df(enrichment_fn):
    enrichment_df = pd.read_csv(enrichment_fn, sep = "\t", header = 0)
    enrichment_df = enrichment_df.rename(columns = {u'state (Emission order)': 'state', u'Genome %' : 'percent_in_genome'})
    enrichment_df.columns = map(get_rid_of_stupid_file_tail, enrichment_df.columns)
    (num_state, num_enr_cont) = (enrichment_df.shape[0] - 1, enrichment_df.shape[1] - 2) # substract the first two columns: state and percent_in_genome 
    percent_genome_of_cont = enrichment_df.iloc[num_state, 2:]
    enrichment_df = enrichment_df.loc[:(num_state - 1)]
    # num_enr_cont : number of enrichment context (TSS, exon, etc.)
    return enrichment_df, num_state, num_enr_cont, percent_genome_of_cont


def calculate_percent_background_of_foreground(fg_df, percent_genome_of_cont_fg, fg_cont_name, frac_gen_bg):
	# this function will calculate the percentage of the background context that one foregrounc context occupies
	# frac_gen_bg: fraction of the genome that are of the background context
	# fg_cont_name: column name of the foreground context that we are trying to calculate
	calculate_df = (fg_df[['percent_in_genome', fg_cont_name]]).copy() 
	calculate_df['frac_gen_in_state'] = calculate_df['percent_in_genome'] / 100.0 # fraction of the genome that is in the state
	calculate_df['frac_cont_in_state'] = calculate_df[fg_cont_name] * calculate_df['frac_gen_in_state'] # #SM/#M = FE * #S/#G --> fraction of the context that is in the state
	calculate_df['frac_gen_in_state_in_cont'] = calculate_df['frac_cont_in_state'] * percent_genome_of_cont_fg[fg_cont_name]  / 100 # #SM/#G = #SM/#M * #M/#G
	frac_gen_in_fg_cont = np.sum(calculate_df['frac_gen_in_state_in_cont']) # fraction of gene in the fg context 
	percent_bg_in_fg_cont = frac_gen_in_fg_cont / frac_gen_bg * 100 # #Mf/#Mb = #Mf/#G * #G/#Mb * 100
	return percent_bg_in_fg_cont

def calculate_fraction_background_in_state(bg_df, percent_genome_of_cont_bg):
	background_cont_name = bg_df.columns[-1] # the last column is the column of the background fold enrichment
	bg_df['frac_gen_in_state'] = bg_df['percent_in_genome'] / 100.0 
	bg_df['frac_bg_in_state'] = bg_df[background_cont_name] * bg_df['frac_gen_in_state'] # #SM/#M = FE * #S/#G
	return np.array(bg_df['frac_bg_in_state'])

def calculate_FE_against_background(background_fn, foreground_fn, output_fn):
	bg_df, num_state_bg, num_enr_cont_bg, percent_genome_of_cont_bg = get_enrichment_df(background_fn)
	fg_df, num_state_fg, num_enr_cont_fg, percent_genome_of_cont_fg = get_enrichment_df(foreground_fn)
	assert num_state_fg == num_state_bg, 'Number of states between the foreground and background data is not the same. Some thing went wrong befor this step. Check your pipeline'
	assert num_enr_cont_bg == 1, 'Number of background enrichment context should be 1'
	assert np.array_equal(bg_df['percent_in_genome'], fg_df['percent_in_genome']), 'The percentage of each state in the genome in the foreground and background datasets are different. This is not good since we suppose you calculate foreground and background enrichment based on the same segmentation file'
	bg_cont_name = bg_df.columns[-1] # the column name of fold enrichment with the background 

	result_df = pd.DataFrame()
	result_df['state'] = bg_df['state']
	frac_bg_in_state = calculate_fraction_background_in_state(bg_df, percent_genome_of_cont_bg) # --> an array, each entry is each state denoting the fraction of the background that is in each state. 
	result_df['percent_in_genome'] = frac_bg_in_state * 100 # the column name is misleading, it should be percent_in_background. But we keep it this way for the sake of other downstream analysis
	fg_cont_name_list = list(fg_df.columns[2:])
	result_df[fg_cont_name_list] = (fg_df[fg_cont_name_list]).div(bg_df[bg_cont_name], axis = 0) # divide each enrichment context in the foreground by the enrichment of background. This is dividing each column from one dataframe by the same series (in this case, the enrichment from background)
	# get the percentage of the background that are in  each of the foreground context
	percent_bg_in_fg_cont = percent_genome_of_cont_fg / percent_genome_of_cont_bg[bg_cont_name] * 100 
	(num_state, num_col) = result_df.shape # get the current shape, so that we can add one last row to the dataframe
	result_df.loc[num_state] = ['Base', np.sum(result_df['percent_in_genome'])] + list(percent_bg_in_fg_cont)
	result_df.to_csv(output_fn, index = False, header = True, sep = '\t')
	print ("Done calculating all the necessary information and printing the fold enrichment of forground against background after: " )

def main():
	if len(sys.argv) != 4:
		usage()
	background_fn = sys.argv[1]
	cmh.check_file_exist(background_fn)
	foreground_fn = sys.argv[2]
	cmh.check_file_exist(foreground_fn)
	output_fn = sys.argv[3]
	cmh.create_folder_for_file(output_fn)
	print ("Done getting command line argument after: ")
	calculate_FE_against_background(background_fn, foreground_fn, output_fn)


def usage():
	print ("python calculate_fold_enrichment_against_background.py")
	print ("background_fn: fn of background enrichment")
	print ("foreground_fn: fn of foreground enrichment")
	print ("output_fn: where the output of FE against background is stored")
	exit(1)

main()