import numpy as np 
import pandas as pd 
import string
import os
import sys
import time
NUM_BP_PER_BIN = 200
def make_dir(directory):
	try:
		os.makedirs(directory)
	except:
		print ( 'Folder' + directory + ' is already created')



def check_file_exist(fn):
	if not os.path.isfile(fn):
		print ( "File: " + fn + " DOES NOT EXISTS")
		exit(1)
	return 

def check_dir_exist(fn):
	if not os.path.isdir(fn):
		print ( "Directory: " + fn + " DOES NOT EXISTS")
		exit(1)
	return 
	
def create_folder_for_file(fn):
	last_slash_index = fn.rfind('/')
	if last_slash_index != -1: # path contains folder
		make_dir(fn[:last_slash_index])
	return 

def get_command_line_integer(arg):
	try: 
		arg = int(arg)
		return arg
	except:
		print ( "Integer: " + str(arg) + " IS NOT VALID")
		exit(1)

		
def get_enrichment_df (enrichment_fn): # enrichment_fn follows the format of ChromHMM OverlapEnrichment's format
	enrichment_df = pd.read_csv(enrichment_fn, sep = "\t")
	# rename the org_enrichment_df so that it's easier to work with
	enrichment_df =	enrichment_df.rename(columns = {"state (Emission order)": "state", "Genome %": "percent_in_genome"})
	return enrichment_df

def get_non_coding_enrichment_df (non_coding_enrichment_fn):
	nc_enrichment_df = pd.read_csv(non_coding_enrichment_fn, sep = '\t')
	if len(nc_enrichment_df.columns) != 3:
		print ( "Number of columns in a non_coding_enrichment_fn should be 3. The provided file has " + str(len(nc_enrichment_df.columns)) + " columns.")
		print ( "Exiting, from ChromHMM_untilities_common_functions_helper.py")
		exit(1)
	# Now, we know that the nc_enrichment_df has exactly 3 columns
	# change the column names
	nc_enrichment_df.columns = ["state", "percent_in_genome", "non_coding"]
	return (nc_enrichment_df)
