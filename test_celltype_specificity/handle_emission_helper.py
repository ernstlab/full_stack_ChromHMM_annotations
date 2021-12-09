import pandas as pd
import numpy as np 
import os
import sys
MARK_NAME = 'mark_name'
CT = 'ct'
CHROM_MARK = 'chrom_mark'
DEFAULT_METADATA_FN = '/Users/vuthaiha/Desktop/window_hoff/ROADMAP_aligned_reads/ROADMAP_metadata_july2013.csv'

def get_emission_df_right(emission_fn): 
	emission_df = pd.read_csv(emission_fn, sep = "\t", header = 0)
	(num_state, num_marks) = (emission_df.shape[0], emission_df.shape[1] - 1)
	emission_df = emission_df.transpose()
	emission_df.columns = emission_df.iloc[0]
	emission_df = emission_df[1:]
	emission_df[MARK_NAME] = emission_df.index # example: E013-H3K36me3
	emission_df[CT] = emission_df[MARK_NAME].apply(lambda x: x.split("-")[0]) # ex: E013
	emission_df[CHROM_MARK] = emission_df[MARK_NAME].apply(lambda x: x.split("-")[1]) # ex: H3K36me3
	# emission_df: columns: <state index, 1 based>, mark_name: experiment name, ct: the cell type associated with this experiment, chrom_mark: the chromatin mark associated with this experiment 
	return(emission_df, num_state, num_marks)

def get_metada (metadata_fn):
	meta_df = pd.read_csv(metadata_fn, sep = ',', header = 0)
	meta_df = meta_df.rename(columns = {u'Epigenome ID (EID)' : "CT_NAME", u'Epigenome name (from EDACC Release 9 directory)' : 'Epig_name'})
	meta_df = meta_df[["CT_NAME", "Epig_name", "GROUP", "TYPE", "ANATOMY"]]
	return meta_df

def make_dir(directory):
	try:
		os.makedirs(directory)
	except:
		print 'Folder' + directory + ' is already created'

def create_folder_for_file(fn):
	last_slash_index = string.rfind(fn, '/')
	if last_slash_index != -1: # path contains folder
		make_dir(fn[:last_slash_index])
	return 

def open_file(fn): 
	"""
	Open a zipped or unzipped file. Return the open file object
	"""
	if fn[-3:] == ".gz":
		# zip file
		F = gzip.open(fn, 'rb')
	else: 
		# non zip file
		F = open(fn, 'r')
	return F
