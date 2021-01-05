import string
import os
import gzip
WINDOW_SIZE = 50000
HALF_WINDOW_SIZE = int(WINDOW_SIZE / 2)
BIN_SIZE = 200
TOTAL_NUM_BINS_INCLUDING_TSS = int(WINDOW_SIZE / BIN_SIZE + 1)
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

def round_down_to_the_hundred(number):
	'''
	number = 125 --> 000
	number = 399 --> 200
	number = negative --> 0
	'''
	return max(BIN_SIZE * (int(number) / BIN_SIZE), 0)

def round_up_to_the_hundred(number):
	'''
	number = 125 --> 200
	number = 399 --> 400
	number = 200 --> 200
	'''
	return BIN_SIZE * (int(number) / BIN_SIZE + 1)

def create_folder_for_file(fn):
	last_slash_index = fn.rfind('/')
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

