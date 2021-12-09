import os
import sys
import string
def make_dir(directory):
	try:
		os.makedirs(directory)
	except:
		print 'Folder' + directory + ' is already created'



def check_file_exist(fn):
	if not os.path.isfile(fn):
		print "File: " + fn + " DOES NOT EXISTS"
		exit(1)
	return 

def check_dir_exist(fn):
	if not os.path.isdir(fn):
		print "Directory: " + fn + " DOES NOT EXISTS"
		exit(1)
	return 
	
def create_folder_for_file(fn):
	last_slash_index = string.rfind(fn, '/')
	if last_slash_index != -1: # path contains folder
		make_dir(fn[:last_slash_index])
	return 

def get_command_line_integer(arg):
	try: 
		arg = int(arg)
		return arg
	except:
		print "Integer: " + str(arg) + " IS NOT VALID"
		exit(1)