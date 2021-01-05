print "Given the data from roadmap about gene info and gene expression in a bunch of cell types, we want to extract the regions that surround the TSS of each gene. Looking into the gene_info file provided by ROADMAP gene_expression data package, the TSS is the smaller coordinate if on the positive strand and the larger coordinate if it is on the negative strand."
print "This is the first step in getting a figure that looks like the supp_fig 2 in Ernst et al. 2011"
print "Gene info file was downloaded online from https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/ on 05/13/2019"
print ''
print ""
print ''
import sys
import os
import time 
import gzip
import replicate_helper as rhelp
from collections import Counter

INPUT_GENE_ENS_ID_INDEX = 0
INPUT_CHR_INDEX = 1
INPUT_BP1_INDEX = 2
INPUT_BP2_INDEX = 3
INPUT_STRAND_INDEX = 4
INPUT_GENE_TYPE = 5 
def get_TSS_surrounding_region(input_gene_info_fn, output_gene_TSS_region_fn, gene_exp_dict, celltype_list):
	# setting up input and output file stream
	inF = rhelp.open_file(input_gene_info_fn)
	rhelp.create_folder_for_file(output_gene_TSS_region_fn)
	outF = gzip.open(output_gene_TSS_region_fn, 'wb')
	# write the header line for the output file
	header_list = ['chr', 'start_tss_region', 'end_tss_region', 'gene_ENS_id', 'strand', 'gene_type', 'tss_bp'] + celltype_list
	outF.write("\t".join(header_list) + "\n")
	added_gene_counter = Counter([]) # need this to get rid of duplicated genes. There are 7 of them. 
	for line in inF:
		line_data = line.strip().split()
		gene_ens_id = line_data[INPUT_GENE_ENS_ID_INDEX]
		if added_gene_counter[gene_ens_id] == 1:
			continue
		if not gene_ens_id in gene_exp_dict:
			print "Gene id: " + gene_ens_id + " is not present in the expression file"
			continue
		added_gene_counter[gene_ens_id] += 1 # now that it passed all the test, we will add this gene to the list of genes to report!
		if line_data[INPUT_STRAND_INDEX] == '-1': #  this gene is on the negative strand. The TSS is the larger coordinate. 
			tss_index = max(int(line_data[INPUT_BP1_INDEX]), int(line_data[INPUT_BP2_INDEX]))
		if line_data[INPUT_STRAND_INDEX] == '1': # this gene is on the positive strand. The TSS is the smaller coodinate
			tss_index = min(int(line_data[INPUT_BP1_INDEX]), int(line_data[INPUT_BP2_INDEX]))
		# get the starting and ending of the region surrounding the TSS. It should be spanning roughly rhelp.WINDOW_SIZE basepairs around the TSS of each gene
		tss_region_start_index = rhelp.round_down_to_the_hundred(tss_index - rhelp.HALF_WINDOW_SIZE)
		tss_region_end_index = rhelp.round_up_to_the_hundred(tss_index + rhelp.HALF_WINDOW_SIZE)
		write_data = ['chr' + line_data[INPUT_CHR_INDEX], str(tss_region_start_index), str(tss_region_end_index), gene_ens_id, line_data[INPUT_STRAND_INDEX], line_data[INPUT_GENE_TYPE], str(tss_index)]
		this_gene_exp_list = gene_exp_dict[gene_ens_id]
		write_data += this_gene_exp_list 
		outF.write("\t".join(write_data) + "\n")
	inF.close()
	outF.close()

def process_gene_exp_data(gene_exp_fn):
	'''
	On 05/15/2019: I checked /u/home/h/havu73/project-ernst/data/roadmap_epigenome/gene_exp/57epigenomes.RPKM.pc.gz and confirmed that all the genes are unique in this file. there is no gene code that appears more than once.
	'''
	GENE_ID_INDEX = 0
	STARTING_INDEX_TO_COLLECT_EXP_DATA = 2 ## skip the first two entries in this line because it is the gene_id and cell type E000, the rest are the cell type code
	gene_exp_F = rhelp.open_file(gene_exp_fn)
	gene_exp_dict = {}
	header_line = gene_exp_F.readline()
	celltype_list = header_line.strip().split()[STARTING_INDEX_TO_COLLECT_EXP_DATA:] # skip the first two entries in this line because it is the gene_id and cell type E000, the rest are the cell type code
	for line in gene_exp_F:
		line_data = line.strip().split()
		gene_code = line_data[GENE_ID_INDEX]
		gene_exp_dict[gene_code] = line_data[STARTING_INDEX_TO_COLLECT_EXP_DATA:] # the following the gene code contains information about expression of the gene in different cell types. We want to keep a record of such data in a dictionary.
 	gene_exp_F.close()
 	return gene_exp_dict, celltype_list

def main():
	if len(sys.argv) != 4:
		usage()
	input_gene_info_fn = sys.argv[1]
	if not os.path.isfile(input_gene_info_fn):
		print ("input_gene_info_fn: " + input_gene_info_fn + " DOES NOT EXIST")
		usage()
	output_gene_TSS_region_fn = sys.argv[2]
	gene_exp_fn = sys.argv[3]
	if not os.path.isfile(gene_exp_fn):
		print ("gene_exp_fn: " + gene_exp_fn + " DOES NOT EXIST")
		usage()
	print ("Done getting command line ")
	gene_exp_dict, celltype_list = process_gene_exp_data(gene_exp_fn) # gene_exp_dict: keys: genes, values: list of gene expression in all cell types. celltype_list: list of cell types for which we got the data of gene expression
	print ("Done process_gene_exp_data")
	get_TSS_surrounding_region(input_gene_info_fn, output_gene_TSS_region_fn, gene_exp_dict, celltype_list) # wrote into output_fn a table with the following headers: ['chr', 'start_tss_region', 'end_tss_region', 'gene_ENS_id', 'strand', 'tss_bp'] + celltype_list. Strand is -1 --> tss is the max coordinate, strand is 1 --> tss is the min coordinate
	print ("Done get_TSS_surrounding_region " )

def usage():
	print ("pypy get_50kb_around_TSS_roadmap_gene_info.py")
	print ("input_gene_info_fn")
	print ("output_gene_TSS_region_fn")
	print ("gene_exp_fn")
	exit(1)

main()
