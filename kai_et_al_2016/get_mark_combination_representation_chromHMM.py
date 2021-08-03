
import sys
import math 
import numpy as np

def to_str(number):
	return str(number)

def write_graph_representation(inputFName, outputFName, Num_marks_to_take, total_marks):
	results = [0] * total_marks
	inF = open(inputFName, 'r')
	outF = open(outputFName, 'w')
	for i in range(Num_marks_to_take):
		line = inF.readline().strip()
		mark_index = int(line.split("_")[1])
		results[mark_index] = 1

	outF.write("".join(map(to_str, results)))
	inF.close()
	outF.close()
# write_graph_representation("../../assay_selection_method/correlation_binary/top80_fixed_03052018/kai_output_mark_rank.txt", \
# 	'../../assay_selection_method/correlation_binary/top80_fixed_03052018/binary_mark_representation_top80.txt', 80, 4649)

def write_graph_representation_increasing_marks(inputFName, outputFName, Num_marks_to_take, total_marks):
	results = [0] * total_marks
	inF = open(inputFName, 'r')
	outF = open(outputFName, 'w')
	chosen_marks = []
	for i in range(Num_marks_to_take):
		line = inF.readline()
		mark_to_include = int(line.strip())
		chosen_marks.append(mark_to_include)
		results[mark_to_include] = 1
		outF.write("".join(map(to_str, results)) + "\n")
	inF.close()
	outF.close()
"""
write_graph_representation_increasing_marks("../../assay_selection_method/correlation_binary/top80_fixed_03052018/kai_output_mark_rank.txt", \
	"../../assay_selection_method/correlation_binary/top80_fixed_03052018/binary_mark_representation_increasing_marks.txt", 4649)
"""
def get_total_num_marks (inputFName):
	inF = open(inputFName, 'r')
	total_marks = len(inF.readlines())
	inF.close()
	return total_marks

def main():
	if len(sys.argv) != 4:
		print "Wrong input man!"
		print "python rabbitsRecurrence.py \
		<input file> \
		<output file> \
		<Num_marks_to_take> \
		<total_marks>"
		exit(1)
	else:
		inputFName = sys.argv[1]
		outputFName = sys.argv[2]
		Num_marks_to_take = int(sys.argv[3])

		total_marks = get_total_num_marks(inputFName)
		# write_graph_representation(inputFName, outputFName, Num_marks_to_take, total_marks)
		write_graph_representation_increasing_marks(inputFName, outputFName, Num_marks_to_take, total_marks)
main()