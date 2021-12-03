# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import sys
from collections import Counter
def to_str(number):
	return str(number)

def write_graph_representation_increasing_marks(inputFName, outputFName, Num_marks_to_take, total_marks):
	results = [0] * total_marks
	inF = open(inputFName, 'r')
	outF = open(outputFName, 'w')
	chosen_marks = []
	for i in range(Num_marks_to_take):
		line = inF.readline()
		mark_to_include = int(line.strip().split()[0])
		chosen_marks.append(mark_to_include)
		results[mark_to_include] = 1
		outF.write("".join(map(to_str, results)) + "\n")
	inF.close()
	outF.close()

def get_num_marks_to_take (inputFName):
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
		Num_marks_to_take = get_num_marks_to_take(inputFName)
		total_marks = int(sys.argv[3])
		# write_graph_representation(inputFName, outputFName, Num_marks_to_take, total_marks)
		write_graph_representation_increasing_marks(inputFName, outputFName, Num_marks_to_take, total_marks)
main()