import sys
def strip(word):
	return word.strip()
	
def get_mark_names_from_output(mark_rank_fname, mark_names_file, output_fName, num_marks):
	rankF = open(mark_rank_fname, 'r')
	markF = open(mark_names_file, 'r')
	mark_names = markF.readlines()
	mark_names = map(strip, mark_names)
	outF = open(output_fName, 'w')
	for i in range(num_marks):
		line = rankF.readline().strip()
		mark_index = int(line) #int(line.split("_")[1])
		outF.write(mark_names[mark_index] + "\n")
	rankF.close()
	markF.close()
	outF.close()

def main():
	if len(sys.argv) != 5:
		print "Wrong input man!"
		print "python rabbitsRecurrence.py \
		<input file> \
		<output file> \
		<Num_marks_to_take> \
		<total_marks>"
		exit(1)
	else:
		kai_rank_FName = sys.argv[1]
		mark_name_FName = sys.argv[2]
		outputFName = sys.argv[3]
		num_marks_to_take = int(sys.argv[4])
		get_mark_names_from_output(kai_rank_FName, mark_name_FName, outputFName, num_marks_to_take)

main()