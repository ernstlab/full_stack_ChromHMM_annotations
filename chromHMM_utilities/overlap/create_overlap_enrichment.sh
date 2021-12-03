# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
if [[ $# -ne 6 ]]   ; then
    echo 'Usage ./create_overlap_enrichment.sh \
    <folder where we store the output .sh file> \
    <name of the .sh file (without .sh)> \
    <fn of state segmentation file> \
    <fn of the folder/file with the target enrichment annotations> \
    <output_folder> \
    <output_prefix: prefix of the output file name(without .txt), which will be outputted by ChromHMM OverlapEnrichment>'
    exit 0
fi

ChromHMM=/u/home/h/havu73/project-ernst/ChromHMM1_18/ChromHMM.jar
shDir=$(readlink -f $1)
sh_file=$2
input_segment=$(readlink -f $3)
input_coord_dir=$(readlink -f $4)
output_folder=$(readlink -f $5)
output_prefix=$6

mkdir -p $shDir
mkdir -p $output_folder
echo ". /u/local/Modules/default/init/modules.sh
module load java
java -jar $ChromHMM OverlapEnrichment -b 1 -lowmem -noimage $input_segment $input_coord_dir $output_folder/$output_prefix"> ${shDir}/$sh_file.sh
    chmod +x ${shDir}/$sh_file.sh
echo "Done. File can be found at ${shDir}/$sh_file.sh"
