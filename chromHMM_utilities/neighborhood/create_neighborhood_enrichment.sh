# Copyright 2021 Ha Vu (havu73@ucla.edu)

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

if [[ $# -ne 6 ]] ; then
    echo 'Usage ./create_neighborhood_enrichment.sh \
    <output directory to store the .sh files> \
    <prefix of .sh files> 
    <segment_fn> \
    <assembly_name: hg18, hg19, hg38, etc.> \
    <output_prefix>'
    exit 0
fi

ChromHMM_folder=/u/home/h/havu73/project-ernst/ChromHMM1_18/
ChromHMM_1_18=/u/home/h/havu73/project-ernst/ChromHMM1_18/ChromHMM.jar
shDir=$(readlink -f $1)
sh_file=$2
sh_fn=${shDir}/$sh_file.sh
segment_fn=$(readlink -f $3)
assembly_name=$4
anchor_folder=${ChromHMM_folder}/ANCHORFILES/${assembly_name}/
if [ ! -d $anchor_folder ]; then # chekc if the anchor folder exist
    echo 'anchor_folder: ${anchor_folder} DOES NOT EXIST. The output of this file therefore will be trash'
fi 
anchor_tss_fn=${anchor_folder}/RefSeqTSS.${assembly_name}.txt.gz
anchor_tes_fn=${anchor_folder}/RefSeqTES.${assembly_name}.txt.gz
# check if the TSS anchor file exists
if [ ! -f $anchor_tss_fn ]; then
    echo 'anchor_tss_fn: ${anchor_tss_fn} DOES NOT EXIST. The .sh file created from this file will not run properly'
fi 

# check if the TES anchor file exists
if [ ! -f $anchor_tes_fn ]; then
    echo 'anchor_tes_fn: ${anchor_tes_fn} DOES NOT EXIST. The .sh file created from this file will not run properly'
fi 

output_folder=$5
output_prefix=$6
mkdir -p $shDir
mkdir -p $output_folder
# by this point, we are done getting and checking output

echo ". /u/local/Modules/default/init/modules.sh
module load java
java -jar $ChromHMM_1_18 NeighborhoodEnrichment -l 50 -r 50 $segment_fn $anchor_tss_fn $output_folder/${output_prefix}_tss_$assembly_name 

java -jar $ChromHMM_1_18 NeighborhoodEnrichment -l 50 -r 50 $segment_fn $anchor_tes_fn $output_folder/${output_prefix}_tes_$assembly_name 
"> $sh_fn


chmod +x $sh_fn

