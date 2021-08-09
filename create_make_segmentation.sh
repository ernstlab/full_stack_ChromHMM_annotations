if [[ $# -ne 6 ]] ; then
    echo 'Usage ./create_job_array_eval_subset.sh \
    <output directory to store the .sh files> \
    <model file> 
    <binary data folder> \
    <segments folder> \
    <prefix output> \
    <mark comb file> \
    <start_index>\
    <end index> \
    <number of cores> \
    <gigs of RAM on each core> \
    <total hours> \
    <highp 0/1>'
    exit 0
fi

ChromHMM=/u/home/h/havu73/project-ernst/ChromHMM1_18/ChromHMM.jar
ChromHMM_1_18=/u/home/h/havu73/project-ernst/ChromHMM1_18/ChromHMM.jar
get_bed_file_right=/u/home/h/havu73/project-ernst/source/model_evaluation/get_bed_file_right.py
shDir=$(readlink -f $1)
sh_file=$2
model_file=$(readlink -f $3)
input_bin_dir=$(readlink -f $4)
output_folder=$(readlink -f $5)
num_state=$6
mkdir -p $shDir
mkdir -p $output_folder
echo ". /u/local/Modules/default/init/modules.sh
module load java
java -jar $ChromHMM_1_18 MakeSegmentation -lowmem $model_file $input_bin_dir $output_folder
#pypy $get_bed_file_right $output_folder/genome_${num_state}_segments.bed $output_folder/corrected_genome_${num_state}_segments.bed"> ${shDir}/$sh_file.sh
chmod +x ${shDir}/$sh_file.sh
# echo "Done. File can be found at ${shDir}/$sh_file.sh"