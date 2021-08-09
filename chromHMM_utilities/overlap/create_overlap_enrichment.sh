if [[ $# -ne 6 ]]   ; then
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
