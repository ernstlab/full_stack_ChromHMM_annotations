
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
shDir=$(readlink -f $1)
sh_file=$2
binary_folder=$(readlink -f $3)
output_model_folder=$(readlink -f $4)
num_state=$5
assembly=$6

mkdir -p $shDir
mkdir -p $output_model_folder
# echo ". /u/local/Modules/default/init/modules.sh
# module load java
# java -jar $ChromHMM LearnModel -holdcolumnorder -pseudo -many -p 4 -n 150 -d -1 -nobed -nobrowser -noimage  -lowmem $binary_folder $output_model_folder $num_state $assembly"> ${shDir}/$sh_file.sh

echo ". /u/local/Modules/default/init/modules.sh
module load java
java -jar $ChromHMM_1_18 LearnModel -init random -splitrows -holdcolumnorder  -pseudo -many -p 6 -n 300 -d -1 -lowmem -gzip    $binary_folder $output_model_folder $num_state $assembly"> ${shDir}/$sh_file.sh



# echo ". /u/local/Modules/default/init/modules.sh
# module load java
# java -jar $ChromHMM LearnModel  -holdcolumnorder -init random -nobed -nobrowser -noimage  -lowmem $binary_folder $output_model_folder $num_state $assembly"> ${shDir}/$sh_file.sh

# 
    chmod +x ${shDir}/$sh_file.sh
# echo "Done. File can be found at ${shDir}/$sh_file.sh"
