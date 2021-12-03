if [[ $# -ne 8 ]] ; then
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

get_nb_data="/u/home/h/havu73/project-ernst/source/naive_bayes/get_naive_bayes_data.py"
naive_bayes="/u/home/h/havu73/project-ernst/source/naive_bayes/naive_bayes_version3.py"
get_mark_names="/u/home/h/havu73/project-ernst/source/naive_bayes/get_mark_names_from_naive_bayes_output.py"
shDir=$(readlink -f $1)
sh_file=$2
bin_dir=$(readlink -f $3)
binFile_fn=$(readlink -f $4)
output_folder=$(readlink -f $5)
num_process=$6
max_num_marks_to_picks=$7
segmentation_fn=$(readlink -f $8)


mkdir -p $shDir
mkdir -p $output_folder
prior_prob_fn=$output_folder/prior_probabilities.txt
cond_prob_fn=$output_folder/conditional_probabilities.txt
draft_folder=$output_folder/draft/
nb_output_fn=$output_folder/top${max_num_marks_to_picks}_marks
nb_mark_names_fn=$output_folder/top_mark_names.txt
echo "pypy $get_nb_data $segmentation_fn $bin_dir $binFile_fn $output_folder" > ${shDir}/$sh_file.sh
echo "pypy $naive_bayes $cond_prob_fn $prior_prob_fn $bin_dir $binFile_fn $segmentation_fn $nb_output_fn $max_num_marks_to_picks $draft_folder $num_process" >> ${shDir}/$sh_file.sh
echo "python $get_mark_names $nb_output_fn $bin_dir $nb_mark_names_fn" >> ${shDir}/$sh_file.sh
chmod +x ${shDir}/$sh_file.sh
