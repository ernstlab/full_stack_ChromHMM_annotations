input_segment_fn=$(readlink -f $1)
output_folder=$(readlink -f $2)
num_state=$3
mkdir -p $output_folder

for state_index in `seq 1 $num_state`
do

done