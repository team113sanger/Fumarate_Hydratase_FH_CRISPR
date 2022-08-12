#!/bin/bash

#BSUB -q normal
#BSUB -J 5972_merge_lane_counts
#BSUB -oo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_merge_lane_counts.o 
#BSUB -eo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_merge_lane_counts.e
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M 2000

export RUN_DIR=/lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours

# Code for merging
#cd $RUN_DIR
#git clone https://github.com/cancerit/crisprReadCounts.git
#cd crisprReadCounts/
#chmod +x setup.sh
#./setup.sh $PWD
export PATH=$RUN_DIR/crisprReadCounts/bin:$PATH
export PERL5LIB=$RUN_DIR/crisprReadCounts/lib/perl5:$PERL5LIB

sample_names=('Plasmid_Mouse_v2_library' 'flfl_1' 'flfl_2' 'flfl_3' 'cl1_1' 'cl1_2' 'cl1_3' 'pFH_1' 'pFH_2' 'pFH_3')

for sn in ${sample_names[@]}; do
    echo "Sample to merge: $sn"

    lane_counts=$(find $RUN_DIR/counts -name *$sn*)
    files_to_merge=$(echo "${lane_counts[@]}" | paste -s -d',')
    echo "Lane counts to merge: $files_to_merge"

    sample_counts="$RUN_DIR/sample_counts/$sn.counts"
    echo "Sample counts written to: $sample_counts"

    crisprMergeResults.pl -i "$files_to_merge" -o "$sample_counts" -p y
done
