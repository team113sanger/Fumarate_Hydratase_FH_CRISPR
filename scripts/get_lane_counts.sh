#!/bin/bash

#BSUB -q normal
#BSUB -J 5972_get_lane_counts
#BSUB -oo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_get_lane_counts.o 
#BSUB -eo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_get_lane_counts.e
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M 2000

export RUN_DIR=/lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours

# Code for counting
#cd $RUN_DIR
#git clone https://github.com/cancerit/crisprReadCounts.git
#cd crisprReadCounts/
#chmod +x setup.sh
#./setup.sh $PWD
export PATH=$RUN_DIR/crisprReadCounts/bin:$PATH
export PERL5LIB=$RUN_DIR/crisprReadCounts/lib/perl5:$PERL5LIB

# samtools
export MODULEPATH=/software/modules/:/software/CGP/modules/modulefiles
module load samtools/1.9

# Required files and paths
reference="$RUN_DIR/reference/Mus_musculus.GRCm38.68.dna.toplevel.fa"
library="$RUN_DIR/library/yusa-crispr-knockout-mouse-v2-cgp.tsv"
plasmid_counts="$RUN_DIR/library/yusa-crispr-knockout-mouse-v2-plasmid.tsv"

# Do counting
for cram_file in $RUN_DIR/cram/*.cram
do
    count_file=$(echo "$cram_file" | sed 's/cram/counts/g')
    echo "Parsing $cram_file to $count_file..."
    crisprReadCounts.pl -l $library -i $cram_file -p $plasmid_counts -o $count_file -r $reference
done
