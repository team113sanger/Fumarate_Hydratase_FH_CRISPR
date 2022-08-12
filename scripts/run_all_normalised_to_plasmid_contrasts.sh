#BSUB -q normal
#BSUB -J 5972_norm2plasmid_analyses
#BSUB -oo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_all_norm2plasmid_analysis.o
#BSUB -eo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_all_norm2plasmid_analysis.e
#BSUB -R "select[mem>3000] rusage[mem=3000]"
#BSUB -M 3000

RUN_DIR=/lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours

counts="$RUN_DIR/analysis/all_vs_plasmid/03_crispr_cleanr_corrected_counts.tsv"
flfl="ERS3957541.sample,ERS3957543.sample,ERS3957542.sample"
cl1="ERS3957545.sample,ERS3957544.sample,ERS3957548.sample"
pFH="ERS3957547.sample,ERS3957546.sample,ERS3957549.sample"
plasmid="ERS3957550.sample"

SINGULARITY_BINDPATH=$RUN_DIR
MODULEPATH=$MODULEPATH:/software/CGP/modules/modulefiles

module load mageck/0.5.9.3

declare -a contrasts=("cl1_vs_flfl" "cl1_vs_pFH" "pFH_vs_flfl" "flfl_vs_plasmid" "pFH_vs_plasmid" "cl1_vs_plasmid")
for contrast in "${contrasts[@]}"
do
    echo "Running MAGeCK for: $contrast"

    RESULTS_DIR="$RUN_DIR/analysis/normalised_to_plasmid/${contrast}"
    mkdir $RESULTS_DIR
    echo "Results found in: $RESULTS_DIR"

    IFS='_' read -r -a conditions <<< "$contrast"
    control="${conditions[2]}"
    treatment="${conditions[0]}"

    control_samples="${!control}"
    treatment_samples="${!treatment}"
    mageck test -k $counts -c $control_samples -t $treatment_samples --norm-method none -n "${RESULTS_DIR}/${contrast}"
done
