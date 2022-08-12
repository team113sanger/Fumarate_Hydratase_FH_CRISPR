#BSUB -q normal
#BSUB -J 5972_flfl_vs_cl1_analysis
#BSUB -oo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_flfl_vs_cl1_analysis.o
#BSUB -eo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_flfl_vs_cl1_analysis.e
#BSUB -R "select[mem>3000] rusage[mem=3000]"
#BSUB -M 3000

RUN_DIR=/lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours

ANALYSIS_DIR="$RUN_DIR/analysis"
RESULTS_DIR="$ANALYSIS_DIR/flfl_vs_cl1"
COUNTS_DIR="$RUN_DIR/sample_counts"

library="$RUN_DIR/library/yusa-crispr-knockout-mouse-v2-crisprcleanr.tsv"
counts="$RESULTS_DIR/flfl_vs_cl1_combined_counts.tsv"
ref_genes="$RUN_DIR/ref_genes_mouse"
num_control=3

mkdir $RESULTS_DIR
cd $RESULTS_DIR

$RUN_DIR/scripts/run_pyCRISPRcleanR_farm5.sh -c $counts -l $library -o $RESULTS_DIR -n $num_control -r $ref_genes
