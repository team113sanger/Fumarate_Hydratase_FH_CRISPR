#BSUB -q normal
#BSUB -J 5972_plasmid_vs_plasmid
#BSUB -oo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_plasmid_vs_plasmid.o
#BSUB -eo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_plasmid_vs_plasmid.e
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M 2000

RUN_DIR=/lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours

ANALYSIS_DIR="$RUN_DIR/analysis"
RESULTS_DIR="$ANALYSIS_DIR/Plasmid_Mouse_v2_library_vs_plasmid"
COUNTS_DIR="$RUN_DIR/sample_counts"

library="$RUN_DIR/library/yusa-crispr-knockout-mouse-v2-crisprcleanr.tsv"
counts="$RESULTS_DIR/combined_counts.tsv"
ref_genes="$RUN_DIR/ref_genes_mouse"
num_control=1

mkdir $RESULTS_DIR

awk -F"\t" 'BEGIN{OFS="\t"} { if (NR == 1) {print "sgRNA", "gene", "plasmid", "Plasmid_Mouse_v2"} else {gsub(/\r/, "", $4); print $1,$2,$4,$3} }' $COUNTS_DIR/Plasmid_Mouse_v2_library.counts > "$RESULTS_DIR/combined_counts.tsv"

$RUN_DIR/scripts/run_pyCRISPRcleanR_farm5.sh -c $counts -l $library -o $RESULTS_DIR -n $num_control -r $ref_genes

