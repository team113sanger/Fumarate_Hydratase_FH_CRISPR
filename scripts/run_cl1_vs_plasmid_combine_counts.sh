#BSUB -q normal
#BSUB -J 5972_cl1_vs_plasmid_combine_counts
#BSUB -oo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_cl1_vs_plasmid_combine_counts.o
#BSUB -eo /lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours/logs/5972_cl1_vs_plasmid_combine_counts.e
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M 2000

export SINGULARITY_CACHEDIR=/software/team113/singularity_cache
export MODULEPATH=$MODULEPATH:/software/CGP/modules/modulefiles

RUN_DIR=/lustre/scratch117/casm/team113/projects/5972_Identification_of_CoOperative_Oncogenic_Events_in_Fumarate_Hydratase_FH_Deficient_Tumours
ANALYSIS_DIR=$RUN_DIR/analysis
COUNTS_DIR=$RUN_DIR/sample_counts
R_IMG=/software/team113/singularity_images/t113-singularity_r-3.6.0.base-1.0.1.sif

count_files="$COUNTS_DIR/Plasmid_Mouse_v2_library.counts,$COUNTS_DIR/cl1_1.counts,$COUNTS_DIR/cl1_2.counts,$COUNTS_DIR/cl1_3.counts"
count_columns="ERS3957550.sample,ERS3957545.sample,ERS3957544.sample,ERS3957548.sample"
contrast_prefix="$ANALYSIS_DIR/cl1_vs_plasmid/cl1_vs_plasmid"

module load ISG/singularity/03.2.0

singularity exec $R_IMG Rscript $RUN_DIR/scripts/combine_counts_pyCRISPRcleanR.R -f $count_files -c $count_columns -p $contrast_prefix
