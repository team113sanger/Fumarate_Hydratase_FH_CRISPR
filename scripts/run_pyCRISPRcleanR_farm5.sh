#!/bin/bash

counts=""
numcontrol=""
library=""
outdir=""
ref_genes="/usr/local/lib/python3.5/dist-packages/pyCRISPRcleanR/config/ref_genes"

usage() {  # Print a help message.
  echo "Usage: $0 [ -c countfile ] [ -n num of controls ] [ -o outdir ] [ -l library ] [ -r ref_genes ]" 1>&2 
}
exit_error() { 
  usage
  exit 1
}
 
while getopts ":c:n:o:l:r:" option; do
case "${option}"
in
	c) counts=${OPTARG};;
	n) numcontrol=${OPTARG};;
	o) outdir=${OPTARG};;
	l) library=${OPTARG};;
	r) ref_genes=${OPTARG};;
esac
done

if [ "$counts" = "" ] || [ "$outdir" = "" ] || [ "$library" = "" ]; then
	exit_error;
elif [ "$numcontrol" = "" ]; then
	echo "Using default number of controls = 1";
	numcontrol=1;
fi

export MODULEPATH=$MODULEPATH:/software/CASM/modules/modulefiles
module load pyCRISPRcleanR/2.0.15
export SINGULARITY_BINDPATH=/lustre
export SINGULARITY_CACHEDIR=/software/team113/singularity_cache
#export PATH=/software/singularity-v3.2.0/bin:$PATH; unset PYTHONPATH; unset R_LIBS; module load ISG/singularity/03.2.0;
pyCRISPRcleanR -f $counts  -l $library -nc $numcontrol -cc -mk -bl -o $outdir -gs $ref_genes


