# 5972 - Identification of cooperative oncogenic events in fumarate hydratase (FH) deficient tumours

Single guide CRISPR screen using Kosuke/Yusa Mouse V2 library ( ~18k genes / ~90k guides ).

**Experimental design:**  
3 conditions in triplicate:  
        WT – wild type  
        KO – FH knock out  
        pFH – FH reconstituted  

**Method summary:**  
Guides were quantified against the [Yusa Mouse V2 library](https://www.addgene.org/pooled-library/yusa-crispr-knockout-mouse-v2/) using [crisprReadCounts v1.3.1](https://github.com/cancerit/crisprReadCounts). Sample counts were merged into a single matrix before CRISPR-bias correction and normalising all sample counts against the plasmid using [pyCRISPRcleanR v2.0.8](https://github.com/cancerit/pyCRISPRcleanR). Conditions were compared using `mageck test` from [MAGeCK v0.5.9.2](https://github.com/vaofford/docker-mageck) with `--norm-method none`. Quality control and post-analysis tables and plots were generated in RStudio v1.2.1578 (R version 3.6.1).

| sample_id         | sample_name | zero_guides | num_reads | total_counts | mapped      |
|-------------------|-------------|-------------|-----------|--------------|-------------|
| ERS3957541.sample | flfl_1      | 0.011781004 | 61828736  | 58184282     | 0.941055661 |
| ERS3957542.sample | flfl_3      | 0.012966862 | 77901002  | 71246556     | 0.914578172 |
| ERS3957543.sample | flfl_2      | 0.013099856 | 67183821  | 63003694     | 0.937780749 |
| ERS3957544.sample | cl1_2       | 0.021445196 | 64412668  | 60608688     | 0.940943604 |
| ERS3957545.sample | cl1_1       | 0.022276405 | 71274845  | 67024514     | 0.940367026 |
| ERS3957546.sample | pFH_2       | 0.010340242 | 77748042  | 73368672     | 0.943672279 |
| ERS3957547.sample | pFH_1       | 0.009331708 | 65262702  | 61544375     | 0.943025237 |
| ERS3957548.sample | cl1_3       | 0.025745318 | 59490688  | 56017501     | 0.941617972 |
| ERS3957549.sample | pFH_3       | 0.008688906 | 68855406  | 64958735     | 0.943407915 |
| ERS3957550.sample | Plasmid     | 0.003801396 | 64243551  | 60980343     | 0.949205672 |


## Repository structure

### Scripts

* **perl** - converting reference gene lists from human to mouse and parsing library into different formats
* **python** - workflow runner for pulling sequencing data out of iRODS
* **scripts** - jobscripts to run analyses on farm5
* **Rscripts** - scripts and functions for building QC and post-analysis plots and tables, combining counts into single matrix and pulling human/mouse gene symbols from BioMart

### Reference gene lists (pyCRISPRcleanR / BAGEL)

* **ref_genes_human**
* **ref_genes_mouse**

### Analyses

* **analysis/qc_runs** - pyCRISPRcleanR results from CRISPR bias correction comparing each condition to the plasmid (to check for essential gene dropout)
* **analysis/all_vs_plasmid** - pyCRISPRcleanR results from CRISPR bias correction and using single matrix of all samples to normalise against the plasmid
* **analysis/normalised_to_plasmid** - MAGeCK results for contrasts using counts that had been normalised to the plasmid (analysis/all_vs_plasmid/03_crispr_cleanr_corrected_counts.tsv)

### Quality control and post-analysis

* **post_analysis/qc_results** - zero guides, sample mapping, plasmid correlation, MDS
* **post_analysis/normalised_results** - contrast plots, tables and summaries

## Methods overview

### Getting and parsing the CRISPR library into useful formats

Download Yusa Mouse V2 library (`library/yusa-crispr-knockout-mouse-v2.tsv`):  
https://www.addgene.org/static/data/plasmids/67/67988/67988-attachment_CgOMdRaU11JJ.xlsx

Parse library:
```
perl perl/parse_library.pl [path/to/library/folder]
```

Library formats:

* CGP (comma-delimited into sgRNA,gene,sequence) - `library/yusa-crispr-knockout-mouse-v2-cgp.tsv`
* pyCRISPRcleanR (tab-delimited into sgRNA,gene,chr,start,end) - `library/yusa-crispr-knockout-mouse-v2-pycrisprcleanr.tsv`
* complete (all of the above combined) - `library/yusa-crispr-knockout-mouse-v2-complete.tsv`
* plasmid counts - `library/yusa-crispr-knockout-mouse-v2-plasmid.tsv`


### Quantify each lane CRAM against CRISPR library and then merge sample counts

Download crisprReadCounts:
```
git clone https://github.com/cancerit/crisprReadCounts.git
cd crisprReadCounts/
chmod +x setup.sh
./setup.sh $PWD
export PATH=$PWD/crisprReadCounts/bin:$PATH
export PERL5LIB=$PWD/crisprReadCounts/lib/perl5:$PERL5LIB
```

NOTES:

*Requires samtools 1.9 to be available in the PATH.*

*N.B. requires the reference used to generate the CRAM file.*

*Counting was performed on a compute cluster (LSF) - commands and requirements are in the jobscripts but may have hardcoded paths*

Job script to generate count file per lane:
```
bsub < scripts/get_lane_counts.sh
```

Job script to merge lane counts into sample counts:
```
bsub < scripts/merge_lane_counts.sh
```

### Convert human gene symbols to mouse gene symbols using BioMart

Download human reference gene sets (`ref_genes_human`):
https://github.com/cancerit/dockstore-pyCRISPRcleanR/blob/master/examples/data/signatures.tar.gz

Convert files from human (ref_genes_human) to mouse (ref_genes_mouse):
```
perl perl/convert_human_ref_genes_to_mouse.pl [path/to/repository] > logs/convert_human_ref_genes_to_mouse.tsv
```

The human and mouse genes are found in `library/hgnc_to_mgi_symbols.tsv`.

### Run each condition against the plasmid to check that essential genes are dropping out

*Commands were run via LSF jobscripts and may contain hard coded paths*

WT (flfl) vs Plasmid (analysis/flfl_vs_plasmid):
```
bsub < scripts/run_flfl_vs_plasmid_combine_counts.sh
bsub < scripts/run_flfl_vs_plasmid_analysis.sh
```

KO (cl1) vs Plasmid (analysis/cl1_vs_plasmid):
```
bsub < scripts/run_cl1_vs_plasmid_combine_counts.sh
bsub < scripts/run_cl1_vs_plasmid_analysis.sh
```

pFH (pFH) vs Plasmid (analysis/pFH_vs_plasmid):
```
bsub < scripts/run_pFH_vs_plasmid_combine_counts.sh
bsub < scripts/run_pFH_vs_plasmid_analysis.sh
```

Results are in `analysis/qc_runs`.

### Correct CRISPR bias and normalise to plasmid

Scripts to run on farm5:
```
bsub < scripts/run_all_vs_plasmid_combine_counts.sh
bsub < scripts/run_all_vs_plasmid_analysis.sh
```
Results are in `analysis/all_vs_plasmid`.


### Compare conditions with MAGeCK using sample counts that were corrected and normalised to plasmid

Get MAGeCK container (Singularity v3.2.0):
```
singularity pull docker://quay.io/vaofford/mageck:v1.0.0
```

Script to run all contrasts (farm5):
```
bsub < scripts/run_all_normalised_to_plasmid_contrasts.sh
```

Results are in `analysis/normalised_to_plasmid`.

## Post-analysis processing and qc_plots

### Generate QC plots

To generate QC plots:

```
Rscript Rscripts/qc_plots.R -d $PWD
```

### Treatment-Control and Treatment-Plasmid quadrant plots

To generate plots and results tables:

```
Rscript Rscripts/plot_mageck_logfc_gene.R -t ${PWD}/analysis/normalised_to_plasmid/cl1_vs_flfl/cl1_vs_flfl.gene_summary.txt -c ${PWD}/analysis/normalised_to_plasmid/cl1_vs_plasmid/cl1_vs_plasmid.gene_summary.txt -m 'cl1 vs flfl' -p ${PWD}/post_analysis/normalised_results/cl1_vs_flfl/cl1_vs_flfl

Rscript Rscripts/plot_mageck_logfc_gene.R -t ${PWD}/analysis/normalised_to_plasmid/cl1_vs_pFH/cl1_vs_pFH.gene_summary.txt -c ${PWD}/analysis/normalised_to_plasmid/cl1_vs_plasmid/cl1_vs_plasmid.gene_summary.txt -m 'cl1 vs pFH' -p ${PWD}/post_analysis/normalised_results/cl1_vs_pFH/cl1_vs_pFH

Rscript Rscripts/plot_mageck_logfc_gene.R -t ${PWD}/analysis/normalised_to_plasmid/pFH_vs_flfl/pFH_vs_flfl.gene_summary.txt -c ${PWD}/analysis/normalised_to_plasmid/pFH_vs_plasmid/pFH_vs_plasmid.gene_summary.txt -m 'pFH vs flfl' -p ${PWD}/post_analysis/normalised_results/pFH_vs_flfl/pFH_vs_flfl
```

### Plot guide LFC for enriched genes

To generate plots:

```
Rscript Rscripts/overlapping_genes_plots.R -d $PWD
```

### Generate summary tables

To generate summary plots and tables (e.g. number of significantly enriched/depleted genes per contrast):

```
Rscript Rscripts/summary_tables_and_plots.R -d $PWD
```
