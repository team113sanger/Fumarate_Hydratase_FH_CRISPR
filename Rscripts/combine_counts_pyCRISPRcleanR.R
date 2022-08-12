#!/usr/bin/env Rscript

# This script creates a raw count file for pyCRISPRcleanR
# Format is sgRNA, gene, control(s), samples

suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages( library("dplyr") )
suppressPackageStartupMessages( library("tidyr") )
suppressPackageStartupMessages( library("stringr") )

option_list = list(
	make_option( c( "-f", "--files" ), type="character", default=NULL,
    			 help="comma-separated list of count filenames", metavar="character" ),
	make_option( c( "-c", "--columns" ), type="character", default=NULL,
                 help="comma-separated list of sample columns to join", metavar="character" ),
	make_option( c( "-p", "--prefix" ), type="character", default=NULL,
				 help="output file prefix", metavar="character" )
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$files) | length(opt$files) == 0){
  print_help(opt_parser)
  stop( "Please provide a list of count files to combine", call.=FALSE )
}

count_files <- unlist( strsplit( opt$files, "," ) )
if ( length(count_files) < 2) {
  print_help(opt_parser)
  stop( "Please provide at least two count files", call.=FALSE )
}

columns_to_merge <- unlist( strsplit( opt$columns, "," ) )
str(columns_to_merge)
if ( length(columns_to_merge) < 1 ) {
  print_help(opt_parser)
  stop( "Please provide at least one column name", call.=FALSE )
}
keep <- append(c("sgRNA","gene"),c(columns_to_merge))
print("Keeping samples:")
str(keep)

for ( i in 1:length(count_files ) ) {
	counts <- read.table(count_files[i], sep="\t", header=TRUE)
	if ( i == 1 ) {
		all_counts = counts
	} else {
		all_counts = all_counts %>% left_join(counts, by=c("sgRNA","gene","plasmid"))
	}
#	print(i)
#	print(head(all_counts))
}

all_counts <- all_counts[,keep,drop=FALSE]
#str(all_counts)

all_counts <- all_counts %>% filter(sgRNA != 'none')
all_counts <- all_counts %>% replace(., is.na(.), 0)

# for now, remove non-targeting sgRNAs for pyCRISPRcleanR
print("Removing NON-TARGETING sgRNA")
all_counts <- subset(all_counts, gene != "NON-TARGETING")


output_filename <- "combined_counts.tsv"
if ( length(opt$prefix) > 0 ) {
    output_filename <- paste(opt$prefix, output_filename, sep='_')
}

write.table(all_counts, output_filename, row.names=FALSE, sep="\t", quote=FALSE)

# make library file

#libinfo <- str_split_fixed(all_counts$sgRNA, "_",)
#libinfo
