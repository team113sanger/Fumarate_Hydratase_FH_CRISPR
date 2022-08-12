suppressPackageStartupMessages( library( tidyverse ) )
suppressPackageStartupMessages( library( optparse ) )
suppressPackageStartupMessages( library( ggrepel ) )
suppressPackageStartupMessages( library( RColorBrewer ) )
suppressPackageStartupMessages( library( GGally ) )

###############################################################################
#* --                                                                     -- *#
#* --                             OPTIONS                                 -- *#
#* --                                                                     -- *#
###############################################################################

option_list = list(
  make_option( c( "-d", "--dir" ), type = "character",
               help = "project directory", metavar = "character" )
);

opt_parser <- OptionParser( option_list = option_list );
opt <- parse_args( opt_parser );

if ( is.null( opt$dir ) ) {
  stop( "Please provide a project directory." ) 
}

if ( !dir.exists( opt$dir ) ) {
  stop( paste( "Project directory does not exist:", opt$dir ) ) 
}

###############################################################################
#* --                                                                     -- *#
#* --                          MAIN SCRIPT                                -- *#
#* --                                                                     -- *#
###############################################################################

# Set project directory
print( paste( "Setting project directory as:", opt$dir ) )
top_dir <- opt$dir

# Set up output directories
analysis_dir <- file.path( top_dir, 'post_analysis/qc_results' )
plot_dir <- file.path( analysis_dir, 'plots/qc' )
table_dir <- file.path( analysis_dir, 'tables/qc' )

# Read in raw counts
print( "Reading raw counts..." )
raw_counts.fn <- file.path( top_dir, 'analysis/all_vs_plasmid/all_vs_plasmid_combined_counts.tsv' )
raw_counts <- read_delim( raw_counts.fn, delim = "\t" )

# Read in corrected counts
print( "Reading corrected counts..." )
corrected_counts.fn <- file.path( top_dir, 'analysis/all_vs_plasmid/03_crispr_cleanr_corrected_counts.tsv' )
corrected_counts <- read_delim( corrected_counts.fn, delim = "\t" )

# Read in sample metadata
print( "Reading sample metadata..." )
sample_metadata.fn <- file.path( top_dir, '5972_sample_annotation.tsv' )
sample_metadata <- read_delim( sample_metadata.fn, delim = "\t", col_names = c( 'sample', 'lane_cram', 'count_sample_id', 'lane_reads' ) )
sample_metadata <- sample_metadata %>% 
                    mutate( sample = ifelse( sample == 'Plasmid_Mouse_v2_library', 'plasmid', sample ),
                            group = case_when( grepl( 'plasmid', sample ) ~ 'plasmid',
                                               grepl( 'flfl', sample ) ~ 'flfl',
                                               grepl( 'pFH', sample ) ~ 'pFH',
                                               grepl( 'cl1', sample ) ~ 'cl1' ) )

# Read in reference gene sets
print( "Reading reference gene sets.." )
ref_genes.files <- list.files( file.path( top_dir, 'ref_genes_mouse' ), full.names = T )
ref_genes <- data.frame()
for ( fn in ref_genes.files ) {
  tmp <- data.frame( 'gene' = scan( fn, character(), quote = "" ), 'ref_set' = str_replace( basename( fn ), '.txt', '' ) )
  if ( nrow( ref_genes ) == 0 ){
    ref_genes = tmp
  } else {
    ref_genes = rbind( ref_genes, tmp )
  }
}


# Gather counts and replace sample ids with sample names
print( "Gathering raw counts..." )
raw_counts.gathered <- raw_counts %>%
                        gather( count_sample_id, counts, -sgRNA, -gene ) %>%
                        left_join( sample_metadata %>% select( sample, count_sample_id ) %>% unique(), by = 'count_sample_id' ) %>%
                        select( sgRNA, gene, sample, counts ) 

# Gather counts and replace sample ids with sample names
print( "Gathering corrected counts..." )
corrected_counts.gathered <- corrected_counts %>%
                              gather( count_sample_id, counts, -sgRNA, -gene ) %>%
                              left_join( sample_metadata %>% select( sample, count_sample_id ) %>% unique(), by = 'count_sample_id' ) %>%
                              select( sgRNA, gene, sample, counts ) 

# Generate QC stats
print( "Generating sample qc stats..." )
count_stats <- raw_counts.gathered %>%
                  group_by( sample ) %>%
                  summarise( mapped_reads = sum( counts ),
                             million_mapped_reads = mapped_reads / 1e+6, 
                             num_guides = n(),
                             zero_sgrnas = sum( counts == 0 ),
                             prop_zero_sgrnas = round( zero_sgrnas / num_guides, 3 ),
                             low_sgrnas = sum( counts < 30 ),
                             prop_low_sgrnas = round( low_sgrnas / num_guides, 3 )  )

sample_qc_stats <- sample_metadata %>% 
                    select( -lane_cram ) %>%
                    group_by( sample, group ) %>%
                    summarise( total_reads = sum( lane_reads ) ) %>%
                    mutate( million_total_reads = total_reads/1e+6 ) %>% 
                    left_join( count_stats, by = 'sample' ) %>%
                    mutate( prop_mapped_reads = round( mapped_reads / total_reads, 2 ) )
sample_qc_stats.fn <- file.path( table_dir, 'sample_qc_stats.tsv' )
write.table( sample_qc_stats, sample_qc_stats.fn, row.names = F, quote = F, sep = "\t" )
print( paste( "Sample qc stats written to:", sample_qc_stats.fn ) )
  
# Plot mapping stats
print( "Plotting total and mapped reads (millions)..." )
mapped_and_total_reads <- sample_qc_stats %>%
                            select( sample, group, million_total_reads, million_mapped_reads ) %>%
                            gather( category, million_reads, -sample, -group )

total_and_mapped_reads.plot <- ggplot( mapped_and_total_reads, aes( x = sample, y = million_reads, fill = category ) ) +
                                  geom_bar( stat="identity", position ="identity", alpha = 0.8 ) +
                                  theme_classic() +
                                  theme( axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) ) +
                                  scale_fill_manual( values = c( 'gray50', 'gray80' ), labels = c( 'mapped', 'total' ), name = '' ) +
                                  ylab( 'Total reads (million)' ) +
                                  xlab( '' ) +
                                  facet_grid( ~ group, scales = 'free_x' )
total_and_mapped_reads.fn <- file.path( plot_dir, 'total_and_mapped_reads' )
ggsave( paste0( total_and_mapped_reads.fn, '.png'), total_and_mapped_reads.plot, device = 'png', width = 12, height = 8, dpi = 200 )
ggsave( paste0( total_and_mapped_reads.fn, '.pdf'), total_and_mapped_reads.plot, device = 'pdf', width = 12, height = 8, dpi = 200 )
print( paste( "Total and mapped reads barplot written to:", total_and_mapped_reads.fn ) )

# Plot low/zero guides
print( "Plotting low (< 30) and zero guides..." )
low_and_zero_guides <- sample_qc_stats %>%
                        select( sample, group, prop_low_sgrnas, prop_zero_sgrnas ) %>%
                        gather( category, prop_sgrnas, -sample, -group )

low_and_zero_guides.plot <- ggplot( low_and_zero_guides, aes( x = sample, y = prop_sgrnas, fill = category ) ) +
                            geom_bar( stat="identity", position ="identity", alpha = 0.8 ) +
                            theme_classic() +
                            theme( axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) ) +
                            scale_fill_manual( values = c( 'gray50', 'gray80' ), labels = c( 'zero', 'low (<30)' ), name = '' ) +
                            ylab( 'Proportion of total guides' ) +
                            xlab( '' ) +
                            ylim( 0, 0.06 ) + 
                            facet_grid( ~ group, scales = 'free_x' )
low_and_zero_guides.fn <- file.path( plot_dir, 'low_and_zero_reads' )
ggsave( paste0( low_and_zero_guides.fn, '.png'), low_and_zero_guides.plot, device = 'png', width = 12, height = 8, dpi = 200 )
ggsave( paste0( low_and_zero_guides.fn, '.pdf'), low_and_zero_guides.plot, device = 'pdf', width = 12, height = 8, dpi = 200 )
print( paste( "Low (< 30) and zero guide barplot written to:", low_and_zero_guides.fn ) )


# Plot MDS
print( "Plotting MDS..." )
counts.mat <- raw_counts.gathered %>% 
                spread( sample, counts ) %>% 
                select( -gene ) %>% 
                column_to_rownames( 'sgRNA' ) %>% 
                as.matrix()
counts.cor <- cor( counts.mat, method = 'spearman' )
mds.cor <- ( 1 - counts.cor ) %>% cmdscale() %>% as_tibble()
colnames( mds.cor ) <- c( "Dimension.1", "Dimension.2" )
mds.cor$sample <- colnames( counts.mat )
mds.cor <- mds.cor %>% left_join( sample_metadata %>% select( sample, group ) %>% unique(), by = 'sample' )

all_sample_mds <- ggplot( mds.cor, aes( x = Dimension.1, y = Dimension.2, color = group, shape = group, label = sample ) ) + 
                          geom_point( size = 3 ) +
                          theme_classic() +
                          geom_text_repel( point.padding = 0.5, segment.color = 'transparent' )
all_sample_mds.fn <- file.path( plot_dir, 'mds.all_samples' )
ggsave( paste0( all_sample_mds.fn, '.png'), all_sample_mds, device = 'png', width = 10, height = 8, dpi = 200 )
ggsave( paste0( all_sample_mds.fn, '.pdf'), all_sample_mds, device = 'pdf', width = 10, height = 8, dpi = 200 )
print( paste( "MDS plot written to:", all_sample_mds.fn ) )

# Pairs plots for counts
print( "Plotting pair-wise scatter..." )
pairs.scatter.fn <- file.path( plot_dir, 'pair_scatter_max_8000.png' )
png( pairs.scatter.fn, width = 2000, height = 2000, units = 'px' )
ggpairs( raw_counts.gathered %>% filter( counts <= 8000 ) %>% spread( sample, counts ) %>% select( -gene ) %>% drop_na() %>% column_to_rownames( 'sgRNA' ) )
dev.off()
print( paste( "Pairs plots written to:", pairs.scatter.fn ) )

# Plot reference set densities (raw)
print( "Plotting raw reference gene set densities (counts <= 5000)..." )
raw_counts.gene_sets <- raw_counts.gathered %>% 
                        left_join( ref_genes, by = 'gene' ) %>%
                        left_join( sample_metadata %>% select( sample, group ) %>% unique(), by = 'sample' ) %>%
                        mutate( ref_set = str_replace( ref_set, '_', ' ' ),
                                ref_set = ifelse( is.na( ref_set ), 'other', ref_set ),
                                ref_set = ifelse( ref_set == 'non essential', 'non-essential', ref_set ) )
ref_set_density.plot <- ggplot( raw_counts.gene_sets, aes( x = counts, fill = group ) ) +
                          geom_density() +
                          xlim( 0, 3000 ) +
                          theme_classic() +
                          facet_grid( sample ~ ref_set )
ref_set_density.fn <- file.path( plot_dir, 'ref_gene_density_max_3000.raw.png' )
ggsave( ref_set_density.fn, ref_set_density.plot, device = 'png', width = 16, height = 12 )
print( paste( "Reference raw gene set densities plot written to:", ref_set_density.fn ) )

# Plot reference set densities (corrected)
print( "Plotting corrected reference gene set densities (counts <= 5000)..." )
corrected_counts.gene_sets <- corrected_counts.gathered  %>%
                                left_join( ref_genes, by = 'gene' ) %>%
                                left_join( sample_metadata %>% select( sample, group ) %>% unique(), by = 'sample' ) %>%
                                mutate( ref_set = str_replace( ref_set, '_', ' ' ),
                                        ref_set = ifelse( is.na( ref_set ), 'other', ref_set ),
                                        ref_set = ifelse( ref_set == 'non essential', 'non-essential', ref_set ) )
corrected_ref_set_density.plot <- ggplot( corrected_counts.gene_sets, aes( x = counts, fill = group ) ) +
                                    geom_density() +
                                    xlim( 0, 500 ) +
                                    theme_classic() +
                                    facet_grid( sample ~ ref_set )
corrected_ref_set_density.fn <- file.path( plot_dir, 'corrected_ref_set_density_500.raw.png' )
ggsave( corrected_ref_set_density.fn, corrected_ref_set_density.plot, device = 'png', width = 16, height = 12 )
print( paste( "Reference corrected gene set densities plot written to:", corrected_ref_set_density.fn ) )

# Plot raw count boxplot
print( "Plotting raw count boxplot..." )
raw_count_boxplot <- ggplot( raw_counts.gathered %>% left_join( sample_metadata %>% select( sample, group ) %>% unique(), by = 'sample' ), aes( x = sample, y = counts, fill = group ) ) +
                      geom_boxplot() +
                      theme_classic() + 
                      xlab('') +
                      facet_grid( ~ group, scales = 'free_x' )
raw_count_boxplot.fn <- file.path( plot_dir, 'raw_count_boxplot' )
ggsave( paste0( raw_count_boxplot.fn, '.png'), raw_count_boxplot, device = 'png', width = 8, height = 10, dpi = 200 )
ggsave( paste0( raw_count_boxplot.fn, '.pdf'), raw_count_boxplot, device = 'pdf', width = 8, height = 10, dpi = 200 )
print( paste( "Raw count boxplot written to:", raw_count_boxplot.fn ) )

# Plot corrected count boxplot
print( "Plotting corrected count boxplot..." )
corrected_count_boxplot <- ggplot( corrected_counts.gathered %>% left_join( sample_metadata %>% select( sample, group ) %>% unique(), by = 'sample' ), aes( x = sample, y = counts, fill = group ) ) +
                            geom_boxplot() +
                            theme_classic() + 
                            xlab('') +
                            facet_grid( ~ group, scales = 'free_x' )
corrected_count_boxplot.fn <- file.path( plot_dir, 'corrected_count_boxplot' )
ggsave( paste0( corrected_count_boxplot.fn, '.png'), corrected_count_boxplot, device = 'png', width = 8, height = 10, dpi = 200 )
ggsave( paste0( corrected_count_boxplot.fn, '.pdf'), corrected_count_boxplot, device = 'pdf', width = 8, height = 10, dpi = 200 )
print( paste( "Corrected count boxplot written to:", corrected_count_boxplot.fn ) )


# Plot raw count violin
print( "Plotting raw count violin..." )
raw_count_violin <- ggplot( raw_counts.gathered %>% left_join( sample_metadata %>% select( sample, group ) %>% unique(), by = 'sample' ), aes( x = sample, y = counts, fill = group ) ) +
                      geom_violin() +
                      theme_classic() + 
                      xlab('') +
                      facet_grid( ~ group, scales = 'free_x' )
raw_count_violin.fn <- file.path( plot_dir, 'raw_count_violin.png' )
ggsave( raw_count_violin.fn, raw_count_violin, device = 'png', width = 8, height = 10 )
print( paste( "Raw count violin written to:", raw_count_violin.fn ) )

# Plot corrected count violin
print( "Plotting corrected count violin..." )
corrected_count_violin <- ggplot( corrected_counts.gathered %>% left_join( sample_metadata %>% select( sample, group ) %>% unique(), by = 'sample' ), aes( x = sample, y = counts, fill = group ) ) +
                          geom_violin() +
                          theme_classic() + 
                          xlab('') +
                          facet_grid( ~ group, scales = 'free_x' )
corrected_count_violin.fn <- file.path( plot_dir, 'corrected_count_violin.png' )
ggsave( corrected_count_violin.fn, corrected_count_violin, device = 'png', width = 8, height = 10 )
print( paste( "Corrected count violin written to:", corrected_count_violin.fn ) )