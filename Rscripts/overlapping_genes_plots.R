suppressPackageStartupMessages( library( tidyverse ) )
suppressPackageStartupMessages( library( optparse ) )
suppressPackageStartupMessages( library( viridis ) )

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
results_dir <- file.path( top_dir, 'analysis/normalised_to_plasmid' )
analysis_dir <- file.path( top_dir, 'post_analysis/normalised_results' )
plot_dir <- file.path( analysis_dir, 'cl1_pFH_overlapping_genes' )

# Read in enriched gene lists
print( "Reading in enriched gene lists..." )
cl1_vs_flfl.enriched.c1_only <- read.table( file.path( analysis_dir, 'cl1_vs_flfl/cl1_vs_flfl_mageck_c1_only_sig_enriched_genes.tsv' ), header = T,
                                            sep = "\t", col.names = c( 'gene', 'cl1_vs_flfl.LFC', 'cl1_vs_flfl.FDR', 'cl1_vs_plasmid.LFC', 'cl1_vs_plasmid.FDR' ) )
cl1_vs_flfl.enriched.c1c2 <- read.table( file.path( analysis_dir, 'cl1_vs_flfl/cl1_vs_flfl_mageck_c1c2_sig_enriched_genes.tsv' ), header = T, 
                                         sep = "\t", col.names = c( 'gene', 'cl1_vs_flfl.LFC', 'cl1_vs_flfl.FDR', 'cl1_vs_plasmid.LFC', 'cl1_vs_plasmid.FDR' ) )
cl1_vs_pFH.enriched.c1_only <- read.table( file.path( analysis_dir, 'cl1_vs_pFH/cl1_vs_pFH_mageck_c1_only_sig_enriched_genes.tsv' ), header = T, 
                                           sep = "\t", col.names = c( 'gene', 'cl1_vs_pFH.LFC', 'cl1_vs_pFH.FDR', 'cl1_vs_plasmid.LFC', 'cl1_vs_plasmid.FDR' ) )
cl1_vs_pFH.enriched.c1c2 <- read.table( file.path( analysis_dir, 'cl1_vs_pFH/cl1_vs_pFH_mageck_c1c2_sig_enriched_genes.tsv' ), header = T, 
                                        sep = "\t", col.names = c( 'gene', 'cl1_vs_pFH.LFC', 'cl1_vs_pFH.FDR', 'cl1_vs_plasmid.LFC', 'cl1_vs_plasmid.FDR' ) )

enriched_genes <- rbind( cl1_vs_flfl.enriched.c1_only %>% mutate( cl1_vs_flfl.is_hit = 'c1_only' ), cl1_vs_flfl.enriched.c1c2 %>% mutate( cl1_vs_flfl.is_hit = 'c1c2' ) ) %>%
                  left_join( rbind( cl1_vs_pFH.enriched.c1_only %>% mutate( cl1_vs_pFH.is_hit = 'c1_only' ), cl1_vs_pFH.enriched.c1c2 %>% mutate( cl1_vs_pFH.is_hit = 'c1c2' ) ), by = 'gene' )

# Read in guide results
print( "Reading in MAGeCK sgrna results..." )
cl1_vs_flfl.sgrna <- read.table( file.path( results_dir, 'cl1_vs_flfl/cl1_vs_flfl.sgrna_summary.txt' ), header = T, sep = "\t" )
cl1_vs_pFH.sgrna <- read.table( file.path( results_dir, 'cl1_vs_pFH/cl1_vs_pFH.sgrna_summary.txt' ), header = T, sep = "\t" )
cl1_vs_plasmid.sgrna <- read.table( file.path( results_dir, 'cl1_vs_plasmid/cl1_vs_plasmid.sgrna_summary.txt' ), header = T, sep = "\t" )

all_sgrna_lfc <- rbind( cl1_vs_flfl.sgrna %>% select( sgrna, gene = Gene, LFC ) %>% mutate( contrast = 'cl1_vs_flfl' ),
                        cl1_vs_pFH.sgrna %>% select( sgrna, gene = Gene, LFC ) %>% mutate( contrast = 'cl1_vs_pFH' ),
                        cl1_vs_plasmid.sgrna %>% select( sgrna, gene = Gene, LFC ) %>% mutate( contrast = 'cl1_vs_plasmid' ) ) %>%
                  mutate( contrast = factor( contrast, levels = c( 'cl1_vs_plasmid', 'cl1_vs_flfl', 'cl1_vs_pFH' ) ) )

# Plot c1c2 overlaps
print( "Plotting c1c2 overlapping gene barplot..." )
c1c2.genes <- enriched_genes %>% filter( cl1_vs_flfl.is_hit == 'c1c2' & cl1_vs_pFH.is_hit == 'c1c2' )
c1c2.sgrna <- all_sgrna_lfc %>% filter( gene %in% c1c2.genes$gene ) 

c1c2.sgrna.plot <- ggplot( c1c2.sgrna, aes( x = LFC, y = sgrna, fill = contrast ) ) +
                    geom_bar( stat = 'identity' ) +
                    scale_fill_viridis( discrete = T ) +
                    theme_bw() +
                    scale_x_continuous( limits = c( -1, 5 ), breaks = c( -1:5 ) ) +
                    xlab('') +
                    facet_grid( gene ~ contrast, scales = 'free_y' )
c1c2.sgrna.fn <- file.path( plot_dir, 'overlap_c1c2.cl1_vs_flfl.cl1_vs_pFH.png' )
ggsave( c1c2.sgrna.fn, c1c2.sgrna.plot, device = 'png', width = 16, height = 6 )
print( paste( "c1c2 overlapping gene barplot written to:", c1c2.sgrna.fn ) )

# Plot c1_only overlaps
print( "Plotting c1_only overlapping gene barplot..." )
c1_only.genes <- enriched_genes %>% filter( cl1_vs_flfl.is_hit == 'c1_only' & cl1_vs_pFH.is_hit == 'c1_only' )
c1_only.sgrna <- all_sgrna_lfc %>% filter( gene %in% c1_only.genes$gene )

c1_only.sgrna.plot <- ggplot( c1_only.sgrna, aes( x = LFC, y = sgrna, fill = contrast ) ) +
                        geom_bar( stat = 'identity' ) +
                        scale_fill_viridis( discrete = T ) +
                        theme_bw() +
                        scale_x_continuous( limits = c( -1, 5 ), breaks = c( -1:5 ) ) +
                        xlab('') +
                        facet_grid( gene ~ contrast, scales = 'free_y' )
c1_only.sgrna.fn <- file.path( plot_dir, 'overlap_c1_only.cl1_vs_flfl.cl1_vs_pFH.png' )
ggsave( c1_only.sgrna.fn, c1_only.sgrna.plot, device = 'png', width = 16, height = 24 )
print( paste( "c1_only overlapping gene barplot written to:", c1_only.sgrna.fn ) )

