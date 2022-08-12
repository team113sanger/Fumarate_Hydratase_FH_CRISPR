#!/usr/bin/env Rscript

# This script will parse MAGeCK gene summary files from 2 contrasts,
# eg: plasmid vs treatment (contrast 2) and control vs treatment
# (contrast 1) and find significantly enriched and depleted genes 
# by comparing both FDR and the direction of the logFC in the 2 contrasts.

suppressPackageStartupMessages( library( optparse ) )
suppressPackageStartupMessages( library( tidyverse ) )

option_list<-list(
	make_option( c( "-t", "--contrast1" ), type = "character", default = NULL,
               help = "MAGeCK gene summary file for treatment vs control", metavar = "character" ),
	make_option( c( "-c", "--contrast2" ), type = "character", default = NULL,
               help = "MAGeCK gene summary file for treatment vs plasmid", metavar = "character" ),
	make_option( c( "-f", "--fdr" ), type = "numeric", default = 0.1,
               help = "FDR threshold to be applied to neg.fdr and pos.fdr [default = %default]", metavar = "numeric" ),
	make_option( c( "-m", "--main" ), type = "character", default = NULL,
               help = "Comparison for plot title, eg: d30 vs d30C", metavar = "character" ),
	make_option( c( "-p", "--prefix" ), type = "character", default = NULL,
	             help = "Output file prefix (optional)", metavar = "character" )
)

opt_parser<-OptionParser( option_list = option_list) ;
opt<-parse_args( opt_parser );

### Functions
label_sig_genes <- function( c1, c2, fdr = 0.1 ) {
  c1 <- c1 %>% rename_at( vars( -id ), function( x ) paste( 'tvsc' , x, sep = '.' ) )
  c2 <- c2 %>% rename_at( vars( -id ), function( x ) paste( 'tvsp' , x, sep = '.' ) )
  all <- c1 %>% full_join( c2, by = 'id' )
  all <- all %>% 
    mutate( 'tvsc.enriched' =  ifelse( tvsc.pos.fdr <= fdr_cutoff & tvsc.pos.lfc > 0, TRUE, FALSE ) ) %>% 
    mutate( 'tvsc.depleted' =  ifelse( tvsc.neg.fdr <= fdr_cutoff & tvsc.neg.lfc < 0, TRUE, FALSE ) ) %>% 
    mutate( 'tvsp.enriched' =  ifelse( tvsp.pos.fdr <= fdr_cutoff & tvsp.pos.lfc > 0, TRUE, FALSE ) ) %>% 
    mutate( 'tvsp.depleted' =  ifelse( tvsp.neg.fdr <= fdr_cutoff & tvsp.neg.lfc < 0, TRUE, FALSE ) ) %>%
    mutate( 'sig.enriched'  =  ifelse( tvsc.enriched == TRUE & tvsp.enriched == TRUE, TRUE, FALSE ) ) %>%
    mutate( 'sig.depleted'  =  ifelse( tvsc.depleted == TRUE & tvsp.depleted == TRUE, TRUE, FALSE ) ) %>%
    mutate( 'sig'  =  ifelse( sig.enriched == TRUE | sig.depleted == TRUE , TRUE, FALSE ) )
  return( all )
}

write_df_to_file <- function( df, outfile, prefix = NULL ) {
  if ( ! is.null( prefix ) ) {
    outfile <- paste( prefix, outfile, sep='_' )
  }
  write.table( df, outfile, row.names = F, sep = "\t", quote = F )
  return( outfile )
}

get_limits <- function( comparison, which ) {
  # pos.lfc == neg.lfc
  if ( which == 'min' ) {
    if ( comparison == 'tvsp' ) { 
      min <- ifelse( length( which( plot_data$tvsp.depleted ) ) > 0, floor( min( plot_data[ which( plot_data$tvsp.depleted ), 'tvsp.neg.lfc' ] ) ), -2 )
      # Make extra room for labels
      min <- ifelse( min > 5, min - 1, min - 0.5 )
    }
    else if ( comparison == 'tvsc' ) {
      min <- ifelse( length( which( plot_data$tvsc.depleted ) ) > 0, floor( min( plot_data[ which( plot_data$tvsc.depleted ), 'tvsc.neg.lfc' ] ) ), -2 )
      # Make extra room for text
      min <- ifelse( min > 10, min - 1, min - 0.2 )
    }
    return ( min )
  } else if ( which == 'max' ) {
    if ( comparison == 'tvsp' ) {
      max <- ifelse( length( which( plot_data$tvsp.enriched ) ) > 0, ceiling( max( plot_data[ which( plot_data$tvsp.enriched ), 'tvsp.pos.lfc' ] ) ), 2 )
    }
    else if ( comparison == 'tvsc' ) {
      max <- ifelse( length( which( plot_data$tvsc.enriched ) ) > 0, ceiling( max( plot_data[ which( plot_data$tvsc.enriched ), 'tvsc.pos.lfc' ] ) ), 2 )
      # Make extra room for text 
      max <- ifelse( max > 10, max + 1, max + 0.2 )
    }
    return ( max )
  }
}



# Validate input options
if ( is.null( opt$contrast1 ) ) {
  print_help( opt_parser )
  stop( "Please provide a MAGeCK gene summary file for contrast 1 (treatment vs control)", call. = FALSE )
}
if ( is.null( opt$contrast2 ) ) {
  print_help( opt_parser )
  stop( "Please provide a MAGecK gene summary file for contrast 2 (treatment vs plasmid)", call. = FALSE )
}
if ( is.null( opt$main ) ) {
  print_help( opt_parser )
  stop( "Please provide a comparison for the plot title", call. = FALSE )
}

# Rename input variables
tvsc_file <- opt$contrast1
tvsp_file <- opt$contrast2
main <- opt$main
fdr_cutoff <- opt$fdr
output_prefix <- opt$prefix

# Feedback input file info
print( paste( "Contrast file for TvsC (y-axis):", tvsc_file ) )
print( paste( "Contrast file for TvsP (x-axis):", tvsp_file ) )
print( paste( "Using FDR cutoff:", fdr_cutoff ) )

# Read in plasmid vs treatment and control vs treatment gene summary files from MAGeCK
# Note: "|" in MAGeCK headers converted to "." by default since check.names isn't used here
tvsc <- read.table( tvsc_file, header = T, sep = "\t" )
tvsp <- read.table( tvsp_file, header = T, sep = "\t" )

# Merge contrasts, label sig genes
merged_contrast_df <- label_sig_genes( tvsc, tvsp, fdr = fdr_cutoff )

# Write combined contrasts to file
merged_contrast_outfile <- 'mageck_c1c2_combined_results_genes.tsv'
merged_contrast_outfile <- write_df_to_file( merged_contrast_df, merged_contrast_outfile, output_prefix )

print( paste( "Merged gene summary written to:", merged_contrast_outfile ) )

# Get sig depleted genes ( in both tvsc and tvsp: neg.fdr < fdr_cutoff & neg.lfc < 0 )
sig_depleted <- merged_contrast_df %>%
                  filter( sig.depleted == TRUE ) %>%
                  select( id, 
                          'treatment_vs_control.LFC' = tvsc.neg.lfc, 'treatment_vs_control.FDR' = tvsc.neg.fdr,
                          'treatment_vs_plasmid.LFC' = tvsp.neg.lfc, 'treatment_vs_plasmid.FDR' = tvsp.neg.fdr )

n_sig_depleted <- nrow( sig_depleted )

print( paste( "Number of significantly depleted genes found:", n_sig_depleted ) )
print( paste( "Significantly depleted genes:", paste( sig_depleted$id, collapse = "," ) ) )

# Write significantly depleted genes to file
sig_depleted_outfile <- 'mageck_c1c2_sig_depleted_genes.tsv'
sig_depleted_outfile <- write_df_to_file( sig_depleted , sig_depleted_outfile, output_prefix )

print( paste( "Significantly depleted genes written to:", sig_depleted_outfile ) )

# Get sig enriched genes ( in both tvsc and tvsp: pos.fdr < fdr_cutoff & pos.lfc > 0 )
sig_enriched <- merged_contrast_df %>%
                  filter( sig.enriched == TRUE ) %>%
                  select( id, 
                          'treatment_vs_control.LFC' = tvsc.pos.lfc, 'treatment_vs_control.FDR' = tvsc.pos.fdr,
                          'treatment_vs_plasmid.LFC' = tvsp.pos.lfc, 'treatment_vs_plasmid.FDR' = tvsp.pos.fdr )

n_sig_enriched <- nrow( sig_enriched )

print( paste( "Number of significantly enriched genes found:", n_sig_enriched ) )
print( paste( "Significantly enriched genes:", paste( sig_enriched$id, collapse = "," ) ) )

# Write significantly enriched genes to file
sig_enriched_outfile <- 'mageck_c1c2_sig_enriched_genes.tsv'
sig_enriched_outfile <- write_df_to_file( sig_enriched , sig_enriched_outfile, output_prefix )

print( paste( "Significantly enriched genes written to:", sig_enriched_outfile ) )

# Get signicicantly depleted genes ( in tvsc only: neg.fdr < fdr_cutoff & neg.lfc < 0 )
sig_c1_depleted <- merged_contrast_df %>%
                  filter( tvsc.depleted == TRUE & tvsp.depleted == FALSE & tvsp.neg.lfc < 0) %>%
                  select( id, 
                          'treatment_vs_control.LFC' = tvsc.neg.lfc, 'treatment_vs_control.FDR' = tvsc.neg.fdr,
                          'treatment_vs_plasmid.LFC' = tvsp.neg.lfc, 'treatment_vs_plasmid.FDR' = tvsp.neg.fdr )

n_sig_c1_depleted <- nrow( sig_c1_depleted )

print( paste( "Number of significantly depleted genes found in c1 only:", n_sig_c1_depleted ) )
print( paste( "Significantly depleted genes in c1 only:", paste( sig_c1_depleted$id, collapse = "," ) ) )

# Write significantly enriched genes in c1 only to file
sig_c1_depleted_outfile <- 'mageck_c1_only_sig_depleted_genes.tsv'
sig_c1_depleted_outfile <- write_df_to_file( sig_c1_depleted , sig_c1_depleted_outfile, output_prefix )

print( paste( "Significantly depleted genes in c1 only written to:", sig_c1_depleted_outfile ) )

# Get sig enriched genes ( in tvsc only: neg.fdr < fdr_cutoff & neg.lfc < 0 )
sig_c1_enriched <- merged_contrast_df %>%
                  filter( tvsc.enriched == TRUE & tvsp.enriched == FALSE & tvsp.pos.lfc > 0) %>%
                  select( id, 
                          'treatment_vs_control.LFC' = tvsc.pos.lfc, 'treatment_vs_control.FDR' = tvsc.pos.fdr,
                          'treatment_vs_plasmid.LFC' = tvsp.pos.lfc, 'treatment_vs_plasmid.FDR' = tvsp.pos.fdr )

n_sig_c1_enriched <- nrow( sig_c1_enriched )

print( paste( "Number of significantly enriched genes found in c1 only:", n_sig_c1_enriched ) )
print( paste( "Significantly enriched genes in c1 only:", paste( sig_c1_enriched$id, collapse = "," ) ) )

# Write significantly enriched genes in c1 only to file
sig_c1_enriched_outfile <- 'mageck_c1_only_sig_enriched_genes.tsv'
sig_c1_enriched_outfile <- write_df_to_file( sig_c1_enriched , sig_c1_enriched_outfile, output_prefix )
print( paste( "Significantly enriched genes in c1 only written to:", sig_c1_enriched_outfile ) )

# Plot
lfc_plot_title <- paste( "Significantly depleted or enriched genes in", main )
lfc_plot_subtitle <- paste( "Comparison of significant genes in treatment vs. plasmid with treatment vs. control (MAGeCK; FDR <= ", fdr_cutoff, ")", sep='' )

plot_data <- merged_contrast_df %>% filter( tvsc.enriched & ! tvsp.enriched | tvsc.depleted & ! tvsp.depleted |
							   tvsc.enriched & tvsp.enriched | tvsc.depleted & tvsp.depleted )

if( is.null( plot_data ) ) {
  print("There are not significant genes in Treatment vs Control. A plot will not be created.")
  q()
}

# Find the axes min and max
xmax <- get_limits( 'tvsp', 'max' )
xmin <- get_limits( 'tvsp', 'min' )
ymax <- get_limits( 'tvsc', 'max' )
ymin <- get_limits( 'tvsc', 'min' )

# Get the relevant rows
plot_data <- plot_data %>% mutate( direction = ifelse( tvsc.enriched & ! tvsp.enriched | tvsc.enriched & tvsp.enriched , 'enriched',
                       ifelse( tvsc.depleted & ! tvsp.depleted | tvsc.depleted & tvsp.depleted, 'depleted', NA ) ) )

# Specify the levels for legends
plot_data$direction <- factor( plot_data$direction, levels = c( "depleted", "enriched" ) )
plot_data$sig <- factor(plot_data$sig, levels = c("TRUE","FALSE") )
 
# Get colour and alpha based on condition
alpha_value <- ifelse( plot_data$tvsc.enriched & ! plot_data$tvsp.enriched | 
					 plot_data$tvsc.depleted & ! plot_data$tvsp.depleted, 0.5, 1 )

# Adjust x-axis shift depending on range
nudge_x <- ifelse( ymax - ymin > 4, -0.1, -0.05 )

# Plot logFCs - treatment vs plasmid on x-axis and treatment vs control on y-axis
# Color - Significantly enriched or depeleted in treatment vs control
# Shape - significant/not significant in treatment vs plasmid
lfc_plot <- ggplot( plot_data, aes( x = tvsp.neg.lfc, y = tvsc.neg.lfc, shape = sig, col = direction ) ) +
              geom_point( alpha = alpha_value ) +
              ylim( ymin, ymax ) +
              xlim( xmin, xmax ) +
              geom_hline( yintercept = 0 ) +
              geom_vline( xintercept = 0 ) +
              theme_bw() +
              theme( plot.title = element_text( size = 11, hjust = 0.5, face = 'bold' ),
                     plot.subtitle = element_text( size = 8, hjust = 0.5 ),
                     axis.title = element_text( size = 8 ),
                     axis.text = element_text( size = 8 ),
                     legend.position = "bottom", legend.margin = margin (),
                     legend.text = element_text( size = 7 ),
                     legend.title = element_text( size = 8 ) ) +
              ggtitle( label = lfc_plot_title, subtitle = lfc_plot_subtitle ) +
              geom_text( data = plot_data %>% filter( tvsc.enriched & tvsp.enriched ), aes( label = id ), hjust = 1, nudge_x = nudge_x, size = 2 ) +
              geom_text( data = plot_data %>% filter( tvsc.depleted & tvsp.depleted ), aes( label = id ), hjust = 1, nudge_x = nudge_x, size = 2 ) +
              geom_text( x = xmax/2, y = ymax + 0.05, label = "Enriched", color = "red", size = 3 ) +
              geom_text( x = xmin/2, y = ymin - 0.05, label = "Depleted", color = "blue", size = 3 ) +
              xlab( "logFC (Treatment vs Plasmid)" ) +
              ylab( "logFC (Treatment vs Control)" ) +
              scale_color_manual( name = 'Treatment vs Control', values = c( 'depleted' = 'blue', 'enriched' = 'red' ),
                labels = c( 'Significantly depleted', 'Significantly enriched' ) , drop = FALSE ) + 
              scale_shape_manual( name = 'Treatment vs Plasmid', values = c( 15, 16 ),
                labels = c( 'Significant', 'Not significant' ), drop = FALSE ) +
              guides( color = guide_legend( order = 1, ncol = 1 , title.position="top", title.hjust = 0.5 ),
                      shape = guide_legend( order = 2, ncol = 1 , title.position="top", title.hjust = 0.5 ) )
              
              
lfc_plot_outfile <- 'mageck_logFC_plot_genes.png'

if ( ! is.null( output_prefix ) ) {
  lfc_plot_outfile <- paste( output_prefix, lfc_plot_outfile, sep='_' )
}          

ggsave( plot = lfc_plot, filename = lfc_plot_outfile, device = 'png', width = 6, height = 6, dpi = 300 )


