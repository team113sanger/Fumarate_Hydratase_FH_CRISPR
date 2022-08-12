suppressPackageStartupMessages( library( tidyverse ) )
suppressPackageStartupMessages( library( VennDiagram ) )
suppressPackageStartupMessages( library( RColorBrewer ) )
suppressPackageStartupMessages( library( optparse ) )

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
table_dir <- file.path( analysis_dir, 'summary_tables' )
plot_dir <- file.path( analysis_dir, 'summary_plots' )

# Contrast list
experimental_contrasts <- c( 'cl1_vs_plasmid', 'flfl_vs_plasmid', 'pFH_vs_plasmid', 'cl1_vs_flfl', 'cl1_vs_pFH', 'pFH_vs_flfl' )

# Read in contrast gene results
print( "Reading in MAGeCK gene results..." )
gene.results <- list()
gene.filenames <- list.files( file.path( results_dir, experimental_contrasts ), pattern = '*.gene_summary.txt' )
for ( gfn in gene.filenames ) {
  contrast <- str_replace( gfn, '.gene_summary.txt', '' )
  gene.results[[ contrast ]] <- read.table( file.path( results_dir, contrast, gfn ), sep = "\t", header = T )
}

# Combine gene results into a single table
print( "Combining MAGeCK gene results..." )
gene.results.df <- data.frame()
for ( cn in names( gene.results ) ) {
  tmp <- gene.results[[ cn ]]
  tmp <- tmp %>% mutate( is_enriched = ifelse( pos.fdr < 0.1, TRUE, FALSE ) )
  tmp <- tmp %>% mutate( is_depleted = ifelse( neg.fdr < 0.1, TRUE, FALSE ) )
  tmp <- tmp %>% rename_at( vars( -id ), ~ paste( cn, ., sep = "." ) )
  if ( nrow( gene.results.df ) == 0 ) {
    gene.results.df <- tmp
  } else {
    gene.results.df <- gene.results.df %>% left_join( tmp, by = 'id' )
  }
}
gene.results.df.fn <- file.path( table_dir, 'all_contrasts.gene.tsv' )
write.table( gene.results.df, gene.results.df.fn, sep = "\t", row.names = F, quote = F )

# Summarise gene results per contrast/direction
print( "Summarising MAGeCK gene results..." )
gene.results.summary <- gene.results.df %>% 
                          select( id, ends_with( 'is_enriched' ), ends_with( 'is_depleted' ) ) %>%
                          gather( cname, is_significant, -id ) %>%
                          separate( cname, sep = '\\.is_', into = c( 'contrast', 'direction' ) ) %>%
                          filter( is_significant ) %>%
                          group_by( contrast, direction ) %>%
                          summarise( n_genes = n() ) %>%
                          spread( direction, n_genes )
gene.results.summary.fn <- file.path( table_dir, 'number_of_sig_genes_per_contrast.tsv' )
write.table( gene.results.summary, gene.results.summary.fn, sep = "\t", row.names = F, quote = F )

# Gene venn diagrams
venn.palette <- brewer.pal( 8, "Dark2" )
plot.new()
vennd.enriched <- draw.triple.venn( area1 = gene.results.df %>% filter( cl1_vs_flfl.is_enriched ) %>% nrow(),
                                    area2 = gene.results.df %>% filter( cl1_vs_pFH.is_enriched ) %>% nrow(),
                                    area3 = gene.results.df %>% filter( cl1_vs_plasmid.is_enriched ) %>% nrow(),
                                    n12 = gene.results.df %>% filter( cl1_vs_flfl.is_enriched & cl1_vs_pFH.is_enriched ) %>% nrow(),
                                    n13 = gene.results.df %>% filter( cl1_vs_flfl.is_enriched & cl1_vs_plasmid.is_enriched ) %>% nrow(),
                                    n23 = gene.results.df %>% filter( cl1_vs_pFH.is_enriched & cl1_vs_plasmid.is_enriched ) %>% nrow(),
                                    n123 = gene.results.df %>% filter( cl1_vs_flfl.is_enriched & cl1_vs_pFH.is_enriched & cl1_vs_plasmid.is_enriched ) %>% nrow(),
                                    category = c( 'cl1 vs flfl', 'cl1 vs pFH', 'cl1 vs plasmid' ),
                                    fill = venn.palette[1:3], alpha = 0.4,
                                    scaled = FALSE,
                                    cex = .6,
                                    fontface = "bold",fontfamily = "sans",
                                    cat.cex = 0.5, cat.fontface = "bold", cat.dist = c( 0.07, 0.07, 0.07 ), cat.default.pos = "outer", cat.fontfamily = "sans",
                                    output = FALSE )
png( file.path( plot_dir, 'cl1_enriched_genes.venn.png' ), width = 1400, height = 1400, res = 400 )
plot.new()
grid.draw( vennd.enriched )
dev.off()

plot.new()
vennd.depleted <- draw.triple.venn( area1 = gene.results.df %>% filter( cl1_vs_flfl.is_depleted ) %>% nrow(),
                                    area2 = gene.results.df %>% filter( cl1_vs_pFH.is_depleted ) %>% nrow(),
                                    area3 = gene.results.df %>% filter( cl1_vs_plasmid.is_depleted ) %>% nrow(),
                                    n12 = gene.results.df %>% filter( cl1_vs_flfl.is_depleted & cl1_vs_pFH.is_depleted ) %>% nrow(),
                                    n13 = gene.results.df %>% filter( cl1_vs_flfl.is_depleted & cl1_vs_plasmid.is_depleted ) %>% nrow(),
                                    n23 = gene.results.df %>% filter( cl1_vs_pFH.is_depleted & cl1_vs_plasmid.is_depleted ) %>% nrow(),
                                    n123 = gene.results.df %>% filter( cl1_vs_flfl.is_depleted & cl1_vs_pFH.is_depleted & cl1_vs_plasmid.is_depleted ) %>% nrow(),
                                    category = c( 'cl1 vs flfl', 'cl1 vs pFH', 'cl1 vs plasmid' ),
                                    fill = venn.palette[1:3], alpha = 0.4,
                                    scaled = FALSE,
                                    cex = .6,
                                    fontface = "bold",fontfamily = "sans",
                                    cat.cex = 0.5, cat.fontface = "bold", cat.dist = c( 0.04, 0.04, 0.04 ), cat.default.pos = "outer", cat.fontfamily = "sans",
                                    output = FALSE )
png( file.path( plot_dir, 'cl1_depleted_genes.venn.png' ), width = 1400, height = 1400, res = 400 )
plot.new()
grid.draw( vennd.depleted )
dev.off()

# Read in contrast sgRNA results
print( "Reading in MAGeCK sgRNA results..." )
sgrna.results <- list()
sgrna.filenames <- list.files( file.path( results_dir, experimental_contrasts ), pattern = '*.sgrna_summary.txt' )
for ( sfn in sgrna.filenames ) {
  contrast <- str_replace( sfn, '.sgrna_summary.txt', '' )
  sgrna.results[[ contrast ]] <- read_delim( file.path( results_dir, contrast, sfn ), delim = "\t" )
}

# Combine gene results into a single table
print( "Combining MAGeCK sgRNA results..." )
sgrna.results.df <- data.frame()
for ( cn in names( sgrna.results ) ) {
  tmp <- sgrna.results[[ cn ]]
  tmp <- tmp %>% mutate( is_enriched = ifelse( FDR < 0.1 & high_in_treatment, TRUE, FALSE ) )
  tmp <- tmp %>% mutate( is_depleted = ifelse( FDR < 0.1 & ! high_in_treatment, TRUE, FALSE ) )
  tmp <- tmp %>% rename_at( vars( -sgrna, -Gene ), ~ paste( cn, ., sep = "." ) )
  if ( nrow( sgrna.results.df ) == 0 ) {
    sgrna.results.df <- tmp
  } else {
    sgrna.results.df <- sgrna.results.df %>% left_join( tmp, by = c( 'sgrna', 'Gene' ) )
  }
}
sgrna.results.df.fn <- file.path( table_dir, 'all_contrasts.sgrna.tsv' )
write.table( sgrna.results.df, sgrna.results.df.fn, sep = "\t", row.names = F, quote = F )

# Summarise sgrna results per contrast/direction
print( "Summarising MAGeCK sgRNA results..." )
sgrna.results.summary <- sgrna.results.df %>% 
                          select( sgrna, Gene, ends_with( 'is_enriched' ), ends_with( 'is_depleted' ) ) %>%
                          gather( cname, is_significant, -sgrna, -Gene ) %>%
                          separate( cname, sep = '\\.is_', into = c( 'contrast', 'direction' ) ) %>%
                          filter( is_significant ) %>%
                          group_by( contrast, direction ) %>%
                          summarise( n_sgrnas = n() ) %>%
                          spread( direction, n_sgrnas )
sgrna.results.summary.fn <- file.path( table_dir, 'number_of_sig_sgrnas_per_contrast.tsv' )
write.table( sgrna.results.summary, sgrna.results.summary.fn, sep = "\t", row.names = F, quote = F )




