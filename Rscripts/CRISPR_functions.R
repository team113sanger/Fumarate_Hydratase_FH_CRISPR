require( tidyverse )
require( viridis )
require( hrbrthemes )
require( ggrepel )

# Function takes the neg.FDR based on the values of pos.goodsgrna and neg.sgrna
# Adds labels for a gene signature based on options given (i.e. BAGEL essential gene list)
plot_mageck_gene_volcano_with_signature <- function( gene_summary, sig_gene_list, fdr_cutoff = 0.1, lfc_cutoff = 0.5, hit_col = "red", miss_col = "grey", hit_label = "BAGEL essential gene", miss_label = "Not BAGEL essential gene", plot_path = '.', plot_prefix = NULL, save = FALSE, plot_width = 14, plot_height = 12, plot_device = 'png', plot_ext = 'png'  ) {
  
  gene_summary <- gene_summary %>%
    mutate( 'LFC' = neg.lfc,
            'FDR'= ifelse(  ( max( c( neg.goodsgrna, pos.goodsgrna ) ) > 1 ) & 
                              ( neg.goodsgrna > pos.goodsgrna ), neg.fdr, pos.fdr ),
            'is_hit' = ifelse( id %in% sig_gene_list, hit_label, miss_label ) ) 
  
  mageck_gene_volcano <- ggplot( gene_summary, aes( x = LFC, y = -log10( FDR ) ) ) +
    geom_point( data = subset( gene_summary %>% filter( is_hit == miss_label ) ), 
                color = miss_col, alpha = 0.2, size = 0.7 ) +
    geom_point( data = subset( gene_summary %>% filter( is_hit == hit_label ) ), 
                color = hit_col, alpha = 0.8, size = 0.7 ) +
    geom_hline( yintercept = -log10( fdr_cutoff ), linetype = "dotted" ) +
    geom_vline( xintercept = c( -lfc_cutoff, lfc_cutoff ), linetype = "dotted" ) +
    theme_ipsum() +
    ylab( '-log10 FDR' ) +
    xlab( 'log fold change' )
  
  if ( isTRUE( save ) ) {
    plot_filename <- 'mageck_gene_volcano_with_signature'
    if ( !is.null( plot_prefix ) ) {
      plot_filename <- paste( plot_prefix, plot_filename, plot_ext, sep = "." )
    }
    ggsave( plot = mageck_gene_volcano, 
            path = plot_path, filename = plot_filename,
            device = plot_device, width = plot_width, height = plot_height ) 
  }
  return( mageck_gene_volcano )
}

# Function takes the neg.FDR based on the values of pos.goodsgrna and neg.sgrna
# Adds labels for top x enriched (pos.rank) or depleted (neg.rank) genes
plot_mageck_gene_volcano <- function( gene_summary, top_x = 10, fdr_cutoff = 0.1, lfc_cutoff = 0.5, pos_col = "darkseagreen", neg_color = "coral3", nonsig_col = "grey", plot_path = '.', plot_prefix = NULL, save = FALSE, plot_width = 14, plot_height = 12, plot_device = 'png', plot_ext = 'png' ) {
  gene_summary <- gene_summary %>%
                    mutate( 'LFC' = neg.lfc,
                            'FDR'= ifelse(  ( max( c( neg.goodsgrna, pos.goodsgrna ) ) > 1 ) & 
                                            ( neg.goodsgrna > pos.goodsgrna ), neg.fdr, pos.fdr ) ) 
  
  mageck_gene_volcano <- ggplot( gene_summary, aes( x = LFC, y = -log10( FDR ) ) ) +
                            geom_point( data = subset( gene_summary, FDR >= fdr_cutoff ), 
                                        color = nonsig_col, alpha = 0.2, size = 0.7 ) +
                            geom_point( data = subset( gene_summary, FDR < fdr_cutoff & LFC < -lfc_cutoff ), 
                                        color = neg_color, alpha = 0.8, size = 0.7 ) +
                            geom_label_repel( data = subset( gene_summary, 
                                                             FDR < fdr_cutoff & LFC < -lfc_cutoff & neg.rank < ( top_x + 1 ) ), 
                                             aes( label = id ), fontface = 'bold', fill = neg_color, 
                                             force = 2, vjust = -1, label.size = 0.1 ) +
                            geom_point( data = subset( gene_summary, FDR < fdr_cutoff & LFC > lfc_cutoff ), 
                                        color = pos_col, alpha = 0.8, size = 0.7 ) + 
                            geom_label_repel( data = subset( gene_summary, 
                                                             FDR < fdr_cutoff & LFC > lfc_cutoff & pos.rank < ( top_x + 1 ) ), 
                                              aes(label = id), fontface = 'bold', fill = pos_col, 
                                              force = 2, vjust = -1, label.size = 0.1) +
                            geom_hline( yintercept = -log10( fdr_cutoff ), linetype = "dotted" ) +
                            geom_vline( xintercept = c( -lfc_cutoff, lfc_cutoff ), linetype = "dotted" ) +
                            theme_ipsum() +
                            ylab( '-log10 FDR' ) +
                            xlab( 'log fold change' )
  
  if ( isTRUE( save ) ) {
    plot_filename <- 'mageck_gene_volcano'
    if ( !is.null( plot_prefix ) ) {
      plot_filename <- paste( plot_prefix, plot_filename, plot_ext, sep = "." )
    }
    ggsave( plot = mageck_gene_volcano, 
            path = plot_path, filename = plot_filename,
            device = plot_device, width = plot_width, height = plot_height ) 
  }
  
  return( mageck_gene_volcano )
}

# Adds labels for top 20 enriched (pos.rank) or depleted (neg.rank) genes
plot_mageck_sgrna_volcano <- function( sgrna_summary, top_x = 10, fdr_cutoff = 0.1, lfc_cutoff = 0.5, pos_col = "darkseagreen", neg_color = "coral3", nonsig_col = "grey", plot_path = '.', plot_prefix = NULL, save = FALSE, plot_width = 14, plot_height = 12, plot_device = 'png', plot_ext = 'png' ) {
  top_x_pos_sgrna <- sgrna_summary %>% 
                          filter( FDR < fdr_cutoff & LFC > lfc_cutoff ) %>% 
                          top_n( top_x, desc( FDR ) ) %>% 
                          dplyr::select( sgrna )
  top_x_neg_sgrna <- sgrna_summary %>% 
                          filter( FDR < fdr_cutoff & LFC < lfc_cutoff ) %>% 
                          top_n( top_x, desc( FDR ) ) %>% 
                          dplyr::select( sgrna )
  
  mageck_sgrna_volcano <- ggplot( sgrna_summary, aes( x = LFC, y = -log10( FDR ) ) ) +
                            geom_point( data = subset( sgrna_summary, FDR >= fdr_cutoff ), 
                                        color = "grey", alpha = 0.2, size = 0.7 ) +
                            geom_point( data = subset( sgrna_summary, FDR < fdr_cutoff & LFC < -lfc_cutoff ), 
                                        color = pos_col, alpha = 0.8, size = 0.7) +
                            geom_label_repel( data = subset( sgrna_summary, FDR < fdr_cutoff & LFC < -lfc_cutoff & 
                                                              sgrna %in% unlist( top_x_neg_sgrna ) ), 
                                              aes( label = Gene ), fontface = 'bold', fill = pos_col, 
                                              force = 2, vjust = -1, label.size = 0.1 ) +
                            geom_point( data = subset( sgrna_summary, FDR < fdr_cutoff & LFC > lfc_cutoff ), 
                                        color = neg_color, alpha = 0.8, size = 0.7 ) + 
                            geom_label_repel( data = subset( sgrna_summary, FDR < fdr_cutoff & LFC > lfc_cutoff & 
                                                               sgrna %in% unlist( top_x_pos_sgrna ) ), 
                                              aes( label = Gene ), fontface = 'bold', fill = neg_color,
                                              force = 2, vjust = -1, label.size = 0.1 ) +
                            geom_hline( yintercept = -log10( fdr_cutoff ), linetype = "dotted" ) +
                            geom_vline( xintercept = c( -lfc_cutoff, lfc_cutoff ), linetype = "dotted" ) +
                            theme_ipsum() +
                            ylab( '-log10 FDR' ) +
                            xlab( 'log fold change' )
  
  if ( isTRUE( save ) ) {
    plot_filename <- 'mageck_sgrna_volcano'
    if ( !is.null( plot_prefix ) ) {
      plot_filename <- paste( plot_prefix, plot_filename, plot_ext, sep = "." )
    }
    ggsave( plot = mageck_sgrna_volcano, 
            path = plot_path, filename = plot_filename,
            device = plot_device, width = plot_width, height = plot_height ) 
  }
  
  return( mageck_sgrna_volcano )  
}

# Barplots all sgRNA LFC for a vector of gene names
plot_mageck_sgrna_lfc_per_gene <- function( sgrna_summary, gene_summary, enriched = TRUE, top_x = 5, ncols = 1, fdr_cutoff = 0.1, plot_path = '.', plot_prefix = NULL, save = FALSE, plot_width = 10, plot_height = 14, plot_device = 'png', plot_ext = 'png', gene_list = NULL, arrange_by = 'pos.rank' ) {

  if ( !is.null( gene_list )  ) {
      top_x_genes <- gene_summary %>% 
                      filter( id %in% gene_list ) %>%
                      arrange( get( arrange_by ) ) %>%
                      dplyr::select( id ) %>%
                      unlist()
  } else {    
    if ( isTRUE( enriched ) ) {
      top_x_genes <- gene_summary %>% 
                        filter( pos.fdr < fdr_cutoff & pos.rank < ( top_x + 1 ) ) %>%
                        arrange( get( arrange_by ) ) %>%
                        dplyr::select( id ) %>%
                        unlist()
    } else {
      top_x_genes <- gene_summary %>% 
                        filter( neg.fdr < fdr_cutoff & neg.rank < ( top_x + 1 ) ) %>%
                        arrange( get( arrange_by ) ) %>%
                        dplyr::select( id ) %>%
                        unlist()
    }
  } 
  
  sgrna_summary_subset <- sgrna_summary %>% 
                            filter( Gene %in% top_x_genes ) %>% 
                            mutate( Gene = factor( Gene, levels = top_x_genes ) )
  
  
  mageck_sgrna_lfc_per_gene_barplot <-  ggplot( sgrna_summary_subset, 
                                                aes( x = sgrna, y = LFC, fill = Gene ) ) +
                                        geom_bar( stat = 'identity' ) +
                                        facet_wrap( ~ Gene, scales = "free_y", ncol = ncols ) +
                                        coord_flip() +
                                        theme_ipsum() +
                                        xlab( '' ) +
                                        theme( panel.border = element_rect( color = "black", fill = NA, size = 1 ), 
                                               strip.background = element_rect( color = "black", size = 1 ) ) +
                                        scale_fill_viridis( discrete = T, option = 'D', name = '', guide = FALSE )
  if ( isTRUE( save ) ) {
    plot_filename <- 'mageck_sgrna_lfc_per_gene_barplot'
    if ( !is.null( plot_prefix ) ) {
      plot_filename <- paste( plot_prefix, plot_filename, plot_ext, sep = "." )
    }
    ggsave( plot = mageck_sgrna_lfc_per_gene_barplot, 
            path = plot_path, filename = plot_filename,
            device = plot_device, width = plot_width, height = plot_height ) 
  }
  return( mageck_sgrna_lfc_per_gene_barplot )
}

# Merge MAGeCK gene summary results from multiple contrasts and add contrast prefix
merge_mageck_gene_summaries <- function( contrast_list, fdr_cutoff = 0.1, table_path = '.', table_prefix = NULL, save = FALSE ) {
  counter <-  1
  for ( contrast_prefix in names( contrast_list ) ) {
    current_contrast <- contrast_list[[ contrast_prefix ]] 
    current_contrast <- current_contrast %>% mutate( hit = ifelse( neg.fdr < fdr_cutoff | pos.fdr < fdr_cutoff, "TRUE", "FALSE" ) )
    current_contrast <- current_contrast %>% mutate( enriched = ifelse( pos.fdr < fdr_cutoff, "TRUE", "FALSE" ) )
    current_contrast <- current_contrast %>% mutate( depleted = ifelse( neg.fdr < fdr_cutoff, "TRUE", "FALSE" ) )
    current_contrast <- current_contrast %>% rename_at( vars( -id ) , function( x ) paste( contrast_prefix, x, sep = '.' ) ) 
    if ( counter == 1 ) {
      all_contrasts <- current_contrast
    } else {
      all_contrasts <- all_contrasts %>% full_join( current_contrast, by = c( 'id' ) )
    }
    counter <- counter + 1
  }
  
  if ( isTRUE( save ) ) {
    table_filename <- 'mageck_gene_summary.combined'
    if ( !is.null( table_prefix ) ) {
      table_filename <- paste( table_prefix, table_filename, 'tsv', sep = "." )
    }
    
    table_filename = paste( table_path, table_filename, sep = "/")
    write.table( all_contrasts, file = table_filename, row.names = F, quote = F, sep = "\t" )
  }
  
  return( all_contrasts )
}

# Merge MAGeCK sgrna summary results from multiple contrasts and add contrast prefix
merge_mageck_sgrna_summaries <- function( contrast_list, fdr_cutoff = 0.1, table_path = '.', table_prefix = NULL, save = FALSE ) {
  counter <-  1
  for ( contrast_prefix in names( contrast_list ) ) {
    current_contrast <- contrast_list[[ contrast_prefix ]] 
    current_contrast <- current_contrast %>% mutate( hit = ifelse( FDR < fdr_cutoff, "TRUE", "FALSE" ) )
    current_contrast <- current_contrast %>% mutate( enriched = ifelse( FDR < fdr_cutoff & high_in_treatment == "True", "TRUE", "FALSE" ) )
    current_contrast <- current_contrast %>% mutate( depleted = ifelse( FDR < fdr_cutoff & high_in_treatment == "False", "TRUE", "FALSE" ) )
    current_contrast <- current_contrast %>% rename_at( vars( -sgrna, -Gene ) , function( x ) paste( contrast_prefix, x, sep = '.' ) ) 
    if ( counter == 1 ) {
      all_contrasts <- current_contrast
    } else {
      all_contrasts <- all_contrasts %>% full_join( current_contrast, by = c( 'sgrna', 'Gene' ) )
    }
    counter <- counter + 1
  }
  if ( isTRUE( save ) ) {
    table_filename <- 'mageck_sgrna_summary.combined'
    if ( !is.null( table_prefix ) ) {
      table_filename <- paste( table_prefix, table_filename, 'tsv', sep = "." )
    }
    
    table_filename = paste( table_path, table_filename, sep = "/")
    write.table( all_contrasts, file = table_filename, row.names = F, quote = F, sep = "\t" )
  }
  return( all_contrasts )
}

write_mageck_hit_gene_lists <- function( df, contrasts, table_path = '.', table_prefix = NULL, save = FALSE ) {
  for ( contrast in contrasts ) {
    enriched_genes <- df %>% 
                        dplyr::select( id,
                                ( .data[[paste0( contrast, '.enriched')]] ), 
                                ( .data[[paste0( contrast, '.pos.fdr')]] ),
                                ( .data[[paste0( contrast, '.pos.lfc')]] ),
                                ( .data[[paste0( contrast, '.pos.rank')]] ) ) %>%
                        filter( ( .data[[paste0( contrast, '.enriched' )]] ) == "TRUE" ) %>%
                        arrange( .data[[!!paste0( contrast, '.pos.rank')]] )
    
    depleted_genes <- df %>% 
                        dplyr::select( id,
                                ( .data[[paste0( contrast, '.depleted')]] ), 
                                ( .data[[paste0( contrast, '.neg.fdr')]] ),
                                ( .data[[paste0( contrast, '.neg.lfc')]] ),
                                ( .data[[paste0( contrast, '.neg.rank')]] ) ) %>%
                        filter( ( .data[[paste0( contrast, '.depleted' )]] ) == "TRUE" ) %>%
                        arrange( .data[[!!paste0( contrast, '.neg.rank')]] )
    if ( isTRUE( save ) ) {
      enriched_table_filename <- paste( contrast, 'enriched', sep = "_" )
      depleted_table_filename <- paste( contrast, 'depleted', sep = "_" )
      if ( !is.null( table_prefix ) ) {
        enriched_table_filename <- paste( table_prefix, enriched_table_filename, 'tsv', sep = "." )
        depleted_table_filename <- paste( table_prefix, depleted_table_filename, 'tsv', sep = "." )
      }
      
      enriched_table_filename = paste( table_path, enriched_table_filename, sep = "/")
      depleted_table_filename = paste( table_path, depleted_table_filename, sep = "/")
      write.table( enriched_genes, file = enriched_table_filename, row.names = F, quote = F, sep = "\t" )
      write.table( depleted_genes, file = depleted_table_filename, row.names = F, quote = F, sep = "\t" )
    }
  }
}

write_mageck_hit_sgrna_lists <- function( df, contrasts, table_path = '.', table_prefix = NULL, save = FALSE ) {
  for ( contrast in contrasts ) {
    enriched_guides <- df %>% 
      dplyr::select( sgrna, Gene,
              ( .data[[paste0( contrast, '.enriched')]] ), 
              ( .data[[paste0( contrast, '.FDR')]] ),
              ( .data[[paste0( contrast, '.LFC')]] ),
              ( .data[[paste0( contrast, '.high_in_treatment')]] ) ) %>%
      filter( ( .data[[paste0( contrast, '.enriched' )]] ) == "TRUE" ) %>%
      arrange( .data[[!!paste0( contrast, '.FDR')]] )
    
    depleted_guides <- df %>% 
      dplyr::select( sgrna, Gene,
              ( .data[[paste0( contrast, '.depleted')]] ), 
              ( .data[[paste0( contrast, '.FDR')]] ),
              ( .data[[paste0( contrast, '.LFC')]] ),
              ( .data[[paste0( contrast, '.high_in_treatment')]] ) ) %>%
      filter( ( .data[[paste0( contrast, '.depleted' )]] ) == "TRUE" ) %>%
      arrange( .data[[!!paste0( contrast, '.FDR')]] )
    if ( isTRUE( save ) ) {
      enriched_table_filename <- paste( contrast, 'enriched', sep = "_" )
      depleted_table_filename <- paste( contrast, 'depleted', sep = "_" )
      if ( !is.null( table_prefix ) ) {
        enriched_table_filename <- paste( table_prefix, enriched_table_filename, 'tsv', sep = "." )
        depleted_table_filename <- paste( table_prefix, depleted_table_filename, 'tsv', sep = "." )
      }
      
      enriched_table_filename = paste( table_path, enriched_table_filename, sep = "/")
      depleted_table_filename = paste( table_path, depleted_table_filename, sep = "/")
      write.table( enriched_guides, file = enriched_table_filename, row.names = F, quote = F, sep = "\t" )
      write.table( depleted_guides, file = depleted_table_filename, row.names = F, quote = F, sep = "\t" )
    }
  }
}

# Takes a vector of gene symbols and runs Reactome pathway analysis
reactome_pathway_analysis <- function( gene_symbols = NULL, db = 'org.Hs.eg.db', organism = 'human', table_path = '.', table_prefix = NULL, save = FALSE ) {
  geneList <- AnnotationDbi::select( db, 
                                     keys = gene_symbols, 
                                     columns = 'ENTREZID', 
                                     keytype = 'SYMBOL')[,2]
  pathways <- as.data.frame( enrichPathway( gene = geneList , pvalueCutoff = 0.05, readable = T, organism = organism ) )
  
  if ( isTRUE( save ) ) {
    table_filename <- 'reactome_pathway_analysis'
    if ( !is.null( table_prefix ) ) {
      table_filename <- paste( table_prefix, table_filename, 'tsv', sep = "." )
    }
    table_filename = paste( table_path, table_filename, sep = "/")
    write.table( pathways, file = table_filename, row.names = F, quote = F, sep = "\t" )
  }
  
  return( pathways )
}

# Takes a data frame of gene symbols (column 1) and log fold changes (column 2) and runs Reactome GSEA
reactome_gsea <- function( df = NULL, db = 'org.Hs.eg.db', organism = 'human', table_path = '.', table_prefix = NULL, save = FALSE ) {
  geneList = df[,2]
  geneIDs <- as.character( AnnotationDbi::select( db, keys = as.vector( df[,1] ), columns = 'ENTREZID', keytype = 'SYMBOL')[,2] )
  names(geneList ) = geneIDs 
  geneList = sort( geneList, decreasing = TRUE )
  gsea_results <- gsePathway( geneList, nPerm = 10000, pvalueCutoff = 0.05,
                              pAdjustMethod = "BH", verbose = FALSE, organism = organism )
  
  if ( isTRUE( save ) ) {
    table_filename <- 'reactome_gsea'
    if ( !is.null( table_prefix ) ) {
      table_filename <- paste( table_prefix, table_filename, 'tsv', sep = "." )
    }
    table_filename = paste( table_path, table_filename, sep = "/")
    write.table( gsea_results, file = table_filename, row.names = F, quote = F, sep = "\t" )
  }
  
  return( gsea_results )
}

