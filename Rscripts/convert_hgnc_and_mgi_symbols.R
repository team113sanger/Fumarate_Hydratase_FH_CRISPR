require( "biomaRt" )
require( "dplyr" )

getMouseGeneList <- function( x ) {
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mouse_symbols = getBM( attributes = "mgi_symbol", mart = mouse )
  return( mouse_symbols )
}

getHumanGeneList <- function( x ) {
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  human_symbols = getBM( attributes = "hgnc_symbol", mart = human )
  return( human_symbols )
}

convertMouseToHuman <- function( x ) {
  human = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
  mouse = useMart( "ensembl", dataset = "mmusculus_gene_ensembl" )

  genes = getLDS( attributes = c( "mgi_symbol" ), 
                  #filters = "mgi_symbol", 
                  #values = x, 
                  mart = mouse, 
                  attributesL = c( "hgnc_symbol" ), 
                  martL = human, uniqueRows = T )
  return( genes )
}

convertHumanToMouse <- function( x ) {
  human = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
  mouse = useMart( "ensembl", dataset = "mmusculus_gene_ensembl" )
  
  genes  = getLDS( attributes = c( "hgnc_symbol" ), 
                   #filters = "hgnc_symbol", 
                   #values = x , 
                   mart = human, 
                   attributesL = c( "mgi_symbol" ), 
                   martL = mouse, uniqueRows = T )
  return( genes )
}

mouse_genes <- getMouseGeneList()
mouse_to_human_genes <- convertMouseToHuman( mouse_genes )
human_to_mouse_genes <- convertHumanToMouse( human_genes )

full_gene_list <- human_to_mouse_genes %>% full_join( mouse_to_human_genes, by = c( 'HGNC.symbol', 'MGI.symbol' ) )

write.table( mouse_to_human_genes, "mgi_to_hgnc_symbols.tsv", sep = "\t", quote = F, row.names = F )
write.table( human_to_mouse_genes, "hgnc_to_mgi_symbols.tsv", sep = "\t", quote = F, row.names = F )
write.table( full_gene_list, "all_hgnc_and_mgi_symbols.tsv", sep = "\t", quote = F, row.names = F )
