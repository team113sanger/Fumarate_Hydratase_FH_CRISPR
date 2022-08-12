#!/usr/bin/perl

# This script is to convert human gene symbols to mouse gene symbols for pyCRISPRcleanR and downstream analyses

use strict;
use warnings FATAL => 'all';

# Set repository path
my $top_dir = $ARGV[0]; 

# Set reference paths based on repository structure
my $gene_list_file = "$top_dir/library/hgnc_to_mgi_symbols.tsv";
my $human_ref_genes = "$top_dir/ref_genes_human";
my $mouse_ref_genes = "$top_dir/ref_genes_mouse";

my $gene_hashref = read_gene_list( $gene_list_file );
opendir( my $REF_DIR, $human_ref_genes ) or die "Cannot open directory ($human_ref_genes): $!";
my @ref_files = readdir $REF_DIR;
closedir $REF_DIR;

print join( "\t", "reference_file", "num_genes", "num_converted_genes", "num_missing_genes" ) . "\n";
foreach my $ref_file ( @ref_files ) {
	if ( $ref_file ne "." && $ref_file ne ".." ) {
		my $ref_genes = read_reference_file( $human_ref_genes . "/" . $ref_file);
		my $output_ref_file = $mouse_ref_genes . "/" . $ref_file;
		my ( $converted_genes, $missing_genes ) = convert_reference_genes( $ref_genes, $gene_hashref );
		print join( "\t", $ref_file, scalar( @{ $ref_genes } ), scalar( @{ $converted_genes } ), scalar( @{ $missing_genes } ) ) . "\n";
		write_converted_genes( $converted_genes, $output_ref_file );
	}
}

sub read_gene_list {
	my ( $list_file ) = @_;
	my %gene_list;

	open( my $LIST, '<', $list_file) or die "Cannot open gene list ($list_file): $!\n";
	while( <$LIST> ){
		my $line = $_;
		chomp( $line );
		if ( $line !~ m/HGNC.symbol/ ) {
        	my ( $human, $mouse ) = split( /\t/, $line );
			next if $human eq "";
			if ( exists $gene_list{ $human } ) {
				$gene_list{ $human } = join( ",", $gene_list{ $human }, $human );
			} else {
				$gene_list{ $human } = $mouse;
			}
		}
	}
	close( $LIST );
	return( \%gene_list );
}

sub read_reference_file {
	my ( $reference_file ) = @_;
	my @reference_genes;

	open( my $REF, '<', $reference_file ) or die "Cannot open reference file ($reference_file): $!\n";
	while( <$REF> ){
        my $line = $_;
        chomp( $line );
		push( @reference_genes, $line );
	}
	close( $REF ); 
	
	return \@reference_genes;
}

sub convert_reference_genes {
	my ( $ref_genes, $gene_hashref ) = @_;

	my @missing_genes;
	my @converted_genes;
	foreach my $gene_to_convert ( @{ $ref_genes } ) {
		if ( exists $gene_hashref->{ $gene_to_convert } ) {
			my @genes_to_add = split( /,/, $gene_hashref->{ $gene_to_convert } );
			foreach my $converted_gene ( @genes_to_add ) {
				push( @converted_genes, $converted_gene );
			}
		} else {
			push( @missing_genes, $gene_to_convert );
		}
	}
	return( \@converted_genes, \@missing_genes );
}

sub write_converted_genes {
	my ( $converted_genes, $output_ref_file ) = @_;
	
	open( my $OUTFILE, '>', $output_ref_file ) or die "Cannot open output reference file ($output_ref_file): $!\n";
	foreach my $converted_gene ( @{ $converted_genes } ) {
		print $OUTFILE $converted_gene . "\n";
	}
	close( $OUTFILE );
}
