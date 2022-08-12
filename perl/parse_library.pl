#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

# Set repository path
my $top_dir = $ARGV[0]; 

# Set reference paths based on repository structure
my $lib_directory = $top_dir . "/library/";
my $lib_to_parse = $lib_directory . "yusa-crispr-knockout-mouse-v2.tsv";
my $plasmid_count_file = $lib_directory . "yusa-crispr-knockout-mouse-v2-plasmid.tsv";
my $cgp_lib_file = $lib_directory . "yusa-crispr-knockout-mouse-v2-cgp.tsv";
my $crisprcleanr_lib_file = $lib_directory . "yusa-crispr-knockout-mouse-v2-crisprcleanr.tsv";
my $complete_lib_file = $lib_directory . "yusa-crispr-knockout-mouse-v2-complete.tsv";

my $parsed_lib = parse_library( $lib_to_parse );
my $lib_size = scalar( keys %{ $parsed_lib } );

write_plasmid_counts( $parsed_lib, $plasmid_count_file );
write_cgp_library( $parsed_lib, $cgp_lib_file  );
write_crisprcleanr_library( $parsed_lib, $crisprcleanr_lib_file );
write_complete_library( $parsed_lib, $complete_lib_file  );

print "Library size: " . $lib_size . "\n";

sub parse_library {
    my ($lib_file) = @_;

    my %lib_data;

    if($lib_file){
        open ( my $LIB, '<', $lib_file ) or die 'Failed to open library ('. $lib_file . '):' . $! . "\n";

          while(<$LIB>){
            my $line = $_;
            chomp( $line );
			if ( $line !~ m/gRNA_ID/ ) {
	            my @data = split( /\t/, $line );
    	        my $id = $data[0];
        	    die "Duplicate IDs in library\n" if ( exists $lib_data{ $id } );
	
    	        $lib_data{ $id }{ 'gene' } = $data[1];
        	    $lib_data{ $id }{ 'seq' } = $data[2];
            	$lib_data{ $id }{ 'plasmid_count' } = $data[3];

            	my @annotations = split( /_/, $id );
				my @coords = split( /[:-]/, $annotations[3] );
        	   	$lib_data{ $id }{ 'chr' } = $coords[0];
    	       	$lib_data{ $id }{ 'start' } = $coords[1];
	           	$lib_data{ $id }{ 'end' } = $coords[2];
			}
          }

          close( $LIB );
    }
    return \%lib_data;
}

sub write_plasmid_counts {
    my ( $parsed_lib, $plasmid_count_file ) = @_;

    open ( my $PLASMID, '>', $plasmid_count_file ) or die 'Failed to open plasmid count file ('. $plasmid_count_file . '):' . $! . "\n";
    print $PLASMID join( "\t", 'sgRNA', 'gene', 'plasmid' ) . "\n";
    foreach my $id ( keys %{ $parsed_lib } ) {
        print $PLASMID join( "\t",     $id, 
                                    $parsed_lib->{ $id }->{ 'gene' },
                                    $parsed_lib->{ $id }->{ 'plasmid_count' } ) . "\n";
    }
    close( $PLASMID );
}

sub write_cgp_library {
    my ( $parsed_lib, $cgp_lib_file ) = @_;

    open ( my $CGP, '>', $cgp_lib_file ) or die 'Failed to open CGP library file ('. $cgp_lib_file . '):' . $! . "\n";
#    print $CGP join( ",", 'sgRNA', 'gene', 'sequence' ) . "\n";
    foreach my $id ( keys %{ $parsed_lib } ) {
        print $CGP join( ",",     $id, 
                                    $parsed_lib->{ $id }->{ 'gene' },
                                    $parsed_lib->{ $id }->{ 'seq' } ) . "\n";
    }
    close( $CGP );
}

sub write_crisprcleanr_library {
    my ( $parsed_lib, $crisprcleanr_lib_file ) = @_;

    open ( my $CCR, '>', $crisprcleanr_lib_file ) or die 'Failed to open CRISPRcleanR library file ('. $crisprcleanr_lib_file . '):' . $! . "\n";
    print $CCR join( "\t", 'sgRNA', 'gene', 'chr', 'start', 'end' ) . "\n";
    foreach my $id ( keys %{ $parsed_lib } ) {
        print $CCR join( "\t",     $id, 
                                    $parsed_lib->{ $id }->{ 'gene' },
                                    $parsed_lib->{ $id }->{ 'chr' },
                                    $parsed_lib->{ $id }->{ 'start' },
                                    $parsed_lib->{ $id }->{ 'end' } ) . "\n";
    }
    close( $CCR );
}

sub write_complete_library {
    my ( $parsed_lib, $complete_lib_file ) = @_;

    open ( my $COMPLETE, '>', $complete_lib_file ) or die 'Failed to open complete library file ('. $complete_lib_file . '):' . $! . "\n";
    print $COMPLETE join( "\t", 'sgRNA', 'gene', 'chr', 'start', 'end', 'sequence', 'plasmid' ) . "\n";
    foreach my $id ( keys %{ $parsed_lib } ) {
        print $COMPLETE join( "\t",     $id,
                                    $parsed_lib->{ $id }->{ 'gene' },
                                    $parsed_lib->{ $id }->{ 'chr' },
                                    $parsed_lib->{ $id }->{ 'start' },
                                    $parsed_lib->{ $id }->{ 'end' },
									$parsed_lib->{ $id }->{ 'seq' },
									$parsed_lib->{ $id }->{ 'plasmid_count' }	) . "\n";
    }
    close( $COMPLETE );
}
