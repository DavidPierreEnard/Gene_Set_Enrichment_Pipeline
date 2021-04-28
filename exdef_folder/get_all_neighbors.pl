#!/usr/bin/perl

#use warnings;

use strict;

my ($coord_file,$out_file,$dist) = @ARGV;

my %prot;

open(DATA,"valid_file.txt");
while(<DATA>){
    chomp $_;
    $prot{$_} = "yes";

}
close DATA;

my %coords;
my %rev_coords;

open(DATA,$coord_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    my $gene = $splitter_line[0];
    if($prot{$gene} eq "yes"){
	my $chrom = $splitter_line[1];
	my $center = int(($splitter_line[2]+$splitter_line[3])/2);
	my $new_center = int($center/10000)*10000;
	$coords{$gene} .= $chrom."_".$new_center;
	$rev_coords{$chrom."_".$new_center}.= $gene." ";
    }
}
close DATA;

open OUT, ">".$out_file;

my $key_01;

foreach $key_01 (sort keys %coords){

    #print $key_01."\n";

    my $coords = $coords{$key_01};
    
    my @splitter_coords = split("_",$coords);

    my $chrom = $splitter_coords[0];
    my $coordinate = $splitter_coords[1];

    my $neighbors_chain = "";
    
    for(my $i=0;$i<100000000;$i++){

	my $left_coord = $coordinate-10000*$i;
	my $right_coord = $coordinate+10000*$i;

	if($rev_coords{$chrom."_".$left_coord} ne ""){
	    my $gene_chain = $rev_coords{$chrom."_".$left_coord};
	    my @splitter_gc = split(" ",$gene_chain);
	    my $sc_gc = scalar(@splitter_gc);
	    for(my $j=0;$j<=$sc_gc-1;$j++){
		my $current_gene = $splitter_gc[$j];
		if($current_gene ne $key_01){
		    $neighbors_chain .= $current_gene." ";
		}
	    }
	}

	if($rev_coords{$chrom."_".$right_coord} ne ""){
            my $gene_chain = $rev_coords{$chrom."_".$right_coord};
            my @splitter_gc = split(" ",$gene_chain);
            my $sc_gc = scalar(@splitter_gc);
            for(my $j=0;$j<=$sc_gc-1;$j++){
                my $current_gene = $splitter_gc[$j];
                if($current_gene ne $key_01){
		    $neighbors_chain .= $current_gene." ";
                }
            }
        }

	if($right_coord-$coordinate>=$dist){
	    last;
	}
    }
    
    if($neighbors_chain ne ""){
	chop $neighbors_chain;
    }

    print OUT $key_01."\t".$neighbors_chain."\n";

}

close OUT;
