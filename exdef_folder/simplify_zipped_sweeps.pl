#!/usr/bin/perl

#use warnings;

use strict;

use Compress::Zlib;

my ($in_file,$out_file,$threshold_chain) = @ARGV;

#print $in_file."\n";

#my @cut_list = (2000,1500,1000,900,800,700,600,500,450,400,350,300,250,200,150,100,80,60,50,40,30,25,20,15,10);

my @cut_list = split(/\,/,$threshold_chain);

my $sc_cl = scalar(@cut_list);

open OUT, ">".$out_file;

my $file = gzopen($in_file, "rb");

while ($file->gzreadline($_)>0){
    chomp $_;
    my @splitter_line = split(" ",$_);
    my $gene = $splitter_line[0];
    my $sc_sl = scalar(@splitter_line);
    my $chain = "";
    my %chain_parts;
    for(my $i=1;$i<=$sc_sl-1;$i++){

	my $value = $splitter_line[$i];
	
	my $new_value = 10000;

	for(my $j=0;$j<=$sc_cl-1;$j++){
	    my $cutoff = $cut_list[$j];				   
	    if($splitter_line[$i]<=$cutoff){
		$new_value = $cutoff;
	    }
	}

	if($new_value<10000){
	    $chain_parts{$new_value} .= $i."_";
	}
    }

    my $key_01;
    foreach $key_01 (sort keys %chain_parts){
	chop $chain_parts{$key_01};
	$chain_parts{$key_01} = $key_01."_".$chain_parts{$key_01};
	$chain .= $chain_parts{$key_01}." ";
    }

    if($chain ne ""){

	print OUT $gene." ".$chain."\n";
    }
}
$file->gzclose();

close OUT;
