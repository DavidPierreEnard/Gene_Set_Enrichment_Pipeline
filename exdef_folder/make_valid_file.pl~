#!/usr/bin/perl

#use warnings;

use strict;

my ($input_files,$out_file) = @ARGV;

my %input_files;

my $file_number = 0;

open(DATA,$input_files);
while(<DATA>){
    chomp $_;
    if(not($_ eq "")){
	my @splitter_line = split(" ",$_);
	$input_files{$splitter_line[0]} = $splitter_line[1]." ".$splitter_line[2]." ".$splitter_line[3];
	$file_number += 1;
    }
}
close DATA;

my %gene_number;

my $key_01;
foreach $key_01 (sort keys %input_files){
    my $chain = $input_files{$key_01};
    my @splitter_chain = split(" ",$chain);
    my $bounds = $splitter_chain[0];
    my $inf = $splitter_chain[1];
    my $sup = $splitter_chain[2];

    open(DATA,$key_01);
    while(<DATA>){
	chomp $_;
	if($bounds eq "no"){
	    my @splitter_line = split(" ",$_);
	    $gene_number{$splitter_line[0]} += 1;
	}
	if($bounds eq "yes"){
            my @splitter_line = split(" ",$_);
	    my $value = $splitter_line[1];
	    if(($value>=$inf)and($value<=$sup)){
		$gene_number{$splitter_line[0]} += 1;
	    }
        }
    }
    close DATA;
}

open OUT, ">".$out_file;
my $counter = 0;
my $key_02;
foreach $key_02 (sort keys %gene_number){
    #print $gene_number{$key_02}." ".$file_number."\n";
    if($gene_number{$key_02}==$file_number){
	print OUT $key_02."\n";
	$counter += 1;
    }
}
close OUT;

print "There are ".$counter." valid genes."."\n";
