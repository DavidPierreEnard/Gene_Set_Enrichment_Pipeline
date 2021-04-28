#!/usr/bin/perl

#use warnings;

use strict;

my ($indir,$keyword) = @ARGV;

while(<$indir/fake*$keyword*>){

    my $file = $_;

    my $line = "rm -f ".$file."\n";

    #print $file."\n";

    system($line);

}
