#!/usr/bin/perl

use warnings;

use strict;

my ($indir) = @ARGV;

while(<$indir/fake*>){

    my $file = $_;

    my $line = "rm -f ".$file."\n";

    #print $file."\n";

    system($line);

}
