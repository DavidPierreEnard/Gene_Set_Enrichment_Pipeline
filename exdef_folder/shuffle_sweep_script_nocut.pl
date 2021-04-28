#!/usr/bin/perl

#use warnings;

use strict;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my ($indir,$outdir,$valid_file,$pair,$coord_file,$prefix,$sweep_size_chain,$segcut,$maxi_rank) = @ARGV;

my %used_genes;

my $used_file = $indir."used_genes.txt";

#print $used_file."\n";

#sleep 20;

open(DATA,$used_file);
while(<DATA>){
    chomp $_;
    $used_genes{$_} = "yes";
}
close DATA;

#for(my $f=1;$f<=2;$f++){

my %go;

#open(DATA,"manual_go_october2015");
#while(<DATA>){
    #chomp $_;
    #my @splitter_line = split(" ",$_);
    #$go{$splitter_line[0]} .= $splitter_line[1];
#}
#close DATA;

my %valid;

open(DATA,$valid_file);

while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    if($used_genes{$splitter_line[0]} eq "yes"){

	#print $splitter_line[0]."\n";

	#sleep 1;

	$valid{$splitter_line[0]} = "yes";
    }
}
close DATA;

my %max;

for(my $i=1;$i<22;$i++){
    $max{$i} = 0;
}

my %coords;

open(DATA,$coord_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    my $gene = $splitter_line[0];
    my $chrom = $splitter_line[1];
    my $center = ($splitter_line[2]+$splitter_line[3])/2;
    my $int_center = int($center/10000)*10000;
    if($chrom =~ /^[123456789]/){
	$coords{$chrom."_".$int_center} .= $gene." ";
	if($int_center>$max{$chrom}){
	    $max{$chrom} = $int_center;
	}
    }
}
close DATA;


my %vips;

my %present_vips;

my %genes;

my @inds = ();

my $counter = 0;

for(my $i=1;$i<=22;$i++){

    my $max_ready = $max{$i}/10000;
    
    for(my $j=0;$j<=$max_ready;$j++){
	my $coord = $j*10000;
	my $coord_info = $i."_".$coord;
	if($coords{$coord_info} ne ""){
	    my $gene_chain = $coords{$coord_info};
	    chop $gene_chain;
	    my @splitter_chain = split(" ",$gene_chain);
	    my $sc_sc = scalar(@splitter_chain);
	    for(my $k=0;$k<=$sc_sc-1;$k++){
		my $gene = $splitter_chain[$k];
		if(($valid{$gene} eq "yes")){ #and(not($go{$gene} =~ /002376/))and(not($go{$gene} =~ /00695[25]/))){
		    $counter += 1;
		    $vips{$counter} = $gene;
		    $genes{$gene} = $counter;
		    push(@inds,$counter);
		}
	    }
	}
    }
}
#close DATA;

my %shuffled_genes;
my %shuffled_genes_2;
my %shuffled_genes_3;
my %shuffled_genes_4;
my %shuffled_genes_5;
my %shuffled_genes_6;
my %shuffled_genes_7;
my %shuffled_genes_8;



my @filelist_cut = ();

my @splitter_size_cut = split("_",$sweep_size_chain);

my $sc_ssc = scalar(@splitter_size_cut);

for(my $h=0;$h<=$sc_ssc-1;$h++){
    my $sweep_file = $prefix."_".$splitter_size_cut[$h];
    push(@filelist_cut,$sweep_file);
}

#my @shuffled_inds = shuffle(@inds);
my $dice=rand(10);
my $sc_ov = scalar(@inds);

#my $seg_size = int($sc_ov/$segcut);
my $seg_size = $segcut;

my %ok_cuts;
my %cut_genes;
my %genes_at_cut;

my $counter_cut = 0;
for(my $i=0;$i<=$sc_ov-1;$i++){

    $counter_cut += 1;

    if($counter_cut==$seg_size){

	my $first_ind = $inds[$i];
	my $second_ind = $inds[$i+1];

	my $first_gene = $vips{$first_ind};
	my $second_gene = $vips{$second_ind};

	$ok_cuts{$i} = "yes";
	$cut_genes{$i} = $first_gene." ".$second_gene;
	$genes_at_cut{$first_gene} = "good";
	$genes_at_cut{$second_gene} = "good";

	$counter_cut = 0;
    }
}

for(my $i=0;$i<=$sc_ssc-1;$i++){
    my $sweep_file = $filelist_cut[$i];

    #print $sweep_file."\n";

    open(DATA,$sweep_file);
    while(<DATA>){
	chomp $_;
	my @splitter_line = split(" ",$_);
	if($genes_at_cut{$splitter_line[0]} eq "good"){
	    my $sc_sl = scalar(@splitter_line);
	    for(my $j=1;$j<=$sc_sl-1;$j++){
		if($splitter_line[$j]<=$maxi_rank){
		    $genes_at_cut{$splitter_line[0]} = "bad";
		}
	    }
	}
    }
    close DATA;
}

my $key_cuts;

foreach $key_cuts (sort keys %cut_genes){

    my $gene_pair = $cut_genes{$key_cuts};
    my @splitter_pair = split(" ",$gene_pair);
    my $gene_1 = $splitter_pair[0];
    my $gene_2 = $splitter_pair[1];

    if(($genes_at_cut{$gene_1} eq "bad")and($genes_at_cut{$gene_2} eq "bad")){
	$ok_cuts{$key_cuts} = "no";
    }
}



my @segments = ();
my $seg = "";
my $counter_001 = 0;
for(my $i=0;$i<=$sc_ov-1;$i++){
    
    $counter_001 += 1;
    if($counter_001<$seg_size){
	#$seg .= $inds[$i]." ";
	if($dice<=5){
	    $seg .= $inds[$i]." ";
	}
	if($dice>5){
	    $seg .= $inds[$sc_ov-1-$i]." ";
	}
    }
    if($counter_001==$seg_size){

	if($ok_cuts{$i} eq "yes"){
	    #print $i." ".$ok_cuts{$i}."\n";
	    #$seg .= $inds[$i]." ";
	    if($dice<=5){
		$seg .= $inds[$i]." ";
	    }
	    if($dice>5){
		$seg .= $inds[$sc_ov-1-$i]." ";
	    }
	    push(@segments,$seg);
	    $seg = "";
	    $counter_001 = 0;
	    $dice=rand(10);
	}

	if($ok_cuts{$i} eq "no"){
	    if($dice<=5){
		$seg .= $inds[$i]." ";
	    }
	    if($dice>5){
		$seg .= $inds[$sc_ov-1-$i]." ";
	    }
	    
	    $counter_001 = 0;
	}
    }
}
if($counter_001<$seg_size){
    push(@segments,$seg);
}

my $sc_segg = scalar(@segments);

#print "Number of segments: ".$sc_segg."\n";

if($sc_segg<10){
    die;
}

my @shuffled_inds = ();
my @shuffled_segments = shuffle(@segments);
my $sc_seg = scalar(@segments);
for(my $i=0;$i<=$sc_seg-1;$i++){
    my $chain = $shuffled_segments[$i];
    chop $chain;
    my @splitter_chain = split(" ",$chain);
    my $sc_sc = scalar(@splitter_chain);
    for(my $k=0;$k<=$sc_sc-1;$k++){
	push(@shuffled_inds,$splitter_chain[$k]);
    }
}
#sleep 1;
#my @shuffled_inds2 = shuffle(@inds);
my @shuffled_inds2 = ();
@shuffled_segments = shuffle(@segments);
my $sc_seg = scalar(@segments);
for(my $i=0;$i<=$sc_seg-1;$i++){
    my $chain = $shuffled_segments[$i];
    chop $chain;
    my @splitter_chain = split(" ",$chain);
    my $sc_sc = scalar(@splitter_chain);
    for(my $k=0;$k<=$sc_sc-1;$k++){
	push(@shuffled_inds2,$splitter_chain[$k]);
    }
}
#sleep 1;
#my @shuffled_inds3 = shuffle(@inds);
my @shuffled_inds3 = ();
@shuffled_segments = shuffle(@segments);
my $sc_seg = scalar(@segments);
for(my $i=0;$i<=$sc_seg-1;$i++){
    my $chain = $shuffled_segments[$i];
    chop $chain;
    my @splitter_chain = split(" ",$chain);
    my $sc_sc = scalar(@splitter_chain);
    for(my $k=0;$k<=$sc_sc-1;$k++){
	push(@shuffled_inds3,$splitter_chain[$k]);
    }
}
#sleep 1;
#my @shuffled_inds4 = shuffle(@inds);
my @shuffled_inds4 = ();
@shuffled_segments = shuffle(@segments);
my $sc_seg = scalar(@segments);
for(my $i=0;$i<=$sc_seg-1;$i++){
    my $chain = $shuffled_segments[$i];
    chop $chain;
    my @splitter_chain = split(" ",$chain);
    my $sc_sc = scalar(@splitter_chain);
    for(my $k=0;$k<=$sc_sc-1;$k++){
	push(@shuffled_inds4,$splitter_chain[$k]);
    }
}
#sleep 1;
#my @shuffled_inds5 = shuffle(@inds);
my @shuffled_inds5 = ();
@shuffled_segments = shuffle(@segments);
my $sc_seg = scalar(@segments);
for(my $i=0;$i<=$sc_seg-1;$i++){
    my $chain = $shuffled_segments[$i];
    chop $chain;
    my @splitter_chain = split(" ",$chain);
    my $sc_sc = scalar(@splitter_chain);
    for(my $k=0;$k<=$sc_sc-1;$k++){
	push(@shuffled_inds5,$splitter_chain[$k]);
    }
}
#sleep 1;
#my @shuffled_inds6 = shuffle(@inds);
my @shuffled_inds6 = ();
@shuffled_segments = shuffle(@segments);
my $sc_seg = scalar(@segments);
for(my $i=0;$i<=$sc_seg-1;$i++){
    my $chain = $shuffled_segments[$i];
    chop $chain;
    my @splitter_chain = split(" ",$chain);
    my $sc_sc = scalar(@splitter_chain);
    for(my $k=0;$k<=$sc_sc-1;$k++){
	push(@shuffled_inds6,$splitter_chain[$k]);
    }
}
#sleep 1;
#my @shuffled_inds7 = shuffle(@inds);
my @shuffled_inds7 = ();
@shuffled_segments = shuffle(@segments);
my $sc_seg = scalar(@segments);
for(my $i=0;$i<=$sc_seg-1;$i++){
    my $chain = $shuffled_segments[$i];
    chop $chain;
    my @splitter_chain = split(" ",$chain);
    my $sc_sc = scalar(@splitter_chain);
    for(my $k=0;$k<=$sc_sc-1;$k++){
	push(@shuffled_inds7,$splitter_chain[$k]);
    }
}
#sleep 1;
#my @shuffled_inds8 = shuffle(@inds);
my @shuffled_inds8 = ();
@shuffled_segments = shuffle(@segments);
my $sc_seg = scalar(@segments);
for(my $i=0;$i<=$sc_seg-1;$i++){
    my $chain = $shuffled_segments[$i];
    chop $chain;
    my @splitter_chain = split(" ",$chain);
    my $sc_sc = scalar(@splitter_chain);
    for(my $k=0;$k<=$sc_sc-1;$k++){
	push(@shuffled_inds8,$splitter_chain[$k]);
    }
}
#sleep 1;


my $sc_ind = scalar(@inds);

for(my $i=1;$i<=$sc_ind;$i++){
    my $original_gene = $vips{$i};
    my $shuffled_ind = $shuffled_inds[$i-1];
    my $new_gene = $vips{$shuffled_ind};
    #print $original_gene."\t".$new_gene."\n";
    $shuffled_genes{$original_gene} = $new_gene;
}

for(my $i=1;$i<=$sc_ind;$i++){
    my $original_gene = $vips{$i};
    my $shuffled_ind = $shuffled_inds2[$i-1];
    my $new_gene = $vips{$shuffled_ind};
    #print $original_gene."\t".$new_gene."\n";
    $shuffled_genes_2{$original_gene} = $new_gene;
}

for(my $i=1;$i<=$sc_ind;$i++){
    my $original_gene = $vips{$i};
    my $shuffled_ind = $shuffled_inds3[$i-1];
    my $new_gene = $vips{$shuffled_ind};
    #print $original_gene."\t".$new_gene."\n";
    $shuffled_genes_3{$original_gene} = $new_gene;
}

for(my $i=1;$i<=$sc_ind;$i++){
    my $original_gene = $vips{$i};
    my $shuffled_ind = $shuffled_inds4[$i-1];
    my $new_gene = $vips{$shuffled_ind};
    #print $original_gene."\t".$new_gene."\n";
    $shuffled_genes_4{$original_gene} = $new_gene;
}

for(my $i=1;$i<=$sc_ind;$i++){
    my $original_gene = $vips{$i};
    my $shuffled_ind = $shuffled_inds5[$i-1];
    my $new_gene = $vips{$shuffled_ind};
    #print $original_gene."\t".$new_gene."\n";
    $shuffled_genes_5{$original_gene} = $new_gene;
}

for(my $i=1;$i<=$sc_ind;$i++){
    my $original_gene = $vips{$i};
    my $shuffled_ind = $shuffled_inds6[$i-1];
    my $new_gene = $vips{$shuffled_ind};
    #print $original_gene."\t".$new_gene."\n";
    $shuffled_genes_6{$original_gene} = $new_gene;
}

for(my $i=1;$i<=$sc_ind;$i++){
    my $original_gene = $vips{$i};
    my $shuffled_ind = $shuffled_inds7[$i-1];
    my $new_gene = $vips{$shuffled_ind};
    #print $original_gene."\t".$new_gene."\n";
    $shuffled_genes_7{$original_gene} = $new_gene;
}

for(my $i=1;$i<=$sc_ind;$i++){
    my $original_gene = $vips{$i};
    my $shuffled_ind = $shuffled_inds8[$i-1];
    my $new_gene = $vips{$shuffled_ind};
    #print $original_gene."\t".$new_gene."\n";
    $shuffled_genes_8{$original_gene} = $new_gene;
}


my @filelist = ();

my @splitter_size = split("_",$sweep_size_chain);

my $sc_ss = scalar(@splitter_size);

for(my $h=0;$h<=$sc_ss-1;$h++){
    my $sweep_file = $prefix."_".$splitter_size[$h];
    push(@filelist,$sweep_file);
}

#my @filelist = ("all_ihsfreqabs_ranks_50kb","all_ihsfreqabs_ranks_100kb","all_ihsfreqabs_ranks_200kb","all_ihsfreqabs_ranks_500kb","all_ihsfreqabs_ranks_1000kb");

#my @filelist = ("all_ihsfreqabs_ranks_50kb","all_ihsfreqabs_ranks_100kb","all_ihsfreqabs_ranks_200kb","all_ihsfreqabs_ranks_500kb","all_ihsfreqabs_ranks_1000kb","all_ihsfreqvar_ranks_50kb",
#    "all_ihsfreqvar_ranks_100kb","all_ihsfreqvar_ranks_200kb","all_ihsfreqvar_ranks_500kb","all_ihsfreqvar_ranks_1000kb","all_ihsfreq_ranks_50kb","all_ihsfreq_ranks_100kb",
 #   "all_ihsfreq_ranks_200kb","all_ihsfreq_ranks_500kb","all_ihsfreq_ranks_1000kb","all_SDS1KG_ranks_50kb","all_SDS1KG_ranks_100kb","all_SDS1KG_ranks_200kb",
  #  "all_SDS1KG_ranks_500kb","all_SDS1KG_ranks_1000kb","all_SDSabs1KG_ranks_50kb","all_SDSabs1KG_ranks_100kb","all_SDSabs1KG_ranks_200kb","all_SDSabs1KG_ranks_500kb","all_SDSabs1KG_ranks_1000kb",
   # "all_SDSvar1KG_ranks_50kb","all_SDSvar1KG_ranks_100kb","all_SDSvar1KG_ranks_200kb","all_SDSvar1KG_ranks_500kb","all_SDSvar1KG_ranks_1000kb","all_SFselect_ranks_50kb","all_SFselect_ranks_100kb",
    #"all_SFselect_ranks_200kb","all_SFselect_ranks_500kb","all_SFselect_ranks_1000kb","all_SDSUK10K_ranks_50kb","all_SDSUK10K_ranks_100kb","all_SDSUK10K_ranks_200kb","all_SDSUK10K_ranks_500kb",
    #"all_SDSUK10K_ranks_1000kb","all_SDSabsUK10K_ranks_50kb","all_SDSabsUK10K_ranks_100kb","all_SDSabsUK10K_ranks_200kb","all_SDSabsUK10K_ranks_500kb","all_SDSabsUK10K_ranks_1000kb",
    #"all_SDSvarUK10K_ranks_50kb","all_SDSvarUK10K_ranks_100kb","all_SDSvarUK10K_ranks_200kb","all_SDSvarUK10K_ranks_500kb","all_SDSvarUK10K_ranks_1000kb","all_SFselect_ranks_50kb",
    #"all_SFselect_ranks_100kb","all_SFselect_ranks_200kb","all_SFselect_ranks_500kb","all_SFselect_ranks_1000kb");

my $sc_file = scalar(@filelist);

for(my $f=0;$f<=$sc_file-1;$f++){

    my $file = $filelist[$f];
    
    my $new_file_1 = $outdir."fake".$file."_1";
    my $new_file_2 = $outdir."fake".$file."_2";
    my $new_file_3 = $outdir."fake".$file."_3";
    my $new_file_4 = $outdir."fake".$file."_4";
    my $new_file_5 = $outdir."fake".$file."_5";
    my $new_file_6 = $outdir."fake".$file."_6";
    my $new_file_7 = $outdir."fake".$file."_7";
    my $new_file_8 = $outdir."fake".$file."_8";
    
    #print $new_file_1."\n";
    
    my %sweeps;
    
    #open OUT, ">".$new_file;
    
    open(DATA,$file);
    while(<DATA>){
	chomp $_;
	my @splitter_line = split(" ",$_);
	my $gene = $splitter_line[0];
	shift @splitter_line;
	$sweeps{$gene} = "@splitter_line";
    }
    close DATA;

    open OUT, ">".$new_file_1."_".$pair;
    my $key_01;
    foreach $key_01 (sort keys %sweeps){
	my $gene = $key_01;
	my $shuffle_gene = $shuffled_genes{$gene};
	my $sweep_line = $sweeps{$shuffle_gene};
	if($sweep_line ne ""){
	    print OUT $gene." ".$sweep_line."\n";
	}
    }
    close OUT;

    open OUT, ">".$new_file_2."_".$pair;
    my $key_01;
    foreach $key_01 (sort keys %sweeps){
	my $gene = $key_01;
	my $shuffle_gene = $shuffled_genes_2{$gene};
	my $sweep_line = $sweeps{$shuffle_gene};
	if($sweep_line ne ""){
	    print OUT $gene." ".$sweep_line."\n";
	}
    }
    close OUT;

    open OUT, ">".$new_file_3."_".$pair;
    my $key_01;
    foreach $key_01 (sort keys %sweeps){
	my $gene = $key_01;
	my $shuffle_gene = $shuffled_genes_3{$gene};
	my $sweep_line = $sweeps{$shuffle_gene};
	if($sweep_line ne ""){
	    print OUT $gene." ".$sweep_line."\n";
	}
    }
    close OUT;
    
    open OUT, ">".$new_file_4."_".$pair;
    my $key_01;
    foreach $key_01 (sort keys %sweeps){
	my $gene = $key_01;
	my $shuffle_gene = $shuffled_genes_4{$gene};
	my $sweep_line = $sweeps{$shuffle_gene};
	if($sweep_line ne ""){
	    print OUT $gene." ".$sweep_line."\n";
	}
    }
    close OUT;

    open OUT, ">".$new_file_5."_".$pair;
    my $key_01;
    foreach $key_01 (sort keys %sweeps){
	my $gene = $key_01;
	my $shuffle_gene = $shuffled_genes_5{$gene};
	my $sweep_line = $sweeps{$shuffle_gene};
	if($sweep_line ne ""){
	    print OUT $gene." ".$sweep_line."\n";
	}
    }
    close OUT;

    open OUT, ">".$new_file_6."_".$pair;
    my $key_01;
    foreach $key_01 (sort keys %sweeps){
	my $gene = $key_01;
	my $shuffle_gene = $shuffled_genes_6{$gene};
	my $sweep_line = $sweeps{$shuffle_gene};
	if($sweep_line ne ""){
	    print OUT $gene." ".$sweep_line."\n";
	}
    }
    close OUT;

    open OUT, ">".$new_file_7."_".$pair;
    my $key_01;
    foreach $key_01 (sort keys %sweeps){
	my $gene = $key_01;
	my $shuffle_gene = $shuffled_genes_7{$gene};
	my $sweep_line = $sweeps{$shuffle_gene};
	if($sweep_line ne ""){
	    print OUT $gene." ".$sweep_line."\n";
	}
    }
    close OUT;

    open OUT, ">".$new_file_8."_".$pair;
    my $key_01;
    foreach $key_01 (sort keys %sweeps){
	my $gene = $key_01;
	my $shuffle_gene = $shuffled_genes_8{$gene};
	my $sweep_line = $sweeps{$shuffle_gene};
	if($sweep_line ne ""){
	    print OUT $gene." ".$sweep_line."\n";
	}
    }
    close OUT;
}

#}
