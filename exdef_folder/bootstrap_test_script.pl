#!/usr/bin/perl

#use warnings;

use strict;

my ($indir,$valid_file,$sweeps_file,$hgnc_file,$vip_file,$factors_table_file,$distance_file,$iterations,$range,$out_file_1,$out_file_2,$dist,$flip,$repok) = @ARGV;

#This script does boostrap with multiple target averages and is useful in the following situation. Imagine you have a group of genes of interest, for example 200 genes involved in autism, or 500 genes involved in interactions with a specific virus, etc... For this group of genes of interest, you are for example interested in their speed of evolution, and you would like to know if your genes of interest evolve faster or slower than the rest of all the other genes in the human genome. The problem is that many other factors in the genome influence the speed of evolution other than the category of the group of interest (being an autism gene, or a gene that interacts with a specific type of virus, etc...). In other words, there are many correlated confounding factors that need to be controlled for before being able to conclude whther or not the group of genes of interest evolves faster or slower. Ideally, you would be able to create a predictive generalized linear or semi-parametric additive model from all the genes outside of the group of interest, and then compare the predictions of the model with the actual observations of the speed of evolution within the group of interest. But sometimes the data you are dealing with is just very badly behaved from a statistical point of view and no generalized model is flexible enough to fit the data properly. Alternatively you did build a model but you want to be able to validate it with a fully independent method.

#Then what the script does is that it creates sets of control genes that have:
#1) the same size as the group of genes of interest.
#2) the same average value for multiple confounding factors as the group of genes of interest.

#By doing this the script makes it possible to compare a group of genes of interest with the rest of the genome for a paramaeter of interest while simultaneously controlling for multiple confounding factors. It is coded in way that can accomodate many (I tested up to 10) confounding factors while maximising the number of genes form the rest of the genome that can be used as control. In summary, it is a fully non-parametric test that can be used to test the effect of a parameter of interest when other parametric or semi-parametric models need to be validated or cannot properly fit the data.

#The script uses a number of tricks to build as many control sets of genes as needed while maximising the number of genes that can be used in the boostrap, thus increasing power:

#1) a "fake seed" of control genes that just make it possible to start build the real control.
#2) addition of pairs of genes to the control instead of genes one by one. It can be shown empirically that adding pairs of genes to build the controls maximizes the number of genes that can be used as controls.
#3) a simple algorithm that makes it possible to add genes to the control while always staying close to the target averages for the different confounding factors, even when the group of genes of interest is markedly different from the rest of the genome in terms of confounding factors.

#######################
#######################
#######################
#######################

#Gathering of the data required

my %used_genes;

my %printed;

my %ppi;

my %hgnc; #checks names of genes to later exclude HLA and HISTone genes if not already excluded from the valid_file. Remove this part if you don't care about removing HLA genes.
open(DATA,$indir.$hgnc_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    $hgnc{$splitter_line[0]} = $splitter_line[1];
}
close DATA;

#my %score;
#open(DATA,$indir.$disease_score_file);
#while(<DATA>){
    #chomp $_;
    #if($_ =~ /yes/){
	#my @splitter_line = split(" ",$_);
	#my $score = $splitter_line[2];
	#$score{$splitter_line[0]} = $score;
    #}
#}
#close DATA;

my %valid; #hash of Ensembl gene ID to include in the test.
open(DATA,$indir.$valid_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    if((not($hgnc{$splitter_line[0]} =~ /HIST/))and(not($hgnc{$splitter_line[0]} =~ /HLA/))){
	#and($splitter_line[2]/$splitter_line[3]<=0.08258)){ #lactase, HLA and HISTone genes excluded here.
	$valid{$splitter_line[0]} = "yes";
    }
}
close DATA;

my %sweeps; # gets sweeps information.
open(DATA,$indir.$sweeps_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    $sweeps{$splitter_line[0]} = "yes";
}
close DATA;

my %vir_numb;

my %factors; # gets confounding factors information.
my $factor_number = 0;
open(DATA,$indir.$factors_table_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);

    #$splitter_line[8] = 0;
    #$splitter_line[9] =	0;
    #$splitter_line[10] = 0;
    #$splitter_line[11] = 0;
    #$splitter_line[12] = 0;
    #$splitter_line[13] = 0;
    #$splitter_line[15] = 0;
    #$splitter_line[16] = 0;
    #$splitter_line[17] = 0;
    #$splitter_line[6] = 0;
    #$splitter_line[7] = 0;


    
    #$splitter_line[14] = 0;
    #$splitter_line[7] = 0;
    if($valid{$splitter_line[0]} eq "yes"){
	my $sc_sl = scalar(@splitter_line);
     
	$vir_numb{$splitter_line[0]} = 1;
	my $gene = $splitter_line[0];
	shift @splitter_line;
	#pop @splitter_line;
	$factor_number = scalar(@splitter_line);
	$factors{$gene} = "@splitter_line";
	#$factor_number = 1;
	#$factors{$gene} = $splitter_line[5];
    }
}
close DATA;

my %include; #defines the intersection of sweeps, valid genes, and genes with defined confounding factors, in case it was not already done to generate the valid_file. The intersection is the list of genes that will be used.
my $key_01;
foreach $key_01 (sort keys %valid){
    if(($sweeps{$key_01} eq "yes")and($factors{$key_01} ne "")){
	$include{$key_01} = "yes";
    }
}

my %distance; #gets distance of control genes from genes in group of interest in the case genes close to group of interest need to be excluded from control.

open(DATA,$indir.$distance_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    $distance{$splitter_line[0]} = $splitter_line[1];
}
close DATA;

my %vips;
my %non_vips;
my $vip_number = 0;
my @nonvips = ();

my $nvip_number = 0;

my @included = ();

if($flip eq "no"){

    open(DATA,$indir.$vip_file); # gets information on genes of interest.
    while(<DATA>){
	chomp $_;
	my @splitter_line = split(" ",$_);
	if($include{$splitter_line[0]} eq "yes"){
	    push(@included,$splitter_line[0]);
	    if(($splitter_line[1] eq "yes")){
		$vips{$splitter_line[0]} = $factors{$splitter_line[0]};
		$vip_number += 1;
		$used_genes{$splitter_line[0]} = "yes";
	    }
	    if(($splitter_line[1] eq "no")and($distance{$splitter_line[0]}>=$dist)){
		$non_vips{$splitter_line[0]} = $factors{$splitter_line[0]};
		push(@nonvips,$splitter_line[0]);
		$nvip_number += 1;
	    }
	}
    }
    close DATA;
}

if($flip eq "yes"){

    open(DATA,$indir.$vip_file); # gets information on genes of interest.                                                                                                                                                                                                                                                                                                 
    while(<DATA>){
        chomp $_;
        my @splitter_line = split(" ",$_);
        if($include{$splitter_line[0]} eq "yes"){
            push(@included,$splitter_line[0]);
            if(($splitter_line[1] eq "yes")){
                #$vips{$splitter_line[0]} = $factors{$splitter_line[0]};
                #$vip_number += 1;
		$non_vips{$splitter_line[0]} = $factors{$splitter_line[0]};
                push(@nonvips,$splitter_line[0]);
                $nvip_number += 1;

            }
            if(($splitter_line[1] eq "no")and($distance{$splitter_line[0]}>=$dist)){
                $vips{$splitter_line[0]} = $factors{$splitter_line[0]};
                $vip_number += 1;
            }
        }
    }
    close DATA;


}

print "There are ".$vip_number." genes of interest and ".$nvip_number." potential control genes at distance of at least ".$dist." bases."."\n";

if($nvip_number<=1.5*$vip_number){

    print "Warning! The number of control genes is less than 1.5 times the number of genes of interest. FDR may be high as a result."."\n";

}

#open OUT, ">"."all_nonVIPs_file.txt";
#my $sc_sn = scalar(@nonvips);
#for(my $d=0;$d<=$sc_sn-1;$d++){
    #print OUT $nonvips[$d]."\n";
#}
#close OUT;


my %averages; #measure averages of confounding factors in the group of interest, and calculate tolerable range of the controls as a function of argument $range.
my $den = 0;
open OUT, ">".$out_file_1; #prints file with Ensembl IDs of the genes in the group of interest.
my $key_001;
foreach $key_001 (sort keys %vips){
    print OUT $key_001."\n";
    my @splitter_fact = split(" ",$factors{$key_001});
    $den+=1;
    my $sc_sf = scalar(@splitter_fact);
    for(my $i=0;$i<=$sc_sf-1;$i++){
	my $numb = $i+1;
	$averages{$numb} += $splitter_fact[$i];
    }
}
close OUT;

my %nvips_av;
my $n_den = 0;
my $key_n001;
foreach $key_n001 (sort keys %non_vips){
    my @splitter_fact = split(" ",$factors{$key_n001});
    $n_den+=1;
    my $sc_sf = scalar(@splitter_fact);
    for(my $i=0;$i<=$sc_sf-1;$i++){
        my $numb = $i+1;
        $nvips_av{$numb} += $splitter_fact[$i];
    }
}
my %highlow;
my $key_n002;
foreach $key_n002 (sort keys %nvips_av){
    $averages{$key_n002} = $averages{$key_n002}/$den;
    $nvips_av{$key_n002} = $nvips_av{$key_n002}/$n_den;
    if($averages{$key_n002}>$nvips_av{$key_n002}){
        $highlow{$key_n002} = "higher";
    }
    else{
        $highlow{$key_n002} = "lower";
    }
}
my %inf;
my %sup;
my $key_002;
foreach $key_002 (sort keys %averages){
    if($highlow{$key_002} eq "higher"){
        $inf{$key_002} = (1-$range)*$averages{$key_002};
        $sup{$key_002} = (1+1*$range)*$averages{$key_002};
    }
    if($highlow{$key_002} eq "lower"){
        $inf{$key_002} = (1-1*$range)*$averages{$key_002};
        $sup{$key_002} = (1+$range)*$averages{$key_002};
    }
}

#my %inf;
#my %sup;
#my $key_002;
#foreach $key_002 (sort keys %averages){
    #$averages{$key_002} = $averages{$key_002}/$den;

    #$inf{$key_002} = (1-$range)*$averages{$key_002};
    #$sup{$key_002} = (1+$range)*$averages{$key_002};

    #print $key_002."\t".$inf{$key_002}."\t".$sup{$key_002}."\n";
#}

#print $vip_number."\t".$nvip_number."\n";

#sleep 10;

my $previous_gn = 0;

my $same_counter = 0;

open OUT, ">".$out_file_2; # opens the file to print with the control sets of genes. One line is one control set of genes represented by their Ensembl gene ID separated by spaces.

my $new_iter = int($iterations/100); #Control sets are created by packs of 100 to speed up control creation while keeping memory usage in check.

for(my $r=1;$r<=$new_iter;$r++){ # main iterations loop.

    #print $r."\n";

    my %used_number;
    my $key_003;
    foreach $key_003 (sort keys %non_vips){
	$used_number{$key_003} = 0;
    }

    my $fake_seed_num = 100; #fake seed. From there control genes are added to control sets as if there were already 100 genes in the control set that fit the group of interest perfectly for all the confounding factors. We later apply a large burn-in so that this initial trick does not affect the final matching of confounding factors between the group of interest and the control sets.

    my $init_fake = $fake_seed_num;

    my %current_factors;
    my %direction;

    my $key_avdir; # to keep the control sets close the averages for confounding factors in the group of interest, we later require that there some level of alternance with the newly added control genes swithcing between increasing the overall average and decreasing it. 
    foreach $key_avdir (sort keys %averages){
	$direction{$key_avdir} = "start";
	$current_factors{$key_avdir} = $averages{$key_avdir};
    }

    my $added_nonvips = 0;

    my @good_nonvips = ();

    my $av_check = 0;
    my $av_den = 0;

    for(my $k=1;$k<=100000000000000000000;$k++){

	my $lim = $factor_number;
	my $sc_nv = scalar(@nonvips);
	my $rand_1 = int(rand($sc_nv)); # the two control genes to choose to add to the current control set.
	my $rand_2 = int(rand($sc_nv));
	for(my $e=1;$e<=100000000000;$e++){
	    if($rand_2==$rand_1){
		$rand_2 = int(rand($sc_nv));
	    }
	    if($rand_2!=$rand_1){
		last;
	    }
	}

	my $gene_1 = $nonvips[$rand_1];
	my $gene_2 = $nonvips[$rand_2];

	my $factor_chain_1 = $non_vips{$gene_1};
	my $factor_chain_2 = $non_vips{$gene_2};
	my @splitter_fc1 = split(" ",$factor_chain_1);
	my @splitter_fc2 = split(" ",$factor_chain_2);
	my $sc_fc = scalar(@splitter_fc1);
	
	my %test_values;
	my %test_values2;
	my $matching = 0;

	my $key_av001; #check if the pair of potential control genes do not bring the average of confounding factors outside of tolerable range.
	foreach $key_av001 (sort keys %averages){
	    my $inf = $inf{$key_av001};
	    my $sup = $sup{$key_av001};
	    #print $inf."\t".$sup."\n";
	    my $ind = $key_av001-1;
	    my $factor_value_1 = $splitter_fc1[$ind];
	    my $factor_value_2 = $splitter_fc2[$ind];
	    my $sc_t  = scalar(@good_nonvips);

	    my $test_value = ($fake_seed_num*$averages{$key_av001}+$current_factors{$key_av001}*$sc_t+$factor_value_1+$factor_value_2)/($sc_t+2+$fake_seed_num); 
	    my $test_value2 = ($current_factors{$key_av001}*$sc_t+$factor_value_1+$factor_value_2)/($sc_t+2);
	    #if($key_av001==11){
		#print $factor_value_1."\t".$factor_value_2."\n";
		#print $key_av001."\t".$test_value."\t".$current_factors{$key_av001}."\t".$averages{$key_av001}."\t".$sc_t."\n";
	    #}
	    $test_values{$key_av001} = $test_value;
	    $test_values2{$key_av001} = $test_value2;
	    if(($test_value>=$inf)and($test_value<=$sup)){
		$matching += 1;
	    }
	}
	
	if($matching>=$lim){ # for pairs of genes that don't bring averages outside tolerable ranges, check if the direction of change in average was opposite to previous one for a sufficient number of of confouding factors. Not doing this usually results in the algorithm becoming "trapped at the border" of the permissive ranges of confouding factors. 
	    my %test_dir;
	    my $key_cc;
	    foreach $key_cc (sort keys %direction){
		if($test_values{$key_cc}>=$current_factors{$key_cc}){
		    $test_dir{$key_cc} = "more";
		}
		else{
		    $test_dir{$key_cc} = "less";
		}
	    }

	    my $other_dir = 0;
	    my $key_dir;
	    foreach $key_dir (sort keys %direction){
		if($test_dir{$key_dir} ne $direction{$key_dir}){
		    $other_dir += 1;
		}
	    }
	    
	    if(($other_dir>=0)and($used_number{$gene_1}/(scalar(@good_nonvips)+1)<=$repok/$vip_number)and($used_number{$gene_2}/(scalar(@good_nonvips)+1)<=$repok/$vip_number)){ # specifically the direction of changes has to switch every time for at least the specified proportion of confounding factors. If is the case, then the two potential control genes are added to the control set. Only use in difficult cases when the scritp gets trapped at a boundary. Otherwise set to zero, or it will slow down the script.
		
		#sleep 1;

		$av_check += $splitter_fc1[10]+$splitter_fc2[10]; 
		$av_den += 2;
		my $ratio = $av_check/$av_den;

		if($fake_seed_num>0){
		    $fake_seed_num = $fake_seed_num-1;
		}

		#print "check: ".$ratio."\t".$av_check."\t".$av_den."\t".$inf{"11"}."\t".$sup{"11"}."\t".$test_values{"11"}."\t".$fake_seed_num."\n";


		#print $k."\t".$matching."\t".$other_dir."\n";
		push(@good_nonvips,$nonvips[$rand_1]);
		push(@good_nonvips,$nonvips[$rand_2]);

		my $nused_1 = $nonvips[$rand_1];
		my $nused_2 = $nonvips[$rand_2];

		$used_genes{$nused_1} = "yes";
		$used_genes{$nused_2} = "yes";

		$used_number{$nonvips[$rand_1]} += 1;
		$used_number{$nonvips[$rand_2]} += 1;
		my $key_curr;
		foreach $key_curr (sort keys %test_values){
		    $current_factors{$key_curr} = $test_values{$key_curr};
		    #print $current_factors{$key_curr}."\n";
		    #sleep 1;

		}
		my $key_dir2;
		foreach $key_dir2 (sort keys %direction){
		    $direction{$key_dir2} = $test_dir{$key_dir2};
		}
	    }
	}
    
	my $sc_gn = scalar(@good_nonvips);

	if($sc_gn==$previous_gn){
	    $same_counter += 1;
	}

	if($sc_gn>$previous_gn){
            $same_counter = 0;
        }


	if($same_counter>=5000000){
	    print "The script is struggling too much to find matching controls. Stopping bootstrap test here. Stop and restart the entire pipeline for safety"."\n";
	    die;
	}

	if($sc_gn/100==int($sc_gn/100)){

	    if($printed{$sc_gn} ne "yes"){
		print "good non-vips"."\t".$sc_gn."\n";
	    }
	    $printed{$sc_gn} = "yes";
	}

	if($sc_gn>=$init_fake+(10+100)*$vip_number){
	    last;
	}

	$previous_gn = $sc_gn;
    }

    my $counter = 0;
    my $key_used;
    foreach $key_used (sort keys %used_number){
	if($used_number{$key_used}>=1){
	    $counter += 1;
	}
    }

    my $sc_t = scalar(@good_nonvips);

    ###########################
    ###########################

    #What the script did is to actually build a very large list of control genes that can now be cut to create the control sets "per se". Creating a big list speeds up the process because it avoids repeating the slow start of adding control genes to the "fake" seed.

    for(my $p=0;$p<=100-1;$p++){

	my $sup_ind = $sc_t-1-$p*$vip_number;
	my $inf_ind = $sc_t-1-$p*$vip_number-$vip_number+1;
	my @iter = @good_nonvips[$inf_ind..$sup_ind];
	my $z=($r-1)*100+$p+1;
	print OUT "sample_".$z." "."@iter"."\n"; #prints the control set number z in the output file of control sets.
	my %iters;

	my $sc_i = scalar(@iter);

	#print $sc_i."\n";

	my %verify;
        
        my $sc_it = scalar(@iter);
        
        for(my $q=0;$q<=$sc_it-1;$q++){
            my $gene = $iter[$q];
            my @splitter_fact = split(" ",$factors{$gene});
            my $sc_sf = scalar(@splitter_fact);
            for(my $it=0;$it<=$sc_sf-1;$it++){
                my $numb = $it+1;
                $verify{$numb} += $splitter_fact[$it];
            }
        }
        
        my $key_verif;
        foreach $key_verif (sort keys %verify){
            $verify{$key_verif} = $verify{$key_verif}/$sc_it;
            #print $averages{$key_verif}."\t".$verify{$key_verif}."\n";
        }


	for(my $q=0;$q<=$sc_i-1;$q++){
	    $iters{$iter[$q]} = "";
	}
	my $counter_iter = 0;
	my $key_iter;
	foreach $key_iter (sort keys %iters){
	    $counter_iter += 1;
	}
	#print "sample_".$z." ".$vip_number." ".$counter_iter."\n";
    }
}

close OUT;

open OUT, ">".$indir."used_genes.txt"."\n";

my $key_used;
foreach $key_used (sort keys %used_genes){
    print OUT $key_used."\n";
}
close OUT;
