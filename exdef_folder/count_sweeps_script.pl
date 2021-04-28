#!/usr/bin/perl

#use warnings;

#use File::Slurp;

use List::Util qw[min max];

use strict;

my ($set_file,$control_file,$sweep_file,$cluster_file,$setpop,$outdir,$use_clust,$threshold_chain,$pop_chain,$group_chain) = @ARGV;

print "Counting with file: ".$sweep_file."\n";

my %iters;

$iters{"1"} = 100;
$iters{"2"} = 1000;
$iters{"3"} = 1000000000000000;

my $direction = "";

for(my $r=1;$r<=3;$r++){

my @splitter_sf = split(/\//,$sweep_file);

my $sc_sf = scalar(@splitter_sf);

#my $suffix = "_001";

my $file_name = $splitter_sf[$sc_sf-1];

    open OUT, ">".$outdir.$file_name;

#print $outdir.$file_name."\n";

#my @popnames = ("ACB","ASW","BEB","CDX","CEU","CLM","CHB","CHS","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI");

my @popnames = split(/\,/,$pop_chain);

#my @continen = ("MIX","MIX","SAS","EAS","EUR","AMR","EAS","EAS","AFR","EUR","EUR","SAS","AFR","EUR","SAS","EAS","EAS","AFR","AFR","AMR","AMR","SAS","AMR","SAS","EUR","AFR");

my @continen = split(/\,/,$group_chain);

my @group_num = ($setpop);
my @group_den = ("AFR","EUR","SAS","EAS","AMR");


#                132    7    40    221    3   1255   80     2     8     5     3   2129   48     2   1138  1041   141    44    26    3     1     56    9    102    15   43

my @cutoffs = split(/\,/,$threshold_chain);

#my @cutoffs = (2000,1500,1000,900,800,700,600,500,450,400,350,300,250,200,150,100,80,60,50,40,30,25,20,15,10);

#my @cutoffs = (2000,1500,1000,900,800,700,600,500,450,400,350,300,250,200,150,100,50,25);
#my @cutoffs = (500,450,400,350,300,250,200,150,100,50);
#my @cutoffs = (500);

my %cutoff_indices;

my $sc_ct = scalar(@cutoffs);

for(my $i=0;$i<=$sc_ct-1;$i++){
    $cutoff_indices{$cutoffs[$i]} = $i;
}

my $sc_cutfirst = scalar(@cutoffs);
my $above = 0;
my %groups;
my %pops;
my %inds;

my $sc_c = scalar(@continen);

for(my $i=0;$i<=$sc_c-1;$i++){
    $groups{$continen[$i]} = "yes";
    $pops{$popnames[$i]} = "yes";
    my $indice = "p".$i."p";
    $inds{$continen[$i]} .= $indice." ";
}

my %set_score;
my %set_recurrent_score;
my %group_score;
my %pop_number;
my %iter_number;
my %control_scores;
my %control_recurrent_scores;

my %set_all_score;
my %control_all_scores;

my %sweeps;
my %minimum;
my %sweep_recurrence;

open(DATA,$sweep_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    my @sweeps = @splitter_line;
    shift @sweeps;
    $sweeps{$splitter_line[0]} = "@sweeps";

    if($sweeps{$splitter_line[0]} ne ""){
	my $min = "yes";
	$minimum{$splitter_line[0]} = $min;
    
	my $gene = $splitter_line[0];
	my $sc_sw = scalar(@sweeps);
	for(my $i=0;$i<=$sc_sw-1;$i++){
	    my $chain = $sweeps[$i];
	    my @splitter_chain = split("_",$chain);
	    my $sc_ch = scalar(@splitter_chain);

	    my $cutoff = $splitter_chain[0];

	    for(my $j=1;$j<=$sc_ch-1;$j++){

		my $indice = $splitter_chain[$j]-1;
		$sweep_recurrence{"All"." ".$gene." ".$cutoff} +=1;
		my $cutoff_ind = $cutoff_indices{$cutoff};
		for(my $k=0;$k<=$cutoff_ind-1;$k++){
		    $sweep_recurrence{"All"." ".$gene." ".$cutoffs[$k]} +=1;
		}

		my $key_01;
		foreach $key_01 (sort keys %inds){
		    my $ind_chain = $inds{$key_01};
		    my $here = "no";
		    my @splitter_ic = split(" ",$ind_chain);
		    my $sc_ic = scalar(@splitter_ic);
		    for(my $l=0;$l<=$sc_ic-1;$l++){
			my $cur_ind = $splitter_ic[$l];
			#print $cur_ind."\n";
			if($cur_ind eq "p".$indice."p"){
			    $here = "yes";
			}
		    }

		    #print $here."\n";

		    if($here eq "yes"){
			
			#print $key_01."\n";
			#print $gene."\n";
			#print $cutoff."\n";

			$sweep_recurrence{$key_01." ".$gene." ".$cutoff} += 1;
			#print $key_01." ".$gene." ".$cutoff."\n";
			my $cutoff_ind = $cutoff_indices{$cutoff};
			for(my $k=0;$k<=$cutoff_ind-1;$k++){
			    $sweep_recurrence{$key_01." ".$gene." ".$cutoffs[$k]} += 1;
			}
		    }
		}
	    }
	}
    }
}
close DATA;

my %clusters;

open(DATA,$cluster_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    my $gene = $splitter_line[0];
    if($minimum{$gene} eq "yes"){
        my $new_chain = "";
        my $sc_sl = scalar(@splitter_line);
        for(my $i=1;$i<=$sc_sl-1;$i++){
            my $ng_gene = $splitter_line[$i];
            if($minimum{$gene} eq "yes"){
                $new_chain .= $ng_gene." ";
            }
        }
        chop $new_chain;

	if($use_clust eq "yes"){
	    $clusters{$gene} = $new_chain;
	}
    }
}
close DATA;


open(DATA,$set_file);
while(<DATA>){
    chomp $_;
    my $gene = $_;
    for(my $i=0;$i<=$sc_c-1;$i++){
	for(my $p=0;$p<=$sc_cutfirst-1;$p++){
	    $set_score{$i." ".$p} = 0;
	    $pop_number{$i} = "yes";
	}
    }
}
close DATA;

my %done;

open(DATA,$set_file);
while(<DATA>){
    chomp $_;
    my $gene = $_;
    my $sweep_chain = $sweeps{$gene};

    if($sweep_chain ne ""){

	my @splitter_chain = split(" ",$sweep_chain);
	my $sc_sc = scalar(@splitter_chain);
	for(my $s=0;$s<=$sc_sc-1;$s++){
	    my $current_cut = $splitter_chain[$s];
	    my @splitter_cut = split("_",$current_cut);
	    my $sc_cut = scalar(@splitter_cut);
	    my $cutoff = $splitter_cut[0];
	    my $cutoff_indice = $cutoff_indices{$cutoff};
	    if($cutoff<=$cutoffs[0]){
		for(my $p=1;$p<=$sc_cut-1;$p++){
		    my $current_ind = $splitter_cut[$p]-1;
		    if($done{$gene." ".$current_ind." ".$cutoff_indice} ne "yes"){
			$set_score{$current_ind." ".$cutoff_indice} += 1;

			my $continent = $continen[$current_ind];
			my $recurr_line = $continent." ".$gene." ".$cutoff;
			my $all_line = "All"." ".$gene." ".$cutoff;

			$set_recurrent_score{$current_ind." ".$cutoff_indice} += 1/$sweep_recurrence{$recurr_line};
			$set_all_score{$current_ind." ".$cutoff_indice} += 1/$sweep_recurrence{$all_line};

                        $done{$gene." ".$current_ind." ".$cutoff_indice} = "yes";
			for(my $h=0;$h<=$cutoff_indice-1;$h++){
                            if($done{$gene." ".$current_ind." ".$h} ne "yes"){
                                $set_score{$current_ind." ".$h} += 1;

				$recurr_line = $continent." ".$gene." ".$cutoffs[$h];
				$all_line = "All"." ".$gene." ".$cutoffs[$h];
				$set_recurrent_score{$current_ind." ".$h} += 1/$sweep_recurrence{$recurr_line};
				$set_all_score{$current_ind." ".$h} += 1/$sweep_recurrence{$all_line};
                                $done{$gene." ".$current_ind." ".$h} = "yes";
                            }
                        }
			my $cluster_genes = $clusters{$gene};
			my @splitter_clust = split(" ",$cluster_genes);
			my $sc_clust = scalar(@splitter_clust);
			for(my $w=0;$w<=$sc_clust-1;$w++){
			    my $clust_gene = $splitter_clust[$w];
			    my $clust_sweep_chain = $sweeps{$clust_gene};
			    my @splitter_sweep_chain = split(" ",$clust_sweep_chain);
			    my $sc_swc = scalar(@splitter_sweep_chain);
			    for(my $t=0;$t<=$sc_swc-1;$t++){
				my $current_cutclust = $splitter_sweep_chain[$t];
				my @splitter_cutclust = split("_",$current_cutclust);
				my $sc_cutclust = scalar(@splitter_cutclust);
				my $cutoffclust = $splitter_cutclust[0];
				my $cutoff_indice_clust = $cutoff_indices{$cutoffclust};
				if($cutoffclust<=$cutoffs[0]){
				    for(my $x=1;$x<=$sc_cutclust-1;$x++){
					my $current_ind_clust = $splitter_cutclust[$x]-1;
					if($current_ind_clust==$current_ind){
					    if($done{$clust_gene." ".$current_ind_clust." ".$cutoff_indice_clust} ne "yes"){
						#$set_score{$current_ind_clust." ".$cutoff_indice_clust} += 1;

						$done{$clust_gene." ".$current_ind_clust." ".$cutoff_indice_clust} = "yes";
						for(my $h=0;$h<=$cutoff_indice_clust-1;$h++){
						    if($done{$clust_gene." ".$current_ind_clust." ".$h} ne "yes"){
							#$set_score{$current_ind_clust." ".$h} += 1;
							$done{$clust_gene." ".$current_ind_clust." ".$h} = "yes";
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}
close DATA;

my $counter = 0;

#my $file_content = read_file($control_file);

#my @splitter_all = split("\n",$file_content);

#my $sc_sa = scalar(@splitter_all);

open(DATA,$control_file);

while(<DATA>){
  
    chomp $_;

    $counter += 1;

    if($counter<=$iters{$r}){

	
	my %iterated;
	my %rand_done;

	$iter_number{$counter} = "yes";
	my @splitter_line = split(" ",$_);
	#print $counter."\n";
	shift @splitter_line;
	my $sc_sl = scalar(@splitter_line);
	for(my $j=0;$j<=$sc_sl-1;$j++){
	    my $gene = $splitter_line[$j];
	    $iterated{$gene} += 1;
	}

	my $sc_sl = scalar(@splitter_line);
        for(my $j=0;$j<=$sc_sl-1;$j++){
            my $gene = $splitter_line[$j];
            my $sweep_chain = $sweeps{$gene};
            if($sweep_chain ne ""){

		my %num;
		my %den;

		my %num_recurr;
		my %den_recurr;

		my %num_all;
		my %den_all;

		my @splitter_chain = split(" ",$sweep_chain);
		my $sc_sc = scalar(@splitter_chain);
		for(my $s=0;$s<=$sc_sc-1;$s++){
		    my $current_cut = $splitter_chain[$s];
		    my @splitter_cut = split("_",$current_cut);
		    my $sc_cut = scalar(@splitter_cut);
		    my $cutoff = $splitter_cut[0];
		    my $cutoff_indice = $cutoff_indices{$cutoff};
		    if($cutoff<=$cutoffs[0]){
			for(my $p=1;$p<=$sc_cut-1;$p++){
			    my $current_ind = $splitter_cut[$p]-1;
			    if($rand_done{$gene." ".$current_ind." ".$cutoff_indice} ne "yes"){
				#$control_scores{$counter." ".$current_ind." ".$cutoff_indice} += 1;
			
				$num{$current_ind." ".$cutoff_indice} += $iterated{$gene};
				$den{$current_ind." ".$cutoff_indice} += 1;

				my $continent = $continen[$current_ind];
				my $recurr_line = $continent." ".$gene." ".$cutoff;
				my $all_line = "All"." ".$gene." ".$cutoff;
				$num_recurr{$current_ind." ".$cutoff_indice} += $iterated{$gene}/$sweep_recurrence{$recurr_line};
				$den_recurr{$current_ind." ".$cutoff_indice} += 1;
				$num_all{$current_ind." ".$cutoff_indice} += $iterated{$gene}/$sweep_recurrence{$all_line};
				$den_all{$current_ind." ".$cutoff_indice} += 1;


				$rand_done{$gene." ".$current_ind." ".$cutoff_indice} = "yes";
				for(my $h=0;$h<=$cutoff_indice-1;$h++){
				    if($rand_done{$gene." ".$current_ind." ".$h} ne "yes"){
					#$control_scores{$counter." ".$current_ind." ".$h} += 1;
					$num{$current_ind." ".$h} += $iterated{$gene};
					$den{$current_ind." ".$h} += 1;

					$recurr_line = $continent." ".$gene." ".$cutoffs[$h];
					$all_line = "All"." ".$gene." ".$cutoffs[$h];
					$num_recurr{$current_ind." ".$h} += $iterated{$gene}/$sweep_recurrence{$recurr_line};
					$den_recurr{$current_ind." ".$h} += 1;
					$num_all{$current_ind." ".$h} += $iterated{$gene}/$sweep_recurrence{$all_line};
					$den_all{$current_ind." ".$h} += 1;

					$rand_done{$gene." ".$current_ind." ".$h} = "yes";
				    }
				}
				my $cluster_genes = $clusters{$gene};
				my @splitter_clust = split(" ",$cluster_genes);
				my $sc_clust = scalar(@splitter_clust);
				for(my $w=0;$w<=$sc_clust-1;$w++){
				    my $clust_gene = $splitter_clust[$w];

				    if($iterated{$clust_gene}>=1){

					my $clust_sweep_chain = $sweeps{$clust_gene};
					my @splitter_sweep_chain = split(" ",$clust_sweep_chain);
					my $sc_swc = scalar(@splitter_sweep_chain);
					for(my $t=0;$t<=$sc_swc-1;$t++){
					    my $current_cutclust = $splitter_sweep_chain[$t];
					    my @splitter_cutclust = split("_",$current_cutclust);
					    my $sc_cutclust = scalar(@splitter_cutclust);
					    my $cutoffclust = $splitter_cutclust[0];
					    my $cutoff_indice_clust = $cutoff_indices{$cutoffclust};
					    if(($cutoffclust<=$cutoffs[0])and($iterated{$clust_gene}>=1)){
						for(my $x=1;$x<=$sc_cutclust-1;$x++){
						    my $current_ind_clust = $splitter_cutclust[$x]-1;
						    if($current_ind_clust==$current_ind){
							if($rand_done{$clust_gene." ".$current_ind_clust." ".$cutoff_indice_clust} ne "yes"){
							    #$control_scores{$counter." ".$current_ind_clust." ".$cutoff_indice_clust} += 1;
							    $num{$current_ind_clust." ".$cutoff_indice_clust} += $iterated{$clust_gene};
							    $den{$current_ind_clust." ".$cutoff_indice_clust} += 1;

							    my $clust_recurr_line = $continent." ".$clust_gene." ".$cutoffclust;
							    #print $clust_recurr_line."\n";
							    #print $sweep_recurrence{$clust_recurr_line}."\n";
							    my $clust_all_line = "All"." ".$clust_gene." ".$cutoffclust;
							    $num_recurr{$current_ind_clust." ".$cutoff_indice_clust} += $iterated{$clust_gene}/$sweep_recurrence{$clust_recurr_line};
							    $den_recurr{$current_ind_clust." ".$cutoff_indice_clust} += 1;
							    $num_all{$current_ind_clust." ".$cutoff_indice_clust} += $iterated{$clust_gene}/$sweep_recurrence{$clust_all_line};
							    $den_all{$current_ind_clust." ".$cutoff_indice_clust} += 1;

							    $rand_done{$clust_gene." ".$current_ind_clust." ".$cutoff_indice_clust} = "yes";
							    for(my $h=0;$h<=$cutoff_indice_clust-1;$h++){
								if($rand_done{$clust_gene." ".$current_ind_clust." ".$h} ne "yes"){
								    #$control_scores{$counter." ".$current_ind_clust." ".$h} += 1;
								    $num{$current_ind_clust." ".$h} += $iterated{$clust_gene};
								    $den{$current_ind_clust." ".$h} += 1;

								    $clust_recurr_line = $continent." ".$clust_gene." ".$cutoffs[$h];
								    $clust_all_line = "All"." ".$clust_gene." ".$cutoffs[$h];
								    $num_recurr{$current_ind_clust." ".$h} += $iterated{$clust_gene}/$sweep_recurrence{$clust_recurr_line};
								    $den_recurr{$current_ind_clust." ".$h} += 1;
								    $num_all{$current_ind_clust." ".$h} += $iterated{$clust_gene}/$sweep_recurrence{$clust_all_line};
								    $den_all{$current_ind_clust." ".$h} += 1;
								  
								    $rand_done{$clust_gene." ".$current_ind_clust." ".$h} = "yes";
								}
							    }
							}
						    }
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}

		my $key_01;

		foreach $key_01 (sort keys %num){
		    $num{$key_01} = $num{$key_01}/$den{$key_01};
		    $control_scores{$counter." ".$key_01} += $num{$key_01};
		    
		    $num_recurr{$key_01} = $num_recurr{$key_01}/$den_recurr{$key_01};
                    $control_recurrent_scores{$counter." ".$key_01} += $num_recurr{$key_01};

		    $num_all{$key_01} = $num_all{$key_01}/$den_all{$key_01};
                    $control_all_scores{$counter." ".$key_01} += $num_all{$key_01};

		}
	    }
	}
    }
}

close DATA;

#open OUT, ">"."/Users/davidenard/Dropbox/Anisa_project/test_outputs/".$file_name;

for(my $p=0;$p<=$sc_cutfirst-1;$p++){

my $total_score = 0;

my $total_pval = 0;
my $total_den = 0;
my $total_av = 0;
my $total_lowCI = 0;
my $total_highCI = 0;
my @total = ();

my $simple_total_pval = 0;
my $simple_total_av = 0;
my $simple_total_lowCI = 0;
my $simple_total_highCI = 0;
my @simple_total = ();

my $outgroup_pval = 0;
my $outgroup_den = 0;
my $outgroup_av = 0;
my $outgroup_lowCI = 0;
my $outgroup_highCI = 0;
my @outgroup = ();

my %group_pval;
my %group_av;
my %group_lowCI;
my %group_highCI;
my %group_all;
my %group_scores;
my $counter_001;

my %pop_pval;
my %pop_av;
my %pop_all;
my %pop_scores;

my $group_ratio_av = 0;
my $group_ratio_pval = 0;

my @group_ratios = ();


my $key_group3;

foreach $key_group3 (sort keys %groups){
    $group_pval{$key_group3} = 0;
    $group_av{$key_group3}  = 0;
    $group_scores{$key_group3}  = 0;
} 

my $key_pop3;

foreach $key_pop3 (sort keys %pops){
    $pop_pval{$key_pop3} = 0;
    $pop_av{$key_pop3}  = 0;
    $pop_scores{$key_pop3}  = 0;
}


my $group_num_score = 1;
my $group_den_score = 1;

my $set_group_ratio = 0;

my $key_03;

foreach $key_03 (sort keys %iter_number){

    #print $key_03."\n";

    my $total_ratio = 0;

    my $random_group_ratio = 0;

    my $random_group_num_score = 0;
    my $random_group_den_score = 0;
    

    my $total_rscore = 0;

    my %current_groups;
    my %current_pops;

    $counter_001 += 1;

    my $den  = 0;
    my $key_04;
    foreach $key_04 (sort{$a<=>$b} keys %pop_number){
	my $ratio = ($set_score{$key_04." ".$p}+0.1)/($control_scores{$key_03." ".$key_04." ".$p}+0.1);
	$total_ratio += $ratio;
	$den += 1;
	$total_rscore += $control_all_scores{$key_03." ".$key_04." ".$p};

	#print $control_all_scores{$key_03." ".$key_04." ".$p}."\n";
	my $key_group = $continen[$key_04];
	$current_groups{$key_group} += $control_recurrent_scores{$key_03." ".$key_04." ".$p};

	if("@group_num" =~ $key_group){
	    $random_group_num_score += $control_recurrent_scores{$key_03." ".$key_04." ".$p};
	}
	if("@group_den" =~ $key_group){

	    if(not("@group_num" =~ $key_group)){
		$random_group_den_score += $control_recurrent_scores{$key_03." ".$key_04." ".$p};
	    }
	}



	if($counter_001==1){
	    $total_score += $set_all_score{$key_04." ".$p};
	    $group_scores{$key_group} += $set_recurrent_score{$key_04." ".$p};
	    
	    if("@group_num" =~ $key_group){
		$group_num_score += $set_recurrent_score{$key_04." ".$p};
	    }
	    if("@group_den" =~ $key_group){
		if(not("@group_num" =~ $key_group)){
		    $group_den_score += $set_recurrent_score{$key_04." ".$p};
		}
	    }


	}

	#$set_group_ratio = $group_num_score/$group_den_score;
	#$random_group_ratio = $random_group_num_score/$random_group_den_score;

	my $key_pop = $popnames[$key_04];
	$current_pops{$key_pop} += $control_scores{$key_03." ".$key_04." ".$p};
        if($counter_001==1){
            $pop_scores{$key_pop} += $set_score{$key_04." ".$p};
        }
    }

    $set_group_ratio = ($group_num_score+1)/($group_den_score+1);
    $random_group_ratio = ($random_group_num_score+1)/($random_group_den_score+1);

    push(@group_ratios,$random_group_ratio);

    if($random_group_ratio>$set_group_ratio){
	$group_ratio_pval += 1;
    }
    $group_ratio_av += $random_group_ratio;


    my $key_group2;
    foreach $key_group2 (sort keys %groups){
	$group_av{$key_group2} += $current_groups{$key_group2};
	$group_all{$key_group2} .= $current_groups{$key_group2}." ";
	if($current_groups{$key_group2}>$group_scores{$key_group2}){
	    $group_pval{$key_group2} += 1;
	}
    }

    my $key_pop2;
    foreach $key_pop2 (sort keys %pops){
        $pop_av{$key_pop2} += $current_pops{$key_pop2};
        $pop_all{$key_pop2} .= $current_pops{$key_pop2}." ";
        if($current_pops{$key_pop2}>$pop_scores{$key_pop2}){
            $pop_pval{$key_pop2} += 1;
        }
    }

    $total_ratio = $total_ratio/$den;
    if($total_ratio<=1){
	$total_pval += 1;
    }
    $total_den += 1;
    $total_av += $total_ratio;
    $simple_total_av += $total_rscore;

    if($total_rscore>$total_score){
	$simple_total_pval += 1;
    }

    push(@total,$total_ratio);
    push(@simple_total,$total_rscore);

    $outgroup_den += 1;
    $outgroup_av += $random_group_den_score;

    if($random_group_den_score>$group_den_score){
        $outgroup_pval += 1;
    }
    push(@outgroup,$random_group_den_score);

}

$total_av = $total_av/$total_den;
$total_pval = $total_pval/$total_den;
my @sorted_total = sort{$a<=>$b}@total;
my $sc_tot = scalar(@sorted_total);
my $lowind = int(0.025*$sc_tot)-1;
my $supind = int(0.975*$sc_tot)-1;
my $low_ratio = $sorted_total[$lowind];
my $high_ratio = $sorted_total[$supind];
#print "All: ".$total_score." ".$total_av." ".$low_ratio." ".$high_ratio." ".$total_pval."\n";

my $key_pop4;
foreach $key_pop4 (sort keys %pops){
    $pop_av{$key_pop4} = $pop_av{$key_pop4}/$total_den;
    $pop_pval{$key_pop4} = $pop_pval{$key_pop4}/$total_den;
    my @simple_pop = split(" ",$pop_all{$key_pop4});
    my @simple_sorted_pop = sort{$a<=>$b}@simple_pop;
    my $sc_sg = scalar(@simple_sorted_pop);
    my $slow = int(0.025*$sc_sg);
    my $ssup = int(0.975*$sc_sg);
    my $s_low = $simple_sorted_pop[$slow];
    my $s_high = $simple_sorted_pop[$ssup];
    my $pop_ratio = ($pop_scores{$key_pop4}+0.1)/($pop_av{$key_pop4}+0.1);
    print OUT $cutoffs[$p]." ".$key_pop4.": ".$pop_ratio." ".$pop_scores{$key_pop4}." ".$pop_av{$key_pop4}." ".$s_low." ".$s_high." ".$pop_pval{$key_pop4}."\n";
}

my $key_group4;
foreach $key_group4 (sort keys %groups){
    $group_av{$key_group4} = $group_av{$key_group4}/$total_den;
    $group_pval{$key_group4} = $group_pval{$key_group4}/$total_den;
    my @simple_group = split(" ",$group_all{$key_group4});
    my @simple_sorted_group = sort{$a<=>$b}@simple_group;
    my $sc_sg = scalar(@simple_sorted_group);
    my $slow = int(0.025*$sc_sg);
    my $ssup = int(0.975*$sc_sg);
    my $s_low = $simple_sorted_group[$slow];
    my $s_high = $simple_sorted_group[$ssup];
    my $group_ratio = ($group_scores{$key_group4}+0.1)/($group_av{$key_group4}+0.1);
    print OUT $cutoffs[$p]." ".$key_group4.": ".$group_ratio." ".$group_scores{$key_group4}." ".$group_av{$key_group4}." ".$s_low." ".$s_high." ".$group_pval{$key_group4}."\n";
}

$simple_total_av = $simple_total_av/$total_den;
$simple_total_pval = $simple_total_pval/$total_den;
my @simple_sorted_total = sort{$a<=>$b}@simple_total;
my $sc_stot = scalar(@simple_sorted_total);
my $slowind = int(0.025*$sc_stot);
my $ssupind = int(0.975*$sc_stot);
my $simple_low = $simple_sorted_total[$slowind];
my $simple_high = $simple_sorted_total[$ssupind];
my $simple_ratio = ($total_score+0.1)/($simple_total_av+0.1);
print OUT $cutoffs[$p]." "."All: ".$simple_ratio." ".$total_score." ".$simple_total_av." ".$simple_low." ".$simple_high." ".$simple_total_pval."\n";

$outgroup_av = $outgroup_av/$outgroup_den;
$outgroup_pval = $outgroup_pval/$outgroup_den;
my @sorted_outgroup = sort{$a<=>$b}@outgroup;
my $sc_out = scalar(@sorted_outgroup);
my $low_out = int(0.025*$sc_out);
my $sup_out = int(0.975*$sc_out);
my $o_low = $sorted_outgroup[$low_out];
my $o_high = $sorted_outgroup[$sup_out];
my $outgroup_ratio = ($group_den_score+0.1)/($outgroup_av+0.1);
print OUT $cutoffs[$p]." "."OUT: ".$outgroup_ratio." ".$group_den_score." ".$outgroup_av." ".$o_low." ".$o_high." ".$outgroup_pval."\n";


$group_ratio_av = $group_ratio_av/$total_den;
$group_ratio_pval = $group_ratio_pval/$total_den;

my @sorted_group_ratios = sort{$a<=>$b}@group_ratios;

my $sc_gr = scalar(@sorted_group_ratios);
my $grind_inf = int(0.025*$sc_gr);
my $grind_sup = int(0.975*$sc_gr);
my $gr_low = $sorted_group_ratios[$grind_inf];
my $gr_high = $sorted_group_ratios[$grind_sup];

my $ratt = ($set_group_ratio+0.1)/($group_ratio_av+0.1);

print OUT $cutoffs[$p]." Group_ratio: ".$ratt." ".$set_group_ratio." ".$group_ratio_av." ".$gr_low." ".$gr_high." ".$group_ratio_pval."\n";

}

    my $higher_95 = 0;
    my $higher_99 = 0;
    my $higher_999 = 0;

    my $stop = "yes";

    open(DATA,$outdir.$file_name);
    while(<DATA>){
	chomp $_;
	if($_ =~ $setpop){
	    my @splitter_line = split(" ",$_);
	    my $pval = $splitter_line[7];

	    #if($pval>=0.95){
	    if(($pval<=0.05)or($pval>=0.95)){
		if($r==1){
		    $stop = "no"; #set back to no
		}
	    }
	   
	    #if($pval>=0.998){
	    if(($pval<=0.002)or($pval>=0.998)){
		if($r==2){
		    $stop = "no"; #set back to no
		}
	    }
	}
    }
    close DATA;

if($r==3){
    $stop = "yes";
}

print OUT $stop."\n";

close OUT;

#print $r." ".$stop."\n";

    if($stop eq "yes"){
	last;
    }
}


