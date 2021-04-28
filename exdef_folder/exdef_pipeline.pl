#!/usr/bin/perl

use warnings;

use strict;

my ($input_file) = @ARGV;

my %parameters;

open(DATA,$input_file);
while(<DATA>){
    chomp $_;
    my @splitter_line = split(" ",$_);
    $parameters{$splitter_line[0]} = $splitter_line[1];
}
close DATA;

my @sweep_files = ();

my $sweep_size_chain = $parameters{"Sweep_sizes"};

my @splitter_ssc = split("_",$sweep_size_chain);

my $sc_ssc = scalar(@splitter_ssc);

my $file_prefix = $parameters{"Sweep_files_prefix"};

for(my $i=0;$i<=$sc_ssc-1;$i++){
    my $sweep_file = $file_prefix."_".$splitter_ssc[$i];
    push(@sweep_files,$sweep_file);
}

my $pipeline_dir = $parameters{"Pipeline_dir"};

my $bootstrap_sweep_file = $sweep_files[0];
my $bootstrap_iterations = $parameters{"Iterations_number"};
my $valid_file = $parameters{"Valid_genes_file"};
my $gene_set_file = $parameters{"Genes_set_file"};
my $control_distance = $parameters{"Min_distance"};
my $control_distance_file = $parameters{"Distance_file"};
my $tolerance_range = $parameters{"Tolerance_range"};
my $factors_table = $parameters{"Factors_table"};
my $hgnc_file = $parameters{"HGNC_file"};
my $flip = $parameters{"Flip"};
my $max_rep = $parameters{"Max_rep"};
my $bootstrap_runs_number = $parameters{"Runs_number"};
my $simult_runs = $parameters{"Simult_runs"};

my $bootstrap_shell_file = $pipeline_dir."bootstrap_rafale.sh";

open OUT, ">".$bootstrap_shell_file;
my $counter_01 = 0;
for(my $i=1;$i<=$bootstrap_runs_number;$i++){
    $counter_01 += 1;
    my $add = " & \n";
    if($counter_01==$simult_runs){   
        $add = "\n"."sleep 5"."\n";
        $counter_01 = 0;
    }

    if($i==$bootstrap_runs_number){
	$add = "\n"."sleep 5"."\n";
    }
    
    my $gene_set_output_file = $pipeline_dir."VIPs_parts/file_".$i;
    my $control_output_file = $pipeline_dir."nonVIPs_parts/file_".$i;

    my $bootstrap_command_line  = $pipeline_dir."bootstrap_test_script.pl ".$pipeline_dir." ".$valid_file." ".$bootstrap_sweep_file." ".$hgnc_file." ".$gene_set_file." ".$factors_table." ".$control_distance_file." ".$bootstrap_iterations." ".$tolerance_range." ".$gene_set_output_file." ".$control_output_file." ".$control_distance." ".$flip." ".$max_rep.$add;

    print OUT $bootstrap_command_line;
}
close OUT;

my $run_bootstrap = $parameters{"Run_bootstrap"};

if($run_bootstrap eq "yes"){
    print "Starting bootstrap test"."\n";
    my $remove_line_1 = "rm -f ".$pipeline_dir."VIPs/*"."\n";
    my $remove_line_2 = "rm -f ".$pipeline_dir."nonVIPs/*"."\n";
    my $remove_line_3 = "rm -f ".$pipeline_dir."VIPs_parts/*"."\n";
    my $remove_line_4 = "rm -f ".$pipeline_dir."nonVIPs_parts/*"."\n";
    system($remove_line_1);
    system($remove_line_2);
    system($remove_line_3);
    system($remove_line_4);
    my $command_line = "chmod +x ".$bootstrap_shell_file." && ".$bootstrap_shell_file."\n";
    system($command_line);
    
    for(my $e=1;$e<=100000000000000;$e++){
	sleep 5;
	my $f_done = 0;
	
	my $pipeline_dirchopped = $pipeline_dir;
	chop $pipeline_dirchopped;

	while(<$pipeline_dirchopped/nonVIPs_parts/*>){
	    my $file = $_;
	    my $reco = "sample_".$bootstrap_iterations;
	    open(DATA,$file);
	    while(<DATA>){
		my @splitter_line = split(" ",$_);
		if($splitter_line[0] eq $reco){
		    $f_done += 1;
		}
	    }
	}
	if($f_done==$bootstrap_runs_number){
	    last;
	}
    }




    system("cp ".$pipeline_dir."VIPs_parts/file_1 ".$pipeline_dir."VIPs/file_1"."\n");
    system("cat ".$pipeline_dir."nonVIPs_parts/file_* > ".$pipeline_dir."nonVIPs/file_1"."\n");
    print "Boostrap test done"."\n";
}

if($run_bootstrap eq "no"){
    print "Skipping bootstrap test. Make sure it was run before. Otherwise the rest of the pipeline will fail"."\n";
}

my $run_sweep_counting = $parameters{"Run_sweeps_count"};

my $pop = $parameters{"Pop_interest"};

my $cluster_distance = $parameters{"Cluster_distance"};

my $gene_coords_file = $parameters{"Gene_coords_file"};

my $threshold_chain = $parameters{"Threshold_list"};

my @splitter_thc = split(/\,/,$threshold_chain);

my $maxi_rank = $splitter_thc[0];

my $pop_chain = $parameters{"Populations_list"};

my $group_chain = $parameters{"Groups_list"};

my $count_sweeps = $parameters{"Count_sweeps"};

if($run_sweep_counting eq "yes"){
    print "Starting sweep or gene counting."."\n";

    if($count_sweeps eq "yes"){
	print "Counting the number of sweeps overlapping genes matched by the bootstrap test"."\n";
    }

    if($count_sweeps eq "no"){
	print "Counting the number of genes matched by the bootstrap test within sweeps"."\n";
	print "Warning: this means that clustring of genes in sweeps will be taken into account only if you run a full FDR analysis."."\n";
    }

    my $sc_sf = scalar(@sweep_files);

    for(my $i=0;$i<=$sc_sf-1;$i++){
	my $current_sweep_file = $sweep_files[$i];
	
	#make gene neighbors file here.

	my $ng_command_line = $pipeline_dir."get_all_neighbors.pl ".$pipeline_dir.$gene_coords_file." ".$pipeline_dir."ensembl_gene_neighbors.txt ".$cluster_distance."\n";

	system("$ng_command_line");

	#simplify sweep file here.

	system("gzip ".$pipeline_dir.$current_sweep_file);

	my $simplify_line = $pipeline_dir."simplify_zipped_sweeps.pl ".$pipeline_dir.$current_sweep_file.".gz ".$pipeline_dir."real_simple/".$current_sweep_file." ".$threshold_chain."\n";
	
	system($simplify_line);

	system("gunzip ".$pipeline_dir.$current_sweep_file.".gz");

	my $command_line = $pipeline_dir."count_sweeps_script.pl ".$pipeline_dir."VIPs/file_1 ".$pipeline_dir."nonVIPs/file_1 ".$pipeline_dir."real_simple/".$current_sweep_file." ".$pipeline_dir."ensembl_gene_neighbors.txt ".$pop." ".$pipeline_dir."/test_outputs/ ".$count_sweeps." ".$threshold_chain." ".$pop_chain." ".$group_chain."\n";

	system($command_line);
    }

}

if($run_sweep_counting eq "no"){
    print "Skipping sweep or gene counting. Make sure it was run before. Otherwise the rest of the pipeline will fail"."\n";
}

my $run_FDR = $parameters{"Run_FDR"};

if($run_FDR eq "no"){
    print "Not running FDR analysis. Make sure you used number of sweeps not number genes in sweeps."."\n"; 
}

if($run_FDR eq "yes"){

    print "Attempting to run FDR analysis. First checking if it is required."."\n";

    my $significant = 0;

    my $sc_sf = scalar(@sweep_files);

    for(my $i=0;$i<=$sc_sf-1;$i++){
        my $current_sweep_file = $sweep_files[$i];
	my $sign_file = $pipeline_dir."test_outputs/".$current_sweep_file."\n";
	open(DATA,$sign_file);
	while(<DATA>){
	    chomp $_;
	    if($_ =~ $pop){
		my @splitter_line = split(" ",$_);
		if(($splitter_line[7]<=0.05)or($splitter_line[7]>=0.95)){
		    $significant += 1;
		}
	    }
	}
	close DATA;
    }

    if($significant<2){
	print "FDR analysis not required. Stopping here"."\n";
    }

    if($significant>=2){

	print "Starting FDR analysis."."\n";
	
	my $shuffle = $parameters{"Run_shuffling"};

	if($shuffle eq "yes"){
	   
	    print "Starting genome shuffling."."\n";

	    #first remove whatever is in ./shuffled_sweep_files/ and in ./shuffled_sweep_files_simplified/ (see below).

	    my $remove_line = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files _1_ &"."\n";
	    system($remove_line);
	    $remove_line = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files _2_ &"."\n";
	    system($remove_line);
	    $remove_line = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files _3_ &"."\n";
	    system($remove_line);
	    $remove_line = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files _4_ &"."\n";
	    system($remove_line);
	    $remove_line = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files _5_ &"."\n";
	    system($remove_line);
	    $remove_line = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files _6_ &"."\n";
	    system($remove_line);
	    $remove_line = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files _7_ &"."\n";
	    system($remove_line);
	    $remove_line = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files _8_"."\n";
	    system($remove_line);

	    my $remove_line_2 = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files_simplified _1_ &"."\n";
    	    system($remove_line_2);
	    $remove_line_2 = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files_simplified _2_ &"."\n";
	    system($remove_line_2);
	    $remove_line_2 = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files_simplified _3_ &"."\n";
	    system($remove_line_2);
	    $remove_line_2 = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files_simplified _4_ &"."\n";
	    system($remove_line_2);
	    $remove_line_2 = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files_simplified _5_ &"."\n";
	    system($remove_line_2);
	    $remove_line_2 = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files_simplified _6_ &"."\n";
	    system($remove_line_2);
	    $remove_line_2 = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files_simplified _7_ &"."\n";
	    system($remove_line_2);
	    $remove_line_2 = $pipeline_dir."remove.pl ".$pipeline_dir."shuffled_sweep_files_simplified _8_"."\n";
	    system($remove_line_2);

	    my $remove_line_3 = $pipeline_dir."remove.pl ".$pipeline_dir."test_outputs _1_ &"."\n";
            system($remove_line_3);
	    $remove_line_3 = $pipeline_dir."remove.pl ".$pipeline_dir."test_outputs _2_ &"."\n";
            system($remove_line_3);
	    $remove_line_3 = $pipeline_dir."remove.pl ".$pipeline_dir."test_outputs _3_ &"."\n";
            system($remove_line_3);
	    $remove_line_3 = $pipeline_dir."remove.pl ".$pipeline_dir."test_outputs _4_ &"."\n";
            system($remove_line_3);
	    $remove_line_3 = $pipeline_dir."remove.pl ".$pipeline_dir."test_outputs _5_ &"."\n";
            system($remove_line_3);
	    $remove_line_3 = $pipeline_dir."remove.pl ".$pipeline_dir."test_outputs _6_ &"."\n";
            system($remove_line_3);
	    $remove_line_3 = $pipeline_dir."remove.pl ".$pipeline_dir."test_outputs _7_ &"."\n";
            system($remove_line_3);
	    $remove_line_3 = $pipeline_dir."remove.pl ".$pipeline_dir."test_outputs _8_ "."\n";
            system($remove_line_3);

	    sleep 60;

	    my $shuffle_number = int($parameters{"FDR_number"}/8);
	    open OUT, ">".$pipeline_dir."shuffle_sweeps.sh";
	    my $counter = 0;
	    for(my $i=1;$i<=$shuffle_number;$i++){
		my $add = " & "."\n";
		$counter += 1;
		if($counter==$simult_runs){
		    $add = "\n"."sleep 0.5"."\n";
		    $counter = 0;
		}
		
		my $segcut = $parameters{"Shuffling_segments_number"};

		print OUT $pipeline_dir."shuffle_sweep_script_nocut.pl ".$pipeline_dir." ".$pipeline_dir."shuffled_sweep_files/ ".$pipeline_dir.$valid_file." ".$i." ".$pipeline_dir.$gene_coords_file." ".$file_prefix." ".$sweep_size_chain." ".$segcut." ".$maxi_rank." && gzip -f ".$pipeline_dir."shuffled_sweep_files/*_".$i.$add;
	    }
	    close OUT;
	    
	    my $command_line = "chmod +x ".$pipeline_dir."shuffle_sweeps.sh && ".$pipeline_dir."shuffle_sweeps.sh"."\n";

	    system($command_line);

	    my $pipeline_dirchopped = $pipeline_dir;
            chop $pipeline_dirchopped;

	    
	    for(my $e=1;$e<=100000000000000;$e++){
		sleep 5;
		my $fake_done = 0;
		while(<$pipeline_dirchopped/shuffled_sweep_files/*.gz>){
		    $fake_done += 1;
		}
		if($fake_done==scalar(@splitter_ssc)*$parameters{"FDR_number"}){
		    last;
		}
	    }

	    #simplify here

	    open OUT, ">".$pipeline_dir."simplify_rafale.sh";
	    my $counter_02 = 0;
	    while(<$pipeline_dirchopped/shuffled_sweep_files/*>){
		$counter_02 += 1;
		my $add = " & \n";
		if($counter_02==$simult_runs){
		    $add = "\n sleep 0.5"; 
		    $counter_02 = 0;
		}
		my $file1 = $_;
		my @splitter_file = split(/\//,$file1);
		my $sc_sff = scalar(@splitter_file);
		my $file_name = $splitter_file[$sc_sff-1];
		chop $file_name;
		chop $file_name;
		chop $file_name;
		my $simplify_line = $pipeline_dir."simplify_zipped_sweeps.pl ".$file1." ".$pipeline_dir."shuffled_sweep_files_simplified/".$file_name." ".$threshold_chain.$add;
		print OUT $simplify_line."\n";
	    }
	    close OUT;

	    my $run_simple_line = "chmod +x ".$pipeline_dir."simplify_rafale.sh && ".$pipeline_dir."simplify_rafale.sh"."\n";

	    system($run_simple_line);

	    for(my $e=1;$e<=100000000000000;$e++){
                sleep 5;
                my $fake_done = 0;
                while(<$pipeline_dirchopped/shuffled_sweep_files_simplified/*>){
		    my $file = $_;
		    my $ok_done = "no";
		    open(DATA,$file);
                    while(<DATA>){
			if($_ =~ /ENSG/){
                            $ok_done = "yes";
                        }
                    }
		    close DATA;
		    if($ok_done eq "yes"){
			$fake_done += 1;
		    }
                }
		
		my $fnum = scalar(@splitter_ssc)*$parameters{"FDR_number"};

		#print $fake_done."\t".$fnum."\n";

                if($fake_done==$fnum){
                    last;
		}
            }

	}

	if($shuffle eq "no"){
	    print "Skipping shuffling. Make sure it was done before"."\n";
	}

	#count shuffled sweeps here.

	my $interrupted = $parameters{"Interrupted_FDR_run"};

	if($interrupted eq "no"){
	    my $remove_line = $pipeline_dir."remove.pl ".$pipeline_dir."test_outputs"."\n";
	    system($remove_line);
	}

	my %already_done;

	my $pipeline_dirchopped = $pipeline_dir;
	chop $pipeline_dirchopped;

	while(<$pipeline_dirchopped/test_outputs/fake*>){
	    my $file = $_;
	    my @splitter_dir = split(/\//,$file);
	    my $sc_sd = scalar(@splitter_dir);
	    my $file_name = $splitter_dir[$sc_sd-1];

	    open(DATA,$file);
	    while(<DATA>){
		chomp $_;
		if($_ =~ /yes/){
		    $already_done{$file_name} = "yes";
		}
	    }
	    close DATA;
	}

	for(my $j=1;$j<=100000000000000000;$j++){
	    sleep 1;
	    my %all;
	    my %all_nosuffix;

	    while(<$pipeline_dirchopped/shuffled_sweep_files_simplified/*>){
		my $file = $_;
		my @splitter_file = split(/\//,$file);
		my $sc_sl = scalar(@splitter_file);
		my $num = $splitter_file[$sc_sl-1];
		$all{$num} = "no";
	    }
	    my $done_counter = 0;
	    my $run_counter = 0;
	    my %done;
	    my @sizes = @splitter_ssc;

	    #my %already_done;

	    while(<$pipeline_dirchopped/test_outputs/fake*>){
		my $file = $_;
		my @splitter_dir = split(/\//,$file);
		my $sc_sd = scalar(@splitter_dir);
		my $file_name = $splitter_dir[$sc_sd-1];
		#$already_done{$file_name} = "yes";
		my @splitter_file = split(/\./,$file_name);
		my $part_file = $splitter_file[0];
		my @splitter_part = split("_",$part_file);
		my $num = $splitter_dir[$sc_sd-1];
		$all{$num} = "yes";
		$run_counter += 1;

		open(DATA,$file);
		while(<DATA>){
		    chomp $_;
		    if(($_ =~ /yes/)){
			$done_counter += 1;
			$done{$num} = "yes";
			$all{$num} = "yes";
		    }
		}
		close DATA;
	    }

	    my $tot = int($parameters{"FDR_number"}/8);

	    my $sc_sizes = scalar(@sizes);

	    if($done_counter>=$parameters{"FDR_number"}*$sc_sizes){
		last;
	    }

	    my $running_number = $run_counter - $done_counter;
	    my $free_number = $simult_runs-$running_number;
	    #print $free_number."\t".$running_number."\t".$done_counter."\n";
	    if($free_number>0){
		my $counter = 0;
		for(my $p=1;$p<=8;$p++){	  
		    for(my $q=1;$q<=$tot;$q++){
			for(my $d=0;$d<=$sc_sizes-1;$d++){
				my $size = $sizes[$d];
				my $file = "fake".$file_prefix."_".$size."_".$p."_".$q;
				if($all{$file} eq "no"){
				    $counter += 1;
				    #print $file."\n";

				    my $com_line = $pipeline_dir."count_sweeps_script.pl ".$pipeline_dir."VIPs/file_1 ".$pipeline_dir."nonVIPs/file_1 ".$pipeline_dir."shuffled_sweep_files_simplified/".$file." ".$pipeline_dir."ensembl_gene_neighbors.txt ".$pop." ".$pipeline_dir."test_outputs/ ".$count_sweeps." ".$threshold_chain." ".$pop_chain." ".$group_chain." &"."\n";

				    #if($already_done{$file} ne "yes"){
					system($com_line);
				    #}

				    if($counter==$free_number){
					last;
				    }
				}
				if($counter==$free_number){
				    last;
				}
			}
			if($counter==$free_number){
			    last;
			}
		    }
		    if($counter==$free_number){
			last;
		    }
		}
	    }
	}

	print "FDR analysis done. You can now estimate the FDR with the files in the test_outputs folder."."\n";

    }

}
