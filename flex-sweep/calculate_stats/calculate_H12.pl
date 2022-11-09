#!/usr/bin/perl

#use warnings;

use strict;

use Compress::Zlib;

my ($hap_file,$map_file,$out_file,$size) = @ARGV;
	# $size = windowsize
my @pops = ("CLMD");
#my @pops = ("ACBD","ASWD","ESND","GWDD","LWKD","YRID","MSLD","CEUD","GBRD","FIND","TSID","IBSD","CHBD","CHSD","CDXD","KHVD","JPTD","GIHD","BEBD","PJLD","ITUD","STUD","MXLD","CLM","PURD","PELD");
#my @pops = ("YRID","MSLD","CEUD","GBRD","FIND","TSID","IBSD","CHBD","CHSD");
#my @pops = ("JPTD","GIHD","BEBD","PJLD","ITUD");
my $sc_p = scalar(@pops);
for(my $z=0;$z<=$sc_p-1;$z++){

    my %okfreq;

    my $counter_03 = 0;
    my $file3 = gzopen($hap_file, "rb");
    while ($file3->gzreadline($_)>0){
        #open(DATA,$hap_file);                                                                                                                                                                                                               
        #while (<DATA>){                                                                                                                                                                                                                     
        chomp $_;
        $counter_03 += 1;
        my @splitter_line = split("|",$_);
        my @grepper = grep(/1/,@splitter_line);
        my $derived_freq = scalar(@grepper)/scalar(@splitter_line);                                                                                                                                                                                 
	if(($derived_freq>=0.05)and($derived_freq<=1)){
            $okfreq{$counter_03} = "yes";
        }
    }
    $file3->gzclose();

    
    my %coords;

    my %true_coords;

    my %count_coords;

    my $counter_01 = 0;

    my $file = gzopen($map_file, "rb");
    while ($file->gzreadline($_)>0){
	#open(DATA,$map_file);
	#while (<DATA>){
	chomp $_;
	my @splitter_line = split(" ",$_);
	$counter_01 += 1;
	if($okfreq{$counter_01} eq "yes"){
	    my $coord = $splitter_line[3];
	    my $int_coord = int($splitter_line[3]/100)*100;
	    $coords{$int_coord} .= $coord." ";
	    $true_coords{$coord} = $counter_01;
	    $count_coords{$counter_01} = $coord;
	}
    }
    #close DATA;
    $file->gzclose();
    
    
    my %haplos;

    my %derived;

    my $counter_02 = 0;
    my $file2 = gzopen($hap_file, "rb");
    while ($file2->gzreadline($_)>0){
	#open(DATA,$hap_file);
	#while (<DATA>){
	chomp $_;
	$counter_02 += 1;
	if($okfreq{$counter_02} eq "yes"){
	    my $coord = $count_coords{$counter_02};
	    $haplos{$coord} = $_;
	    my @splitter_line = split("|",$_);
	    my @grepper = grep(/1/,@splitter_line);
	    my $derived_freq = scalar(@grepper)/scalar(@splitter_line);
	    #print $derived_freq."\n";
	    if(($derived_freq>=0.25)and($derived_freq<=0.95)){
		#$derived{$coord} = $derived_freq;
	    }
	}
    }
    $file2->gzclose();

    #close DATA;

    $derived{"600000"} = 1;
    $coords{"600000"} = 600000;

    my %maxhaplos;

    my %secondhaplos;

    my %thirdhaplos;

    my %keep_haplo_freq;

    my $key_001;

    foreach $key_001 (sort{$a<=>$b} keys %derived){

	#print $key_001."\n";

	my $pi_d = 0;
        my $pi_a = 0;
        my $coord = $key_001;
        #print $key_01."\n";                                                                                                                                                                                                                                                                                                                                               
        my $int_coord = int($coord/100)*100;
        my $inf = $int_coord-$size/2;
        my $sup = $int_coord+$size/2;
        my $hap_line = "1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111";
        my @hap = split("|",$hap_line);
        my @grepper = grep(/1/,@hap);
        my $n_d = scalar(@grepper);
        my $n_a = scalar(@hap)-$n_d;
        my $s_a = 0;
        my $s_d = 0;
        my $diff_num = 0;
        my $diff_den = 0;
        my $diff_num2 = 0;
        my $diff_den2 = 0;
        my $low_bound = 0;
        my $up_bound = 0;

        my $snp_left = 0;
        my $snp_right = 0;

	my %ongoing_haplos;

        for(my $i=1;$i<=$size/200;$i++){

	    my $coo = $i*100;

            my $inf_i = $int_coord-$i*100;

            $low_bound = $inf_i;
            my $finish = "no";

            if($inf_i<=0){
                $finish = "yes";
            }

	    if($coords{$inf_i} ne ""){
                my $chain = $coords{$inf_i};

		#print $chain."\n";
		
                my @splitter_chain = split(" ",$chain);
                my $sc_sc = scalar(@splitter_chain);
                for(my $j=0;$j<=$sc_sc-1;$j++){
                    my $true_coord = $splitter_chain[$j];
                    if($true_coord ne $coord){

                        my $haplotype = $haplos{$true_coord};
                        my @current_haplo = split("|",$haplotype);
                        my $sc_h = scalar(@current_haplo);
                        for(my $k=0;$k<=$sc_h-1;$k++){

                            if($hap[$k]==1){
                                $ongoing_haplos{$k} .= $current_haplo[$k]." ";
                            }
			}
		    }

		}
	    }

	    if(($coo>=$size/2)){
		$finish = "yes";
            }

	    if($finish eq "yes"){
		last;
	    }
	}

	for(my $i=1;$i<=$size/200;$i++){

            my $sup_i = $int_coord+$i*100;

	    my $coo = $i*100;

	    $up_bound = $sup_i;
	    my $finish = "no";

	    if($sup_i>=1200000){
		$finish = "yes";
	    }

            if($coords{$sup_i} ne ""){
                my $chain = $coords{$sup_i};
                my @splitter_chain = split(" ",$chain);
                my $sc_sc = scalar(@splitter_chain);
                for(my $j=0;$j<=$sc_sc-1;$j++){
                    my $true_coord = $splitter_chain[$j];
                    if($true_coord ne $coord){

                        my $haplotype = $haplos{$true_coord};
                        my @current_haplo = split("|",$haplotype);
                        my $sc_h = scalar(@current_haplo);
                    
                        for(my $k=0;$k<=$sc_h-1;$k++){

                            if($hap[$k]==1){
                                $ongoing_haplos{$k} .= $current_haplo[$k]." ";
                            }
			}
                    }

		}
            }

	    if(($coo>=$size/2)){
		$finish = "yes";
	    }

            if($finish eq "yes"){
		last;
            }
        }
	
	my %haplos_number;
	my $key_ongo;
	foreach $key_ongo (sort keys %ongoing_haplos){
	    my $haplo = $ongoing_haplos{$key_ongo};
	    $haplos_number{$haplo} += 1;
	}

	my $max_haplo = "";
	my $second_haplo = "";
	my $third_haplo = "";

	my %best_haplos;

	my %revert_number;

	my $key_numb;
	foreach $key_numb (sort keys %haplos_number){
	    my $number = $haplos_number{$key_numb};

	    #if($revert_number{$number} eq ""){
		$revert_number{$number} .= $key_numb."_";
	    #}
	}

	my $counter_rev = 0;

	my $done_rev = 0;

	my $key_rev;

	#my %keep_haplo_freq;
	
	foreach $key_rev (sort{$b<=>$a} keys %revert_number){
	
	    #print $key_001."\t".$key_rev."\n";

	    my $chain = $revert_number{$key_rev};

	    my @splitter_chain = split("_",$chain);

	    my $sc_ch = scalar(@splitter_chain);

	    for(my $f=0;$f<=$sc_ch-1;$f++){
		$done_rev += 1;
		
		$best_haplos{$done_rev} = $splitter_chain[$f];

		$keep_haplo_freq{$done_rev} = $key_rev;

		if($done_rev==1){
		    #$max_haplo = $splitter_chain[$f];
		}
		if($done_rev==2){
		    #$second_haplo = $splitter_chain[$f];
		}

		if($done_rev==3){
                    #$third_haplo = $splitter_chain[$f];
                }

	    }

	    $counter_rev += $done_rev;

	    if($counter_rev>=10){
		last;
	    }
	}

	my %similar_pairs;

	my %done;

	my $key_compf;

	foreach $key_compf (sort{$a<=>$b} keys %best_haplos){
	    $similar_pairs{$key_compf} = "";
	}

	my $key_comp;

	foreach $key_comp (sort{$a<=>$b} keys %best_haplos){

	    my $key_comp2;

	    foreach $key_comp2 (sort{$a<=>$b} keys %best_haplos){

		if(($key_comp2 ne $key_comp)and($done{$key_comp." ".$key_comp2} ne "yes")){

		    my $haplo_1 = $best_haplos{$key_comp};
		    my $haplo_2 = $best_haplos{$key_comp2};

		    my @splitter_haplo_1 = split(" ",$haplo_1);
		    my @splitter_haplo_2 = split(" ",$haplo_2);

		    my $sc_sh = scalar(@splitter_haplo_1);

		    my $identical = 0;
		    my $different = 0;

		    my $total = 0;

		    for(my $f=0;$f<=$sc_sh-1;$f++){
			if(($splitter_haplo_1[$f]==1)and($splitter_haplo_2[$f]==1)){
			    $identical += 1;
			    $total += 1;
			}
			if(($splitter_haplo_1[$f]==0)and($splitter_haplo_2[$f]==1)){
                            $different += 1;
                            $total += 1;
                        }
			if(($splitter_haplo_1[$f]==1)and($splitter_haplo_2[$f]==0)){
                            $different += 1;
                            $total += 1;
                        }
		    }

		    if(($different/$total<=0.2)){

			$similar_pairs{$key_comp} .= $key_comp2." "; 
			$done{$key_comp2." ".$key_comp} = "yes";
			$done{$key_comp." ".$key_comp2}= "yes";
		    }

		}

	    }

	}

	my $counter_rev2 = 0;
        my $done_rev2 = 0;
        my $key_rev2;

	my %exclude;

        foreach $key_rev2 (sort{$a<=>$b} keys %similar_pairs){

	    #print $key_rev2." --- ".$similar_pairs{$key_rev2}."\n";

            if($exclude{$key_rev2} ne "yes"){

		my $chain = $best_haplos{$key_rev2};

		my $similar = $similar_pairs{$key_rev2};

		#print $key_rev2." - ".$similar."\n";

		if($similar ne ""){

		    my @splitter_similar = split(" ",$similar);

		    my $sc_si = scalar(@splitter_similar);

		    for(my $f=0;$f<=$sc_si-1;$f++){
			my $cur_rev = $splitter_similar[$f];
			$exclude{$cur_rev} = "yes";
			$chain .= "_".$best_haplos{$cur_rev};
		    }
		}

		$counter_rev2 += 1;

		if($counter_rev2==1){
		    $max_haplo = $chain;
		}
		if($counter_rev2==2){
                    $second_haplo = $chain;
                }
		if($counter_rev2==3){
                    $third_haplo = $chain;
                }

		if($counter_rev2==3){
		    #print "done"."\n";
		    last;
		}
	    }
        }


	my $freq_1 = 0;

	my $freq_2 = 0;

	my $freq_3 = 0;

	my $toto = 0;

	my $key_ongo2;
        foreach $key_ongo2 (sort keys %ongoing_haplos){

	    my $ongoing = $ongoing_haplos{$key_ongo2};

	    $toto += 1;

	    if($max_haplo =~ $ongoing){
		$maxhaplos{$key_001} .= "_".$key_ongo2."_";
		$freq_1 += 1;
	    }

	    if($second_haplo =~ $ongoing){
		$secondhaplos{$key_001} .= "_".$key_ongo2."_";
		$freq_2 += 1;
            }

	    if($third_haplo =~ $ongoing){
                $thirdhaplos{$key_001} .= "_".$key_ongo2."_";
		$freq_3 += 1;
            }
	}

	my $H12 = ($freq_1/$toto+$freq_2/$toto)*($freq_1/$toto+$freq_2/$toto);

	open OUT, ">".$out_file;

	print OUT $out_file."\t".$H12."\n";
	#print $out_file."\t".$H12."\n";

	close OUT;

    }
    
}

    #my $key_01;

    #foreach $key_01 (sort{$a<=>$b} keys %derived){
	#my $pi_d = 0;
	#my $pi_a = 0;
	#my $coord = $key_01;
	#print $key_01."\n";
	#my $int_coord = int($coord/100)*100;
	#my $inf = $int_coord-$size/2;
	#my $sup = $int_coord+$size/2;
	#my $hap_line = $haplos{$coord};
	#my @hap = split("|",$hap_line);
	#my @grepper = grep(/1/,@hap);
	#my $n_d = scalar(@grepper);
	#my $n_a = scalar(@hap)-$n_d;
	#my $s_a = 0;
	#my $s_d = 0;
	#my $diff_num = 0;
	#my $diff_den = 0;
	#my $diff_num2 = 0;
	#my $diff_den2 = 0;
	#my $low_bound = 0;
	#my $up_bound = 0;

	#my $snp_left = 0;
	#my $snp_right = 0;
	
	#for(my $i=1;$i<=$size/200;$i++){
	    
	    #my $coo = $i*100;

	    #print $coo."\n";
	    
	    #my $inf_i = $int_coord-$i*100;

	    #$low_bound = $inf_i;
	    #my $finish = "no";

	    #if($inf_i<=0){
		#$finish = "yes";
	    #}
	    
	    #if($coords{$inf_i} ne ""){
		#my $chain = $coords{$inf_i};
		#my @splitter_chain = split(" ",$chain);
		#my $sc_sc = scalar(@splitter_chain);
		#for(my $j=0;$j<=$sc_sc-1;$j++){
		    #my $true_coord = $splitter_chain[$j];
		    #if($true_coord ne $coord){

			#$snp_left += 1;
			
			#my $haplotype = $haplos{$true_coord};
			#my @current_haplo = split("|",$haplotype);
			#my $sc_h = scalar(@current_haplo);
			#my $freq_a = 0;
			#my $freq_d = 0;
			#my $den_d = 0;
			#my $freq_d2 = 0;
			#my $den_d2 = 0;
			#my $freq_d3 = 0;
                        #my $den_d3 = 0;
			#my $freq_dall = 0;


			#my $tot_freq = 0;
			#for(my $k=0;$k<=$sc_h-1;$k++){
			    #if($hap[$k]==0){
				#$freq_a += $current_haplo[$k];
			    #}

			    #my $stress1 = "_".$k."_";

			    #if(($hap[$k]==1)and($maxhaplos{$key_01} =~ $stress1)){
				#$freq_d += $current_haplo[$k];
				#$den_d += 1;
			    #}

			    #if(($hap[$k]==1)and($secondhaplos{$key_01} =~ $stress1)){
                                #$freq_d2 += $current_haplo[$k];
                                #$den_d2 += 1;
                            #}

			    #if(($hap[$k]==1)and($thirdhaplos{$key_01} =~ $stress1)){
                                #$freq_d3 += $current_haplo[$k];
                                #$den_d3 += 1;
                            #}

			    #if(($hap[$k]==1)){
				#$freq_dall += $current_haplo[$k];
			    #}

			#}

			#print $den_d."\t".$den_d2."\t".$den_d3."\t".$n_d."\n";

			#$tot_freq = ($freq_a+$freq_dall)/($n_a+$n_d);

			#$freq_a = $freq_a/$n_a;
			#if($den_d>0){
			    #$freq_d = $freq_d/$den_d;
			#}
			#if($den_d2>0){
			    #$freq_d2 = $freq_d2/$den_d2;
			#}
			#if($den_d3>0){
			    #$freq_d3 = $freq_d3/$den_d3;
			#}

			#if($freq_a>0){
			    #$s_a += 1;
			#}
			#if($freq_d>0){
			    #$s_d += 1;
			#}

			#$diff_num += ($freq_d+$freq_d2)*($freq_d+$freq_d2);
			#$diff_den += 1;


			#$freq_dall=$freq_dall/$n_d;
			#$freq_dall=($freq_d+$freq_d2+$freq_d3)/($den_d+$den_d2+$den_d3);

			#if(($freq_d==$freq_d2)and($freq_d3==$freq_d2)){

			    #if(($freq_dall>$freq_a)and($freq_dall<=1)and($freq_a<=0.25)and($tot_freq>=0.25)){
				#$diff_num += (($freq_dall)*($freq_dall))-(($freq_a)*($freq_a));
				#$diff_den += 1;
			    #}
			#}

			#else{

			    #$diff_den += 1;

			    #if(($freq_d>$freq_a)and($freq_d<=1)and($freq_a<=0.25)and($tot_freq>=0.25)and($den_d>=0.25*$n_d)and($den_d>=10)){
				#my $oppfreq_d = 1-$freq_d;
				#my $oppfreq_a = 1-$freq_a;
				#3if($tot_freq>=0){
				
				    #if($den_d>=0*$n_d){
					#$diff_num += (($freq_d)*($freq_d))-(($freq_a)*($freq_a));
				    #}
				    #$diff_den += 1;
				    #$diff_num2 += ($freq_a)*($freq_a);
				#}
				#if(($coo>=$size/2)or($inf_i<=0)){
				    #$finish = "yes";
				#}
			    #}


			    #if(($freq_d2>$freq_a)and($freq_d2<=1)and($freq_a<=0.25)and($tot_freq>=0.25)and($den_d2>=0.25*$n_d)and($den_d2>=10)){
				#if($tot_freq>=0){
				 #   if($den_d2>=0*$n_d){
					#$diff_num += (($freq_d2)*($freq_d2))-(($freq_a)*($freq_a));
				  #  }
				    #$diff_den += 1;
				#}
				#if(($coo>=$size/2)or($inf_i<=0)){
				 #   $finish = "yes";
				#}
			    #}

			    #if(($freq_d3>$freq_a)and($freq_d3<=1)and($freq_a<=0.25)and($tot_freq>=0.25)and($den_d3>=0.25*$n_d)and($den_d3>=10)){
				#if($tot_freq>=0){
				 #   if($den_d3>=0*$n_d){
					#$diff_num += (($freq_d3)*($freq_d3))-(($freq_a)*($freq_a));
				  #  }   
				    #$diff_den += 1;
				#}
				#if(($coo>=$size/2)or($inf_i<=0)){
				 #   $finish = "yes";
				#}
			    #}
			#}


			#if(($freq_a>0)and($freq_a<=0.5)and($freq_d>=0)){
			    #my $oppfreq_d =1-$freq_d;
			    #my $oppfreq_a =1-$freq_a;
			    #if($tot_freq>=0){
				#$diff_num2 += (((1-$freq_a)*(1-$freq_a))); #-((1-$freq_a)*(1-$freq_a)));
				#$diff_den2 += 1;
			    #}
			    #if(($coo>=$size/2)or($inf_i<=0)){
				#$finish = "yes";
			    #}
			#}
			

			#if(($n_a>=2)and($n_d>=2)){
			 #   $pi_a += $freq_a*(1-$freq_a)*$n_a/($n_a-1);
			  #  $pi_d += $freq_d*(1-$freq_d)*$n_d/($n_d-1);
			#}
		    #}
		#}
	    #}
	    #if($finish eq "yes"){
		#last;
	    #}
	#}

	#for(my $i=1;$i<=$size/200;$i++){
	    
        #my $sup_i = $int_coord+$i*100;

	#my $coo = $i*100;
	
	#$up_bound = $sup_i;
        #my $finish = "no";

	#if($sup_i>=1200000){
	 #   $finish = "yes";
	#}
       
	#if($coords{$sup_i} ne ""){
         #   my $chain = $coords{$sup_i};
          #  my @splitter_chain = split(" ",$chain);
          #  my $sc_sc = scalar(@splitter_chain);
           # for(my $j=0;$j<=$sc_sc-1;$j++){
            #    my $true_coord = $splitter_chain[$j];
             #   if($true_coord ne $coord){

		#    $snp_right += 1;
		    
                 #   my $haplotype = $haplos{$true_coord};
                  #  my @current_haplo = split("|",$haplotype);
                   # my $sc_h = scalar(@current_haplo);
                   # my $freq_a = 0;
                   # my $freq_d = 0;
		   # my $den_d = 0;
		   # my $freq_d2 = 0;
                   # my $den_d2 = 0;
		   # my $freq_d3 = 0;
                   # my $den_d3 = 0;
                   # my $tot_freq = 0;
		   # my $freq_dall = 0;

                    #for(my $k=0;$k<=$sc_h-1;$k++){
                        #if($hap[$k]==0){
                         #   $freq_a += $current_haplo[$k];
                        #}

			#my $stress1 = "_".$k."_";

                        #if(($hap[$k]==1)and($maxhaplos{$key_01} =~ $stress1)){
                         #   $freq_d += $current_haplo[$k];
			  #  $den_d += 1;
                        #}

			#if(($hap[$k]==1)and($secondhaplos{$key_01} =~ $stress1)){
			 #   $freq_d2 += $current_haplo[$k];
			  #  $den_d2 += 1;
			#}

			#if(($hap[$k]==1)and($thirdhaplos{$key_01} =~ $stress1)){
                         #   $freq_d3 += $current_haplo[$k];
                          #  $den_d3 += 1;
                        #}
			
			#if(($hap[$k]==1)){
			 #   $freq_dall += $current_haplo[$k];
			#}


                    #}
                    #$tot_freq = ($freq_a+$freq_dall)/($n_a+$n_d);

		    #print $freq_d."-".$den_d."\t".$freq_d2."-".$den_d2."\t".$freq_d3."-".$den_d3."\t".$freq_dall."-".$n_d."\n";

		    #if($tot_freq>0.5){

			#$tot_freq = 1-$tot_freq;

		   # }
		    
                   # $freq_a = $freq_a/$n_a;
		    #if($den_d>0){
			#$freq_d = $freq_d/$den_d;
		    #}
		    #if($den_d2>0){
			#$freq_d2 = $freq_d2/$den_d2;
		    #}
		    #if($den_d3>0){
			#$freq_d3 = $freq_d3/$den_d3;
		    #}
                    #if($freq_a>0){
                     #   $s_a += 1;
                    #}
                    #if($freq_d>0){
                     #   $s_d += 1;
                    #}


		    #$diff_num += ($freq_d+$freq_d2)*($freq_d+$freq_d2);
		    #$diff_den += 1;

		    #$freq_dall=$freq_dall/$n_d;
		    #if(($freq_d==$freq_d2)and($freq_d3==$freq_d2)){

			#if(($freq_dall>$freq_a)and($freq_dall<=1)and($freq_a<=0.25)and($tot_freq>=0.25)){
			    #$diff_num += (($freq_dall)*($freq_dall))-(($freq_a)*($freq_a));
			    #$diff_den += 1;
			#}
		    #}

		    #else{

			#$diff_den += 1;

			#print $freq_d."-".$den_d."\t".$freq_d2."-".$den_d2."\t".$freq_d3."-".$den_d3."\t".$freq_dall."-".$n_d."\n";

			#if(($freq_d>$freq_a)and($freq_d<=1)and($freq_a<=0.25)and($tot_freq>=0.25)and($den_d>=0.25*$n_d)and($den_d>=10)){

			 #   my $oppfreq_d = 1-$freq_d;
			  #  my $oppfreq_a = 1-$freq_a;
			
			   # if($tot_freq>=0){
				#$diff_num += (($freq_d*$freq_d)-($freq_a*$freq_a));

				#if($den_d>=0*$n_d){
				    #$diff_num += (($freq_d)*($freq_d))-(($freq_a)*($freq_a));
				#}
				#$diff_den += 1;
				#$diff_num2 += ($freq_a)*($freq_a);
			    #}
			
			    #if(($coo>=$size/2)or($sup_i>=260000000)){
				#$finish = "yes";
			    #}
			#}

			#if(($freq_d2>$freq_a)and($freq_d2<=1)and($freq_a<=0.25)and($tot_freq>=0.25)and($den_d2>=0.25*$n_d)and($den_d2>=10)){

			 #   if($tot_freq>=0){
				#if($den_d2>=0*$n_d){
				    #$diff_num += (($freq_d2)*($freq_d2))-(($freq_a)*($freq_a));
				#}
				#$diff_den += 1;
				#$diff_num2 += ($freq_a)*($freq_a);
			    #}

			    #if(($coo>=$size/2)or($sup_i>=260000000)){
				#$finish = "yes";
			    #}
			#}

			#if(($freq_d3>$freq_a)and($freq_d3<=1)and($freq_a<=0.25)and($tot_freq>=0.25)and($den_d3>=0.25*$n_d)and($den_d3>=10)){

			    #if($tot_freq>=0){
				#if($den_d3>=0*$n_d){
				    #$diff_num += (($freq_d3)*($freq_d3))-(($freq_a)*($freq_a));
				#}
				#$diff_den += 1;
				#$diff_num2 += ($freq_a)*($freq_a);
			    #}

			    #if(($coo>=$size/2)or($sup_i>=260000000)){
				#$finish = "yes";
			    #}
			#}
		    #}


		    #if(($freq_a>$freq_d)and($freq_a<=1)and($freq_d>=0)){
			#my $oppfreq_d =1-$freq_d;
			#my $oppfreq_a =1-$freq_a;
			#if($tot_freq>=0){
			    #$diff_num2 += (((1-$freq_a)*(1-$freq_a))); #-((1-$freq_a)*(1-$freq_a)));
			    #$diff_den2 += 1;
			#}
			#if(($coo>=$size/2)or($sup_i>=260000000)){
			    #$finish = "yes";
			#}
		    #}
		    
		    #if(($n_a>=2)and($n_d>=2)){
		    
			#$pi_a += $freq_a*(1-$freq_a)*$n_a/($n_a-1);
			#$pi_d += $freq_d*(1-$freq_d)*$n_d/($n_d-1);
		    #}
                #}
            #}
        #}
        #if($finish eq "yes"){
            #last;
        #}
    #}

	#if($pi_a==0){
	 #   $pi_a = 0.001;
	#}

    #if($pi_d==0){

	#if($pi_a==0.001){
	 #   $pi_d = 1*$pi_a;
	#}
	#if($pi_a>0.001){
	 #   $pi_d = 0.1*$pi_a;
	#}

    #}

	#print $diff_den."\n";
	
    #if(($diff_den>0)){


	#my $iDAF = ($s_a/($s_d+0.001));

	#my $iDAF = $diff_num/$diff_den;

	#my $rat = $iDAF/$iDAF2;
	
	#my $iDAF= $pi_a/$pi_d;
	
	#my $dind = ($pi_a)/($pi_d);

	#print OUT $coord."\t".$dind."\t".$pi_a."\t".$pi_d."\t".$n_a."\t".$n_d."\t".$s_a."\t".$s_d."\n";

	#print OUT $coord."\t".$diff_num."\t".$diff_den."\t".$derived{$key_01}."\t".$low_bound."\t".$up_bound."\n";

	#print OUT $coord."\t".$iDAF."\t".$derived{$key_01}."\t".$diff_den."\n";

	#print $coord."\t".$iDAF."\t".$derived{$key_01}."\n";
	
	#print $coord."\t".$diff_num."\t".$diff_den."\t".$derived{$key_01}."\t".$low_bound."\t".$up_bound."\n";
	
	#print $coord."\t".$diff_num."\t".$diff_den."\n";

	#print $coord."\t".$dind."\t".$pi_a."\t".$pi_d."\t".$n_a."\t".$n_d."\t".$s_a."\t".$s_d."\n";

    #}

    #}

#close OUT;

#}
