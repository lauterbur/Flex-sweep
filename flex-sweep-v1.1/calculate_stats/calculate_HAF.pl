#!/usr/bin/perl

#use warnings;

use strict;

use Compress::Zlib;

my ($hap_file,$map_file,$out_file,$size) = @ARGV;

my @pops = ("ACBD");

my $sc_p = scalar(@pops);
for(my $z=0;$z<=$sc_p-1;$z++){
    my @splitter_hap = split("_",$hap_file);
    my $sc_sh = scalar(@splitter_hap);
    $splitter_hap[$sc_sh-2] = $pops[$z];
    my @splitter_map = split("_",$map_file);
    my $sc_sm = scalar(@splitter_map);
    $splitter_map[$sc_sm-2] = $pops[$z];
    my @splitter_out = split("_",$out_file);
    my $sc_so = scalar(@splitter_out);
    $splitter_out[$sc_so-3] = $pops[$z];
    my $new_hap = "";
    my $new_map = "";
    my $new_out = "";
    for(my $i=0;$i<=$sc_sh-1;$i++){
	$new_hap .= $splitter_hap[$i]."_";
    }
    chop $new_hap;
    for(my $i=0;$i<=$sc_sm-1;$i++){
	$new_map .= $splitter_map[$i]."_";
    }
    chop $new_map;
    for(my $i=0;$i<=$sc_so-1;$i++){
	$new_out .= $splitter_out[$i]."_";
    }
    chop $new_out;
    $hap_file=$hap_file;
    $map_file=$map_file;
    $out_file=$out_file;
    #print $hap_file."\n";
    #print $map_file."\n";
    #print $out_file."\n";
    

    my %coords;
    my %true_coords;
    my %count_coords;
    my $counter_01 = 0;
    my $max_coord = 0;
    my $file = gzopen($map_file, "rb");
    while ($file->gzreadline($_)>0){
	chomp $_;
	my @splitter_line = split(" ",$_);
	$counter_01 += 1;
	my $coord = $splitter_line[3];
	$max_coord = $coord;
	my $int_coord = int($splitter_line[3]/100)*100;
	$coords{$int_coord} .= $coord." ";
	$true_coords{$coord} = $counter_01;
	$count_coords{$counter_01} = $coord;
    }
    $file->gzclose();

    my %haplos;
    my %derived;
    my $counter_02 = 0;
    my $file2 = gzopen($hap_file, "rb");
    while ($file2->gzreadline($_)>0){
	chomp $_;
	$counter_02 += 1;
	my $coord = $count_coords{$counter_02};
	$haplos{$coord} = $_;
	my @splitter_line = split("|",$_);
	my @grepper = grep(/1/,@splitter_line);
	my $derived_freq = scalar(@grepper)/scalar(@splitter_line);
	if(($derived_freq>=0.0000001)and($derived_freq<=0.9999999)){
	    $derived{$coord} = $derived_freq;
	}
    }
    $file2->gzclose();

    open OUT, ">".$out_file;

    for(my $i=0;$i<=0;$i++){

	my $inf=600000-$size/2;
	my $sup=600000+$size/2;
	
	#if($sup>=$max_coord){
	    #last;
	#}

	my %HAF_num;
	my $pol_num=0;

	my %HAF_den;
    
	for(my $j=$inf;$j<=$sup;$j+=100){
	    my $coords_chain = $coords{$j};
	    chop $coords_chain;
	    if($coords_chain ne ""){
		my @splitter_chain = split(" ",$coords_chain);
		my $sc_sc = scalar(@splitter_chain);
		for(my $k=0;$k<=$sc_sc-1;$k++){
		    $pol_num += 1;
		    my $coord = $splitter_chain[$k];
		    my $haplos = $haplos{$coord};
		    my $freq = $derived{$coord};
		    my @splitter_haplos = split("|",$haplos);
		    my $sc_haplo = scalar(@splitter_haplos);
		    for(my $l=0;$l<=$sc_haplo-1;$l++){
			my $state=$splitter_haplos[$l];
			if($state eq "1"){
			    $HAF_num{$l} += $freq;
			    $HAF_den{$l} += 1;
			}
		    }
		}
	    }
	}

	my @HAF = ();

	my $key_01;

	foreach $key_01 (sort keys %HAF_num){
	    if($HAF_den{$key_01}>0){
		my $value = $HAF_num{$key_01}/$HAF_den{$key_01};
		push(@HAF,$value);
	    }
	}

	my @sorted_HAF = sort{$a<=>$b}@HAF;

	my $sc_HAF = scalar(@HAF);
    
	my $ind1 = int(0.1*$sc_HAF);
	my $ind2 = int(0.9*$sc_HAF);

	my $high_HAF = 0;
	my $low_HAF = 0;

	for(my $m=0;$m<=$ind1;$m++){
	    $low_HAF += $sorted_HAF[$m];
	}

	for(my $m=$ind2;$m<=$sc_HAF-1;$m++){
	    $high_HAF += $sorted_HAF[$m];
	}
    
	%HAF_num = ();
	%HAF_den = ();

	print OUT "600000"." ".$high_HAF."\n";
	#print "600000"." ".$high_HAF."\n";
    }
    close OUT;
}
