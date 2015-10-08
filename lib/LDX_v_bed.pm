package LDX_v_bed;
use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);

our @EXPORT = qw(LDX_v_bed heatmap_intersect);

	# This module intersects an LDX or LDXI file with a BED file.

sub LDX_v_bed
{
    my $LDXinput = shift;
    my $targetBED = shift;

    open LDX, $LDXinput or die $!;
    my $raBED =  "$$"."$^T".".txt";
    open( OUT, ">$raBED") or die "$raBED could not be opened: $!\n";
    my $raSNP = "";
    while(<LDX>){
	chomp;
	if( /raSNP/ or eof(LDX) ){
	# The end of an raSNP/ldSNP cluster
	# is marked by either another raSNP or the EOF.
	    if( eof(LDX) and !/raSNP/ ){
		print OUT $_."\t".$raSNP."\n";
	    }
	    if( !eof(LDX) ){ # new raSNP
		my @line  = split/ /;
		$raSNP = $line[2];
	    }
	} else {
	    print OUT $_."\t".$raSNP."\n";
	}
    }
    close OUT;
    close LDX;
    my $tally = `intersectBed -a $raBED -b $targetBED |cut -f 4|uniq |wc -l`;
    unlink $raBED;
    chomp $tally;
    return $tally;
}

sub heatmap_intersect
{
    my $LDXinput = shift;
    my $targetBED = shift;

    open LDX, $LDXinput or die $!;
    my $raBED =  "$$"."$^T".".txt";
    open( OUT, ">$raBED") or die "$raBED could not be opened: $!\n";
    my $raSNP = "";
    while(<LDX>){
        chomp;
        if( /raSNP/ or eof(LDX) ){
        # The end of an raSNP/ldSNP cluster                                                                                                                                                               
        # is marked by either another raSNP or the EOF.                                                                                                                                                   
            if( eof(LDX) and !/raSNP/ ){
                print OUT $_."\t".$raSNP."\n";
            }
            if( !eof(LDX) ){ # new raSNP                                                                                                                                                                  
                my @line  = split/ /;
                $raSNP = $line[2];
            }
        } else {
            print OUT $_."\t".$raSNP."\n";
        }
    }
    close OUT;
    close LDX;
    open(my $int, '-|', "bedtools intersect -a $raBED -b $targetBED |cut -f 4|uniq") or die "Could not intersect. Quitting.\n";
    my %tally;
    while (<$int>){
	chomp;
	$tally{$_}=1;
    }
    close $int;
    unlink $raBED;
#    chomp $tally;
    return \%tally;
}
