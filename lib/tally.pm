package tally;
use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);

our @EXPORT = qw(tally);

sub tally
{
    my $infile = shift;
    my $outfile = shift;
    my $tally = 0;
    my $inside = 0;
    my $raSNP = "";
    my $highestTally = 0;
    open (LDX, "<$infile") or die $!;
    open (OUT, ">$outfile") or die $!;
    while (<LDX>){
	chomp;
	if( /raSNP/ or eof ){
	    # The end of an raSNP/ldSNP cluster
	    # is marked by either another raSNP or the EOF.
	    if( eof and !/raSNP/ ){
		# The last ldSNP in the file.
		$tally++;
	    }
	    if( $inside ){
		print OUT "$raSNP $tally\n";
	    }
	    if( !eof ){ # new raSNP
		my @line  = split/ /;
		$raSNP  = $line[2];
		$tally = 0;
		$inside = 1;
	    }
	} else {
	    #a ldSNP
	    $tally++;
	    $highestTally = $tally if $tally > $highestTally;
	}
    }
    close OUT;
    close LDX;
    return $highestTally;
}
