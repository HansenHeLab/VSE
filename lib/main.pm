package main;
use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);

our @EXPORT = qw(LDX_to_LDXI);

sub LDX_to_LDXI
{
    my @line;
    my $raSNP;
    my $ldSNP;
    my $ldx;
    my $ldxi;
    my $file = shift;
    my $outfile = shift;
    open LDX, $file or die $!;
    while(<LDX>){
	chomp;
	if( /raSNP/ ) { # New raSNP.
	    @line  = split/ /;
	    $raSNP = $line[2]; # raSNP rsid
	} else {
	    # Add raSNPs to ldSNP hash (usually it's the other way around).
	    push @{ $ldx->{ $_ } }, $raSNP;
	}
    }
    close LDX;
    open (LDXI, ">$outfile") or die $!;
    for my $k ( keys %{ $ldx } ){
	push @{ $ldxi->{ join( "_", @{ $ldx->{ $k } }) }}, $k;
    }
    for my $k ( keys %{ $ldxi } ){
	print LDXI "# raSNP ".$k."\n".join( "\n", @{ $ldxi->{ $k } })."\n";
    }
    close LDXI;
}

sub tally
{
    my $infile = shift;
    my $outfile = shift;
    my $tally = 0;
    my $inside = 0;
    my $raSNP = "";
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
	}
    }
}
