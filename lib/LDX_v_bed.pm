package LDX_v_bed;
use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);

our @EXPORT = qw(LDX_v_bed);

        # This module intersects an LDX or LDXI file with a BED file.

sub LDX_v_bed
{
    my $LDXinput = shift;
    my $POS2ref = shift;
    my %POS2 = %$POS2ref;
    my %POS1;

    foreach my $n (1..22,"X","Y"){
	my $chr ="chr".$n; my @TMP; 
	@{$POS1{$chr}} = @TMP;
    }
    my $raSNP = "";
    open (LDX, "<$LDXinput") or die $!;
    
    while(<LDX>){
	chomp;
	my ($chr,$s,$e) = split (/\t/, $_) if ($_ !~ m/raSNP/);
	if( /raSNP/ or eof(LDX) ){ # The end of an raSNP/ldSNP cluster is marked by either another raSNP or the EOF.
            if( eof(LDX) and !/raSNP/ ){
		my $rec = {S=>$s, N=>$raSNP};
		push @{$POS1{$chr}},$rec;
            }
            if( !eof(LDX) ){ # new raSNP
                my @line  = split/ /;
                $raSNP = $line[2];
            }
        } else {
	    my $rec = {S=>$s, N=>$raSNP};
	    push @{$POS1{$chr}},$rec;
        }
    }
    close LDX;

    foreach my $chr (keys %POS1){
	@{$POS1{$chr}} = sort {$a->{S} <=> $b->{S}} @{$POS1{$chr}};
    }

    my %seenENTRY;
    my $tallyCount = 0;

    foreach my $chr (keys %POS1){
	my @P1 = @{$POS1{$chr}};
	next if !$POS2{$chr};
	my @P2=@{$POS2{$chr}};
	my $n=0; #keep track of the ending position of last entry in the same chromosome
	foreach my $p1 (@P1){ #each SNP
	    my $s =$p1->{S}; 
	    my $rsID = $p1->{N};
	    if ($seenENTRY{$rsID}){next} #skip identical entries
	    my ($fr,$lr)=($n,$#P2);
	    my $i=$n;
	    while ($fr <= $lr){
		$i= int(($fr+$lr)/2);
		if ($P2[$i]->{S} > $s){
		    $lr =$i-1;
		} elsif ($P2[$i]->{E} < $s){
		    $fr =$i+1;
		} elsif ($s>=$P2[$i]->{S} && $s<=$P2[$i]->{E}){
		    $tallyCount++;
		    $seenENTRY{$rsID}=1;
		    $n=$i; 
		    last;
		}
	    }
	}
    }
    return ($tallyCount,\%seenENTRY);
}
