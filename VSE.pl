#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Carp;
use File::Slurp;
use File::Basename;
use DBI;
use lib 'lib';
use LDX_to_LDXI;
use tally;
use LDX_v_bed;
my %opt;
getopts('r:hvYf:l:s:d:p:n:A:', \%opt);

sub usage
{
    print "This tool will generate MRVs\n";
    print "usage: vse.sh -f snpListBed -s suffix -d dirLocation [-i shuffleNo] [-r r2Value] [-v y/n] [-Y y/n] | [-h]\n";
    print "Options:\n";
    print "-i[int]      No of shuffles; default: 100\n";
    print "-r[float]    R2 value to find SNPs in LD; default: 0.8\n";
    print "-v[y/n]      Verbose; default: y\n";
    print "-Y[y/n]      To analyze chrY or not; default: n\n";
    print "-f[path]     Location of tagSNP list. The file must be a bed file\n";
    print "-l[path]     Location of LD snp list. Must be in this format: chr\ts\te\tldSNP\ttagSNP\n";
    print "-s[char]     Suffix for filenames\n";
    print "-d[path]     Path to feature directory\n";
    print "-A[char]     Suffix for AVS/MRV files. Only used when -p is xml in order to avoid having to generate AVS/MRV again.\n";
    print "-p           Which part to run; default: all\n";
    print "-n           No. of iteration; default: 100\n";
    print "-h           This message\n";
}

sub printLog
{
    my $line = shift;
    if ($opt{v}){
	print localtime()." $line\n";
    }
}

my $dbfile = "null/snpinfo.db";
my $dsn      = "dbi:SQLite:dbname=$dbfile";
my $user     = "";
my $password = "";

die usage if ($opt{h});
$opt{n} = 100 if !$opt{n};
die "parameter missing\n" if (!$opt{f} || !$opt{l} || !$opt{s} || !$opt{d});
$opt{p} = "all" if !$opt{p};
$opt{i} = 100 if !$opt{i};
$opt{n} = 100 if !$opt{n};
$opt{r} = 0.8 if !$opt{r};
die "-r must be between 0.6 and 0.8\n" if $opt{r} <0.6 || $opt{r} > 0.8;
#INPUT QC
die "tagSNP file not found\n" if (! -e $opt{f});
#checking feature locations
if ($opt{d} =~ m/(bed|peak|gz)$/i){
    printLog("single bed file provided. Output will be printed on screen");
    die "$opt{d} could not be located\n" if (! -e $opt{d});
} else {
    if (! -d $opt{d}){
	print "Directory path not found\n";
    }
}
#-------------------#

my $totalTagSnp=`wc -l $opt{f} | cut -d " " -f 1`;

print "$opt{f} does not have any snp\n" if ($totalTagSnp == 0);

if ($opt{p} eq "AVS" || $opt{p} eq "all"){
    printLog("Preparing AVS file");
    mkdir $opt{s}.".AVS" if (! -d $opt{s}.".AVS" );
    my $AVS_out_file = $opt{s}.".AVS/".$opt{s}.".LDX.bed";
    if (-e $AVS_out_file){ unlink $AVS_out_file; }
    my @lines = read_file($opt{f}, chomp => 1);
    die "$opt{f} looks empty\n" if scalar @lines == 0;
    foreach my $line (@lines){
	my @f = split /\t/, $line;
	my $header= "# raSNP $f[3]\n";
	write_file($AVS_out_file, {append=>1}, $header);
	open (IN, "awk -v tag=$f[3] '\$5==tag' $opt{l} |") or die;
	while (<IN>){
	    chomp;
	    my @g = split /\t/;
	    if ($g[4] ne $g[3]){
		my $printline = "$g[0]\t$g[1]\t$g[2]\n";
		write_file($AVS_out_file, {append=>1}, $printline);
	    }
	}
	close IN;
    }
    printLog("$AVS_out_file created");
#Fixing tagSNP to LDXI and making tally
    printLog("Generating LDXI for AVS/$opt{s}.LDX.bed");
    LDX_to_LDXI($opt{s}.".AVS/$opt{s}.LDX.bed",$opt{s}.".AVS/$opt{s}.LDXI.bed");
    printLog($opt{s}.".AVS/$opt{s}.LDXI.bed created.");
    printLog("Tallying ".$opt{s}.".AVS/$opt{s}.LDXI.bed");
    tally($opt{s}.".AVS/$opt{s}.LDXI.bed",$opt{s}.".AVS/$opt{s}.LDXI.tally.txt");
    printLog($opt{s}.".AVS/$opt{s}.LDXI.tally.txt created.");
    printLog("Generating LOD MRVs");
    if ( ! -d "MRVs" ){
	mkdir "MRVs";
    }
}

if ($opt{p} eq "MRV" || $opt{p} eq "all"){
    if (! -d $opt{s}.".MRVs"){
	mkdir $opt{s}.".MRVs";
    }
    my $prefix = $opt{s}.".MRVs/";	
    my @SNPS;
    my $null_tally_file = "null/ld_tally_r".$opt{r}.".txt";
    open (SNPS, "<$null_tally_file") or die "$!\n";
    while(<SNPS>){ 
	chomp;
	my @snp = split /\t/;
	next if /SNP/;
	push @{ $SNPS[$snp[1]] }, $snp[0];
    }
    close SNPS;
    my @AV;
    my @RA;
    open (AV, "<", $opt{s}.".AVS/".$opt{s}.".LDXI.tally.txt") or die "$!\n";
    while (<AV>){
	chomp;
	print "$_\n" if eof(AV);
	my @snp = split/ /;
	push @AV, $snp[1];
	push @RA, $snp[0];
	foreach my $i (0..$#{ $SNPS[ $snp[1] ] }){
	    if( $snp[0] =~ ${ $SNPS[ $snp[1] ] }[$i] ){
		splice( @{ $SNPS[ $snp[1] ] }, $i, 1 );
	    }
	}
    }
    close AV;
    my %blocks;
    printLog("Saving LD block into memory");
    foreach my $i (1..22,"X","Y"){
	printLog("chr$i");
	open (IN, "<null/chr${i}.ld") or die;
	my $header=<IN>;
	while (<IN>){
	    chomp;
	    my @f = split /\t/;
	    next if $f[$#f] < $opt{r};
	    my $chr = "chr".$i;
	    my $tag = "chr".$i.":".$f[1];
	    my $ld = "chr".$i.":".$f[4];
	    if (exists $blocks{$chr}->{$tag}){
		$blocks{$chr}->{$tag} .= ",".$ld;
	    } else {
		$blocks{$chr}->{$tag} = $ld;
	    }
	}
	close IN;
    }
    printLog("Writing MRVs in $prefix");
    foreach my $n (0..($opt{n}-1)){
	my $file = $prefix.$opt{s}."_".sprintf("%04d",$n).".MRVS.txt";
	printLog("Generating $file");
	open (OUT, ">$file") or die;
	foreach my $i (0..$#AV){
	    if ($#{$SNPS[$AV[$i]]} < 1){
		$AV[$i]-=1 until $#{$SNPS[$AV[$i]]} > 1;
	    }
#	    print "i: $i\t$AV[$i]\t\$\#SNPS[AV[i]]:".$#{$SNPS[$AV[$i]]}."\n";
	    my $r = int(rand( $#{ $SNPS[ $AV[$i]] } + 1));
	    my ($LDchr,$LDpos) = ${ $SNPS[ $AV[$i] ] }[$r] =~ m/:/ ? split /:/, ${ $SNPS[ $AV[$i] ] }[$r] : die ${ $SNPS[ $AV[$i] ] }[$r];
	    if (!exists $blocks{$LDchr}->{${ $SNPS[ $AV[$i]] }[$r]}){
		printLog(${$SNPS[$AV[$i]]}[$r]." does not exists in block\n");
		$r = int(rand( $#{ $SNPS[ $AV[$i]] } + 1)) until (exists $blocks{$LDchr}->{${ $SNPS[ $AV[$i]] }[$r]});
#		print "r changed to $r now\n";
	    }
#	    print "r: $r\t".${$SNPS[$AV[$i]]}[$r]."\n";
	    
#	    print ${$SNPS[$AV[$i]]}[$r]."\n";		    
	    print OUT "# raSNP ".$RA[ $i ]."\n";
	    
#	    print STDERR $RA[$i]."\t".${ $SNPS[ $AV[$i] ] }[$r]."\n" if (!exists $blocks{${ $SNPS[ $AV[$i] ] }[$r]});
	    if ($blocks{$LDchr}->{${$SNPS[$AV[$i]]}[$r]} =~ m/,/){
		my @LDsnps = split /,/,$blocks{$LDchr}->{${ $SNPS[ $AV[$i]] }[$r]};
		foreach my $LDsnp (@LDsnps){
		    my ($c,$s) = split /:/, $LDsnp;
		    my $e = $s+1;
		    print OUT "$c\t$s\t$e\n";
		}
	    } else {
		my ($c,$s) = split /:/,$blocks{$LDchr}->{${ $SNPS[ $AV[$i] ] }[$r]};
		my $e =$s+1;
		print OUT "$c\t$s\t$e\n";
	    }

#	    $block->execute(${$SNPS[ $AV[$i] ] }[$r]);
#	    while (my $ldsnps = $block->fetchrow_hashref){
#		print "$ldsnps->{ld}\n";
#		$sth->execute($ldsnps->{ld});
#		while (my $row = $sth->fetchrow_hashref) {
#		    my $e = $row->{pos}+1;
#		    print OUT $row->{chr}."\t".$row->{pos}."\t$e\n";
#		}
#	    }
	}
	close OUT;
	printLog("Done");
    }
#    perl MRVSs.pl dat/ld_tally2.txt AVS/$opt{s}.LDXI.tally.txt $noOfShuffles
#-----------------#
}
if ($opt{p} eq "xml" || $opt{p} eq "all" ){
    my %tally;
#     mkdir $opt{s}/outputs/
#     mkdir $opt{s}/jobscripts/
#     mkdir $opt{s}/sh/
    my $AVSsuffix = $opt{A} ? $opt{A} : $opt{s}; #user can choose to use AVS/MRV files previously generated for new sets of beds	
    if ($opt{d} =~ m/(bed|peak|gz)$/i){
	my $bedfile = $opt{d};
	die "Is the bed file zipped?\n" if $bedfile =~ m/(gz|zip|tar)$/;
	printLog("Reading $bedfile");
        my $tally = LDX_v_bed($AVSsuffix.".AVS/".$AVSsuffix.".LDXI.bed", $bedfile   );
        $tally{$bedfile}->{'AVS'} = $tally;
        opendir (MRVdir, $AVSsuffix.".MRVs") or die "Error in opening dir $AVSsuffix.MRVs\n";
        while (my $mrvfile = readdir(MRVdir)){
            next if ($mrvfile !~ m/txt$/i);
            printLog("Reading $mrvfile");
            my $tally = LDX_v_bed($AVSsuffix.".MRVs/".$mrvfile, $bedfile    );
            $tally{$bedfile}->{$mrvfile} = $tally;
        }
        closedir (MRVdir);
	foreach my $bed ( keys %tally){
	    print $tally{$bed}->{'AVS'};
	    foreach my $mrv (keys %{$tally{$bed}}){
		next if $mrv eq "AVS";
		print "\t".$tally{$bed}->{$mrv};
	    }
	    my @bedfileName_arr = split /\./, $bed;
	    my $bedfileName = $bedfileName_arr[0];
	    print "\t".$bedfileName."\n";
	}
    } else {
	opendir (DIR, $opt{d}) or die "Error in opening dir $opt{d}\n";
	while (my $bedfile= readdir(DIR)){
	    next if ($bedfile !~ m/bed$/i && $bedfile !~ m/peak$/i);
	    printLog("Reading $bedfile");
	    my $tally = LDX_v_bed($AVSsuffix.".AVS/".$AVSsuffix.".LDXI.bed", $opt{d}."/".$bedfile	);
	    $tally{$bedfile}->{'AVS'} = $tally;
	    opendir (MRVdir, $AVSsuffix.".MRVs") or die "Error in opening dir $AVSsuffix.MRVs\n";
	    while (my $mrvfile = readdir(MRVdir)){
		next if ($mrvfile !~ m/txt$/i);
		printLog("Reading $mrvfile");
		my $tally = LDX_v_bed($AVSsuffix.".MRVs/".$mrvfile, $opt{d}."/".$bedfile	);
		$tally{$bedfile}->{$mrvfile} = $tally;
	    }
	    closedir (MRVdir);
	}
	closedir (DIR);
	if (! -d $opt{s}.".output"){
	    mkdir $opt{s}.".output";
	}
	open (OUT, ">", $opt{s}.".output/".$opt{s}.".VSE.txt") or die;
	foreach my $bed ( keys %tally){
	    print OUT $tally{$bed}->{'AVS'};
	    foreach my $mrv (keys %{$tally{$bed}}){
		next if $mrv eq "AVS";
		print OUT "\t".$tally{$bed}->{$mrv};
	    }
	    print STDERR $bed."\t";
	    my @bedfileName_arr = split /\./, $bed;
	    my $bedfileName = $bedfileName_arr[0];
	    print OUT "\t".$bedfileName."\n";
	    print STDERR "\t".$bedfileName."\n";
	}
	close OUT;
    }
}
#     for bed in $opt{d}/*.bed; do
# 	bedid=$(basename $bed ".bed");
# 	echo -n '<VSE>' >$opt{s}/outputs/$opt{s}.$bedid.o
# 	echo -n '<AVS>' >>$opt{s}/outputs/${s}.$bedid.o
# 	tally=`perl LDX_v_BED.pl AVS/$opt{s}.LDXI.bed $bed | bash VSE_score.sh | wc -l`
# 	printf " %2d" $tally >>$opt{s}/outputs/${s}.$bedid.o
# 	echo "AVS done. $tally found."
# 	echo -n '<\AVS>' >>$opt{s}/outputs/${s}.$bedid.o
# 	echo -n '<MRVS>' >>$opt{s}/outputs/${s}.$bedid.o
# 	for i in MRVs/*.txt; do 
# 	    tally=`perl LDX_v_BED.pl $i $bed | bash VSE_score.sh | wc -l`
# 	    printf " %2d" $tally >>$opt{s}/outputs/${s}.$bedid.o
# 	    echo "$i done for $bedid. $tally found."
# 	done
# 	echo -n "<\MRVS>" >>$opt{s}/outputs/${s}.$bedid.o
# 	echo -n "<BED>$bedid<\BED>" >>$opt{s}/outputs/${s}.$bedid.o
# 	echo -n "<\VSE>" >>$opt{s}/outputs/${s}.$bedid.o
#     done 
#     echo "$bed finished"
# fi

if ($opt{p} eq "R" || $opt{p} eq "all"){
    open (OUT, ">", $opt{s}.".output/".$opt{s}.".VSE.stat.txt") or die;
    printLog("Generating boxplot");
#     ls $opt{s}/outputs/${s}*.o | while read o; do
# 	    cat $o | grep '<VSE>.*<\\VSE>'; # changed from grep -p to just grep
#     done  > $opt{s}/outputs/xmloutput.xml

#     cat $opt{s}/outputs/xmloutput.xml | \
# 	perl -nle 's/<\\*(VSE|AVS|BED|MRVS)>/ /g; print;' | \
# 	perl -nle 's/ +/\t/g; print;' > $opt{s}/outputs/${s}.VSE.txt
#     module load R/3.1.1

#     echo -n "Generating boxplot"
     my $Routput = `Rscript VSE.tmp.R $opt{s}.output/$opt{s}.VSE.txt $opt{s}.output/$opt{s}.final_boxplot.pdf | grep \"^\\[1\\]\" | sort -k8rn`;
    $Routput =~ s/\[1\]//g;
    print OUT $Routput;
# 	perl -nle \'s/\/.*\// /g; s/\[1\]//g; print;\' > 
# 	".$opt{s}.".output/".$opt{s}.".VSE.stats.txt");
#     echo "...Summarized in $opt{s}/outputs/$opt{s}.VSE.stats.txt";
#     echo "Done"
# fi
    close OUT;
}
