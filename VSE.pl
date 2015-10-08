#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Carp;
use Cwd 'abs_path';
use File::Slurp;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/lib";
use LDX_to_LDXI;
use tally;
use LDX_v_bed;
my %opt;
getopts('r:hvYf:l:s:d:p:n:A:', \%opt);

#----------subroutines---------------
sub usage
{
    print "usage: vse.sh -f snpListBed -s suffix -d dirLocation [-r r2Value] [-v y/n] [-Y y/n] | [-h]\n";
    print "Options:\n";
    print "-r[0.6/0.7/0.8/0.9/1]    R2 value to find SNPs in LD. default: 0.8\n";
    print "-v[y/n]      Verbose; default: y\n";
    print "-Y[y/n]      To analyze chrY or not; default: n\n";
    print "-f[path]     Location of tagSNP list. The file must be a bed file\n";
    print "-l[path]     Location of LD snp list. Must be in this format: chr\ts\te\tldSNP\ttagSNP\n";
    print "-s[char]     Suffix for filenames\n";
    print "-d[path]     Path to feature directory\n";
    print "-A[char]     Suffix for AVS/MRV files. Only used when -p is xml in order to avoid having to generate AVS/MRV again.\n";
    print "-p           [all/AVS/MRV/xml/R] Modular run; default: all\n";
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

sub removeExtension
{
    my $name = shift;
    my @bedfileName_arr = split /\./, $name;
    my $bedfileName = $bedfileName_arr[0];
    return $bedfileName;
}

sub loadTagSNPs
{
    my $AVSsuffix = shift;
    my $AVSLDX = $AVSsuffix.".AVS/".$AVSsuffix.".LDX.bed";
    my @tagSNPs;
    open (HIN, "grep \"#\" $AVSLDX | cut -d \" \" -f 3 |") or die "Could not open $AVSLDX file: $!\n";
    while (<HIN>){
        chomp;
        push (@tagSNPs, $_);
    }
    close HIN;
    return \@tagSNPs;
}

sub loadFileIntoArray
{
    my $file = shift;
    my @lines;
    open (INF, "<$file") or die "Could not open $file: $!\n";
    while (<INF>){
	chomp;
	push @lines, $_;
    }
    close INF;
    return \@lines;
}
#-------------------------

#----PARAMETERS-----------
die usage if ($opt{h});
$opt{n} = 100 if !$opt{n};
$opt{p} = "all" if !$opt{p};
$opt{i} = 100 if !$opt{i};
$opt{n} = 100 if !$opt{n};
$opt{r} = 0.8 if !$opt{r};
my $AVSsuffix = $opt{A} ? $opt{A} : $opt{s}; #user can choose to use AVS/MRV files previously generated for new sets of beds
my $scriptPath = abs_path($0);
my $scriptDir = dirname($scriptPath);
printLog("Script location: $scriptDir");
#------------------------

#-------INPUT QC---------
die "-r must be between 0.6 and 0.8\n" if $opt{r} <0.6 || $opt{r} > 1;
die "parameter missing\n" if (!$opt{f} || !$opt{l} || !$opt{s} || !$opt{d});
die "$opt{f} file not found\n" if (! -e $opt{f});
die "$opt{l} file not found\n" if (! -e $opt{l});

if ($opt{d} =~ m/(bed|peak|gz)$/i){
    printLog("single bed file provided. Output will be printed on screen");
    die "$opt{d} could not be located\n" if (! -e $opt{d});
} else {
    if (! -d $opt{d}){
	print "Directory path not found\n";
    }
}
my $totalTagSnp=`wc -l $opt{f} | cut -d " " -f 1`;
print "$opt{f} does not have any snp\n" if ($totalTagSnp == 0);
die "Unknown -p\n" if $opt{p} !~ m/(all|AVS|MRV|xml|R)/;
my $totalColumnInLD = `awk '{print NF}' $opt{l} | sort -u | wc -l`;
die "Column number is not 6 for all lines in $opt{l}\n" if $totalColumnInLD != 1;

if ($opt{A} && $opt{p} ne "R"){ #if predefined AVS/MRVs set is provided
    my $LDXfile = $AVSsuffix.".AVS/".$AVSsuffix.".LDX.bed";
    die "$LDXfile not found. Please check if -A parameter is correct. If problem persists, run the -p AVS again. Exiting.\n" if (! -e $LDXfile);
    my $LDXIfile = $AVSsuffix.".AVS/".$AVSsuffix.".LDXI.bed";
    die "$LDXIfile not found. Please check if -A parameter is correct. If problem persists, run the -p AVS again. Exiting.\n" if (! -e $LDXIfile);
    my $LDXItallyfile = $AVSsuffix.".AVS/".$AVSsuffix.".LDXI.tally.txt";
    die "$LDXItallyfile not found. Please check if -A parameter is correct. If problem persists, run the -p AVS again. Exiting.\n" if (! -e $LDXItallyfile);
    if ($opt{p} eq "xml"){
	my $MRVSdir = $AVSsuffix.".MRVs/";
	die "$MRVSdir not found. Please check if -A parameter is correct. If problem persists, run the -p MRV again. Exiting.\n" if (! -d $MRVSdir);
	opendir my $dh, $MRVSdir or die "Can't opendir '$MRVSdir': $!";
	my $file_count = scalar grep { -f "$MRVSdir/$_" } readdir $dh;
	closedir $dh or die "Can't closedir: $!";
	die "Number of MRVS found $file_count, but expected 100\n" if ($file_count !=100);
    }
}
#------------------------

#--------AVS-------------
if ($opt{p} eq "AVS" || $opt{p} eq "all"){
    printLog("Preparing AVS file");
    mkdir $opt{s}.".AVS" if (! -d $opt{s}.".AVS" );
    my $AVS_out_file = $opt{s}.".AVS/".$opt{s}.".LDX.bed";
    if (-e $AVS_out_file){ unlink $AVS_out_file; }
    my $lines = loadFileIntoArray($opt{f});
    die "$opt{f} looks empty\n" if scalar @$lines == 0;
    my $lineCounter = 1;
    open (OUT, ">$AVS_out_file") or die "Count not open $AVS_out_file for writing: $!\n";
    foreach my $line (@$lines){
	my @f = split /\t/, $line;
	my $totalColumns = $#f + 1;
	die "Line $lineCounter has $totalColumns columns instead of 4 in $opt{f}\n" if $totalColumns < 3;
	my $header= "# raSNP $f[3]\n";
	print OUT $header;
	open (IN, "awk -v tag=$f[3] '\$5==tag' $opt{l} |") or die;
	while (<IN>){
	    chomp;
	    my @g = split /\t/;
	    if ($g[4] ne $g[3]){
		my $printline = "$g[0]\t$g[1]\t$g[2]\n";
		print OUT $printline;
	    }
	}
	close IN;
	$lineCounter++;
    }
    close OUT;
    printLog("$AVS_out_file created");
#Fixing tagSNP to LDXI and making tally
    printLog("Generating LDXI for AVS/$opt{s}.LDX.bed");
    LDX_to_LDXI($opt{s}.".AVS/$opt{s}.LDX.bed",$opt{s}.".AVS/$opt{s}.LDXI.bed");
    printLog($opt{s}.".AVS/$opt{s}.LDXI.bed created.");
    printLog("Tallying ".$opt{s}.".AVS/$opt{s}.LDXI.bed");
    tally($opt{s}.".AVS/$opt{s}.LDXI.bed",$opt{s}.".AVS/$opt{s}.LDXI.tally.txt");
    printLog($opt{s}.".AVS/$opt{s}.LDXI.tally.txt created.");
}
#-------------------------------

#------------MRV----------------
if ($opt{p} eq "MRV" || $opt{p} eq "all"){
    my $AVSsuffix = $opt{A} ? $opt{A} : $opt{s}; #user can choose to use AVS/MRV files previously generated for new sets of beds
    if (! -d $AVSsuffix.".MRVs"){
	mkdir $AVSsuffix.".MRVs";
    }
    my $prefix = $AVSsuffix.".MRVs/";	
    my @SNPS;
    my $null_tally_file = $scriptDir."/data/ld_tally_r".$opt{r}.".txt.gz";
    die "-r seems wrong!\n" if (!-e $null_tally_file);
    open (SNPS, "gunzip -c $null_tally_file |") or die "$!\n";
    while(<SNPS>){ 
	chomp;
	my @snp = split /\t/;
	next if /SNP/;
	push @{ $SNPS[$snp[1]] }, $snp[0];
    }
    close SNPS;
    my @AV;
    my @RA;
    open (AV, "<", $AVSsuffix.".AVS/".$AVSsuffix.".LDXI.tally.txt") or die "$!\n";
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
	open (IN, "gunzip -c $scriptDir/data/chr${i}.ld.gz |") or die;
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
	my $file = $prefix.$AVSsuffix."_".sprintf("%04d",$n).".MRVS.txt";
	printLog("Generating $file");
	open (OUT, ">$file") or die;
	foreach my $i (0..$#AV){
	    if ($#{$SNPS[$AV[$i]]} < 1){
		$AV[$i]-=1 until $#{$SNPS[$AV[$i]]} > 1;
	    }
	    my $r = int(rand( $#{ $SNPS[ $AV[$i]] } + 1));
	    my ($LDchr,$LDpos) = ${ $SNPS[ $AV[$i] ] }[$r] =~ m/:/ ? split /:/, ${ $SNPS[ $AV[$i] ] }[$r] : die ${ $SNPS[ $AV[$i] ] }[$r];
	    if (!exists $blocks{$LDchr}->{${ $SNPS[ $AV[$i]] }[$r]}){
		printLog(${$SNPS[$AV[$i]]}[$r]." does not exists in block\n");
		$r = int(rand( $#{ $SNPS[ $AV[$i]] } + 1)) until (exists $blocks{$LDchr}->{${ $SNPS[ $AV[$i]] }[$r]});
	    }
	    print OUT "# raSNP ".$RA[ $i ]."\n";
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
	}
	close OUT;
	printLog("Done");
    }
}
#-----------------------------

#-----------xml---------------
if ($opt{p} eq "xml" || $opt{p} eq "all" ){
    my %tally;
    my %heatmap;
    my $AVSsuffix = $opt{A} ? $opt{A} : $opt{s}; #user can choose to use AVS/MRV files previously generated for new sets of beds
    my $LDXfile = $AVSsuffix.".AVS/".$AVSsuffix.".LDX.bed";
    my $LDXIfile = $AVSsuffix.".AVS/".$AVSsuffix.".LDXI.bed";
    my $tagSNPs = loadTagSNPs($AVSsuffix);
    if ($opt{d} =~ m/(bed|peak|gz)$/i){ #single annotation provided
	my $bedfile = $opt{d};
	die "Is the bed file zipped?\n" if $bedfile =~ m/(gz|zip|tar)$/;
	printLog("Reading $bedfile");
        my $tally = LDX_v_bed($LDXIfile, $bedfile); #count tally for enrichment and density plot
        $tally{$bedfile}->{'AVS'} = $tally;
	my $LDX_overlap = heatmap_intersect($LDXfile, $bedfile); #check overlapping status for each risk locus
	foreach my $tagSNP (@$tagSNPs){
	    $heatmap{$bedfile}->{$tagSNP} = exists $LDX_overlap->{$tagSNP} ? 1 : 0;
	}
        opendir (MRVdir, $AVSsuffix.".MRVs") or die "Error in opening dir $AVSsuffix.MRVs\n";
        while (my $mrvfile = readdir(MRVdir)){
            next if ($mrvfile !~ m/txt$/i);
  #          printLog("Reading $mrvfile");
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
	    my $bedfileName = removeExtension($bedfile);
	    print "\t".$bedfileName."\n";
	}
    } else {
	opendir (DIR, $opt{d}) or die "Error in opening dir $opt{d}\n";
	while (my $bedfile= readdir(DIR)){
	    next if ($bedfile !~ m/bed$/i && $bedfile !~ m/peak$/i);
	    printLog("Reading $opt{d}/$bedfile");
	    my $tally = LDX_v_bed($LDXIfile, $opt{d}."/".$bedfile);
	    $tally{$bedfile}->{'AVS'} = $tally;
	    my $LDX_overlap = heatmap_intersect($LDXfile, $opt{d}."/".$bedfile); #check overlapping status for each risk locus
	    foreach my $tagSNP (@$tagSNPs){
		$heatmap{$bedfile}->{$tagSNP} = exists $LDX_overlap->{$tagSNP} ? 1 : 0;
	    }
	    opendir (MRVdir, $AVSsuffix.".MRVs") or die "Error in opening dir $AVSsuffix.MRVs\n";
	    while (my $mrvfile = readdir(MRVdir)){
		next if ($mrvfile !~ m/txt$/i);
#		printLog("Reading $mrvfile");
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
	    my $bedfileName = removeExtension($bed);
	    print OUT "\t".$bedfileName."\n";
	}
	close OUT;
    }
    open (OUT, ">", $opt{s}.".output/".$opt{s}.".matrix.txt") or die;
    print OUT "Annotation";
    for my $tagSNP (@$tagSNPs){
	print OUT "\t$tagSNP";
    }
    print OUT "\n";
    foreach my $bedName (keys %heatmap){
	my $bedfileName = removeExtension($bedName);
	print OUT $bedfileName;
	for my $tagSNP (@$tagSNPs){
	    print OUT "\t".$heatmap{$bedName}->{$tagSNP};
	}
	print OUT "\n";
    }
    close OUT;
}
#------------------------------

#--------------R---------------
if ($opt{p} eq "R" || $opt{p} eq "all"){
    my $AVSsuffix = $opt{A} ? $opt{A} : $opt{s}; #user can choose to use AVS/MRV files previously generated for new sets of beds
    open (OUT, ">", $AVSsuffix.".output/".$AVSsuffix.".VSE.stat.txt") or die;
    printLog("Generating boxplot");
     my $Routput = `Rscript $scriptDir/lib/stat.r $AVSsuffix | grep \"^\\[1\\]\" | sort -k8rn`;
    $Routput =~ s/\[1\]//g;
    print OUT $Routput;
    close OUT;
}
#-----------------------------
