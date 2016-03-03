# VSE
VSE is a Perl command line tool to calculate the enrichment of associated variant set (AVS) for an array of genomic regions. The AVS is the collection of disjoint LD blocks computed from a list of disease associated SNPs and their linked (LD) SNPs. VSE generates a null distribution of matched random variant sets (MRVSs) from 1000 Genome Project Phase III data that are identical to AVS, LD block by block. It then computes the enrichment of AVS intersecting with user provided genomic features (e.g., histone marks or transcription factor binding sites) compared with the null distribution.

**Program requirement:**
```
- R (3.1.1) with package car
- Perl (5.14 or higher) with package Statistics::R
```

####Installation:

**Step 1: Install R**

R is required for calculating the statistics and generating the plots. You can download R at [www.r-project.org](http://www.r-project.org). Version 3.1.1 is preferred, the later versions are not tested but should work. The R package ```car``` is required. The package will be installed automatically if not installed. If it is not automatically installed or you get an error, see at the end of this document for more help.

**Step 2: Install Perl**

You can get the Perl software from [their website](https://www.perl.org/get.html) and install. Version 5.14+ is required. Perl module Statistics::R must be installed. In typical linux environment, you can install by typing ```sudo cpanm Statistics::R``` in your command prompt. Please see at the end of this document for more help.

**Step 3: Download VSE**

You must download the file ``vse`` and the supporting 1000 Genome Project LD data from http://helab.uhnresearch.ca/?q=node/26. Please untar the downloaded file. Optionally, you can make the file ```vse``` excutable by the command ```chmod +x vse``` and add it to your ```$PATH``` to make it executable from any location.


###Using VSE

VSE can be run minimally by the following command:

```
perl vse --tag example.SNPs/NHGRI-BCa.bed --ld example.SNPs/ld_BCa.bed --beds example.beds/* --dataDir VSE_data/ --output vseOutput --visualize
```

#####Required parameters:
```
--tag	  Location to tag SNPs' list. 
--ld  	  Location to LD SNPs' list.
--output  Output file prefix.
--beds    Path to bed files.
--dataDir Location to the VSE supporting data directory. 
	  The data can be dowloaded from 
	  http://helab.uhnresearch.ca/?q=node/26
```

#####Optional parameters:
```
--labels    Labels for bed files. Must be in the same order as the bed files. 
	    If not provided, the filenames of the bed files will be used.
--r 	    The r value used for determining LD SNPs. Must be between 0.6 and 
	    0.9. Default: 0.8
--bgSize    The number of MRVS to compute. Default: 500
--n 	    The number of threads to use. Default: 10
--normality The p-value cutoff to be used for Kolmogorov-Smirnov test. A higher 
	    p-value represents more normality of the null distribution. Anything 
	    lower than 0.05 means the null distribution significantly differs 
	    from normality, which will render the enrichment result invalid. 
	    Default: 0.9. Can be lowered to make the program faster.
--visualize To generate figures from the result.
--p 	    P-value cutoff. Only useful if --visualize switch is on.
--padjust   Cutoff for adjust P-value using Bonferroni correction. Only useful 
	    if --visualize switch is on. Default: 0.01
--seed 	    Predefined seed value for reproducibility of the result.
--keepTmp   Save all the temporary files. Useful for advanced users.
--quiet     Remove verbosity.
```

####Input File Formats
VSE requires four input files.
######The tagSNPs:
It should be a standard tab delimited bed file and provided with ```--tag``` parameter. There should be no header line.
Example:
```
chr1  1000  1001  rs00001
```
######The LD SNPs
The list of LD SNPs should be a **tab delimited** bed file in the following format and should be provided with ```--ld``` parameter. 
```
chr1  900 901 LDSNPid tagSNPid other_optional_columns
```
There should be no header line.

######The genomic ranges
The genomic regions should be in standard bed format.

####Output
VSE produces multiple output files:

```prefix.avs``` The disjoint LD Loci that, when considered together, creates the AVS

```prefix.stat.txt``` Contains the final enrichment analysis output table

```prefix.log``` The log file. Can be useful to find out the seed of the run for reproducibility.

```prefix.boxplot.pdf``` visualizes the null distribution and enrichment of AVS

```prefix.matrix.pdf``` binary representation of overlapping between each locus and annotation. Overlapping is defined as at least one SNP (associated or linked) overlaps the corresponding genomic feature. The dendrogram are drawn using hierarchical structuring.

```prefix.matrix.txt``` The text matrix file used to produce prefix.matrix.pdf. A value of 1 represents at least one SNP (associated or linked) overlaps the corresponding genomic feature.

#####There are several factors to consider:
1. VSE is sensitive to the number of tagSNPs. From our trial and error tests, too low number of tagSNPs (below 15) provide imprecise result.
2. The quality of ChIP-seq data is very important. We recommend users to confirm the quality of the ChIP-seq data and to only use data that are of good quality to avoid false enrichment. There are tools like ChIPQC or Chillin for quality control of ChIP-seq data.
3. Make sure that you use the same r2 cutoff that you used to determine your LD SNPs.
4. Make sure that the input files are from the same genomic build (e.g., both SNPs and genomic features are in Hg19).
4. Null size is a critical factor for reproducible results. Higher the null size better the normalcy of the distribution. Default of 500 is sufficiently large, but you may try 1000 in needed.

####Further help for installing Perl module and R package:
If the R package ```car``` can not be installed automatically be VSE, you can install it manually from CRAN.

In most cases, you can install the packages by:

1. Go to R command line by typing ```R``` in your terminal
2. type ```install.packages("car")```
3. Exit R environment by typing ```q()```

The Perl module ```Statistics::R``` can be installed by typing ```sudo cpanm Statistics::R``` in the command line prompt. If ```cpanm``` command is not found, you can install ```cpanm``` by typing ```sudo cpan App::cpanminus``` in the command line prompt, and then type ```sudo cpanm Statistics::R```. You must have root access.

If you have issues, you can refer to [CPAN help page](http://www.cpan.org/modules/INSTALL.html) for more details.
