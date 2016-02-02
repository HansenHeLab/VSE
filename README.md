# VSE
VSE is a Perl/Rscript command line tool to calculate the enrichment of associated variant set (AVS) for an array of genomic regions.

**Program requirement:**
```
- Perl (5.18 or higher)
- R (3.1.1)

The program further requires R packages ```ggplot2```,```reshape``` and ```car```. These packages are installed automatically if not installed in the system. The system must have access to the internet.
```

####Installation:
**Step 1: Install Perl**

  You can get the Perl software from [their website](https://www.perl.org/get.html) and install. Version 5.x is required.

**Step 2: Install R**

R is also required for calculating the statistics and generating the plots. You can download R at [www.r-project.org](http://www.r-project.org). Version 3.1.1 is preferred, the later versions may work.

The following R packages are required: ```ggplot2```,```reshape```,```car```. These packages are installed automatically from the internet. If connection to internet is not possible, the packages can be downloaded and installed from CRAN.

In most cases, you can install the packages by:

1. Go to R command line by typing ```R``` in your terminal
2. type

  ```
  install.packages("ggplot2")  
  install.packages("reshape")  
  install.packages("car")
  ```

3. Exit R environment by typing ```q()```

```Rscript``` should be an executable command from any location.

**Step 4: Download VSE**

You must download the file ``vse`` and the directories ```lib``` and ```data```. The directory structure must be intact; i.e., ```lib``` and ```data``` directories must reside in the same directory as ```VSE.pl```.


###Using VSE

VSE has three modules that can be run independantly. The details of each module is provided below.
#####Module "preprocess":
The first module, ```preprocess```, computes the AVS and generates the matching random variant sets (MRVSs) from either phase 1 or phase 3 genotyping data from 1000 Genome Project.

```
perl vse preprocess --tag example.SNPs/NHGRI-BCa.bed \
            --ld example.SNPs/ld_BCa.bed \
            -out run1
```

#####Options:
```
--tag Location to tag SNPs' list. 
--ld Location to LD SNPs' list.
--out Output file prefix.
--r The r value used for determining LD SNPs. Must be between 0.6 and 1. Default: 0.8
--phase The 1000 Genome Project data used to calculate LD blocks. Must be 1 or 3. Default: 1
--bgSize The number of MRVS to compute. The default of 100 is sufficient for normal uses. Higher number is computionally heavier and is not recommended for typical uses.
--n The number of threads to use. Default: 4

#####Module "vse":
The second module, ```vse```, takes the outputs from the first module and computes the overlapping of AVS and MRVSs with the provided genomic features.

```
perl vse vse --avs run1.avs \
     	     --tally run1.avs.tally \
	     --bedDir beds/ \
	     --mrvsDir vse.MRVS/ \
	     --out run1
```
#####Options:
```
--avs The AVS file outputted from the first module.
--tally The tally file outputted from the first module.
--mrvsDir Location to the directory where MRVS files are as outputted from the first module.
--bedDir Location to directory that contains all genomic features. The files in the directory must be bed or narrowPeak files.
--out The prefix for the output file. Can be same as the prefix for the first module.
```

#####Module "stat":
The last module, ```stat```, takes the outputs from the second module and calculates the enrichment. It also generates visualization of the analysis.

```
perl vse stat --vse run1.VSE.txt --matrix run1.matrix.txt --out run1-result
```

#####Options:
```
--vse The VSE output file from the second module.
--matrix The matrix output file from the second module.
--out The prefix for the output file.
```
####Input File Formats
VSE requires three input files.
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
There should be no header line. The LD SNPs must be calculated based on EUR population of 1000 Genome Project Phase III genotyping data (Release 2013/05) or Phase I genotyping data (Release 2010/11/23). It's recommended to use the phase III genotyping information for calculating LD because the newer releases are more accurant than the older ones because of inclusion of more individuals.

######The genomic ranges
The genomic regions should be a single directory and the directory path should be given in ```--bedDir``` parameter in the ```vse``` module. The genomics regions should be bed files. The name of the file is used as labels for the regions. For examples, ```DHS.bed``` will be denoted as DHS in final figures and table.
```--bedDir``` can also be a path to a bed file in ```vse``` module. In that case, the tallies will be printed on screen.

####Output
VSE produces multiple output files:

```prefix.density.pdf``` shows the density of overlapping tallies from the null

```prefix.VSE.stat.txt``` contains the statistics table

```prefix.boxplot.pdf``` visualizes the null distribution and enrichment of AVS

```prefix.matrix.pdf``` binary representation of overlapping between each locus and annotation. Overlapping is defined as at least one SNP (associated or linked) is within the annotation.

```prefix.VSE.txt``` is a matrix of all overlapping tallies by AVS and MRVS. The first column is the AVS tally and the rest are MRVS.

#####There are several factors to consider:
1. VSE is sensitive to the number of tagSNPs. From our trial and error tests, too low number of tagSNPs (below 15) provide imprecise result.
2. The quality of ChIP-seq data is very important. We recommend users to confirm the quality of the ChIP-seq data and to only use data that are of good quality to avoid false enrichment. There are tools like ChIPQC or Chillin for quality control of ChIP-seq data.
3. Make sure that you use the same r2 cutoff that you used to determine your LD SNPs. Also, the LD SNPs must be calculated using 1000 Genome Project Phase III (May 2013) release.
