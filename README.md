# VSE
VSE is a Perl/Rscript command line tool to calculate the enrichment of associated variant set (AVS) for an array of genomic regions.

**Program requirement:**
```
- Perl (5.14 or higher)
- R (3.1.1)

The program further requires R package ```car```. This package will be installed automatically if not installed in the system (internet required).
```

####Installation:
**Step 1: Install Perl**

  You can get the Perl software from [their website](https://www.perl.org/get.html) and install. Version 5.14+ is required.

**Step 2: Install R**

R is also required for calculating the statistics and generating the plots. You can download R at [www.r-project.org](http://www.r-project.org). Version 3.1.1 is preferred, the later versions may work.

The following R package is required: ```car```. The package will be installed automatically from the internet. If connection to internet is not possible, the packages can be downloaded and installed from CRAN.

In most cases, you can install the packages by:

1. Go to R command line by typing ```R``` in your terminal
2. type

  ```
  install.packages("car")
  ```

3. Exit R environment by typing ```q()```


**Step 4: Download VSE**

You must download the file ``vse`` and the supporting 1000 Genome Project LD data from http://helab.uhnresearch.ca/?q=node/26. You can make the file ```vse``` excutable by the command ```chmod +x vse``` and move to your ```bin``` directory to make it executable from any location.


###Using VSE

VSE can be run just by one command:

```
vse --tag example.SNPs/NHGRI-BCa.bed \
    --ld example.SNPs/ld_BCa.bed \
    --beds example.beds/* \
    --dataDir VSE_data/ \
    --output vseOutput
```

#####Options:
```
--tag [required] Location to tag SNPs' list. 
--ld [required] Location to LD SNPs' list.
--output [required] Output file prefix.
--beds [required] List of bed files.
--dataDir [required] Location to the VSE supporting data directory. The data can be dowloaded from http://helab.uhnresearch.ca/?q=node/26
--labels Labels for bed files. Must be in the same order as the bed files. If not provided, the filenames of the bed files will be used.
--r The r value used for determining LD SNPs. Must be between 0.6 and 0.9. Default: 0.8
--bgSize The number of MRVS to compute. Default: 500
--n The number of threads to use. Default: 10
--normality The p-value cutoff to be used for Kolmogorov-Smirnov test. A higher p-value indicates more normality of the null distribution. Anything lower than 0.05 means the null distribution significantly differs from normality, which will render the enrichment result invalid. Default: 0.9. Can be lowered to make the program faster.
--visualize To generate figures from the result.
--p P-value cutoff. Only useful if --visualize switch is on.
--padjust Cutoff for adjust P-value using Bonferroni correction. Only useful if --visualize switch is on. Default: 0.01
--seed Predefined seed value for reproducibility of the result.
--keepTmp Save all the temporary files. Useful for advanced users. Check below.
--quiet Remove verbosity.
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
There should be no header line. The LD SNPs must be calculated based on EUR population of 1000 Genome Project Phase III genotyping data (Release 2013/05) or Phase I genotyping data (Release 2010/11/23). It's recommended to use the phase III genotyping information for calculating LD because the newer releases are more accurant than the older ones because of inclusion of more individuals.

######The genomic ranges
The genomic regions should be in standard bed format.

####Output
VSE produces multiple output files:

```prefix.avs``` The disjoint LD Loci that, when considered together, creates the AVS

```prefix.density.pdf``` shows the density of overlapping tallies from the null

```prefix.VSE.stat.txt``` contains the statistics table

```prefix.boxplot.pdf``` visualizes the null distribution and enrichment of AVS

```prefix.matrix.pdf``` binary representation of overlapping between each locus and annotation. Overlapping is defined as at least one SNP (associated or linked) is within the annotation.

```prefix.matrix.txt``` is a matrix of all overlapping tallies by AVS and MRVS. The first column is the AVS tally and the rest are MRVS.

#####There are several factors to consider:
1. VSE is sensitive to the number of tagSNPs. From our trial and error tests, too low number of tagSNPs (below 15) provide imprecise result.
2. The quality of ChIP-seq data is very important. We recommend users to confirm the quality of the ChIP-seq data and to only use data that are of good quality to avoid false enrichment. There are tools like ChIPQC or Chillin for quality control of ChIP-seq data.
3. Make sure that you use the same r2 cutoff that you used to determine your LD SNPs. Also, the LD SNPs must be calculated using 1000 Genome Project Phase III (May 2013) release.



#####Advanced options:

VSE computationally comprises of three modules - ```preprocess```, ```vse``` and ```stat```. These modules can be run independantly if needed. For typical uses of VSE, these advanced options are not necessary; but for cases where a large number of MRVSs are needed, it is time efficient to run ```preprocess``` ones and use these MRVSs for different batches of bed files.

```preprocess``` generates the AVS and MRVS files. Example run:

```
perl vse --program preprocess --tag example.SNPs//NHGRI-BCa.bed --ld example.SNPs/ld_BCa.bed --output vseOutput
```

The following additional arguments can be passed: ```--r```, ```--bgSize```, ```--n``` and ```--phase```.

The second module, ```vse```, takes the outputs from the first module and computes the overlapping of AVS and MRVSs with the provided genomic features.

```
perl vse vse --avs vseOutput.avs \
     	     --tally vseOutput.avs.tally \
	     --bedDir example.beds/ \
	     --mrvsDir vseOutput.MRVS/ \
	     --output vseOutput
```
Arguments:
```
--avs The AVS file outputted from the first module.
--tally The tally file outputted from the first module.
--mrvsDir Location to the directory where MRVS files are as outputted from the first module.
--bedDir Location to directory that contains all genomic features. The files in the directory must be bed or narrowPeak files.
--output The prefix for the output file. Can be same as the prefix for the first module.
```

The last module, ```stat```, takes the outputs from the second module and calculates the enrichment. It also generates visualization of the analysis.

```
perl vse stat --vse vseOutput.VSE.txt --matrix vseOutput.matrix.txt --output vseOutput
```

#####Options:
```
--vse The VSE output file from the second module.
--matrix The matrix output file from the second module.
--output The prefix for the output file.
```

