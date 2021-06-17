# iRNA-seq Analysis

Description: 

iRNA-seq utilizes intron coverage to determine genome-wide transcriptional activity. This pipeline has been adapted for use with genomes not available through the UCSC Genome Browser. For genomes included in the UCSC Genome Browser see (Madsen et al, 2015). 

Madsen, J. G., Schmidt, S. F., Larsen, B. D., Loft, A., Nielsen, R., & Mandrup, S. (2015). iRNA-seq: computational method for genome-wide assessment of acute transcriptional regulation from total RNA-seq data. Nucleic acids research, 43(6), e40. https://doi.org/10.1093/nar/gku1365

Installation:

To install the iRNA-seq pipeline you must download all available files.

The following softwares need to be installed and executable:

	R software: Please see http://www.r-project.org/ for instructions on installation.
	edgeR package: Please see http://www.bioconductor.org/ for instructions on edgeR installation within R.
	featureCounts: Part of teh Subread package, please see http://subread.sourceforge.net/ for installation.
	bedToGenePred and genePredToGtf: Only needed if adding a new genome, included in iRNA folder or see http://hgdownload.cse.ucsc.edu/admin/exe/ for installation. 

This pipeline utilizes Perl. 

Usage:

For the best results, iRNA-seq analysis should only be performed with total RNA-seq samples that have undergone rRNA depletion. These samples will have a higher quantity of intronic reads than samples that have undergone poly(A) selection for mRNA. 

Work within the iRNA directory. Download the SAM or BAM files to be analyzed to the iRNA directory (SAM or BAM files should be aligned to the genome you will be using). Files for alignments are available in the "Genomes" folder.  

$ perl iRNA_###.pl -a [ Condition A SAM/BAM files ] -b [ Condition B SAM/BAM files ] -group [ Grouping levels ] -g [Genome] -count [ exon, intron, gene or pol ] -p [ # processes ] -s [ 0 = unstranded (Default), 1 = stranded, 2 = reversely stranded ] -n [ name ] -dir [ directory to output results to ]

	Necessary parameters:
		iRNA_###.pl	: iRNA script for the genome used [iRNA_AGPv4_ctg.pl, iRNA_AGPv4.pl, iRNA_AGPv3.pl, iRNA_TAIR10.pl].
		-a 		: Path to input file(s) in either SAM or BAM format for conditions A
		-b 		: Path to input file(s) in either SAM or BAM format for conditions B
		-g 		: Name of the genome to be used. Supported genomes are
		-count 		: The features to be quantified. Either exon, intron, gene or pol. Using option pol excludes the first 500bp of all genes from counted regions.
		-n 		: The name of this analysis. Result will be outputted in the installation directory under this name.

	Optional parameters:
		-group 		: Use to normalize batch effects. Do not provide for unreplicated or unpaired analysis. The grouping factor provides information to downstream scripts for pairing of samples.
                         	  Example: 2 conditions, 2 sets of paired samples. Run: -a CondA_Set1 CondA_Set2 -b CondB_Set1 CondB_Set2 -group 1 2 1 2
		-p 		: The number of processors to use during tag counting. Default = 1
		-s 		: Use for stranded analysis. Options are 0 = unstranded, 1 = stranded, 2 = reversely stranded. Default = 0
		-dir 		: Path to directory to output results to. Default = Current directory
		-pe 		: Indicate if data is paired-end or not. Options are 'yes' or 'no'. Default = no. For paired-end files please sort them by queryname.

To add new genomes not available through the UCSC Genome Browser, use the "AddNewGenome.sh" script. The "AddNewGenome.sh" script must be edited at Lines 6, 13, 18, 431, 491, 511, 512, and 513. You will need to build the "Gene.Dump" and "mRNA.Dump" files for your genome and add them to the "tmp" folder before use. The "Analyze.R" script must be modified to accommodate genome gene ID lengths. Example scripts and tables are in the "tbls" folder. See our paper for further details on adding genomes (---CITE---). For genomes available through the UCSC Genome Browser, please see (Madsen et al, 2015) to add genomes.

Citiation:

Please cite ----our paper----
