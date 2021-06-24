#!/usr/bin/perl -w

# Define default and necessary variables
use Cwd;
use Cwd 'abs_path';
my $pair = 'no';
my @group = (1, 1);
my $replication = 'yes';
my $commondisp = 0.05;
my $start = time;
my $BIN = '/bin';
my $TMP = '/tmp';
my $DATA = '/data';
my $ANN = '/ann';
my $strand = 0;
my $countstrand = 'unstranded';
my $p = 1;
my $ResultDir = getcwd;
my $Mode = 'DE';
my $PE = 'no';
$|=1;

# SET DIRECTORY VARIABLES
my $INSTALLDIR = abs_path($0);
$INSTALLDIR = substr $INSTALLDIR, 0, -8;

my $BINDIR = $INSTALLDIR . $BIN;
my $TMPDIR = $INSTALLDIR . $TMP;
my $DATADIR = $INSTALLDIR . $DATA;
my $ANNDIR = $INSTALLDIR . $ANN;

# Read the installed genomes
my $filename = $DATADIR . "/Genomes";
open (FILE, $filename);
my @genomesfull = "TAIR10";
close(FILE);
chomp(@genomesfull);
my %unique = map { $_ => 1 } @genomesfull;
my @genomes = keys %unique;

# Detailed help function for Count
sub printCountCMD {
	print STDERR "\n\tScript is running in Count mode! Use -mode DE to change to differential expression mode\n";
	print STDERR "\n\tUsage: iRNA.pl -mode Count -a [ SAM/BAM files ] -count [ exon, intron, gene or pol ] -p [ # processes ] -s [ 0 = unstranded (Default), 1 = stranded, 2 = reversely stranded ] -n [ name ] -dir [ directory to output results to ]\n";
	print STDERR "\n\t\tObligatory parameters\n";
	print STDERR "\t\t\t-a\t:\t Path to input files in either SAM or BAM format - Must supply atleast 2 files\n";
	print STDERR "\t\t\t-g\t:\t Name of the genome to be used. Supported genomes are @genomes \n";
	print STDERR "\t\t\t-count\t:\t The features to be quantified. Either exon, intron, gene or pol. Using option pol excludes the first 500bp of all genes from counted regions. \n";
	print STDERR "\t\t\t-n\t:\t The name of this analysis. Result will be outputted in the installation directory under this name. \n";
	print STDERR "\n";
	print STDERR "\n\t\tOptional parameters\n";
	print STDERR "\t\t\t-p\t:\t The number of processors to use during tag counting. Default = 1\n";
	print STDERR "\t\t\t-s\t:\t Use for stranded analysis. Options are 0 = unstranded, 1 = stranded, 2 = reversely stranded. Default = 0  \n";
	print STDERR "\t\t\t-dir\t:\t Path to directory to output results to. Default = Current directory  \n";
	print STDERR "\t\t\t-pe\t:\t Indicate if data is paired-end or not. Options are 'yes' or 'no'. Default = no. For paired-end files please sort them by queryname. \n";
	print STDERR "\t\t\t-norm\t:\t Indicate if data should be output as normalized or not. Options are 'yes' or 'no'. Default = no.\n";
	print STDERR "\n";
	exit;
}

# Detailed help function for DE 
sub printDECMD {
	print STDERR "\n\tScript is running in DE mode! Use -mode Count to change to counting mode for data with multiple conditions\n";
	print STDERR "\n\tUsage: iRNA.pl -a [ Condition A SAM/BAM files ] -b [ Condition B SAM/BAM files ] -group [ Grouping levels ] -g [Genome] -count [ exon, intron, gene or pol ] -p [ # processes ] -s [ 0 = unstranded (Default), 1 = stranded, 2 = reversely stranded ] -n [ name ] -dir [ directory to output results to ]\n";
	print STDERR "\n\t\tObligatory parameters\n";
	print STDERR "\t\t\t-a\t:\t Path to input file(s) in either SAM or BAM format for conditions A\n";
	print STDERR "\t\t\t-b\t:\t Path to input file(s) in either SAM or BAM format for conditions B\n";
	print STDERR "\t\t\t-g\t:\t Name of the genome to be used. Supported genomes are @genomes \n";
	print STDERR "\t\t\t-count\t:\t The features to be quantified. Either exon, intron, gene or pol. Using option pol excludes the first 500bp of all genes from counted regions. \n";
	print STDERR "\t\t\t-n\t:\t The name of this analysis. Result will be outputted in the installation directory under this name. \n";
	print STDERR "\n";
	print STDERR "\n\t\tOptional parameters\n";
	print STDERR "\t\t\t-group\t:\t Use to normalize batch effects. Do not provide for unreplicated or unpaired analysis.\n";
	print STDERR "\t\t\t\t\t The grouping factor provides information to downstream scripts for pairing of samples.\n";
        print STDERR "\t\t\t\t\t Example: 2 conditions, 2 sets of paired samples. Run: -a CondA_Set1 CondA_Set2 -b CondB_Set1 CondB_Set2 -group 1 2 1 2\n";
	print STDERR "\t\t\t-p\t:\t The number of processors to use during tag counting. Default = 1\n";
	print STDERR "\t\t\t-s\t:\t Use for stranded analysis. Options are 0 = unstranded, 1 = stranded, 2 = reversely stranded. Default = 0  \n";
	print STDERR "\t\t\t-dir\t:\t Path to directory to output results to. Default = Current directory  \n";
	print STDERR "\t\t\t-pe\t:\t Indicate if data is paired-end or not. Options are 'yes' or 'no'. Default = no. For paired-end files please sort them by queryname.\n";
	print STDERR "\n";
	exit;
}

# Get and check which mode user wants
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-mode') {
		$Mode = $ARGV[++$i];
	}	
}

if ($Mode ne 'Count' && $Mode ne 'DE') {
	print STDERR "\n\n\t Mode must be either 'Count' or 'DE'. Default = DE\n\n";
	exit;
}

# Generate random string
my $string;
my @chars = ("A".."Z", "a".."z");
$string .= $chars[rand @chars] for 1..8;

# Parse shared variables
my $check=0;
my $ArraySize=@ARGV;

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-h' || $ARGV[$i] eq '--help') {
	if ($Mode eq 'DE' ) {
	printDECMD()
	}
	if ($Mode eq 'Count' ) {
	printCountCMD()
	}
	}	
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-g') {
		$g = $ARGV[++$i];
		$check = $check + 1;
		last;
	} elsif ( $i == $ArraySize - 1 ) {
	print STDERR "\n\n\tYou must supply a genome (-g).\n";
        if ($Mode eq 'DE' ) {
        printDECMD()
        }
        if ($Mode eq 'Count' ) {
        printCountCMD()
        }
	}	
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-pe') {
		$PE = $ARGV[++$i];
	}	
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-n') {
		$name = $ARGV[++$i];
		$check = $check + 1;
		last;
	} elsif ( $i == $ArraySize - 1 ) {
        print STDERR "\n\n\tYou must supply a name for the analysis (-n).\n";
        if ($Mode eq 'DE' ) {
        printDECMD()
        }
        if ($Mode eq 'Count' ) {
        printCountCMD()
        }
        }

}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-dir') {
		$ResultDir = $ARGV[++$i];
	}	
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-p') {
		$p = $ARGV[++$i];
	}	
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-s') {
		$strand = $ARGV[++$i];
	}	
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-count') {
		$c = $ARGV[++$i];
		$check = $check + 1;
		last;
	} elsif ( $i == $ArraySize - 1 ) {
        print STDERR "\n\n\tYou must a region to count tags in (-count).\n";
        if ($Mode eq 'DE' ) {
        printDECMD()
        }
        if ($Mode eq 'Count' ) {
        printCountCMD()
        }
        }
	
}


for (my $i=0;$i<@ARGV;$i++) {
        if ($ARGV[$i] eq '-a') {
                $check = $check + 1;
                print STDERR "\n\tInput files:";
                print STDERR "\n\t\tCondition A:";
                while ($ARGV[++$i] !~ /^\-/) {
                        push(@AFiles, $ARGV[$i]);
                        print STDERR "\n\t\t\t$ARGV[$i]";
                }
        }
}


# Set stranded parameter

if ($strand != 0) {
	$countstrand = "strand";
}
	

# Check if featureCounts is installed
my $FC_path = `which featureCounts 2>/dev/null`;
print STDERR "\n\nfeatureCounts is not executable, please see: http://subread.sourceforge.net/\n\n" unless ( $FC_path );
exit unless ( $FC_path );

# Get normalization variable

my $norm = 'no';
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-norm') {
		$norm = $ARGV[++$i];
	}	
}

if ($Mode eq 'Count') {
# Check for correct input of parameters
if ($check != 4) {
	print STDERR "\n\n\tYou must input the necessary parameters.\n";
	printCountCMD()
}

if ($c ne 'intron' && $c ne 'exon' && $c ne 'gene' && $c ne 'pol') {
	print STDERR "\n\n\tYou must choose a valid analysis method (-count). \n";
	printCountCMD()
}

if (! grep (/$g/, @genomes) ) {
	print STDERR "\n\n\tYou must choose a supported genome (-g).\n";
	printCountCMD()
}

my $AFilesSize = @AFiles;
if ($AFilesSize == 1) {
        print STDERR "\n\n\tYou must supply more than 1 file (-a).\n";
        printCountCMD()
}

# Count tags in proper file
print STDERR "\n\n\tSTATUS:\n";
print STDERR "\t\tCounting tags \n";

my @Files = (@AFiles);
`mkdir $TMPDIR 2> /dev/null`;

if ($PE eq 'yes') {
`featureCounts --primary -p -B -C -F GTF -T $p -f -s $strand -a $DATADIR/$g.$c.$countstrand.gtf -o $TMPDIR/$string @Files 2> /dev/null`;
`sed '1d' $TMPDIR/$string > $TMPDIR/$string.tmp`;
`mv $TMPDIR/$string.tmp $TMPDIR/$string`;
`rm $TMPDIR/$string.summary`;
}

if ($PE eq 'no') {
`featureCounts --primary -F GTF -T $p -f -s $strand -a $DATADIR/$g.$c.$countstrand.gtf -o $TMPDIR/$string @Files 2> /dev/null`;
`sed '1d' $TMPDIR/$string > $TMPDIR/$string.tmp`;
`mv $TMPDIR/$string.tmp $TMPDIR/$string`;
`rm $TMPDIR/$string.summary`;
}

# Fetch gene annotation and run R scripts

print STDERR "\t\tStitching regions together \n";

`cp  $ANNDIR/$g.annotation.txt $TMPDIR/$string.ann`;
open (R, "Rscript $BINDIR/Stich.R $INSTALLDIR $c $string $norm @AFiles|");
while ( <R> ) {
        print STDOUT;

}
}

if ($Mode eq 'DE') {
# Check that tools are executable
my $R_path = `which R 2>/dev/null`;
print STDERR "\n\nR is not executable, please see: http://www.r-project.org/\n\n" unless ( $R_path );
exit unless ( $R_path );

my $Package = `Rscript $BINDIR/Check.Package.R`;
if ($Package eq 'FALSE') {
print STDERR "\n\nedgeR is not installed under R, please see: http://http://bioconductor.org/\n\n";
exit;
}

# Parse mode-specific arguments
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-group') {
		$pair = 'yes';
		undef @group;
		while ($i < @ARGV && $ARGV[++$i] !~ /^\-/) {
			push(@group, $ARGV[$i]);
		}
	} 
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-b') {
		$check = $check + 1;
		print STDERR "\n\n\t\tCondition B:";
		while ($ARGV[++$i] !~ /^\-/) {
			push(@BFiles, $ARGV[$i]);
			print STDERR "\n\t\t\t$ARGV[$i]";
		}
	}
}

# Check that user has given all necessary arguments and correct format

if ($check != 5) {
	print STDERR "\n\n\tYou must input the necessary parameters.\n";
	printDECMD()
}

if ($c ne 'intron' && $c ne 'exon' && $c ne 'gene' && $c ne 'pol') {
	print STDERR "\n\n\tYou must choose a valid analysis method (-count). \n";
	printDECMD()
}

if (! grep (/$g/, @genomes) ) {
	print STDERR "\n\n\tYou must choose a supported genome (-g).\n";
	printDECMD()
}


# Check for replicated samples 
my $AFilesSize = @AFiles;
my $BFilesSize = @BFiles;

if ($AFilesSize == 1 && $BFilesSize == 1) {
	print STDERR "\n\n\t##################################################################################### \n";
	print STDERR "\t RUNNING UNREPLICATED ANALYSIS!! SETTING COMMON.DISPERSION TO $commondisp.\n";
	print STDERR "\t PLEASE BE AWARE THE RESULTS ARE LESS THAN OPTIMAL. PLEASE REFER TO THE EDGER MANUAL. \n";
	print STDERR "\t##################################################################################### \n";
	$replication = 'no';
}

# Count tags in proper file

print STDERR "\n\n\tSTATUS:\n";
print STDERR "\t\tCounting tags \n";

my @Files = (@AFiles, @BFiles);
`mkdir $TMPDIR 2> /dev/null`;

if ($PE eq 'yes') {
`featureCounts --primary -p -B -C -F GTF -T $p -f -s $strand -a $DATADIR/$g.$c.$countstrand.gtf -o $TMPDIR/$string @Files 2> /dev/null`;
`sed '1d' $TMPDIR/$string > $TMPDIR/$string.tmp`;
`mv $TMPDIR/$string.tmp $TMPDIR/$string`;
`rm $TMPDIR/$string.summary`;
}

if ($PE eq 'no') {
`featureCounts --primary -F GTF -T $p -f -s $strand -a $DATADIR/$g.$c.$countstrand.gtf -o $TMPDIR/$string @Files 2> /dev/null`;
`sed '1d' $TMPDIR/$string > $TMPDIR/$string.tmp`;
`mv $TMPDIR/$string.tmp $TMPDIR/$string`;
`rm $TMPDIR/$string.summary`;
}

# Fetch gene annotation and run R scripts

print STDERR "\t\tRunning analysis in R \n";


`cp  $ANNDIR/$g.annotation.txt $TMPDIR/$string.ann`;
open (R, "Rscript $BINDIR/Analyze.R $pair $INSTALLDIR $replication $commondisp $c @AFiles SEP @BFiles SEP2 @group SEP3 $string|");
while ( <R> ) {
        print STDOUT;
}
}

# Clean up
`mv $TMPDIR/out.txt $ResultDir/$name.txt`;
`rm $TMPDIR/$string*`;

# Output execution time
my $duration = time - $start;
print "\n\n\tExecution time: $duration s\n\n";


