#!/usr/bin/perl
## CONVERTS PHASE (CHROMOPAINTER) FORMAT TO (CHROMOPAINTER)V2 
use strict;
use warnings;
use Getopt::Long;

sub help {
print("CONVERTS FROM CHROMOPAINTER v1 FORMAT TO v2\n");

print("usage: perl chromopainter2chromopainterv2.pl <phasefile> <outputphasefile>\n");

print("with:\n");
print("<phasefile>:		ChromoPainter/PHASE style SNP file\n");
print("<outputphasefile>:       Output phase file\n\n");

print("<options>:\n");
print("-p <val> : Ploidy\n");
print("-v: Verbose mode\n");
die "\n";
}

###############################
## ARGUMENT PROCESSING
my $phasefile="";
my $outfile="";
my $verbose=0;
my $ploidy=2;

GetOptions ('v|verbose' => \$verbose,
    'p|phase=s' => \$ploidy);
if(@ARGV != 2) {help();}

$phasefile=$ARGV[0];
$outfile=$ARGV[1];

####################################
## Define global variables
my @snplocs; # location of the SNPS

my $numsnps=0; # number of SNPS defined in the file
my $numinds=0; # number of individuals defined in the file
my $numhaps=0; # number of haplotypes observed

####################################
## File IO

## Check we can read the input files
open PHASEFILE, $phasefile or die $!;

## Create output files
open OUTFILE, ">", $outfile or die $!;

####################################
## Functions we need
sub trim($){  # remove whitespace from beginning and end of the argument
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

####################################
## Read the phasefile 

## read the PHASEFILE header
my $skip=1;
my @tmarr;
my $detectv2=1;
while ($skip) {
	my $tmp=<PHASEFILE>;
	my @tmpvals = split(/ /, $tmp);
	if($tmpvals[0] eq "P"){ # found the line with all the SNP locations
		@snplocs= split(/ /, $tmp);
		shift @snplocs;
		my $floc=tell(PHASEFILE);
		$tmp=<PHASEFILE>; # read the line of S's, if it exists
		if(substr($tmp, 0, 1) eq "S"){
		    $detectv2=0;
		}
		$numsnps=trim(pop @tmarr);
		$numinds=trim(pop @tmarr);
		$numhaps=$numinds*$ploidy;
		$skip=0;
	}else {
		push @tmarr, $tmpvals[0];
	}
}
if($detectv2==0){
    print "Detected Chromopainter v1 format\n";
    print "Detected $numinds individuals\n";
}else{
    print "Detected Chromopainter v2 format\n";
    print "Detected $numinds haplotypes\n";
    die("Require v1 file to convert to v2!\n");
}
print "And $numsnps SNPs\n";

## Print the new header

print OUTFILE "$numhaps\n$numsnps\nP";
for(my $i=0;$i <$numsnps;++$i){
    print OUTFILE " $snplocs[$i]";
}
#print OUTFILE "\n";

# remaining lines are SNPs
while (my $tmp=<PHASEFILE>) {
    if($verbose){print "Reading haplotype $numhaps\n";}
    ++$numhaps;
    my @tarr=split(//,trim($tmp));
    if(scalar(@tarr)!=$numsnps){
	my $tmp=scalar(@tarr);
	die "Expected $numsnps SNPs on haplotype $numhaps, but received $tmp\n";
    }
    for(my $i=0;$i <$numsnps;++$i){
	print OUTFILE $tarr[$i];
    }
    print OUTFILE "\n";
}
close PHASEFILE;
close OUTFILE;
