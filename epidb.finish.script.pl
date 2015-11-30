#!/usr/bin/perl

use warnings;
use strict;

my $dirname = $ARGV[0];
my $tissuename = $ARGV[1];
my $species = $ARGV[2];
my $vcfmergeArgs = "~/vcftools_0.1.12b/bin/vcf-merge";
my $systemArgs = "";

opendir(DIR, $dirname) || die "Error opening dir $dirname\n";
while((my $filename = readdir(DIR)))
{	
	if($filename =~ /.vcf.gz/ && !($filename =~ /.tbi/))
	{
		$vcfmergeArgs .= " ".$filename
	}
}

$vcfmergeArgs .= " > ".$tissuename.".all.vcf";
print $vcfmergeArgs."\n";
system("$vcfmergeArgs");
my $parseArgs = "grep -v \"\#\" ".$tissuename.".all.vcf | awk '{print\$1\":\"\$2}' > ".$tissuename.".snv.list";
system("$parseArgs");

my $joinCount = 0;
my $joinArgs = "";
my $geneArgs = "";
opendir(DIR, $dirname) || die "Error opening dir $dirname\n";
while((my $filename = readdir(DIR)))
{	
	if($filename =~ /.snp.indel.vcf.gz/ && !($filename =~ /.tbi/))
	{
		my $srabase = $filename;
		$srabase =~ s/.snp.indel.vcf.gz//;
		
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -nt 12 -R /home/epidb/index/".$species."/current_build.fa -T UnifiedGenotyper -I ".$srabase.".recalc.bam -o ".$srabase.".all.vcf -L ".$tissuename.".snv.list -dt NONE -out_mode EMIT_ALL_SITES -allowPotentiallyMisencodedQuals";
		system("$systemArgs");
		$systemArgs = "echo \"SNP.name ".$srabase.".ref ".$srabase.".alt\" > ".$srabase.".parse.header";
		system("$systemArgs");
		$systemArgs = "grep -v \"\#\" ".$srabase.".all.vcf | perl ~/parse.vcf.for.ase.pl > ".$srabase.".parse.body";
		system("$systemArgs");
		$systemArgs = "cat ".$srabase.".parse.header ".$srabase.".parse.body > ".$srabase.".parse.txt";
		system("$systemArgs");
		if($joinCount == 0)
   		{
   			$joinArgs = "join ".$srabase.".parse.txt";
   			$geneArgs = "join ".$srabase.".genes.txt";
   			$joinCount++;
   		}
   		elsif($joinCount == 1)
   		{
   			$joinArgs .= " ".$srabase.".parse.txt";
   			$geneArgs .= " ".$srabase.".genes.txt";
   			$joinCount++;
   		}
   		else
   		{
   			$joinArgs .= " | join - ".$srabase.".parse.txt";
   			$geneArgs .= " | join - ".$srabase.".genes.txt";
   		}
	}
}
$joinArgs .= " > ".$tissuename.".parse.txt";
print $joinArgs."\n";
system($joinArgs);
$geneArgs .= ' | grep -v "gene_id" | perl ~/average.gene.expression.pl > '.$tissuename.".genes.txt";
print $geneArgs."\n";
system($geneArgs);

$systemArgs = "Rscript ~/ASE.analysis_auto_V1.R ~/processing/".$tissuename.".txt/ ".$tissuename.".parse.txt ~/processing/".$tissuename.".txt ".$tissuename.".ase.output.txt";
system("$systemArgs");

#$systemArgs = "rm *.idx *.tbi *.vcf.gz *.parse.header *.parse.body *.sra *single.fastq";
#system("$systemArgs");
#$systemArgs = "../bzip2.all.sh";
#system("$systemArgs");
#$systemArgs = "rm *.fastq";
#system("$systemArgs");