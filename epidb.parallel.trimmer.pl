#!/usr/bin/perl

use warnings;
use strict;
use lib '/home/epidb/perl/lib';
use Parallel::ForkManager;

my $filename = $ARGV[0];
my $sraname = $ARGV[1];
my $library = $ARGV[2];
my $systemArgs = "";

if($library eq "PAIRED")
{
	my $srabase = $sraname;
   	$srabase =~ s/.sra//;
   	my $fastq1 = $sraname;
   	$fastq1 =~ s/.sra/_1.fastq/;
   	my $fastq1Split = $fastq1."_split";
   	my $fastq2 = $sraname;
   	$fastq2 =~ s/.sra/_2.fastq/;
   	my $fastq2Split = $fastq2."_split";
   	my $trimmed1 = $sraname;
   	$trimmed1 =~ s/.sra/.trimmed_1.fastq/;
   	my $trimmed2 = $sraname;
   	$trimmed2 =~ s/.sra/.trimmed_2.fastq/;
   	
   	my $readFile = $fastq1.".top.read";
   	$systemArgs = "head -4 ".$fastq1." > ".$readFile;
   	system("$systemArgs");
   	
   	my $minReadLength = 0;
   	my $testCount = 0;
   	open (RLENGTH, "<$readFile");
   	while(my $finput = <RLENGTH>)
   	{
   		chomp $finput;
   		if($testCount == 0)
   			{$testCount++;}
   		if($testCount == 1)
   		{
   			my @splitLine = split(//, $finput);
   			$minReadLength = @splitLine;
   			$minReadLength = int($minReadLength*0.8);
   		}
   	}
   	
   	my $fastqLength = $fastq1.".length";
   	$systemArgs = "wc -l ".$fastq1." | awk '{print\$1}' > ".$fastqLength;
   	system("$systemArgs");
   	
	open (FLENGTH, "<$fastqLength");
	while(my $finput = <FLENGTH>)
	{
  		chomp $finput;
  		my $divisor = int($finput/12)+1;
  		while(($divisor%4) != 0)
  			{$divisor++;}
  		$systemArgs = "split -l ".$divisor." ".$fastq1." ".$fastq1Split;
  		system("$systemArgs");
  		$systemArgs = "split -l ".$divisor." ".$fastq2." ".$fastq2Split;
  		system("$systemArgs");
	}
	close FLENGTH;
	
	my $pm = new Parallel::ForkManager(12);
	my $dirname = "/home/epidb/processing/".$filename."/";
	opendir(DIR, $dirname) || die "Error opening dir $dirname\n";
	while((my $fastqFile = readdir(DIR)))
	{	
		if($fastqFile =~ /$fastq1Split/ && !($fastqFile =~ /trimmed/))
		{
			$pm->start and next;
			my $fastqFile2 = $fastqFile;
			$fastqFile2 =~ s/1/2/;
			$systemArgs = "~/sickle-master/sickle pe -t sanger -l ".$minReadLength." -f ".$fastqFile." -r ".$fastqFile2." -o ".$fastqFile."_trimmed -p ".$fastqFile2."_trimmed -s ".$fastqFile."_trimmed_single";
			system($systemArgs);
			$pm->finish;
		}
	}
	$pm->wait_all_children;
	
	my $catArgs1 = "cat";
	my $catArgs2 = "cat";
	opendir(DIR, $dirname) || die "Error opening dir $dirname\n";
	while((my $fastqFile = readdir(DIR)))
	{
		if($fastqFile =~ /$fastq1Split/ && !($fastqFile =~ /trimmed/))
		{
			my $fastqFile2 = $fastqFile;
			$fastqFile2 =~ s/1/2/;
			$catArgs1 .= " ".$fastqFile."_trimmed";
			$catArgs2 .= " ".$fastqFile2."_trimmed";
		}
	}
	$catArgs1 .= " > ".$trimmed1;
	$catArgs2 .= " > ".$trimmed2;
	system("$catArgs1");
	system("$catArgs2");
	$systemArgs = "rm ".$fastq1Split."* ".$fastq2Split."*";
	system("$systemArgs");
}	
elsif($library eq "SINGLE")
{
	my $srabase = $sraname;
   	$srabase =~ s/.sra//;
   	my $fastqname = $sraname;
   	$fastqname =~ s/.sra/.fastq/;
   	my $fastqSplit = $fastqname."_split";
   	my $trimmed = $sraname;
   	$trimmed =~ s/.sra/.trimmed.fastq/;
	my $fastqc = $trimmed;
	$fastqc =~ s/.fastq/_fastqc/;
	my $fastqLength = $fastqname.".length";
   	
   	$systemArgs = "wc -l ".$fastqname." | awk '{print\$1}' > ".$fastqLength;
   	system("$systemArgs");
   	
   	my $readFile = $fastqname.".top.read";
   	$systemArgs = "head -4 ".$fastqname." > ".$readFile;
   	system("$systemArgs");
   	
   	my $minReadLength = 0;
   	my $testCount = 0;
   	open (RLENGTH, "<$readFile");
   	while(my $finput = <RLENGTH>)
   	{
   		chomp $finput;
   		if($testCount == 0)
   			{$testCount++;}
   		if($testCount == 1)
   		{
   			my @splitLine = split(//, $finput);
   			$minReadLength = @splitLine;
   			$minReadLength = int($minReadLength*0.8);
   		}
   	}
   	
	open (FLENGTH, "<$fastqLength");
	while(my $finput = <FLENGTH>)
	{
  		chomp $finput;
  		my $divisor = int($finput/12)+1;
  		while(($divisor%4) != 0)
  			{$divisor++;}
  		$systemArgs = "split -l ".$divisor." ".$fastqname." ".$fastqSplit;
  		system("$systemArgs");
	}
	close FLENGTH;
	
	my $pm = new Parallel::ForkManager(12);
	my $dirname = "/home/epidb/processing/".$filename."/";
	opendir(DIR, $dirname) || die "Error opening dir $dirname\n";
	while((my $fastqFile = readdir(DIR)))
	{	
		if($fastqFile =~ /$fastqSplit/ && !($fastqFile =~ /trimmed/))
		{
			$pm->start and next;
			$systemArgs = "~/sickle-master/sickle se -t sanger -l ".$minReadLength." -f ".$fastqFile." -o ".$fastqFile."_trimmed";
			system($systemArgs);
			$pm->finish;
		}
	}
	$pm->wait_all_children;
	
	my $catArgs = "cat";
	opendir(DIR, $dirname) || die "Error opening dir $dirname\n";
	while((my $fastqFile = readdir(DIR)))
	{
		if($fastqFile =~ /$fastqSplit/ && !($fastqFile =~ /trimmed/))
			{$catArgs .= " ".$fastqFile."_trimmed";}
	}
	$catArgs .= " > ".$trimmed;
	system("$catArgs");
	$systemArgs = "rm ".$fastqSplit."*";
	system("$systemArgs");
}