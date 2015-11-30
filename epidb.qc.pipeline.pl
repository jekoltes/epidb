#!/usr/bin/perl

use warnings;
use strict;

my $filename = $ARGV[0];
my $species = $ARGV[1];
my $library = $ARGV[2];
my $sraname = $ARGV[3];
my $qcCheck = 0;
my $avgCheck = 0;
my $duplicateP1 = 0;
my $duplicateP2 = 0;
my $duplicateP = 0;
my $systemArgs = "";
   	
if($library eq "PAIRED")
{
	my $srabase = $sraname;
   	$srabase =~ s/.sra//;
   	my $fastq1 = $sraname;
   	$fastq1 =~ s/.sra/_1.fastq/;
   	my $fastq2 = $sraname;
   	$fastq2 =~ s/.sra/_2.fastq/;
   	my $trimmed1 = $sraname;
   	$trimmed1 =~ s/.sra/.trimmed_1.fastq/;
   	my $trimmed2 = $sraname;
   	$trimmed2 =~ s/.sra/.trimmed_2.fastq/;
   	my $trimmedSingle = $sraname;
   	$trimmedSingle =~ s/.sra/.trimmed_single.fastq/;
	my $fastqc1 = $trimmed1;
	$fastqc1 =~ s/.fastq/_fastqc/;
	my $fastqc2 = $trimmed2;
	$fastqc2 =~ s/.fastq/_fastqc/;
	my $fileCheck1 = "/home/epidb/processing/".$filename."/".$fastqc1."/summary.txt";
	my $fileCheck2 = "/home/epidb/processing/".$filename."/".$fastqc2."/summary.txt";
	open (FCHECK1, "<$fileCheck1");
	while(my $finput = <FCHECK1>)
	{
  		chomp $finput;
  		my @splitLine = split("\t",$finput);
		if($splitLine[1] eq "Per base N content" || $splitLine[1] eq "Kmer Content")
		{
			if($splitLine[0] eq "FAIL")
				{$qcCheck++;}
		}
	}
	close FCHECK1;
	open (FCHECK2, "<$fileCheck2");
	while(my $finput = <FCHECK2>)
	{
		chomp $finput;
  		my @splitLine = split("\t",$finput);
		if($splitLine[1] eq "Per base N content" || $splitLine[1] eq "Kmer Content")
		{
			if($splitLine[0] eq "FAIL")
				{$qcCheck++;}
		}
	}
	close FCHECK2;
	my $statCheck1 = "/home/epidb/processing/".$filename."/".$fastqc1."/fastqc_data.txt";
	my $statCheck2 = "/home/epidb/processing/".$filename."/".$fastqc2."/fastqc_data.txt";
	my $totalSeq = 0;
	my $totalBases = 0;
	my $totalBCount = 0;
	my $lengthCheck1 = 0;
	my $lengthCheck2 = 0;
	open (SCHECK1, "<$statCheck1");
	while(my $finput = <SCHECK1>)
	{
  		chomp $finput;
  		my @splitLine = split("\t",$finput);
		if($splitLine[0] eq "Total Sequences")
			{$totalSeq += $splitLine[1];}
		elsif($splitLine[0] eq ">>Sequence Length Distribution")
			{$lengthCheck1 = 1;}
		elsif(!($splitLine[0] =~ />>/) && $lengthCheck1 == 1 && !($splitLine[0] =~ /#/))
		{
			my @lengthSplit = split(/-/,$splitLine[0]);
			if($splitLine[1] =~ /E/)
			{
				my @splitSciNo = split(/E/,$splitLine[1]);
				my $tempSeqCount = $splitSciNo[0]*(10**$splitSciNo[1]);
				$totalBases += ((($lengthSplit[0] + $lengthSplit[1])/2)*$tempSeqCount);
				$totalBCount += $tempSeqCount;
			}
			else
			{
				$totalBases += ((($lengthSplit[0] + $lengthSplit[1])/2)*$splitLine[1]);
				$totalBCount += $splitLine[1];
			}
		}
		elsif($splitLine[0] =~ />>/)
			{$lengthCheck1 = 0;}
	}
	close SCHECK1;
	open (DPCHECK1, "<$statCheck1");
	while(my $finput = <DPCHECK1>)
	{
  		chomp $finput;
  		if($finput =~ /#Total Duplicate Percentage/)
		{
			my @splitDP = split(/\t/, $finput);
			$duplicateP1 = $splitDP[1];
		}
	}
	close DPHECK1;
	open (SCHECK2, "<$statCheck2");
	while(my $finput = <SCHECK2>)
	{
		chomp $finput;
  		my @splitLine = split("\t",$finput);
		if($splitLine[0] eq "Total Sequences")
			{$totalSeq += $splitLine[1];}
		elsif($splitLine[0] eq ">>Sequence Length Distribution")
			{$lengthCheck2 = 1;}
		elsif(!($splitLine[0] =~ />>/) && $lengthCheck2 == 1 && !($splitLine[0] =~ /#/))
		{
			my @lengthSplit = split(/-/,$splitLine[0]);
			if($splitLine[1] =~ /E/)
			{
				my @splitSciNo = split(/E/,$splitLine[1]);
				my $tempSeqCount = $splitSciNo[0]*(10**$splitSciNo[1]);
				$totalBases += ((($lengthSplit[0] + $lengthSplit[1])/2)*$tempSeqCount);
				$totalBCount += $tempSeqCount;
			}
			else
			{
				$totalBases += ((($lengthSplit[0] + $lengthSplit[1])/2)*$splitLine[1]);
				$totalBCount += $splitLine[1];
			}
		}
		elsif($splitLine[0] =~ />>/)
			{$lengthCheck2 = 0;}
	}
	close SCHECK2;
	open (DPCHECK2, "<$statCheck2");
	while(my $finput = <DPCHECK2>)
	{
  		chomp $finput;
  		if($finput =~ /#Total Duplicate Percentage/)
		{
			my @splitDP = split(/\t/, $finput);
			$duplicateP2 = $splitDP[1];
		}
	}
	close DPHECK2;
	$avgCheck = $totalBases/$totalBCount;
	print $srabase." Total Sequences ".$totalSeq."\n";
	print $srabase." Average Sequence Length ".$avgCheck."\n";
}	
elsif($library eq "SINGLE")
{
	my $srabase = $sraname;
   	$srabase =~ s/.sra//;
   	my $fastqname = $sraname;
   	$fastqname =~ s/.sra/.fastq/;
   	my $trimmed = $sraname;
   	$trimmed =~ s/.sra/.trimmed.fastq/;
	my $fastqc = $trimmed;
	$fastqc =~ s/.fastq/_fastqc/;
	my $fileCheck = "/home/epidb/processing/".$filename."/".$fastqc."/summary.txt";
	open (FCHECK, "<$fileCheck");
	while(my $finput = <FCHECK>)
	{
  		chomp $finput;
  		my @splitLine = split("\t",$finput);
		if($splitLine[1] eq "Per base N content" || $splitLine[1] eq "Kmer Content")
		{
			if($splitLine[0] eq "FAIL")
				{$qcCheck++;}
		}
	}
	my $statCheck = "/home/epidb/processing/".$filename."/".$fastqc."/fastqc_data.txt";
	my $totalSeq = 0;
	my $totalBases = 0;
	my $totalBCount = 0;
	my $lengthCheck = 0;
	open (SCHECK, "<$statCheck");
	while(my $finput = <SCHECK>)
	{
  		chomp $finput;
  		my @splitLine = split("\t",$finput);
		if($splitLine[0] eq "Total Sequences")
			{$totalSeq += $splitLine[1];}
		elsif($splitLine[0] eq ">>Sequence Length Distribution")
			{$lengthCheck = 1;}
		elsif(!($splitLine[0] =~ />>/) && $lengthCheck == 1 && !($splitLine[0] =~ /#/))
		{
			my @lengthSplit = split(/-/,$splitLine[0]);
			if($splitLine[1] =~ /E/)
			{
				my @splitSciNo = split(/E/,$splitLine[1]);
				my $tempSeqCount = $splitSciNo[0]*(10**$splitSciNo[1]);
				$totalBases += ((($lengthSplit[0] + $lengthSplit[1])/2)*$tempSeqCount);
				$totalBCount += $tempSeqCount;
			}
			else
			{
				$totalBases += ((($lengthSplit[0] + $lengthSplit[1])/2)*$splitLine[1]);
				$totalBCount += $splitLine[1];
			}
		}
		elsif($splitLine[0] =~ />>/)
			{$lengthCheck = 0;}
	}
	close SCHECK;
	open (DPCHECK, "<$statCheck");
	while(my $finput = <DPCHECK>)
	{
  		chomp $finput;
  		if($finput =~ /#Total Duplicate Percentage/)
		{
			my @splitDP = split(/\t/, $finput);
			$duplicateP = $splitDP[1];
		}
	}
	close DPHECK;
	$avgCheck = $totalBases/$totalBCount;
	print $srabase." Total Sequences ".$totalSeq."\n";
	print $srabase." Average Sequence Length ".$avgCheck."\n";
}
if($library eq "PAIRED")
{
	if($qcCheck < 2 && $avgCheck > 39 && $duplicateP1 < 75 && $duplicateP2 < 75)
	{
		my $srabase = $sraname;
   		$srabase =~ s/.sra//;
   		my $fastq1 = $sraname;
   		$fastq1 =~ s/.sra/_1.fastq/;
   		my $fastq2 = $sraname;
   		$fastq2 =~ s/.sra/_2.fastq/;
   		my $trimmed1 = $sraname;
   		$trimmed1 =~ s/.sra/.trimmed_1.fastq/;
   		my $trimmed2 = $sraname;
   		$trimmed2 =~ s/.sra/.trimmed_2.fastq/;
   		my $trimmedSingle = $sraname;
   		$trimmedSingle =~ s/.sra/.trimmed_single.fastq/;
		my $fastqc1 = $trimmed1;
		$fastqc1 =~ s/.fastq/_fastqc/;
		my $fastqc2 = $trimmed2;
		$fastqc2 =~ s/.fastq/_fastqc/;
		$systemArgs = "~/tophat2/tophat2 -p 12 -o ".$srabase.".out /home/epidb/index/".$species."/current_build ".$trimmed1." ".$trimmed2;
		system("$systemArgs");
		$systemArgs = "~/tophat2/samtools sort ".$srabase.".out/accepted_hits.bam ".$srabase.".out/accepted_hits.sorted";
		system("$systemArgs");
		$systemArgs = "java -jar ~/picard/AddOrReplaceReadGroups.jar INPUT=".$srabase.".out/accepted_hits.sorted.bam OUTPUT=".$srabase.".fixed.bam RGLB=1 RGPL=illumina RGPU=allruns RGSM=".$srabase." VALIDATION_STRINGENCY=SILENT";
		system("$systemArgs");
		$systemArgs = "~/tophat2/samtools index ".$srabase.".fixed.bam";
		system("$systemArgs");
		$systemArgs = "~/tophat2/samtools view -h ".$srabase.".fixed.bam > ".$srabase.".fixed.sam";
		system("$systemArgs");
		$systemArgs = "python ~/HTSeq-0.6.1/HTSeq/scripts/count.py -m intersection-nonempty -s no -i transcript_id ".$srabase.".fixed.sam /home/epidb/index/bos_taurus/Bos_taurus.UMD3.1.75.fixed.gtf > ".$srabase.".htseq.txt";
		system("$systemArgs");
		$systemArgs = "rm ".$srabase.".fixed.sam";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -R /home/epidb/index/".$species."/current_build.fa -T SplitNCigarReads -I ".$srabase.".fixed.bam -o ".$srabase.".split.bam -allowPotentiallyMisencodedQuals -U ALLOW_N_CIGAR_READS";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -nt 12 -R /home/epidb/index/".$species."/current_build.fa -T RealignerTargetCreator -I ".$srabase.".split.bam -o ".$srabase.".intervals";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -R /home/epidb/index/".$species."/current_build.fa -T IndelRealigner -I ".$srabase.".split.bam -targetIntervals ".$srabase.".intervals -o ".$srabase.".realigned.bam --maxReadsForRealignment 50000";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -nct 12 -R /home/epidb/index/".$species."/current_build.fa -T BaseRecalibrator -I ".$srabase.".realigned.bam -o ".$srabase.".table -knownSites ~/index/".$species."/dbSNP.vcf";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -nct 12 -R /home/epidb/index/".$species."/current_build.fa -T PrintReads -I ".$srabase.".realigned.bam -o ".$srabase.".recalc.bam -BQSR ".$srabase.".table";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -nct 12 -R /home/epidb/index/".$species."/current_build.fa -T HaplotypeCaller -I ".$srabase.".recalc.bam -o ".$srabase.".snp.indel.vcf -dontUseSoftClippedBases -dt NONE -stand_call_conf 20.0 -stand_emit_conf 20.0";
		system("$systemArgs");
		$systemArgs = "~/tabix-0.2.6/bgzip ".$srabase.".snp.indel.vcf";
		print $systemArgs."\n";
		system("$systemArgs");
		$systemArgs = "~/tabix-0.2.6/tabix -p vcf ".$srabase.".snp.indel.vcf.gz";
		print $systemArgs."\n";
		system("$systemArgs");
		$systemArgs = "rm -r ".$srabase.".trimmed_1_fastqc.zip ".$srabase.".trimmed_2_fastqc.zip ".$srabase.".out";
		system("$systemArgs");
		$systemArgs = "~/rsem-1.2.20/rsem-calculate-expression -p 12 --bowtie2 --bowtie2-path ~/tophat2/ --paired-end ".$trimmed1." ".$trimmed2." ~/index/".$species."/rsem/current_build_transcriptome ".$srabase.".rsem.output";
		system("$systemArgs");
		$systemArgs = "rm -r ".$srabase.".rsem.output.transcript.bam ".$srabase.".rsem.output.stat";
		system("$systemArgs");
		$systemArgs = "awk '{print\$1\"\\t\"\$6}' ".$srabase.".rsem.output.genes.results > ".$srabase.".genes.txt";
		system("$systemArgs");
	}
	else
	{
		my $srabase = $sraname;
   		$srabase =~ s/.sra//;
		print $srabase." Failed ".$qcCheck." QC Filters\n";
		print $srabase." Duplicate Percent: ".$duplicateP1." and ".$duplicateP2."\n";
	}
}	
elsif($library eq "SINGLE")
{
	if($qcCheck < 1 && $avgCheck > 39 && $duplicateP < 75)
	{
		my $srabase = $sraname;
   		$srabase =~ s/.sra//;
   		my $fastqname = $sraname;
   		$fastqname =~ s/.sra/.fastq/;
   		my $trimmed = $sraname;
   		$trimmed =~ s/.sra/.trimmed.fastq/;
		my $fastqc = $trimmed;
		$fastqc =~ s/.fastq/_fastqc/;
		$systemArgs = "~/tophat2/tophat2 -p 12 -o ".$srabase.".out /home/epidb/index/".$species."/current_build ".$trimmed;
		system("$systemArgs");
		$systemArgs = "~/tophat2/samtools sort ".$srabase.".out/accepted_hits.bam ".$srabase.".out/accepted_hits.sorted";
		system("$systemArgs");
		$systemArgs = "java -jar ~/picard/AddOrReplaceReadGroups.jar INPUT=".$srabase.".out/accepted_hits.sorted.bam OUTPUT=".$srabase.".fixed.bam RGLB=1 RGPL=illumina RGPU=allruns RGSM=".$srabase." VALIDATION_STRINGENCY=SILENT";
		system("$systemArgs");
		$systemArgs = "~/tophat2/samtools index ".$srabase.".fixed.bam";
		system("$systemArgs");
		$systemArgs = "~/tophat2/samtools view -h ".$srabase.".fixed.bam > ".$srabase.".fixed.sam";
		system("$systemArgs");
		$systemArgs = "python /home/epidb/HTSeq-0.6.1/HTSeq/scripts/count.py -m intersection-nonempty -s no -i transcript_id ".$srabase.".fixed.sam /home/epidb/index/bos_taurus/Bos_taurus.UMD3.1.75.fixed.gtf > ".$srabase.".htseq.txt";
		system("$systemArgs");
		$systemArgs = "rm ".$srabase.".fixed.sam";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -R /home/epidb/index/".$species."/current_build.fa -T SplitNCigarReads -I ".$srabase.".fixed.bam -o ".$srabase.".split.bam -allowPotentiallyMisencodedQuals -U ALLOW_N_CIGAR_READS";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -nt 12 -R /home/epidb/index/".$species."/current_build.fa -T RealignerTargetCreator -I ".$srabase.".split.bam -o ".$srabase.".intervals";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -R /home/epidb/index/".$species."/current_build.fa -T IndelRealigner -I ".$srabase.".split.bam -targetIntervals ".$srabase.".intervals -o ".$srabase.".realigned.bam --maxReadsForRealignment 50000";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -nct 12 -R /home/epidb/index/".$species."/current_build.fa -T BaseRecalibrator -I ".$srabase.".realigned.bam -o ".$srabase.".table -knownSites ~/index/".$species."/dbSNP.vcf";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -nct 12 -R /home/epidb/index/".$species."/current_build.fa -T PrintReads -I ".$srabase.".realigned.bam -o ".$srabase.".recalc.bam -BQSR ".$srabase.".table";
		system("$systemArgs");
		$systemArgs = "~/java1.7/bin/java -jar -Xmx48g ~/GATK3.4/GenomeAnalysisTK.jar -nct 12 -R /home/epidb/index/".$species."/current_build.fa -T HaplotypeCaller -I ".$srabase.".recalc.bam -o ".$srabase.".snp.indel.vcf -dontUseSoftClippedBases -dt NONE -stand_call_conf 20.0 -stand_emit_conf 20.0";
		system("$systemArgs");
		$systemArgs = "~/tabix-0.2.6/bgzip ".$srabase.".snp.indel.vcf";
		print $systemArgs."\n";
		system("$systemArgs");
		$systemArgs = "~/tabix-0.2.6/tabix -p vcf ".$srabase.".snp.indel.vcf.gz";
		print $systemArgs."\n";
		system("$systemArgs");
		$systemArgs = "rm -r ".$srabase.".trimmed_fastqc.zip ".$srabase.".out";
		system("$systemArgs");
		$systemArgs = "~/rsem-1.2.20/rsem-calculate-expression -p 12 --bowtie2 --bowtie2-path ~/tophat2/ ".$trimmed." ~/index/".$species."/rsem/current_build_transcriptome ".$srabase.".rsem.output";
		system("$systemArgs");
		$systemArgs = "rm -r ".$srabase.".rsem.output.transcript.bam ".$srabase.".rsem.output.stat";
		system("$systemArgs");
		$systemArgs = "awk '{print\$1\"\\t\"\$6}' ".$srabase.".rsem.output.genes.results > ".$srabase.".genes.txt";
		system("$systemArgs");
	}
	else
	{
		my $srabase = $sraname;
   		$srabase =~ s/.sra//;
		print $srabase." Failed ".$qcCheck." QC Filters\n";
		print $srabase." Duplicate Percent: ".$duplicateP."\n";
	}
}