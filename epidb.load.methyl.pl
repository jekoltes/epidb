#!/usr/bin/perl

use warnings;
use strict;
use DBI;
 
my $db2 ="EpiDB";
my $user = "eric";
my $pass = "R1ftw4lker";
my $host="localhost";
my $dbh2 = DBI->connect("DBI:mysql:$db2:$host", $user, $pass);
my $filename = $ARGV[0];
my $species = $ARGV[1];
my $header = 0;

my %srsHash = ();
my $tempSRSID = "";
my $checkQuery = "SELECT SRS_ID from EpiDB_Methyl_track";
my $checkSTH = $dbh2->prepare($checkQuery);
$checkSTH->execute();
$checkSTH->bind_columns(\$tempSRSID);
while($checkSTH->fetch())
{
	$srsHash{$tempSRSID} = $tempSRSID;
}
    
open (GENE, "<$filename");
while(my $finput = <GENE>)
{
        chomp $finput;
        	my @lineSplit = split(/\t/, $finput);
        	my $temp_SRS_ID = $lineSplit[0];
        	my $temp_Updated = $lineSplit[1];
        	my $temp_NumberOfBases = $lineSplit[2];
        	my $temp_RunCenter = $lineSplit[3];
        	my $temp_SRAdb_experiment_id = $lineSplit[4];
        	my $temp_GEO_GSM_ID = $lineSplit[5];
        	my $temp_experiment_title = $lineSplit[6];
        	my $temp_GEO_GSE_ID = $lineSplit[7];
        	my $temp_sample_name = $lineSplit[8];
        	my $temp_library_layout = $lineSplit[9];
        	my $temp_instrument_model = $lineSplit[10];
        	my $temp_SRAdb_sample_id = $lineSplit[11];
        	my $temp_GSM_ID = $lineSplit[12];
        	my $temp_sample_attribute = $lineSplit[13];
        	my $temp_SRAdb_study_id = $lineSplit[14];
        	my $temp_study_alias = $lineSplit[15];
        	my $temp_study_title = $lineSplit[16];
        	my $temp_study_abstract = $lineSplit[17];
        	my $temp_study_attribute = $lineSplit[18];
        	my $temp_SRA_ID = $lineSplit[19];
        	my $temp_submission_center = $lineSplit[20];
        	my $temp_submission_lab = $lineSplit[21];
        	my $temp_submission_date = $lineSplit[22];
        	my $temp_SRA_experiment_id = $lineSplit[23];
        	my $temp_SRA_study_id = $lineSplit[24];
        	my $temp_SRA_sample_id = $lineSplit[25];
        	my $temp_SRA_run_id = $lineSplit[26];
        	my $temp_asperaPath = $lineSplit[27];
        	my $temp_ftpPath = $lineSplit[28];
        	my $temp_tissue = $lineSplit[30];
        	
        	if($temp_Updated =~ /\//)
        	{
        		my @splitUpdated = split(/\//, $temp_Updated);
        		if($splitUpdated[0] < 10)
        		{
        			if($splitUpdated[1] < 10)
        			{
        				$temp_Updated = "20".$splitUpdated[2]."-0".$splitUpdated[0]."-0".$splitUpdated[1];
        			}
        			else
        			{
        				$temp_Updated = "20".$splitUpdated[2]."-0".$splitUpdated[0]."-".$splitUpdated[1];
        			}
        		}
        		else
        		{
        			if($splitUpdated[1] < 10)
        			{
        				$temp_Updated = "20".$splitUpdated[2]."-".$splitUpdated[0]."-0".$splitUpdated[1];
        			}
        			else
        			{
        				$temp_Updated = "20".$splitUpdated[2]."-".$splitUpdated[0]."-".$splitUpdated[1];
        			}
        		}
        	}
        	if($temp_submission_date =~ /\//)
        	{
        		my @splitSubDate = split(/\//, $temp_submission_date);
        		if($splitSubDate[0] < 10)
        		{
        			if($splitSubDate[1] < 10)
        			{
        				$temp_submission_date = "20".$splitSubDate[2]."-0".$splitSubDate[0]."-0".$splitSubDate[1];
        			}
        			else
        			{
        				$temp_submission_date = "20".$splitSubDate[2]."-0".$splitSubDate[0]."-".$splitSubDate[1];
        			}
        		}
        		else
        		{
        			if($splitSubDate[1] < 10)
        			{
        				$temp_submission_date = "20".$splitSubDate[2]."-".$splitSubDate[0]."-0".$splitSubDate[1];
        			}
        			else
        			{
        				$temp_submission_date = "20".$splitSubDate[2]."-".$splitSubDate[0]."-".$splitSubDate[1];
        			}
        		}
        	}
        	
        	$temp_tissue =~ s/^\s//;
        	$temp_tissue =~ s/\s/_/g;
        	$temp_sample_attribute =~ s/\"//g;
        	$temp_library_layout = substr($temp_library_layout, 0, 6);
        	
        	my $otherQuery = "SELECT tissue_ID FROM EpiDB_Tissues WHERE tissuetype='".$temp_tissue."'";
        	my $tissueSTH = $dbh2->prepare($otherQuery);
        	$tissueSTH->execute();
			my $tempTissueID = "";
			$tissueSTH->bind_columns(\$tempTissueID);
			$tissueSTH->fetch();
			
			if(!($srsHash{$temp_SRS_ID}))
			{
				my $iquery1 = "INSERT INTO EpiDB_Methyl_track (SRS_ID,Updated,NumberOfBases,RunCenter,SRAdb_experiment_id,GEO_GSM_ID,experiment_title,GEO_GSE_ID,sample_name,library_layout,instrument_model,SRAdb_sample_id,GSM_ID,sample_attribute,SRAdb_study_id,study_alias,study_title,study_abstract,study_attribute,SRA_ID,submission_center,submission_lab,submission_date,SRA_experiment_id,SRA_study_id,SRA_sample_id,SRA_run_id,asperaPath,ftpPath,tissue_ID,species_ID) VALUES ('".$temp_SRS_ID."','".$temp_Updated."','".$temp_NumberOfBases."','".$temp_RunCenter."','".$temp_SRAdb_experiment_id."','".$temp_GEO_GSM_ID."','".$temp_experiment_title."','".$temp_GEO_GSE_ID."','".$temp_sample_name."','".$temp_library_layout."','".$temp_instrument_model."','".$temp_SRAdb_sample_id."','".$temp_GSM_ID."','".$temp_sample_attribute."','".$temp_SRAdb_study_id."','".$temp_study_alias."','".$temp_study_title."','".$temp_study_abstract."','".$temp_study_attribute."','".$temp_SRA_ID."','".$temp_submission_center."','".$temp_submission_lab."','".$temp_submission_date."','".$temp_SRA_experiment_id."','".$temp_SRA_study_id."','".$temp_SRA_sample_id."','".$temp_SRA_run_id."','".$temp_asperaPath."','".$temp_ftpPath."','".$tempTissueID."','".$species."')";
        		my $isth1 = $dbh2->prepare($iquery1);
            	$isth1->execute();
            }
}
close GENE;
