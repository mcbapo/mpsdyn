#!/usr/bin/perl

# Launch a series of jobs that run the imaginary time evolution 
# to find thermal states of the Schwinger model

$inCluster=1;

$EXEC="./therm_eLoLe";

$MSUB="msub_modified_tqo097";


$basicOutputDir="/ptmp/mpq/banulsm/SchwingerThrm/Hana3";

$jobsDir="/afs/ipp/u/banulsm/jobsToDo/";

# Arguments
$mg=0;
$L0=0;
$alpha=0;
$distanceBeta=0.001;

# L  mg  x  L0  alpha  D  M  delta  init_MPS  outfname  bkdrctnm_init  Intrvl_bckup  distanceBeta  app

#Sets I need to launch:
# x N delta
#55  126  3e-5, 1e-5
#55  150  2e-5, 1e-5
#55 170  2e-5, 1e-5
#60 130  3e-5, 1e-5
#60 160  2e-5, 1e-5
#60 180  2e-5, 1e-5
#65  140  3e-5, 1e-5
#65  160  2e-5, 1e-5
#65  190  1e-5, 5e-6


my $scriptfile="ScriptThermHana_ALL.sh";
#if($inCluster){$scriptfile="ScriptLanzadorDeRFMI.sh";};
open MYFILE, '>', $scriptfile or die "Cannot open ${scriptfile} $!";
print MYFILE "#!/bin/bash\n\n";

for (my $case=6;$case>=0;$case=$case-1){

if($case==0){
@Xs=(55);@Ls=(126);@deltas=(1E-5,3E-5);
$parallel="-P 8 "; # parallel
}
if($case==1){
@Xs=(55);@Ls=(150,170);@deltas=(1E-5,2E-5);
$parallel="-P 10 "; # parallel
}
if($case==2){
@Xs=(60);@Ls=(130);@deltas=(1E-5,3E-5);
$parallel="-P 8 "; # parallel
}
if($case==3){
@Xs=(60);@Ls=(160,180);@deltas=(1E-5,2E-5);
$parallel="-P 10 "; # parallel
}
if($case==4){
@Xs=(65);@Ls=(140);@deltas=(1E-5,3E-5);
$parallel="-P 10 "; # parallel
}
if($case==5){
@Xs=(65);@Ls=(160);@deltas=(1E-5,2E-5);
$parallel="-P 10 "; # parallel
}
if($case==6){
@Xs=(65);@Ls=(190);@deltas=(5E-6,1E-5);
$parallel="-P 10 "; # parallel
}

@Ds=(60,80,100);

#@Xs=(5);@Ds=(20);@deltas=(1E-6);$parallel="";
# my $scriptfile="ScriptThermHana_x{@Xs[0]}_L{@Ls[0]}.sh";
# #if($inCluster){$scriptfile="ScriptLanzadorDeRFMI.sh";};
# open MYFILE, '>', $scriptfile or die "Cannot open ${scriptfile} $!";
# print MYFILE "#!/bin/bash\n\n";

$useMPS=0; # not using earlier file
$nrFile=1; # hence, saved to 1 for later
$maxBeta=4.0; # max value of beta we want to reach
$append=0;
$freqBackup=50; # frequency of backups

foreach $X (@Xs){
    foreach $L (@Ls){
	foreach $delta (@deltas){
	    foreach $D (@Ds){
		$outputdir="${basicOutputDir}/x${X}";
		$mpsdir="${basicOutputDir}/x${X}/MPS";
		print MYFILE "mkdir -p $outputdir \n";
		print MYFILE "mkdir -p $mpsdir \n";
		# Nr steps
		$M=int($maxBeta/(2.*sqrt($X)*2*$delta)+1.);
		# output file
		$outFileName="${outputdir}/output.file.N${L}D${D}dlt${delta}_Thrm_${nrFile}";
		$backupFileName="${mpsdir}/backup_N${L}D${D}dlt${delta}_Thrm_${nrFile}";
		## That was wrong, should have been
		## $outFileName="${outputdir}/output.file.N${L}D${D}dlt${delta}_Thrm";
		## $backupFileName="${mpsdir}/backup_N${L}D${D}dlt${delta}_Thrm_";
		my $args="$L $mg $X $L0 $alpha $D $M $delta $useMPS $outFileName $backupFileName $freqBackup $distanceBeta $append";

		print "Launching for x=$X, N=$L, delta=$delta, D=$D \n";				
		my $NAME="x${X}.N${L}.dt${delta}.D${D}";
		print MYFILE "$MSUB ${parallel} -N $NAME -raw $EXEC ${args} \n";
		print MYFILE "sleep 1\n";
	    }
	}
    }
}

#close(MYFILE);
#system("chmod a+x ${scriptfile} ");

}

close(MYFILE);
system("chmod a+x ${scriptfile} ");
