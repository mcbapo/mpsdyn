#!/usr/bin/perl

# Launch a series of jobs that run the imaginary time evolution of the
# Schwinger model, to find the thermal states
# This will be just a few tests to compare Lcut= 5, 10, 15 with the exact,
# for x=55, for which I have the exact value

$inCluster=1;

$EXEC="./thermExL";

$MSUB="msub_modified_tqo097";

#$parallel=""; #"-P 2 "; # parallel
$parallel="-P 4 "; # parallel

$basicOutputDir="/ptmp/mpq/banulsm/SchwingerThrm/exactL";
$configFile="/afs/ipp/u/banulsm/SVN/mpsdyn/src/config/thermExL.conf";

$jobsDir="/afs/ipp/u/banulsm/jobsToDo/";

@Xs=(55);
@Ds=(60,80,100,120,140);
#@Ls=(16,20,24);
@Ns=(126,150,170); # by hand as for Hana's code
@Ks=(1,2,4); # fractions of the maximum delta step
$betaMax=4.;
#$stepBeta=0.01;
$stepBeta=0.1;
@Lmaxs=(5,8,10,12,15);
## relaunch 15 D=140 with more memory
#@Lmaxs=(15);
#@Ds=(140);

#$Lmax=15;
$mg=0;
$tol=1E-8;
$append=1;

$parallel="-P 2 "; # parallel


my $scriptfile="ScriptLaunch_Lcuts.sh";
open MYFILE, '>', $scriptfile or die "Cannot open ${scriptfile} $!";
print MYFILE "#!/bin/bash\n\n";

foreach $x (@Xs){
    # Which system sizes
#    foreach $factor (@Ls){
    foreach $N (@Ns){
	#my $N=int($factor*sqrt($x));
	my $delta=$stepBeta*.25/sqrt($x); # max delta step: stepBeta/(2 sqrt(x))*1/2
        # I take that delta, and then 1/k of it
	foreach $k (@Ks){
	    my $deltak=$delta/$k;
	    my $M=int($betaMax/$stepBeta+1)*$k; # nr of steps
	    my $outputdir="$basicOutputDir/results/x$x";
	    my $mpsdir="$basicOutputDir/mps/x$x";
	    my $savingFreq=2*$k;
	    foreach $D (@Ds){
		foreach $Lmax (@Lmaxs){
		    my $outfile=sprintf("x%d.N%d.dt%.4E.D%d.Lmax%d.p%d",$x,$N,$deltak,$D,$Lmax,int(-log($tol)/log(10)));
		    my $args=" $configFile -L=$N -mg=$mg -x=$x -delta=$deltak -M=$M -outputdir=$outputdir -outfile=$outfile -mpsdir=$mpsdir -append=$append -tol=$tol -D=$D -blockSteps=$k -freqBkup=$savingFreq -Lmax=$Lmax ";

		    my $NAME=sprintf("x%d.N%d.Lc%d.dt%.4E.D%d",$x,$N,$Lmax,$deltak,$D);
		    # Estimate the required memory (in MB): for an MPS is
		    # N*D^2*dphys*16(two doubles) For an optimization we
		    # have two MPs of the same size, plus the temporary
		    # storage, namely 2N terms, each of them size D*D*Xi,
		    # where Xi is the largest bond dimension plus the
		    # MPOs, which are N times d^2*Xi*d^2*Xi (approx) for
		    # each of the L terms and d^4*4^2 for each of the He, Ho (three of them)
		    my $sizeMPS=16*$N*$D*$D*4;
		    my $sizeTmp=16*$N*2*$D*$D*(2*$Lmax+1);
		    my $sizeMPOsL=16*(16*(2*$Lmax+1)*(2*$Lmax+1)*$N);
		    my $sizeMPOsx=16*$N*16*4*4*3;
		    my $sizeAuxMPOs=16*4*(2*$Lmax+1)*(2*$Lmax+1)*$N;
		    my $MEM=int((3*$sizeMPS+$sizeTmp+2*$sizeMPOsL+3*$sizeMPOsx+$sizeAuxMPOs)*1.2*1E-6)+200; # I will add a 20% AND 200M
		    if($MEM<=1000){
			$MEMreq="-l h_vmem=2000M ";
		    }
		    else{
			if($MEM<2500){ # Memory request below 2500
			    $MEMreq="";
			}
			else{
			    if($MEM<=4000){
				$MEMreq="-l h_vmem=5000M ";
				if($Lmax==15&&$D==140){
				    $MEMreq="-l h_vmem=6000M "; # try larger mem 
				}
			    }
			    else{
				my $myMEM=$MEM+1000;
				if($Lmax==15&&$D==140){
				    $myMEM=$myMEM+1000; # try larger mem 
				}
				$MEMreq="-l h_vmem=${myMEM}M "
			    }
			}
		    }
		    print "Job $NAME requires $MEM M, requesting $MEMreq \n";
		    # Check if the job is running!
		    my $isRunning=`qstat -j $NAME |grep job_number`;
		    if($isRunning =~ /(\d+)/){
			$jobNr=$1;
			print MYFILE "##****!!!!*****\n# Job running $jobNr for $NAME \n##****** !!!!!\n";
		    }
		    else{
			# Check if already finished
			my $fileName="$outputdir/$outfile";
			my $lastLine=`tail -n 1 $fileName`;
			my @entries=split(/\s+/,$lastLine); # entries in the last line of the file (3rd one is g*beta)
			my $gbeta=@entries[2];
			if($gbeta<$betaMax){
			    print "Job $NAME not running and apparently not finished either (gbeta=$gbeta)\n";
			    print MYFILE "mkdir -p $outputdir \n";
			    print MYFILE "mkdir -p $mpsdir \n";
			    print MYFILE "$MSUB $parallel ${MEMreq} -N $NAME -raw $EXEC ${args} \n";
			    print MYFILE "sleep 1 \n";
			}
			else{
			    print "Job $NAME was already completed with gbeta=$gbeta !\n";
			    print MYFILE "cleanjobs $NAME\n";
			}
		    }
		} # end for Lmax
	    } # end for D
	}# end for k (delta)
    }# end for factor (N)
}# end for x

close(MYFILE);
system("chmod a+x ${scriptfile} ");

