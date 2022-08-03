#!/usr/bin/perl

# Launch a series of jobs that run the imaginary time evolution of the
# Schwinger model, to find the thermal states
# This one should check, before launching, if existing file contains already the data


$EXEC="./thermExL";

my $MSUB="msub_slurm";


$basicOutputDir="/ptmp/mpq/irenepa/SchwingerThrm/exactL";
$configFile="/afs/ipp/u/irenepa/SVN/mpsdyn/src/config/thermExL.conf";

# Everything!
my @Xs=(1,4,9,16,25,36,49,64,81,100,121,144,169,196,256,400,484);
my @Ds=(80,100,120,140);
my @Ls=(16,20,24);
my @Ks=(1,2,4); # fractions of the maximum delta step
my $betaMax=8;
my $stepBeta=0.1;
my $Lmax=10;
#$Lmax=15;
my $mg=0;
my $tol=1E-8;
my $append=1;
my $savingFreq=20;

my $parallel="-P 4 ";


my $scriptfile="ScriptReLaunchThermExL.sh";




# Check for existing jobs
sub isRunning{
    my ($jobname)=@_;
    #print "Checking whether job $jobname exists\n";
    my $result=`squeue --noheader -n $jobname -o "%o" `;
    if (length $result){
	print "Job $jobname exists: $result \n";
	return 1;
    }
    return 0;
};



open MYFILE, '>', $scriptfile or die "Cannot open ${scriptfile} $!";
print MYFILE "#!/bin/bash\n\n";

foreach $x (@Xs){
    # Which system sizes
    foreach $factor (@Ls){
	my $N=int($factor*sqrt($x));
	my $delta=$stepBeta*.25/sqrt($x); # max delta step: stepBeta/(2 sqrt(x))*1/2
        # I take that delta, and then 1/k of it
	foreach $k (@Ks){
	    my $deltak=$delta/$k;
	    my $M=int($betaMax/$stepBeta+1)*$k; # nr of steps
	    my $outputdir="$basicOutputDir/results/x$x";
	    my $mpsdir="$basicOutputDir/mps/x$x";
	    foreach $D (@Ds){
		my $outfile=sprintf("x%d.N%d.dt%.4E.D%d.Lmax%d.p%d",$x,$N,$deltak,$D,$Lmax,int(-log($tol)/log(10)));
		my $args=" $configFile -L=$N -mg=$mg -x=$x -delta=$deltak -M=$M -outputdir=$outputdir -outfile=$outfile -mpsdir=$mpsdir -append=$append -tol=$tol -D=$D -blockSteps=$k -freqBkup=$savingFreq -Lmax=$Lmax ";

		my $JOBNAME=sprintf("x%d.N%d.Lc%d.dt%.4E.D%d",$x,$N,$Lmax,$deltak,$D);
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
		    $MEMreq="-mem 2000M ";
		}
		else{
		    if($MEM<2500){ # Memory request below 2500
			$MEMreq="";
		    }
		    else{
			if($MEM<=4000){
			    $MEMreq="-mem 5000M "
			}
			else{
			    my $myMEM=$MEM+1000;
			    $MEMreq="-mem ${myMEM}M "
			}
		    }
		}
		# Check if the job is running!
		# Decide if I am still (already) running
		my $alreadyRunning=isRunning($JOBNAME);
		if(! $alreadyRunning ){
		    # Check if the data have already been computed
		    my $fileName="$outputdir/$outfile";
		    my $lastLine=`tail -n 1 $fileName`;
		    my @entries=split(/\s+/,$lastLine); # entries in the last line of the file (3rd one is g*beta)
		    my $gbeta=@entries[2];
		    if($gbeta<$betaMax){
			print "Launching NEW: Job $JOBNAME requires $MEM M, requesting $MEMreq \n";
			print MYFILE "mkdir -p $outputdir \n";
			print MYFILE "mkdir -p $mpsdir \n";
			print MYFILE "$MSUB $parallel ${MEMreq} -N $JOBNAME -- $EXEC ${args} \n";
			print MYFILE "sleep 1 \n";
		    }
		    else{
			print "Job $JOBNAME was already completed with gbeta=$gbeta !\n";
			print MYFILE "cleanjobs $JOBNAME\n";
		    }
		}
		else{
		    print MYFILE "##****!!!!*****\n# Job $JOBNAME already running\n##****** !!!!!";
		}
	    } # end for D
	}# end for k (delta)
    }# end for factor (N)
}# end for x

close(MYFILE);
system("chmod a+x ${scriptfile} ");

