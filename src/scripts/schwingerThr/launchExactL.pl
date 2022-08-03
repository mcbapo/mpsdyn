#!/usr/bin/perl

# Launch a series of jobs that run the imaginary time evolution of the
# Schwinger model, to find the thermal states

$inCluster=1;

$EXEC="./thermExL";

$MSUB="msub_modified_tqo097";

#$parallel=""; #"-P 2 "; # parallel
$parallel="-P 4 "; # parallel

$basicOutputDir="/ptmp/mpq/banulsm/SchwingerThrm/exactL";
$configFile="/afs/ipp/u/banulsm/SVN/mpsdyn/src/config/thermExL.conf";

$jobsDir="/afs/ipp/u/banulsm/jobsToDo/";

# Everything!
@Xs=(1,4,9,16,25,36,49,64,81,100,121,144,169,196,256,400,484);
@Ds=(60,80,100,120,140);
@Ls=(16,20,24);
@Ks=(1,2,4); # fractions of the maximum delta step
$betaMax=8;
#$stepBeta=0.01;
$stepBeta=0.1;
$Lmax=10;
#$Lmax=15;
$mg=0;
$tol=1E-8;
$append=1;
$savingFreq=20;

# Special, larger x, only for small beta
#@Xs=(529,625);
#@Ds=(80,100,120,140);
#$betaMax=.5;
@Xs=(1024,2500); # 32^2, 50^2
@Ds=(100,120,140);
$parallel="-P 8 "; # parallel
#$mem="-l h_vmem=32000M";
$betaMax=.5;

@Xs=(1024); # 32^2, 50^2
@Ds=(140,160);
$parallel="-P 20 "; # parallel
$mem="-l h_vmem=60000M";
$betaMax=.5;
@Ks=(1);
@Ls=(24);



my $scriptfile="ScriptLaunchThermExL.sh";
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
			    $MEMreq="-l h_vmem=5000M "
			}
			else{
			    my $myMEM=$MEM+1000;
			    $MEMreq="-l h_vmem=${myMEM}M "
			}
		    }
		}
		print "Job $NAME requires $MEM M, requesting $MEMreq \n";
		# Check if the job is running!
		my $isRunning=`qstat -j $NAME |grep job_number`;
		if($isRunning =~ /(\d+)/){
		    $jobNr=$1;
		    print MYFILE "##****!!!!*****\n# Job running $jobNr for $NAME \n##****** !!!!!";
		}
		else{
		    print MYFILE "mkdir -p $outputdir \n";
		    print MYFILE "mkdir -p $mpsdir \n";
		    print MYFILE "$MSUB $parallel ${MEMreq} -N $NAME -raw $EXEC ${args} \n";
		    print MYFILE "sleep 1 \n";
		}
	    } # end for D
	}# end for k (delta)
    }# end for factor (N)
}# end for x

close(MYFILE);
system("chmod a+x ${scriptfile} ");

