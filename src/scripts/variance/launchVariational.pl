#!/usr/bin/perl

# Script to prepare a series of jobs using the minH2 program

# Name of the script file to be created
my $scriptFile="lanzaMinVar.sh";

my $EXEC="./minVar";

my $configFile="../config/randomMPO.conf";

##my $basicDir="/ptmp/mpq/banulsm/variance/Ising/ChebyDelta/";
my @models=("Ising","Heisenberg","XY");
#my $model=0; #"Ising";
my $model=1; #"Heisenberg";
#my $model=2; #"XY";

my $modelStr=$models[$model];
my $basicDir="/ptmp/mpq/banulsm/variance/2/${modelStr}/varDelta";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

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

my $E0=0;  
my $penH=100;
@Ns=(20,40,60,50,80,100);
# Fitst bunch
#my $D0=20;
my $D0=40;
my $incrD=20;
my $maxD=100; # then from 100 to 1000
my $par=" -P 2 "; # to start
# Larger Ds
#my $D0=200;
#my $incrD=100;
#my $maxD=500; 
#my $par=" -P 10 "; # to start

# TRick for small Ds
my $D0=2;
my $incrD=2;
my $maxD=10; # then from 100 to 1000
my $par=" -P 2 "; # to start


if($model==0){ # Ising case
    print "Case Ising \n";
    # Hamiltonian parameters
    @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 ");
    @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0"); # and an integrable one
    
    @hamilPars=(" -J=1 -g=-1.05 -h=0.5 ");
    @idStrings=("J1_gm105_h05");


}
else{
    if($model==1){ # Heisenberg
	print "Case Heisenberg \n";
	#    @hamilPars=(" -Jx=1 -Jy=1 -Jz=1 -h=0. "," -Jx=1 -Jy=1 -Jz=1 -h=0.5 ");
	#    @idStrings=("Jx1_Jy1_Jz1_h0","Jx1_Jy1_Jz1_h05");
	#@hamilPars=(" -Jx=1.1 -Jy=-1 -Jz=.9 -h=0. -hx=0. -hy=0. "," -Jx=1 -Jy=1 -Jz=1 -h=0.5 -hx=0. -hy=0. ");
	#@idStrings=("Jx1o1_Jym1_Jz0o9_h0","Jx1_Jy1_Jz1_h05");
	@hamilPars=(" -Jx=1.1 -Jy=-1 -Jz=.9 -h=0.6 -hx=0. -hy=0. ");#," -Jx=1 -Jy=1 -Jz=1 -h=0.5 -hx=0. -hy=0. ");
	@idStrings=("Jx1o1_Jym1_Jz0o9_h0o6"); #,"Jx1_Jy1_Jz1_h05");
	@labParH=("3");
    }
    else{ # XY
	@hamilPars=(" -Jx=1. -Jy=1. -Jz=0. -hx=0.8 -hy=0 -h=0.5 ");
	@idStrings=("Jx1_Jy1_hx0o8_hy0_hz0o5");
    }
}
my $app=0;

my $randSeed=117; # to start always the same
my $avrg=1;
my $Lcmax=6;
my $Lb=6;

$basicDir="${basicDir}/E0${E0}/";

# Directory for results
my $resultsDir="${basicDir}/results";

my $MSUB="msub_slurm";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

print MYFILE "mkdir -p ${resultsDir} \n";
#print MYFILE "mkdir -p ${mpsDir} \n";


my $cnt=0;
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
    print MYFILE "mkdir -p ${resultsDir}/${strL}/results \n";
    # TODO: Set a file for energies (Emin Emax) and check for it
    my $mpsdir="${resultsDir}/${strL}/MPS";
    print MYFILE "mkdir -p ${mpsdir} \n";
    foreach $N (@Ns){
	my $locpar=$par;
	if($N==100&&$D==1000){$locpar=" -P 40 ";}
	else{$locpar=$par;}
	my $suffix="_N${N}E0${E0}";
	# Name of the output file
	my $outfile="${resultsDir}/${strL}/results/data_${suffix}";
	my $outfileS="${resultsDir}/${strL}/results/dataSchmidt_${suffix}";
	my $JOBNAME="minVar.${modelStr}.N${N}.D${D0}.E0${E0}.par${cnt}";
	# TODO Check that the job is not yet running!!
	my $alreadyRunning=isRunning($JOBNAME);
	if(! $alreadyRunning ){
	    my $args=" -L=${N} ${paramSet} -D=${D0} -maxD=${maxD} -incrD=${incrD} -output=${outfile} -outputS=${outfileS} -E=${E0} -penH=${penH} -checkRDM=1 -Lcmax=$Lcmax -average=$avrg -Lb=$Lb -parNr=${cnt} -mpsdir=${mpsdir}  -jobsdir=${jobsdir} ";
	    if($model!=0){$args="$args -ising=0 ";}
	    if(!$alreadyRunning){
		print MYFILE "${MSUB} ${locpar} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
		print MYFILE "sleep 1\n";
	    }
	    else{
		print MYFILE "CHECK: ${MSUB} ${locpar} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
	    }
	}
	else{
	    print "Not launching $JOBNAME for it exists\n";
	}
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


