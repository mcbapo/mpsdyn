#!/usr/bin/perl

# Script to prepare a series of jobs using the twofermhubbard program

# Name of the script file to be created
my $scriptFile="lanzaFH2.sh";

my $EXEC="/u/banulsm/NewSVN/mpsdyn/src/bin/fH2gs";

my $basicDir="/ptmp/mpq/banulsm/fermihubbard2/new";
my $jobsdir="/u/banulsm/jobsToRun";
my $configfile="/u/banulsm/NewSVN/mpsdyn/src/config/twoorbfermihub.conf";


my $MSUB="msub_slurm";

my $te=1.;
my $te0=$te;
my $te1=$te;
my $tg=.5*$te;
my $tg0=$tg;
my $tg1=$tg;
my $U=10.;
my $Ug=$U;
my $Ue=$U;
my $mug0=0.;
my $mug1=0.;
my $mue0=0.;
my $mue1=0.;

my @Vs=(0..19);
#my @Vs=(9.2,9.4,9.6,9.8);
my @Vs=(8);


# dimensions for which I run the searchO program
my $D0=20; #, 60, 80, 100, 120, 140, 160);
my $D=20;
my $incrD=20;
my $MAXD=200;

my $tol=1E-5;

# frequency to save tmp MPS
my $freqSv=20;

# min number of states
#MINNUMBER="25"
my $MINNUMBER=2;
#Ds="40 60 80 100 120"
my $ONLYVEC=" -onlyvec=1 ";

my $app=1;
my $restart=0;
my $Drestart=200;

my $L=(50);

my $PENALTYNTOT=100; #100;
my $PENALTYSZ=100;

my $Ntot=$L; # quarter filling, as in the 2011 paper
my $Sztarget=0;

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


open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

# I need to check each imbalance. Now what I do is using (for each set
# of parameters) a file where to write the lowest energy found and, if
# the run with a certain imbalance finds energy that is higher than
# that by over 20%, it stops. It is nt 100% safe strategy, though
# (race condition for lock file, or differences at different Ds that
# may be larger than the cutoff but reversed later in D)
# (Q: a smarter strategy, e.g. some derivative of energy? to know when I'm done?)
#my @Ns=(0..$Ntot);


foreach my $V (@Vs){
    my $Vex=.5*($U-1.*$V);
    my $SIZEMEM=int(16*$L*$D*$D*(10+4*$MINNUMBER+2)*1.1*1E-6)+200;
	
    my $LONGMEM="-mem ${SIZEMEM}M";
    if ($SIZEMEM<3000) {
	$LONGMEM="";
    }
	
    my $par="-P 10";
    if ($D>=40){
	$par="-P 20";
    }
    if ($D>=80){
	$par="-P 40";
    }
	
	
    my $outdir="${basicDir}/L${L}/results";
	
    print MYFILE "mkdir -p $outdir \n";
	
    my $outfile="L${L}N${Ntot}V${V}D20";
	
    my $dirMPSs="${basicDir}/L${L}/mps";
    print MYFILE "mkdir -p $dirMPSs \n";
	
    my $JOBNAME="D${D}.N${Ntot}.V${V}";
	
    my $args="${configfile} -output=${outdir}/${outfile} -mpsdir=${dirMPSs} -jobsdir=${jobsdir} -tol=${tol} -L=${L} -freqSv=${freqSv}";
    $args="${args} -tg0=${tg0} -tg1=${tg1} -te0=${te0} -te1=${te1} -tge0=0. -tge1=0. ";
    $args="${args} -Ug=${Ug} -Ue=${Ue} -V=${V} -Vex=${Vex} ";
    $args="${args} -mu_g0=${mug0} -mu_g1=${mug1} -mu_e0=${mue0} -mu_e1=${mue1} -append=${app} -restart=${restart} -Drestart=${Drestart} ";
    if($PENALTYNTOT>0){
	$args="${args} -Ntot=${Ntot} -penaltyNtot=${PENALTYNTOT} ";
    }
    if($PENALTYSZ>0){
	$args="${args} -Sz=${Sztarget} -penaltySz=${PENALTYSZ} ";
    }
    $args="${args} -D=${D} -incrD=${incrD} -maxD=${MAXD} ";
	
    
    print MYFILE "${MSUB} -N ${JOBNAME} ${par} ${LONGMEM} -- ${EXEC} ${args}\n"; 
    print MYFILE "sleep 1\n";
#    print MYFILE "sleep 3600\n";
}

