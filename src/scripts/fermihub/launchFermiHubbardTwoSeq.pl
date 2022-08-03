#!/usr/bin/perl

# Script to prepare a series of jobs using the twofermhubbard program

# Name of the script file to be created
my $scriptFile="refineFH2.sh";

my $EXEC="/u/banulsm/NewSVN/mpsdyn/src/bin/fH2gs_v0";

my $basicDir="/ptmp/mpq/banulsm/fermihubbard2/all";
my $jobsdir="/u/banulsm/jobsToRun";
my $configfile="/u/banulsm/NewSVN/mpsdyn/src/config/twoorbfermihub.conf";


my $MSUB="msub_slurm";

my $te=1.;
my $te0=$te;
my $te1=$te;
my $tg=.5*$te;
my $tg0=$tg;
my $tg1=$tg;
my $Ug=10.;
my $Ue=$Ug;
my $mug0=0.;
my $mug1=0.;
my $mue0=0.;
my $mue1=0.;

my @Vs=(0..19);
# Vex=.5*(U-V)

# dimensions for which I run the searchO program
my $D=200; #, 60, 80, 100, 120, 140, 160);
my $incrD=20;
my $MAXD=200;

# try increasing precision (could try to increase D first)
my $tol=1E-5; # first run was with 1E-3
my $refine=1;

# frequency to save tmp MPS
my $freqSv=4;

# min number of states
#MINNUMBER="25"
my $MINNUMBER=2;
#Ds="40 60 80 100 120"
my $ONLYVEC=" -onlyvec=1 ";

my $app=1;

my $L=(50);

my $PENALTYNTOT=0; #100;
my $PENALTYNG=100;
my $PENALTYNE=100;
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
my $MULTIIMB=0;
my @Ns=(0..$Ntot);


foreach my $V (@Vs){
    my $Vex=.5*($U-1.*$V);
    foreach my $Ngtarget (@Ns){
	my $Netarget=$Ntot-$Ngtarget;
	my $SIZEMEM=int(16*$L*$D*$D*(10+4*$MINNUMBER+2)*1.1*1E-6)+200;
	
	my $LONGMEM="-mem ${SIZEMEM}M";
	if ($SIZEMEM<3000) {
	    $LONGMEM="";
	}
	
	my $par="-P 4";
	if ($D>=40){
	    $par="-P 8";
	}
	if ($D>=80){
	    $par="-P 10";
	}
	
	
	my $outdir="${basicDir}/N${L}/results";
	my $lowEdir="${basicDir}/N${L}/multiImbalance";
	
	print MYFILE "mkdir -p $outdir \n";
	print MYFILE "mkdir -p $lowEdir \n";
	
	my $outfile="N${L}Ng${Ngtarget}V${V}D20";
	
	my $dirMPSs="${basicDir}/N${L}/mps";
	print MYFILE "mkdir -p $dirMPSs \n";
	
	my $JOBNAME="D${D}.Ng${Ngtarget}.V${V}";
	
	my $args="${configfile} -output=${outdir}/${outfile} -mpsdir=${dirMPSs} -jobsdir=${jobsdir} -tol=${tol} -L=${L} -freqSv=${freqSv}";
	$args="${args} -tg0=${tg0} -tg1=${tg1} -te0=${te0} -te1=${te1} -tge0=0. -tge1=0. ";
	$args="${args} -Ug=${Ug} -Ue=${Ue} -V=${V} -Vex=${Vex} ";
	$args="${args} -mu_g0=${mug0} -mu_g1=${mug1} -mu_e0=${mue0} -mu_e1=${mue1} ";
	if($PENALTYNTOT>0){
	    $args="${args} -Ntot=${Ntarget} -penaltyNtot=${PENALTYNTOT} ";
	}
	if($PENALTYNG>0){
	    $args="${args} -Ng=${Ngtarget} -penaltyNg=${PENALTYNG} ";
	}
	if($PENALTYNE>0){
	    $args="${args} -Ne=${Netarget} -penaltyNe=${PENALTYNE} ";
	}
	if($PENALTYSZ>0){
	    $args="${args} -Sz=${Sztarget} -penaltySz=${PENALTYSZ} ";
	}
	if($MULTIIMB>0){
	    $args="${args} -multiImbalance=1 -lowEdir=${lowEdir} ";
	}
	$args="${args} -D=${D} -incrD=${incrD} -maxD=${MAXD} -refine=${refine} -append=${app} ";
	
	#-L=${L} -mg=${mg} -x=${x} -D=${D} -outputdir=${outdir} -penalty=${PENALTY} -mpsdir=${dirMPSs} -initmpsdir=${dirMPSs} -initD=${INITD} -tol=1E-8 -minnumber=${MINNUMBER}  ${ONLYVEC} ";

	print MYFILE "${MSUB} -N ${JOBNAME} ${par} ${LONGMEM} -- ${EXEC} ${args}\n"; 
	print MYFILE "sleep 1\n";
    }
}

