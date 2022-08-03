#!/usr/bin/perl

# Script to prepare a series of jobs using the twofermhubbard program.
# This use version v0, which conserves independently Ne and Ng

# Name of the script file to be created
my $scriptFile="lanzaFH2.sh"; # first run

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
my $U=$Ug;
my $mug0=0.;
my $mug1=0.;
my $mue0=0.;
my $mue1=0.;

my @Vs=(0..19);
# Vex=.5*(U-V)

# dimensions for which I run the searchO program
my $D0=20; #, 60, 80, 100, 120, 140, 160);
my $D=20;
my $incrD=40;
my $MAXD=200;

#my $D=250;
#my $incrD=50;
#my $MAXD=300;
my $tol=1E-3;

# frequency to save tmp MPS
my $freqSv=4;

# min number of states
#MINNUMBER="25"
my $MINNUMBER=2;
#Ds="40 60 80 100 120"
my $ONLYVEC=" -onlyvec=1 ";

my $app=1;
my $refine=0;
my $restart=0;
my $Drestart=200;
my $tryVs=1; my $tryVsstep=1;

my $L=(50);my $refine=1;
#my $L=(20);
my $tol=1E-5;
#my @Vs=(0..19);my $tryVsstep=1;
my @Vs=(0..13);my $tryVsstep=1;
my @Vs=(5.2,5.4,5.6,5.8,6.2,6.4,6.6,6.8,7.2,7.4,7.6,7.8,8.2);my $tryVsstep=.2;my $refine=0;
#my @Vs=(8.4,8.6,8.8,9.2,9.4,9.6,9.8);my $tryVsstep=.2;my $refine=0;
my @Vs=(7.2,7.4,7.6,7.8,8.2,8.4,8.6,8.8,9.2,9.4,9.6,9.8);my $tryVsstep=.2;my $refine=0;
#my @Vs=(1..4);my $tryVsstep=1;my $refine=1;
#my @Vs=(5..10);my $tryVsstep=.2;my $refine=1;
#my @Vs=(11..13);my $tryVsstep=1;my $refine=1;
my $D=140; my $MAXD=100;
#my @Vs=(9.2,9.4,9.6,9.8);my $tryVsstep=.2;
my $PENALTYNTOT=0; #100;
my $PENALTYNG=100;
my $PENALTYNE=100;
my $PENALTYSZ=100;

my $Ntot=$L; # quarter filling, as in the 2011 paper
my $Sztarget=0;

# I need to check each imbalance. Now what I do is using (for each set
# of parameters) a file where to write the lowest energy found and, if
# the run with a certain imbalance finds energy that is higher than
# that by over 20%, it stops. It is nt 100% safe strategy, though
# (race condition for lock file, or differences at different Ds that
# may be larger than the cutoff but reversed later in D)
# (Q: a smarter strategy, e.g. some derivative of energy? to know when I'm done?)
# OBSOLETE!!! Not reliable

my $MULTIIMB=0;
my @Ns=(0..$Ntot/2+10);
#my @Ns=(0..$Ntot);

###### refine
#my $scriptFile="lanzaFH2_v0_refine.sh"; # to refine
#my @Vs=(0..5);
##my @Vs=(5..13);
#my $refine=1; # even if it exists, run again
#my $restart=1; # require that it exists, actually
#my $tol=1E-5;
#my $Drestart=180;
#my $D=180;
#my $incrD=40;
#my $maxD=200;


###### refine


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
	if ($D>=180){
	    $par="-P 10";
	}
	
	
	my $outdir="${basicDir}/N${L}/results";
	#my $lowEdir="${basicDir}/N${L}/multiImbalance";
	
	print MYFILE "#mkdir -p $outdir \n";
	#print MYFILE "mkdir -p $lowEdir \n";
	
	#my $outfile="N${L}Ng${Ngtarget}V${V}D${D}";
	my $outfile="N${L}Ng${Ngtarget}V${V}D20";
	
	my $dirMPSs="${basicDir}/N${L}/mps";
	print MYFILE "#mkdir -p $dirMPSs \n";
	
	my $JOBNAME="D${D}.Ng${Ngtarget}.V${V}";
	
	my $args="${configfile} -output=${outdir}/${outfile} -mpsdir=${dirMPSs} -jobsdir=${jobsdir} -tol=${tol} -L=${L} -freqSv=${freqSv}";
	$args="${args} -tg0=${tg0} -tg1=${tg1} -te0=${te0} -te1=${te1} -tge0=0. -tge1=0. ";
	$args="${args} -Ug=${Ug} -Ue=${Ue} -V=${V} -Vex=${Vex} ";
	$args="${args} -mu_g0=${mug0} -mu_g1=${mug1} -mu_e0=${mue0} -mu_e1=${mue1} -append=${app} -restart=${restart} -Drestart=${Drestart} -refine=${refine} -tryVs=${tryVs} -tryVsstep=${tryVsstep} ";
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
	#if($MULTIIMB>0){
	#    $args="${args} -multiImbalance=1 -lowEdir=${lowEdir} ";
	#}
	$args="${args} -D=${D} -incrD=${incrD} -maxD=${MAXD} ";
	
	#-L=${L} -mg=${mg} -x=${x} -D=${D} -outputdir=${outdir} -penalty=${PENALTY} -mpsdir=${dirMPSs} -initmpsdir=${dirMPSs} -initD=${INITD} -tol=1E-8 -minnumber=${MINNUMBER}  ${ONLYVEC} ";
	my $alreadyRunning=isRunning($JOBNAME);
	if(! $alreadyRunning){
	    print MYFILE "${MSUB} -N ${JOBNAME} ${par} ${LONGMEM} -- ${EXEC} ${args}\n"; 
	    print MYFILE "sleep 1\n";
	}
	else{
	    print "Job $JOBNAME already running \n";
	}
    }  # loop Ntarget
#    print MYFILE "sleep 7200\n";
} # loop V

