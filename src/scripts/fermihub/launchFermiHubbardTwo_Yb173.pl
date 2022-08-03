#!/usr/bin/perl

# Script to prepare a series of jobs using the twofermhubbard program.
# This use version v0, which conserves independently Ne and Ng

# Name of the script file to be created
my $scriptFile="lanzaYb173.sh"; # first run

my $EXEC="/u/banulsm/NewSVN/mpsdyn/src/bin/fH2gs_v0";

my $basicDir="/ptmp/mpq/banulsm/fermihubbard2/all/cold_atoms";
my $jobsdir="/u/banulsm/jobsToRun";
my $configfile="/u/banulsm/NewSVN/mpsdyn/src/config/twoorbfermihub.conf";

my $suffix="";

my $block=1;
if($block){
    $EXEC="/u/banulsm/NewSVN/mpsdyn/src/bin/fH2gsB";
    $basicDir="/ptmp/mpq/banulsm/fermihubbard2/all/cold_atoms/block";
    $suffix=".B";
}


my $MSUB="msub_slurm";

my $te=0.2591;
my $te0=$te;
my $te1=$te;
my $tg=1.;
my $tg0=$tg;
my $tg1=$tg;
my $Ug=9.238;
my $Ue=18.13;
my $Vex=25.646;
my $mug0=0.;
my $mug1=0.;
my $mue0=0.;
my $mue1=0.;

my $Vdir=37.031;
# Vex=.5*(U-V)

# dimensions for which I run the searchO program
my $D0=20; #, 60, 80, 100, 120, 140, 160);
my $D=20;
my $incrD=40;
my $MAXD=260;

# For warmup!
#$D=2; $MAXD=2; 
#my $warmup=30; # nr rounds
#my $app=0;
my $app=1;
my $warmup=0; 

# frequency to save tmp MPS
my $freqSv=4;

# min number of states
#MINNUMBER="25"
#my $MINNUMBER=2;
#Ds="40 60 80 100 120"
#my $ONLYVEC=" -onlyvec=1 ";

#my $app=1;
my $restart=0;
my $Drestart=200;
my $tryVs=0; my $tryVsstep=1;

my $L=(20); 
my $refine=0;
my $tol=1E-8;
#my $tol=1E-5;

#my $L=(12); $D=220;$MAXD=260;
#my $refine=0;
#my $tol=1E-8;

my $PENALTYNTOT=0; #100;
my $PENALTYNG=100;
my $PENALTYNE=100;
my $PENALTYSZ=100;

#my $Ntot=$L; # quarter filling, as in the 2011 paper
my $Nmax=2*$L; # max for each orbital
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
#my @Ns=(0..$Ntot/2+10);
my @Ngs=(0..$Nmax);
my @Nes=(0..$L);
#my @Ngs=(0..$L); $D=260;$refine=0;
#my @Nes=($L+1..$Nmax);

#my @Ngs=(26);
#my @Nes=(16);my $D=260;my $MAXD=260;
#my $D=180;

#my $D=220;my $MAXD=260;
#my @Ngs=(17..26);
#my @Nes=(6..12);

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

my $TIME="-t 7:00:00";  # "";

foreach my $Ngtarget (@Ngs){
    foreach my $Netarget (@Nes){	
	my $V=$Vdir;
	#my $Netarget=$Ntot-$Ngtarget;
	my $Ntot=$Netarget+$Ngtarget;
	if( ($Ntot % 2)==0){ # o.w. cannot be Sz=0!
	    $Sztarget=0;
	}
	else{
	    $Sztarget=1;
	}
	my $SIZEMEM=int(16*$L*$D*$D*(10+4*$MINNUMBER+2)*1.1*1E-6)+200;
    
	my $LONGMEM="-mem ${SIZEMEM}M";
	if ($SIZEMEM<1600) {
	    $LONGMEM="";
	}
	
	my $par="-P 4";
	if ($D>=220){
	    $par="-P 8";
	}
	if ($D>=280){
	    $par="-P 10";
	}
	
	
	my $outdir="${basicDir}/N${L}/results";
	#my $lowEdir="${basicDir}/N${L}/multiImbalance";
    
	print MYFILE "mkdir -p $outdir \n";
	#print MYFILE "mkdir -p $lowEdir \n";
	
	#my $outfile="N${L}Ng${Ngtarget}V${V}D${D}";
	my $outfile="N${Ntot}Ng${Ngtarget}V${V}D20";
	
	my $dirMPSs="${basicDir}/N${L}/mpsSym";
	if($block){
	    $dirMPSs="${basicDir}/N${L}/mps";
	}
	print MYFILE "mkdir -p $dirMPSs \n";
	
	my $JOBNAME="D${D}.L${L}.Ng${Ngtarget}.Ne${Netarget}.V${V}${suffix}";
	
	my $args="${configfile} -output=${outdir}/${outfile} -mpsdir=${dirMPSs} -jobsdir=${jobsdir} -tol=${tol} -L=${L} -freqSv=${freqSv}";
	$args="${args} -tg0=${tg0} -tg1=${tg1} -te0=${te0} -te1=${te1} -tge0=0. -tge1=0. ";
	$args="${args} -Ug=${Ug} -Ue=${Ue} -V=${V} -Vex=${Vex} ";
	$args="${args} -mu_g0=${mug0} -mu_g1=${mug1} -mu_e0=${mue0} -mu_e1=${mue1} -append=${app} "; #-restart=${restart} -Drestart=${Drestart} -refine=${refine} -tryVs=${tryVs} -tryVsstep=${tryVsstep} ";

	$args="${args} -warmup=${warmup} ";
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
	$args="${args} -D=${D} -incrD=${incrD} -maxD=${MAXD} ";
	
	my $alreadyRunning=isRunning($JOBNAME);
	if(0&&$alreadyRunning){
	    print "Job $JOBNAME already running \n";
	}
	else{
	    # Check for larger D, just in case
	    my $curD=$D+$incrD;
	    my $runningLargerD=0;
	    my $LATERJOBNAME="";
	    while($curD<$MAXD){
		$curD=$curD+$incrD;				
		$LATERJOBNAME="D${curD}.L${L}.Ng${Ngtarget}.Ne${Netarget}.V${V}${suffix}";
		$runningLargerD=isRunning($LATERJOBNAME);
	    }
	    if($runningLargerD){
		print "Job $LATERJOBNAME already running \n";
	    }
	    else{ # TODO Check for the older version
		print MYFILE "${MSUB} -N ${JOBNAME} ${par} ${LONGMEM} ${TIME} -- ${EXEC} ${args}\n"; 
		print MYFILE "sleep 1\n";
		if($Netarget==$L){
		  #  print MYFILE "sleep 3599\n";
		}
	    }
	} # Running or not
    }  # loop Netarget
    
} # loop Ngtarget
    
