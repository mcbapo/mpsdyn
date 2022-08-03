#!/usr/bin/perl

# Script to prepare a series of jobs using the twofermihubbardGSmod program
# this one launches, ony by one, all the values of D, using as starting point the results for larger tge

# Name of the script file to be created
my $scriptFile="lanzaFH2.sh";

my $EXEC="/u/banulsm/NewSVN/mpsdyn/src/bin/fH2gs";

my $basicDir="/ptmp/mpq/banulsm/fermihubbard2/tge";
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

#my @Vs=(0..19);

my $inittge=0.1;
my $tge=0.01;
my $Din=60; 
my $incrD=10;
my $MAXD=60;
my $tol=1E-5;
my @Vs=(0..10);
my $scriptFile="lanzaFH2_lowV_${MAXD}.sh";


#my $inittge=0.0001;
#my $tge=0.; #0001;
#my $Din=300; 
#my $incrD=20;
#my $MAXD=300;
#my $tol=1E-12;
#my @Vs=(12..19);
#my $scriptFile="lanzaFH2_highV_${MAXD}.sh";

#my $inittge=0.001;
#my $tge=0.0001;
#my $Din=300; 
#my $incrD=20;
#my $MAXD=300;
#my $tol=1E-5;
#my @Vs=(11);
#my $scriptFile="lanzaFH2_recu.sh";




my $tge0=$tge;
my $tge1=$tge;
my $inittge0=$inittge;
my $inittge1=$inittge;

# dimensions for which I run the searchO program
my @Ds;
for (my $d=$Din;$d<=$MAXD; $d+=$incrD){
    push @Ds, $d ;
}

print "Ds is now ", @Ds ;

# ferquency to save tmp MPS
my $freqSv=4;

# min number of states
#MINNUMBER="25"
my $MINNUMBER=2;
#Ds="40 60 80 100 120"
my $ONLYVEC=" -onlyvec=1 ";

my $app=1;

my $L=(50);

my $PENALTYNTOT=100;
my $PENALTYSZ=100;

my $Ntarget=$L; # quarter filling, as in the 2011 paper
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

foreach my $V (@Vs){
    my $Vex=.5*($U-1.*$V);
    foreach my $D (@Ds){
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
	    $par="-P 20";
	}
	
	
	my $outdir="${basicDir}/N${L}/results";
	my $lowEdir="${basicDir}/N${L}/multiImbalance";
	
	print MYFILE "mkdir -p $outdir \n";
	print MYFILE "mkdir -p $lowEdir \n";
	
	my $outfile="N${L}Ng${Ngtarget}V${V}";
	
	my $dirMPSs="${basicDir}/N${L}/mps";
	print MYFILE "mkdir -p $dirMPSs \n";
	
	my $JOBNAME="D${D}.tge${tge}.V${V}";
	
	my $args="${configfile} -output=${outdir}/${outfile} -mpsdir=${dirMPSs} -jobsdir=${jobsdir} -tol=${tol} -L=${L} -freqSv=${freqSv}";
	$args="${args} -tg0=${tg0} -tg1=${tg1} -te0=${te0} -te1=${te1} -tge0=${tge0} -tge1=${tge1} ";
	$args="${args} -Ug=${Ug} -Ue=${Ue} -V=${V} -Vex=${Vex} ";
	$args="${args} -mu_g0=${mug0} -mu_g1=${mug1} -mu_e0=${mue0} -mu_e1=${mue1} ";
	if($PENALTYNTOT>0){
	    $args="${args} -Ntot=${Ntarget} -penaltyNtot=${PENALTYNTOT} ";
	}
	if($PENALTYSZ>0){
	    $args="${args} -Sz=${Sztarget} -penaltySz=${PENALTYSZ} ";
	}
	$args="${args} -D=${D} ";#-incrD=${incrD} -maxD=${MAXD} ";
	$args="${args} -relaunch=0 -initWithTgeResults=1 -initTge0=${inittge0} -initTge1=${inittge1} -append=${app} ";
	
	#-L=${L} -mg=${mg} -x=${x} -D=${D} -outputdir=${outdir} -penalty=${PENALTY} -mpsdir=${dirMPSs} -initmpsdir=${dirMPSs} -initD=${INITD} -tol=1E-8 -minnumber=${MINNUMBER}  ${ONLYVEC} ";

	print MYFILE "${MSUB} -N ${JOBNAME} ${par} ${LONGMEM} -- ${EXEC} ${args}\n"; 
	print MYFILE "sleep 1\n";
    }
}


