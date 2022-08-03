#!/usr/bin/perl

# Script to prepare a series of jobs using the lowH2post program
# to process the results of lowH2

# Name of the script file to be created
my $scriptFile="lanzaVarThr2.sh";

my $EXEC="./lowH2post";

my $configFile="../config/randomMPO.conf";

my $basicDir="/ptmp/mpq/banulsm/variance/Ising/VarThr";
# Directory for results
my $resultsDir="${basicDir}/results";
#my $mpsDir="${basicDir}/MPS";

# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 ");
my @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0"); # and an integrable one

my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 ");
my @idStrings=("J1_gm105_h05");

my $D1=2;
my $D2=2000;
my $incrD=2;
my $tolC=1E-2; # tolreance for convergence
my @deltaTh=(1.,0.1,0.01);
my @deltaTh=(0.01); # this should contain all!

# lengths of the system
#my @Ns=(20,40,60,80,100);
#my @Ns=(140,200,260,400);
my @Ns=(20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,400,500,600,700,800,1000,1500);
#my @Ns=(400,500,600,700,800,1000,1500);
#my $cont=" -continue=1 -append=1";
my $app=" -append=0 "; # Redo all

my $MSUB="msub_slurm";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

print MYFILE "mkdir -p ${resultsDir} \n";
#print MYFILE "mkdir -p ${mpsDir} \n";


my $cnt=0;
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
#    print MYFILE "mkdir -p ${resultsDir}/${strL}/results \n";
#    print MYFILE "mkdir -p ${resultsDir}/${strL}/MPS \n";

    foreach $varT (@deltaTh){
	foreach $N (@Ns){

#	    if($N>=120){
	    $par=" -P 8 ";
#	    }
#	    else{
#		if($N>=80){
#		    $par=" -P 4 ";
#		}
#		else{
#		    $par="";
#		}
#	    }
	    my $suffix="N${N}_thr${varT}";
	    my $outfile="${resultsDir}/${strL}/results/distanceAvrg_${suffix}";
	    my $mpsdir="${resultsDir}/${strL}/MPS";
	    my $mpsfile="mps_${suffix}";
	    my $JOBNAME="lowH2post.N${N}.d${varT}.${cnt}";
	    # Decide if I am still (already) running
	    # Check that the job is not yet running!!
	    my $isRunning=`squeue -h -n $JOBNAME `;
	    if($isRunning =~ /(\d+)/){ # a job nr exists
		$jobNr=$1;
		print MYFILE "### Job already/still running $jobNr for $JOBNAME ## \n";
	    }
	    else{
		# Check how far the first attempt went
		my $lastLine=`tail -n 1 ${outfile} `;
		my @entries=split(/\s+/,$lastLine);	
		#      print "\t which has entries @entries \n";
		my $lastLc=@entries[7];
		my $lastD=@entries[4];
		printf "In file ${outfile} last D used ${lastD}, up to Lc=${lastLc}\n";
		if($lastLc==2){
		    $lastD+=$incrD;
		}		
		if($lastD==""){$lastD=$D1 ;}
		#my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
		$lastD=$D1;
		my $args=" -L=${N} ${paramSet} -D1=${lastD} -D2=${D2} -incrD=${incrD} -output=${outfile} -mpsdir=${mpsdir} -mpsfile=${mpsfile} ${app}";
		print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
		print MYFILE "sleep 1\n";
	    }
	}
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


