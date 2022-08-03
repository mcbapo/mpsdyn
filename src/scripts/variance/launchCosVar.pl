#!/usr/bin/perl

# Script to prepare a series of jobs using the lowH2 program

# Name of the script file to be created
my $scriptFile="lanzaCosVar.sh";

my $EXEC="./cosH2";

my $configFile="../config/randomMPO.conf";

my $basicDir="/ptmp/mpq/banulsm/variance/Ising/CosVar";
# Directory for results
my $resultsDir="${basicDir}/results";
#my $mpsDir="${basicDir}/MPS";

# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 ");
my @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0"); # and an integrable one

my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 ");
my @idStrings=("J1_gm105_h05");

my @Ds=(40,60,80,100,120,140,160,180,200,240,300);
my $app=1;
my $tolC=1E-6; # tolreance for convergence of the variance
my $stepTol=10; # if for these nr of steps does not chg more than tolC, stop
my $saveFreq=5;

my $expon=0;
my $deltaX=0.1;
my $dt=0.01;
my $M=50;

# lengths of the system
my @Ns=(20,40,60,80,100);
#my @Ns=(140,200,260,400);
#my @Ns=(40,60,80,100,120,140,160,180,200,220,240,260,280,300,400);
my $cont=" -continue=1 -append=1";

my $MSUB="msub_modified_tqo097";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

print MYFILE "mkdir -p ${resultsDir} \n";
#print MYFILE "mkdir -p ${mpsDir} \n";


my $cnt=0;
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
    print MYFILE "mkdir -p ${resultsDir}/${strL}/results \n";
    print MYFILE "mkdir -p ${resultsDir}/${strL}/MPS \n";
    
    foreach $N (@Ns){
	foreach $D (@Ds){       
	    if($D>=120){
		$par=" -P 8 ";
	    }
	    else{
		if($D>=80){
		    $par=" -P 4 ";
		}
		else{
		    $par="";
		}
	    }
	    my $suffix="N${N}D${D}_x${deltaX}";
	    if($expon){$suffix="${suffix}_exp";}
	    my $outfile="${resultsDir}/${strL}/results/data_${suffix}";
	    my $mpsdir="${resultsDir}/${strL}/MPS";
	    my $mpsfile="mps_${suffix}";
	    my $JOBNAME="cosH2.N${N}.D${D}";
	    # Decide if I am still (already) running
	    my $isRunning=`qstat -j $JOBNAME |grep job_number`;
	    if($isRunning =~ /(\d+)/){
		$jobNr=$1;
		print MYFILE "### Job already/still running $jobNr for $JOBNAME ## \n";
	    }
	    else{
		#my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
		my $args=" -L=${N} ${paramSet} -D=${D} -M=${M} -x=${deltaX} -dt=${dt} -exponential=${expon} -output=${outfile} -mpsdir=${mpsdir} -mpsfile=${mpsfile} -convTol=${tolC} -stepsConvTol=${stepTol} -saveFreq=${saveFreq} -append=${app} ";
		print MYFILE "msub_modified_tqo097 ${par} -N ${JOBNAME} -raw ${EXEC} ${configFile} ${args} \n";
		print MYFILE "sleep 1\n";

	    }
	}
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


