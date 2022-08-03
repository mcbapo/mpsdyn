#!/usr/bin/perl

# Script to prepare a series of jobs using the minH2 program

# Name of the script file to be created
my $scriptFile="lanzaVarCos.sh";

my $EXEC="./minH2cL";

my $configFile="../config/randomMPO.conf";

my $basicDir="/ptmp/mpq/banulsm/variance/Ising/Cos";

# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 ");
my @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0"); # and an integrable one

#my @hamilPars=(" 1 0.8 1"," 1 -.5 .5");
#my @idStrings=("J1_g08_h1","J1_gm05_h05");

my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "); #," -J=1 -g=0.905 -h=0.809 ");
my @idStrings=("J1_gm105_h05"); #,"J1_g0905_h0809");
#my @hamilPars=(" 1 0.8 0");
#my @idStrings=("J1_g08_h0"); # an integrable one


my @Ds=(20,40,60,80,100,120,140,200);
my @pars=("-P 4","-P 4","-P 4","-P 4","-P 8","-P 8","-P 10","-P 12"); # use parallel environment

my $app=0;
my $initStaggered=0; # start with staggered X+- (o.w. X+)
my $initRandom=0;
my $randSeed=117; # to start always the same
my $initY=1;
my $avrg=1;
my $Lcmax=6;
my $Lb=6;

if(${initStaggered}){
    $basicDir="${basicDir}/Xst/";
}
else{
    if($initRandom){
#	$basicDir="${basicDir}/rnd/";
	$basicDir="${basicDir}/rndTI/";
    }
    else{
	if($initY){
	    $basicDir="${basicDir}/Yplus/";
	}
	else{
	    $basicDir="${basicDir}/Xplus/";
	}
    }
}
# Directory for results
my $resultsDir="${basicDir}/results";
#my $mpsDir="${basicDir}/MPS";

my $dt=0.01;
my $var2=0.01; # target delta^2
my $rate=100;
my $saveRate=500;

# lengths of the system
my @Ns=(20,40,60,80,100);

my $MSUB="msub_slurm";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

print MYFILE "mkdir -p ${resultsDir} \n";
#print MYFILE "mkdir -p ${mpsDir} \n";


my $cnt=0;
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
    print MYFILE "mkdir -p ${resultsDir}/${strL}/results \n";
    print MYFILE "mkdir -p ${resultsDir}/${strL}/MPS \n";
    my $cntD=0;
    foreach $D (@Ds){
	my $par=$pars[$cntD];
	foreach $N (@Ns){
	    my $suffix="dX1N_N${N}D${D}";
	    my $outfile="${resultsDir}/${strL}/results/data_${suffix}";
	    my $mpsdir="${resultsDir}/${strL}/MPS";
	    my $mpsfile="mps_${suffix}";
	    my $deltaX=1./${N}; #.1/${N};
	    my $Nsteps=(1./${var2})*${N}*${N}*2;
	    # Name of the output file
	    my $JOBNAME="varC.dX1_N.N${N}.D${D}.${cnt}";
	    # Check that the job is not yet running!!
	    #my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
	    my $args=" -L=${N} ${paramSet} -D=${D} -output=${outfile} -mpsdir=${mpsdir} -mpsfile=${mpsfile} -x=${deltaX} -dt=${dt} -M=${Nsteps} -rate=${rate} -saveFreq=${saveRate} -staggered=${initStaggered} -checkRDM=1 -Lcmax=$Lcmax -average=$avrg -Lb=$Lb -random=${initRandom} -seed=$randSeed -initY=${initY} ";
	    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
	    print MYFILE "sleep 1\n";
	}
	$cntD=$cntD+1;
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


