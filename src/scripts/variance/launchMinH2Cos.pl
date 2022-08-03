#!/usr/bin/perl

# Script to prepare a series of jobs using the minH2 program

# Name of the script file to be created
my $scriptFile="lanzaVarCos.sh";

my $EXEC="./minH2c";

my $configFile="../config/randomMPO.conf";

my $basicDir="/ptmp/mpq/banulsm/variance/Ising/CosShort";
# Directory for results
my $resultsDir="${basicDir}/results";
#my $mpsDir="${basicDir}/MPS";

# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 ");
my @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0"); # and an integrable one

#my @hamilPars=(" 1 0.8 1"," 1 -.5 .5");
#my @idStrings=("J1_g08_h1","J1_gm05_h05");

#my @hamilPars=(" 1 0.8 0");
#my @idStrings=("J1_g08_h0"); # an integrable one

my $initSt="+1"; #could be +2 (Y) +3 (Z) or negative or staggered (4,5,6)

my @Ds=(20,40,60,80,100,120,140,200);
my @pars=("","","-P 2","-P 4","-P 8","-P 8","-P 8","-P 10"); # use parallel environment
my $app=0;

my $dt=0.01;
my $var2=0.05; # target delta^2
my $rate=5;
my $eps=.1; # argument in cos is H*eps/N
my $a=2; # a*sqrt(M) terms kept

# lengths of the system
my @Ns=(20,40,60,80,100);
my @Ns=(20);
my $deltaX=0.1;

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
    my $cntD=0;
    foreach $D (@Ds){
	my $par=$pars[$cntD];
	foreach $N (@Ns){
	    #my $deltaX=1.*${eps}/${N};
	    my $Nsteps=(1./${var2})*${N}*${N}; # M
	    my $suffix="N${N}D${D}x${deltaX}";
	    my $outfile="${resultsDir}/${strL}/results/data_${suffix}";
	    my $mpsdir="${resultsDir}/${strL}/MPS";
	    my $mpsfile="mps_${suffix}";
	    # Name of the output file
	    my $JOBNAME="varCs.N${N}.D${D}.${cnt}";
	    # Check that the job is not yet running!!
	    #my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
	    my $args=" -L=${N} ${paramSet} -D=${D} -output=${outfile} -x=${deltaX} -dt=${dt} -M=${Nsteps} -a=${a} -rate=${rate} -initSt=${initSt} ";
	    print MYFILE "msub_modified_tqo097 ${par} -N ${JOBNAME} -raw ${EXEC} ${configFile} ${args} \n";
	    print MYFILE "sleep 1\n";
	}
	$cntD=$cntD+1;
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


