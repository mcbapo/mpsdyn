#!/usr/bin/perl

# Script to prepare a series of jobs using the minH2 program

# Name of the script file to be created
my $scriptFile="lanzaVarIsingX.sh";

my $EXEC="./minH2x";
my $par="-P 8"; # use parallel environment

my $configFile="../config/randomMPO.conf";

my $basicDir="/ptmp/mpq/banulsm/variance/Ising/space";
# Directory for results
my $resultsDir="${basicDir}/results";
my $mpsDir="${basicDir}/MPS";

# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 ");
my @idStrings=("J1_gm105_h05","J1_g0905_h0809");

#my @hamilPars=(" 1 0.8 1"," 1 -.5 .5");
#my @idStrings=("J1_g08_h1","J1_gm05_h05");

#my @hamilPars=(" 1 0.8 0");
#my @idStrings=("J1_g08_h0"); # an integrable one

my $D1=20;
my $D2=160;
my $incrD=10;
my $app=0;
my $noise=1E-3;

# lengths of the system
my @Ns=(20,40,80,100);

my $MSUB="msub_modified_tqo097";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

print MYFILE "mkdir -p ${resultsDir} \n";
print MYFILE "mkdir -p ${mpsDir} \n";

#$D1=1;$D2=1; # try just products

my $cnt=0;
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
    print MYFILE "mkdir -p ${resultsDir}/${strL} \n";
    print MYFILE "mkdir -p ${mpsDir}/${strL} \n";
    foreach $N (@Ns){
	my $suffix="N${N}";
	my $outfile="${resultsDir}/${strL}/data_${suffix}";
	my $outMPS="${mpsDir}/${strL}/mps_${suffix}";
	# Name of the output file
	my $JOBNAME="varX.N${N}.${cnt}";
	# Check that the job is not yet running!!
	#my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
	my $args=" -L=${N} ${paramSet} -D1=${D1} -D2=${D2} -incrD=${incrD} -noise=${noise} -output=${outfile} -append=${app} -outputMPS=${outMPS} ";
	print MYFILE "msub_modified_tqo097 ${par} -N ${JOBNAME} -raw ${EXEC} ${configFile} ${args} \n";
	print MYFILE "sleep 1\n";
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


