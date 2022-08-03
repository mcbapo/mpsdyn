#!/usr/bin/perl
use POSIX;
# Script to prepare a series of jobs using the lowH2 program

# Name of the script file to be created
my $scriptFile="lanzaLowVar.sh";

my $EXEC="./lowVar";

my $MSUB="msub_slurm";

my $configFile="../config/randomMPO.conf";

my $mpsdirBas="/ptmp/mpq/banulsm/variance/Ising/variational";
my $jobsdir="/afs/ipp/u/banulsm/jobsVarH2/";

# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 "," -J=1 -g=1 -h=0 ");
my @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0","J1_g1_h0"); # and an integrable one

my $tol=1E-5; # tolerance for convergence

#my @Ns=(12,20,30,40);
my @Ns=(20);
#my @Ns=(20,40,50,80,100,200,300);
#my $Dmin=20; 
#my $Dmax=40;
my $Dmin=20; 
my $Dmax=300;
my $stepD=20; 
my $E0=0; # only want one energy
my $penH=20;
my $nrStates=0;

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

my $cnt=0;
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
    foreach $N (@Ns){
	my $gendir="${mpsdirBas}/${strL}/N${N}";
	my $mpsdir="${gendir}/MPS";
	print MYFILE "mkdir -p ${gendir}\n";
	print MYFILE "mkdir -p ${mpsdir}\n";
	
	my $suffix="N${N}_${strL}_E${strE0}_D${Dn}";
	my $JOBNAME="lowV0.${cnt}.N${N}.D${Dmin}";
	my $args=" -L=${N} ${paramSet} -D=${Dmin} -mpsdir=${mpsdir} -tol=${tol} -penH=${penH} -Dmin=${Dmin} -Dmax=${Dmax} -nrStates=${nrStates} -jobsdir=${jobsdir} -tol=${tol} -E=${E0} ";
	# Check that the job is not yet running!! Add -l h_vmem=${MEM}M 
	#my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
	# Now write the istruction to launch the job to the script
	print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} -greedy=0 \n";
	print MYFILE "sleep 1\n";
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


