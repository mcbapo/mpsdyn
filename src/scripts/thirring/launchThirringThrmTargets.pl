#!/usr/bin/perl
use POSIX;
# Script to prepare a series of jobs using the thThirringP program for thermal states of Thirring model, using projector onto S=0, and targetting (and saving) a certain number of values of beta

# Name of the script file to be created
my $scriptFile="lanzaThThirring.sh";

#my $EXEC="./thThirring";
my $EXEC="./thThirringP";

my $MSUB="msub_slurm";


#my $mpsdirBas="/ptmp/mpq/banulsm/thirring/thermal";
my $mpsdirBas="/ptmp/mpq/banulsm/thirring/thermalProjTarget";
#my $jobsdir="/u/banulsm/jobsToRun";

# Hamiltonian parameters
# In this case a particular one, for which we target some betas
#my $Delta=-0.9;
#my $ma=0.1;
#my $configFile="../config/thermalThirring_Delta-0.9.conf";

my $Delta=0.5;
my $ma=0.5;
my $configFile="../config/thermalThirring_Delta0.5_ma0.5.conf";

# Does not work like this: needs to be in config file!!
#my $targetBs=""

my $backupFreq=4;

#my @Ds=(40,60,80);
#my @Ns=(100,200); 
my @Ds=(120,100); #60,80); 
my @Ns=(300,200); 

my $tol=1E-6; # tolerance for convergence
my $penH=20;
my $nrStates=0;
my $mu=0; # No chemical potential
my $Starget=0; # sector of Sz=0

my $beta=3;
my $dt=0.01; # Trotter step
my $rate=2; # frequency to record data
my $app=1; # append
my $findGS=0; # not looking for GS
my $par=" -P 20 ";#" -P 8 ";
my $mem=" -mem 12G ";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

my $cnt=0;
foreach $D (@Ds){
    if ($D>=60){ $par=" -P 40 ";}
    foreach $N (@Ns){
	if ($N>=200){$mem=" -mem 20G ";}
	else{$mem=" -mem 12G ";}
	my $gendir="${mpsdirBas}/N${N}";
	my $mpsdir="${gendir}/MPS";
	print MYFILE "mkdir -p ${gendir}\n";
	print MYFILE "mkdir -p ${mpsdir}\n";
	
	my $suffix="N${N}_ma${ma}_Delta${Delta}_D${D}";
	my $outputFile="${gendir}/thermalEvol_${suffix}";
	
	my $JOBNAME="thThP.N${N}.ma${ma}.Delta${Delta}.D${D}";
	my $args=" -L=${N} -ma=${ma} -gtilde2=0 -g2=${Delta} -lambda=${penH} -mu=${mu} -tol=${tol} -Starget=${Starget} -D=${D} -beta=${beta} -delta=${dt} -stepBlock=${rate}  -output=${outputFile} -append=${app} -findGS=${findGS} -mpsdir=${mpsdir} -backupFreq=${backupFreq} "; #-targets=\"${targetBs}\" ";
	# Check that the job is not yet running!! Add -l h_vmem=${MEM}M 
	#my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
	# Now write the istruction to launch the job to the script
	print MYFILE "${MSUB} ${par} ${mem} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args}\n";
	print MYFILE "sleep 1\n";
	
    }
}


close(MYFILE);

system("chmod a+x $scriptFile ");


