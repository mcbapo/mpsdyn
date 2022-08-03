#!/usr/bin/perl
use POSIX;
# Script to prepare a series of jobs using the thThirringP program for thermal states of Thirring model, using projector onto S=0, and targetting (and saving) a certain number of values of beta

# Name of the script file to be created
my $scriptFile="lanzaThThirringRead.sh";

#my $EXEC="./thThirring";
my $EXEC="./thThirringPread";

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

my $backupFreq=50;

#my @Ds=(40,60,80);
#my @Ns=(100,200); 
#my @Ds=(60,80); 
my @Ds=(100,120); 
my @Ns=(200,300); 

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
my $par=" -P 10 ";#" -P 8 ";
my $mem=" -mem 20G ";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

my $cnt=0;
foreach $D (@Ds){
    if ($D>=80){ $par=" -P 10 ";}
    foreach $N (@Ns){
	if ($N>=200){$mem=" -mem 20G ";}
	else{$mem=" -mem 12G ";}
	my $gendir="${mpsdirBas}/N${N}";
	my $mpsdir="${gendir}/MPS";
	print MYFILE "mkdir -p ${gendir}\n";
	print MYFILE "mkdir -p ${mpsdir}\n";

	my $minX=floor(.25*$N);
	my $maxX=floor(.75*$N);
	my $suffix="N${N}_ma${ma}_Delta${Delta}_D${D}";
	my $outputFile="${gendir}/thermalCorr_${suffix}";
	
	my $JOBNAME="corrThP.N${N}.ma${ma}.Delta${Delta}.D${D}";
	my $args=" -L=${N} -ma=${ma} -gtilde2=0 -g2=${Delta} -lambda=${penH} -mu=${mu} -tol=${tol} -Starget=${Starget} -D=${D} -beta=${beta} -delta=${dt} -stepBlock=${rate}  -output=${outputFile} -append=${app} -findGS=${findGS} -mpsdir=${mpsdir} -backupFreq=${backupFreq} -minX=${minX} -maxX=${maxX} ";
	# Check that the job is not yet running!! Add -l h_vmem=${MEM}M 
	#my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
	# Now write the istruction to launch the job to the script
	print MYFILE "${MSUB} ${par} ${mem} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args}\n";
	print MYFILE "sleep 1\n";
	
    }
}


close(MYFILE);

system("chmod a+x $scriptFile ");


