#!/usr/bin/perl
use POSIX;
# Script to prepare a series of jobs using the lowH2 program

# Name of the script file to be created
my $scriptFile="lanzaThThirring.sh";

#my $EXEC="./thThirring";
my $EXEC="./thThirringP";

my $MSUB="msub_slurm";

my $configFile="../config/thermalThirring.conf";
#my $configFile="../config/thermalThirring_Delta-0.9.conf";

#my $mpsdirBas="/ptmp/mpq/banulsm/thirring/thermal";
my $mpsdirBas="/ptmp/mpq/banulsm/thirring/thermalProj";
#my $jobsdir="/u/banulsm/jobsToRun";

# Hamiltonian parameters
#my @Deltas=(-0.9,0.9);
#my @Ms=(0.,0.1,0.3);
#my @Deltas= (-0.3,0.3,-0.2,0.2,-0.1,0.1,0.);
my @Deltas= (-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7);
my @Deltas=(-0.9);
my @Deltas=(0.5);
my @Ms=(0.1,0.5);
my @Deltas=(-0.8); #(-0.6,-0.8);
my @Ms=(0.5,1.0,1.5);
print "Deltas: @Deltas \n" ;

#my @Ds=(40,60,80);
#my @Ns=(100,200); 
#my @Ds=(60,80); 
#my @Ns=(200,300); 
my @Deltas= (-0.7,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9);
my @Ds=(60,80);my @Ns=(100);

my $tol=1E-6; # tolerance for convergence
my $penH=20;
my $nrStates=0;
my $mu=0; # No chemical potential
my $Starget=0; # sector of Sz=0

my $beta=10;
my $dt=0.01; # Trotter step
my $rate=2; # frequency to record data
my $freqBackup=10; # frequency to store in disk
my $app=0; # no append
my $findGS=0; # not looking for GS
my $par=" -P 10 ";#" -P 8 ";
my $par=" -P 8 ";#" -P 8 ";
my $mem=" -mem 12G ";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

my $cnt=0;
foreach $D (@Ds){
    if ($D>=80){ $par=" -P 12 ";}
    foreach $N (@Ns){
	if ($N>=200){$mem=" -mem 20G ";}
	else{$mem=" -mem 12G ";}
	foreach $Delta (@Deltas){
	    #my $g2=pi-2*acos($Delta);
	    foreach $ma (@Ms){
		my $gendir="${mpsdirBas}/N${N}";
		my $mpsdir="${gendir}/MPS";
		print MYFILE "mkdir -p ${gendir}\n";
		print MYFILE "mkdir -p ${mpsdir}\n";
	
		my $suffix="N${N}_ma${ma}_Delta${Delta}_D${D}";
		my $outputFile="${gendir}/thermalEvol_${suffix}";
		
		my $JOBNAME="thThP.N${N}.ma${ma}.Delta${Delta}.D${D}";
		my $args=" -L=${N} -ma=${ma} -gtilde2=0 -g2=${Delta} -lambda=${penH} -mu=${mu} -tol=${tol} -Starget=${Starget} -D=${D} -beta=${beta} -delta=${dt} -stepBlock=${rate}  -output=${outputFile} -append=${app} -findGS=${findGS} -backupFreq=${freqBackup} -mpsdir=${mpsdir} ";
		# Check that the job is not yet running!! Add -l h_vmem=${MEM}M 
		#my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
		# Now write the istruction to launch the job to the script
		print MYFILE "${MSUB} ${par} ${mem} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args}\n";
		print MYFILE "sleep 1\n";

	    }
	}
    }
}


close(MYFILE);

system("chmod a+x $scriptFile ");


