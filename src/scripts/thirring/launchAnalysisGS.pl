#!/usr/bin/perl
use POSIX;
# Script to prepare a series of jobs that explore GS of Thirring model, computing expectation vales of individual Hamiltonian terms

# Name of the script file to be created
my $scriptFile="lanzaGSAnalysisThirring.sh";

#my $EXEC="./thThirring";
my $EXEC="./thirringGS";

my $MSUB="msub_slurm";

my $configFile="../config/thermalThirring.conf";
#my $configFile="../config/thermalThirring_Delta-0.9.conf";

#my $mpsdirBas="/ptmp/mpq/banulsm/thirring/thermal";
my $mpsdirBas="/ptmp/mpq/banulsm/thirring/GS";
#my $jobsdir="/u/banulsm/jobsToRun";

# Hamiltonian parameters
#my @Deltas=(-0.9,0.9);
#my @Ms=(0.,0.1,0.3);
#my @Deltas= (-0.3,0.3,-0.2,0.2,-0.1,0.1,0.);

my @Deltas=(-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9); #-0.8
my @Ms=(0.1,0.5,1.0,1.5);
my @Ms=(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5);
print "Deltas: @Deltas \n" ;

#my $D0=20; # in the name of file, to append
my $app=0; # no append
my @Ns=(300);
my @Ds=(180); # Data with up to 180

my $tol=1E-8; # tolerance for convergence
my $penH=20;
my $nrStates=0;
my $mu=0; # No chemical potential
my $Starget=0; # sector of Sz=0


my $par=" -P 4 ";#" -P 8 ";
my $mem=" -mem 12G ";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

my $cnt=0;
foreach $D (@Ds){
    foreach $N (@Ns){
	if ($N>=300){$mem=" -mem 20G ";}
	else{$mem=" -mem 12G ";}
	foreach $Delta (@Deltas){
	    #my $g2=pi-2*acos($Delta);
	    foreach $ma (@Ms){
		my $gendir="${mpsdirBas}/N${N}";
		my $mpsdir="${gendir}/MPS";
		my $outputdir="${gendir}/EVs";
		
		print MYFILE "mkdir -p ${outputdir}\n";
	
		my $suffix="N${N}_ma${ma}_Delta${Delta}_D${D}";
		my $outputFile="${outputdir}/evGS_${suffix}";
		
		my $JOBNAME="thGScorr.N${N}.ma${ma}.Delta${Delta}.D${D}";
		my $args=" -L=${N} -ma=${ma} -gtilde2=0 -g2=${Delta} -lambda=${penH} -mu=${mu} -tol=${tol} -Starget=${Starget} -D=${D} -output=${outputFile} -append=${app} -mpsdir=${mpsdir} ";
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


