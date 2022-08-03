#!/usr/bin/perl

# Script to prepare a series of jobs using the dos program

# Name of the script file to be created
my $scriptFile="lanzadorDoSrnd.sh";

my $EXEC="./dosHeisRnd";
my $par=" -P 2 "; # use parallel environment
#my $par=" "; # do not use parallel environment

$basicdir="/ptmp/mpq/banulsm/DoS/HeisRandom/E-Var";
$configFile="/afs/ipp/u/banulsm/SVN/mpsdyn/src/config/dos/dos.conf";

# where are the random coefficients
$sampleFile="../config/randomDataLoc_31012014.dat";

###### PARAMETERS TO TUNE ##########

## Instances
my $nrI1=1;
my $maxI=10;
my $nrI2=$nrI1+$maxI;

my @Ns=(4,6,8,10); # system sizes

# Parameters for the MPS simulation (again, they could be ranges of
# values or just one)
my @Ds=(1,2,4,8,16);
my $nrBins=30;
#my $nrBinsV=1; # no variance binning
my $nrBinsV=30; # variance binning
my $maxVar=-1; #default
my $J=1;
#my $g=1.;
my @Hs=(0.,10.);
my $f0=exp(1);
my @Mlogs=(4);
my $maxF=30000; # max nr of steps in f
my $isXY=1;


my $MSUB="msub_slurm";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";
print MYFILE "#!/bin/bash\n\n";


foreach $Mlog (@Mlogs){
    my $M=10**$Mlog;
    foreach $N (@Ns){
	$maxVar=$N;
	my $exactD=2**(.5*$N);
	foreach $h (@Hs){
	    my $outputdir="${basicdir}/Jz0/h${h}/N${N}";
	    print MYFILE "mkdir -p ${outputdir} \n";
	    for (my $nrI=${nrI1};$nrI<=${nrI2};$nrI=$nrI+1){
		foreach $D (@Ds){
		    if($D<=$exactD){
			my $suffix="N${N}_J${J}_h${h}_D${D}_M1e${Mlog}";
			$outdataFile="${outputdir}/results_${suffix}";
			#$outhistoFile="${outputdir}/histos_${suffix}.m";
			print "Now g=$g N=$N D=$D M=$M \n";
			
			my $JOBNAME="dosRnd.N${N}.g${g}.D${D}.M1e${Mlog}";
			# Check that the job is not yet running!! Add -l h_vmem=${MEM}M 
			#my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
			# Now write the istruction to launch the job to the script
			print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} -L=${N} -J=${J} -h=${h} -isXY=${isXY} -paramFile=${sampleFile} -nrInst=${nrI} -D=$D -Nbins=${nrBins} -NbinsV=${nrBinsV} -f0=${f0} -outputFile=${outdataFile} -maxCnt=${M} -maxCntF=${maxF} -BPcrit=1 -maxVar=${maxVar} \n";
			
			print MYFILE "sleep 1\n";
		    }
		}	 
	    }   
	}
    }
}

close(MYFILE);

system("chmod a+x $scriptFile ");


