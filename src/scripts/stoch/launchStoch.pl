#!/usr/bin/perl

# Script to prepare a series of jobs using the specSto program

# Name of the script file to be created
my $scriptFile="lanzaStoch.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/specSto";

my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $basicDir="/ptmp/mpq/banulsm/stoch/basic";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

my $D0=20; # initial D
my $par="-P 10";
my $incrD=20;
my $maxD=100;
my $app=0;
my $tol=1E-10; # tolerance 
my $eigtol=1E-12; # tolerance for the eigensolver 
my $nlev=2; # how many levels to extract
my $saveFreq=10;
my $refine=""; #" -refine=1 "; # to rerun with increased tolerance
my $knr=1; # fro PRIMME

# lengths of the system
my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(20,40,60,79,80,81,99,100,101,199,200,201,300,400);

# To explore small values of c
#my @Cs=(.5,.1,.05,.02);
my @Cs=(0.2,0.3,0.03); #(.1,.5);
my @Cs=(0.05); # left out before!

# Values of s
my $smin = -10.;
my $smax = -1.;
my $sstep = 1.;
my @sValues0 = map $smin+$sstep*$_, 0..($smax-$smin)/$sstep;

## REFINE
$smin = -.9;
$smax = -.5;
$sstep = .1;
my @sValues1 = map $smin+$sstep*$_, 0..($smax-$smin)/$sstep;

$smin = -.45;
$smax=-0.15;
$sstep=0.05;
my @sValues2 = map $smin+$sstep*$_, 0..($smax-$smin)/$sstep;

$smin = -.1;
$smax=-0.01;
$sstep=0.01;
my @sValues3 = map $smin+$sstep*$_, 0..($smax-$smin)/$sstep;

$smin = -2.5;
$smax=-1.;
$sstep=0.1;
my @sValues4 = map $smin+$sstep*$_, 0..($smax-$smin)/$sstep;

$smin = -0.019;
$smax=-0.001;
$sstep=0.001;
my @sValues5 = map $smin+$sstep*$_, 0..($smax-$smin)/$sstep;

# Extra sets for some cases
$smin =-0.045;
$smax=-0.02;
$sstep=0.005;
my @sValues6 = map $smin+$sstep*$_, 0..($smax-$smin)/$sstep;

$smin = -.095;
$smax=-0.055;
$sstep=0.005;
my @sValues6 = map $smin+$sstep*$_, 0..($smax-$smin)/$sstep;

$smin =-0.0119;
$smax=-0.0111;
$sstep=0.0001;
my @sValues7 = map $smin+$sstep*$_, 0..($smax-$smin)/$sstep;


#my @sValues=(@sValues6); #-0.003); #(-0.02),@sValues5); #
my @sValues=(@sValues0,@sValues1,@sValues2,@sValues3,@sValues4,@sValues5,@sValues6,@sValues7); #

#print @sValues;exit 1;
my $reorder = " -reorder=1 "; # to force reordering of the levels
#print @sValues ;exit 1;
#my $reorder="";
my $noise= 0.01; # to increase the bond dimension

my $MSUB="msub_slurm";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

print MYFILE "killall reLaunchJobs.pl \n";

foreach $c (@Cs){
    foreach $N (@Ns){
# Directory for results
	my $resultsDir="${basicDir}/c${c}/N${N}/results";
	my $mpsDir="${basicDir}/c${c}/N${N}/MPS";
	if($c==0.02){
	    #$resultsDir="${basicDir}/c${c}/N${N}/results2";
	    $mpsDir="${basicDir}/c${c}/N${N}/MPS2";
	}
	print MYFILE "mkdir -p ${resultsDir} \n";
	print MYFILE "mkdir -p ${mpsDir} \n";
## LOOp over s values!
	foreach $s (@sValues){
	    #if($s!=0.01){
	    my $suffix="s${s}_c${c}_N${N}_D${D0}";

	    my $outfile="${resultsDir}/data_${suffix}";
	    my $polfile="${resultsDir}/polX_${suffix}";

	    my $JOBNAME="sto.c${c}.N${N}.D${D0}.s${s}";
    
	# TODO!!! Decide if I am still (already) running
	#my $isRunning=`qstat -j $JOBNAME |grep job_number`;
	#if($isRunning =~ /(\d+)/){
	#$jobNr=$1;
	#print MYFILE "### Job already/still running $jobNr for $JOBNAME ## \n";
#}
#	else{
	    #my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
	#my $initmpsDir="."; # trick to avoid using D=40!
	    my $initmpsDir="${mpsDir}"; 
	    my $args=" -L=${N} -D=${D0} -s=${s} -c=${c} -nlevel=${nlev} -outputfile=${outfile} -polfile=${polfile} -mpsdir=${mpsDir} -initmpsdir=${initmpsDir} -tol=${tol} -tmpSavingFreq=${saveFreq} -firstN=1 ${reorder} -noise=${noise} ${reorder} ${refine} -jobsdir=${jobsdir} -maxD=${maxD} -incrD=${incrD} -knr=${knr} -eigtol=${eigtol}";
	    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
	    print MYFILE "sleep 1\n";
	    #}
	}
    }
    print MYFILE "reLaunchJobs.pl ${jobsdir} 1 1\n";	
}
print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 500 300 &\n";



close(MYFILE);

system("chmod a+x $scriptFile ");


