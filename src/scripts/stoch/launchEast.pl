#!/usr/bin/perl

# Script to prepare a series of jobs using the specSto program

# Name of the script file to be created
my $scriptFile="lanzaEast.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/specSto";

my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $basicDir="/ptmp/mpq/banulsm/stoch/basic";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

my $D0=20; # initial D
my $par="-P 2";
my $incrD=20;
my $maxD=100;
my $app=0;
my $tol=1E-8; # tolerance 
my $eigtol=1E-10; # tolerance for the eigensolver 
my $nlev=2; # how many levels to extract
my $saveFreq=10;
my $refine=" -refine=1 "; # ""; # to rerun with increased tolerance
my $knr=1; # for PRIMME
my $maxOverlap=1E-5;

# lengths of the system
#my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(20,40,60,79,80,81,99,100,101,199,200,201,300,400);

# To explore small values of c
my @Cs=(.5,.1,.2,.05,.3);
#my @Cs=(0.5);my @Ns=(20); #(.1,.5);
my @Cs=(.05);my @Ns=(400);
#my @Ns=(20,40,60,80,100,200);
#my @sValues=map $_*1E-4, (1..30); # neet for N=60 c=0.5!

my $use_sN=0; # If 1, the values of s are used as s*N

my $pmin =-50; #139; #127;
my $pmax=-30; #194; #223;
my $prop=-1.; # values are $prop times the sequence of powers
my $base=10;
my $dec=10.; # powers of base are from $pmin/$dec to $pmax/$dec
my @sValuesNegSm = map $prop*(exp(log($base)*$_/$dec)), $pmin..$pmax;


#my @sValues=(@sValuesNegSm);

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

$smin = -0.0009;
$smax=-0.00005;
$sstep=0.00005;
my @sValues5min = map $smin+$sstep*$_, 0..($smax-$smin)/$sstep;
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

# Some positive values
$smin = 2;
$smax = 10;
$sstep = 1;
my @sValuesPos = map {($smin+$sstep*$_)*1E-3} 0..($smax-$smin)/$sstep;

#$smin = 6;
#$smax = 20;
#$sstep = 1;
#my @sValuesPos2 = map {($smin+$sstep*$_)/10} 0..($smax-$smin)/$sstep;


# Equally spaced in log scale
$pmin =-5.52;
$pmax=0.7;
$nrPts=30;
$sstep=($pmax-$pmin)/($nrPts-1);
my @sValues03 = map(exp($pmin+$sstep*$_),0..$nrPts-1) ;


$pmin =log(5E-8);
$pmax=log(1E-3);
$nrPts=50;
$sstep=($pmax-$pmin)/($nrPts-1);
my @sValues005 = map(exp($pmin+$sstep*$_),0..$nrPts-1) ;

my $smin = -95;
my $smax = -5;
my $sstep = 5;
my @sValuesNeg1 = map { ($smin+$sstep*$_)/100 } 0..($smax-$smin)/$sstep;

my $smin = -50;
my $smax = -5;
my $sstep = 5;
my @sValuesNeg2 = map { ($smin+$sstep*$_)/1000 } 0..($smax-$smin)/$sstep;

my @sValuesTot=(@sValuesPos,@sValuesPos2);


##my @sValues=(@sValues6); #-0.003); #(-0.02),@sValues5); #
##my @sValues=(@sValues0,@sValues1,@sValues2,@sValues3,@sValues4,@sValues5,@sValues6,@sValues7); #
##my @sValues=(0.,@sValues2,@sValues3,@sValues4,@sValues5,@sValues6,@sValues7); #
#my @sValues=(0.,@sValues2,@sValues3);#,@sValues4,@sValues5,@sValues6,@sValues7); #
my @sValues=(@sValuesTot);
my @sValues=(1.12221092868506e-07,1.37356942837137e-07,1.68122848060939e-07,2.05779856891804e-07,2.51871473692038e-07, 3.08286924765208e-07);

#print @sValues;exit 1;
my $reorder = " -reorder=1 "; # to force reordering of the levels
#print @sValues ;exit 1;
#my $reorder="";
my $noise= 0.01; # to increase the bond dimension

my $MSUB="msub_slurm";



# Check for existing jobs
sub isRunning{
    my ($jobname)=@_;
    #print "Checking whether job $jobname exists\n";
    my $result=`squeue --noheader -n $jobname -o "%o" `;
    if (length $result){
	print "Job $jobname exists: $result \n";
	return 1;
    }
    return 0;
};


open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

print MYFILE "killall reLaunchJobs.pl \n";

#foreach $s (@sValues){$s=-$s;} # for Pos!!

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
	foreach $s_ (@sValues){
	    if($use_sN){$s=$s_/$N;}
	    else{$s=$s_;}
	    #print "For N=${N} for sN=${s_} using s=$s \n";
	    #if($s!=0.01){
	    my $suffix="s${s}_c${c}_N${N}_D${D0}";

	    my $outfile="${resultsDir}/data_${suffix}";
	    my $polfile="${resultsDir}/polX_${suffix}";

	    my $JOBNAME="E.c${c}.N${N}.D${D0}.s${s}";
    
	    # Decide if I am still (already) running
	    my $alreadyRunning=isRunning($JOBNAME);
	    if(! $alreadyRunning ){
		my $initmpsDir="${mpsDir}"; 
		my $args=" -L=${N} -D=${D0} -s=${s} -c=${c} -nlevel=${nlev} -outputfile=${outfile} -polfile=${polfile} -mpsdir=${mpsDir} -initmpsdir=${initmpsDir} -tol=${tol} -tmpSavingFreq=${saveFreq} -firstN=1 ${reorder} -noise=${noise} ${reorder} ${refine} -jobsdir=${jobsdir} -maxD=${maxD} -incrD=${incrD} -knr=${knr} -eigtol=${eigtol}";
		print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
		#print MYFILE "sleep 1\n";
	    }
	    else{
#		print "${JOBNAME} already running!\n";
	    }
	}
	print MYFILE "cleanjobs \n";
	print MYFILE "reLaunchJobs.pl ${jobsdir} 1 1\n";
    }
#    print MYFILE "sleep 30\n";
}
print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 200 1800 &\n";


close(MYFILE);

system("chmod a+x $scriptFile ");


