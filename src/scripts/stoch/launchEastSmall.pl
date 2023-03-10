#!/usr/bin/perl

# Script to prepare a series of jobs using the specSto program

# Name of the script file to be created
my $scriptFile="lanzaEast.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/specSto";

my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $basicDir="/ptmp/mpq/banulsm/stoch/smallD/basic";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

my $D0=2; # initial D
my $par="-P 1";
my $incrD=2;
my $maxD=20;
my $app=0;
my $tol=1E-12; # tolerance 
my $eigtol=1E-12; # tolerance for the eigensolver 
my $nlev=4; # how many levels to extract
my $saveFreq=10;
my $refine=" -refine=0 "; # ""; # to rerun with increased tolerance
my $knr=1; # for PRIMME
my $maxOverlap=1E-5;

# lengths of the system
my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(20,40,60,79,80,81,99,100,101,199,200,201,300,400);

# To explore small values of c
my @Cs=(.5,.1,.2,.3,.05);
my @Cs=(.5);@Ns=(100,200,300,400);$par="-P 2";
#my @Cs=(0.5);my @Ns=(20); #(.1,.5);
#my @Cs=(.05);my @Ns=(400);
#my @Ns=(20,40,60,80,100,200);
#my @sValues=map $_*1E-4, (1..30); # neet for N=60 c=0.5!

my $use_sN=0; # If 1, the values of s are used as s*N


# Equally spaced in log scale
$pmin =-6;
$pmax=0;
$nrPts=200;
$sstep=($pmax-$pmin)/($nrPts);
my @sValuesTot = map(10**($pmin+$sstep*$_),0..$nrPts) ;



my $reorder = " -reorder=1 "; # to force reordering of the levels
#print @sValuesTot ;exit 1;
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
	foreach $s_ (@sValuesTot){
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
	print MYFILE "sleep 1\n";
	print MYFILE "cleanjobs \n";
	print MYFILE "reLaunchJobs.pl ${jobsdir} 1 1\n";
    }
#    print MYFILE "sleep 30\n";
}
print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 2000 30 &\n";


close(MYFILE);

system("chmod a+x $scriptFile ");


