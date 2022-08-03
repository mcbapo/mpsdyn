#!/usr/bin/perl

# Script to prepare a series of jobs using the specSto program. This time, for s<0, but very large (s=10), for all sizes 

# Name of the script file to be created
my $scriptFile="lanzaEastlim.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/specSto";

my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $basicDir="/ptmp/mpq/banulsm/stoch/basic";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

my $D0=2; # initial D
my $par="-P 2";
my $incrD=2;
my $maxD=20;
my $app=0;
my $tol=1E-8; # tolerance 
my $eigtol=1E-12; # tolerance for the eigensolver 
my $nlev=2; # how many levels to extract
my $saveFreq=10;
my $refine=" -refine=1 "; # to rerun with increased tolerance
my $knr=1; # for PRIMME: probably useless
my $penZero=10; # penalty for state will all zeros!
my $ignoretmp=" -ignoretmp=1 ";

my $pbc=1;

if($pbc==1){
    $basicDir="${basicDir}/PBC";
}

# lengths of the system
#my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(20,40,60,80,100);
#my @Ns=(20,40,60,79,80,81,99,100,101,199,200,201,300,400);

# To explore small values of c
#my @Cs=(.5,.2,.1,.05,.02);
my $c=.5;

#my @Cs=(0.5); #(.1,.5);
my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(200,300,400); $tol=1E-6; 
#my @Ns=(200,300,400); #(20,40,60,80,100); #300,400); #60,80,100,200);



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


my @sValues=(-1111);
foreach $N (@Ns){
# Directory for results
    my $resultsDir="${basicDir}/c${c}/N${N}/results";
    my $mpsDir="${basicDir}/c${c}/N${N}/MPS";
    print MYFILE "mkdir -p ${resultsDir} \n";
    print MYFILE "mkdir -p ${mpsDir} \n";
    my $mytol=$tol;
    if($N>=200){ $mytol=1E-6; }
## LOOp over s values!
    foreach $sV (@sValues){
	my $s=$sV;
	my $sStr=sprintf("%g",$s);
	#my $suffix="s${s}_c${c}_$N${N}_D${D0}";
	my $suffix="s${sStr}_c${c}_N${N}_D${D0}";

	my $outfile="${resultsDir}/data_${suffix}";
	my $polfile="${resultsDir}/polX_${suffix}";

	my $JOBNAME="E.c${c}.N${N}.D${D0}.s${sStr}";
	if($pbc){
	    $JOBNAME="EPBC.c${c}.N${N}.D${D0}.s${sStr}";
	}
    
	# Decide if I am still (already) running
	my $alreadyRunning=isRunning($JOBNAME);
	if(! $alreadyRunning ){
	    
	    my $initmpsDir="${mpsDir}"; 
	    my $args=" -L=${N} -D=${D0} -s=${s} -c=${c} -nlevel=${nlev} -outputfile=${outfile} -polfile=${polfile} -mpsdir=${mpsDir} -initmpsdir=${initmpsDir} -tol=${mytol} -tmpSavingFreq=${saveFreq} -firstN=1 ${reorder} -noise=${noise}  ${refine} -jobsdir=${jobsdir} -maxD=${maxD} -incrD=${incrD} -knr=${knr} -eigtol=${eigtol} -penaltyZero=$penZero -pbc=${pbc} ${ignoretmp} ";
	    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
	    print MYFILE "sleep 1\n";
	}
	else{
	    print "${JOBNAME} already running!\n";
	}
    }
    print MYFILE "reLaunchJobs.pl ${jobsdir} 1 1\n";
}
print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 200 1800 &\n";


close(MYFILE);

system("chmod a+x $scriptFile ");


