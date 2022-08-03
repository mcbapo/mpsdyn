#!/usr/bin/perl

# Script to prepare a series of jobs using the FASpec program.
# In this case, fix the values of s, then pass the product s*N to the program, as it expects

# Name of the script file to be created
my $scriptFile="lanzaFA.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/FASpec";

my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $basicDir="/ptmp/mpq/banulsm/stoch/smallD/FAmodel";
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
my $refine=" -refine=1 "; # to rerun with increased tolerance
my $knr=1; # for PRIMME: probably useless
my $penZero=10; # penalty for state will all zeros!
my $maxOverlap=1E-6;

my $pbc=0;

if($pbc==1){
    $basicDir="${basicDir}/PBC";
}

# lengths of the system
my @Ns=(20,40,60,80,100,200,300,400);
#@Ns=(20,40,60,80,100,200);
#@Ns=(300);


# Values of c
my @Cs=(.5,.3,.2,.1,.05);
#@Cs=(.1);@Ns=(10);
#my @Ns=(20,40,60,80,100,200);
#@Cs=(0.2,0.5,0.02);@Ns=(20,40); # just for density plot!
#@Cs=(0.02);####@Ns=(20,40); # just for density plot!
#@Ns=(300,400);
#@Cs=(0.1,0.3);

my $use_sN=0; # if not, the values of s are multiplied by N, if 1, used directly as s*N

# Values of s


# Equally spaced in log scale
$pmin =-6;
$pmax=0;
$nrPts=200;
$sstep=($pmax-$pmin)/($nrPts);
my @sValuesTot = map(10**($pmin+$sstep*$_),0..$nrPts) ;

# For pbc it s has to be larger/smaller 
#$pmin =-8;
#$pmax=-6;
#$nrPts=int(($pmax-$pmin)/($sstep));
###$sstep=($pmax-$pmin)/($nrPts);
#my @sValuesTot = map(10**($pmin+$sstep*$_),0..$nrPts) ;



#foreach $s (@sValuesTot){
#    print "$s ";
#}
#exit 1;
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

foreach $c (@Cs){
    foreach $N (@Ns){
# Directory for results
	my $resultsDir="${basicDir}/c${c}/N${N}/results";
	my $mpsDir="${basicDir}/c${c}/N${N}/MPS";
	print MYFILE "mkdir -p ${resultsDir} \n";
	print MYFILE "mkdir -p ${mpsDir} \n";
	my $mytol=$tol;
	#if($N>=200){ $mytol=1E-6; }
	my $refine_=$refine;
	if($N>=200){$refine_=" -refine=1 ";} # to rerun with increased tolerance (before it was too small)
## LOOp over s values!
	foreach $s_ (@sValuesTot){
	    my $s=$s_*$N;
	    if($use_sN){$s=$s_;}
	    #if($s!=0.01){
	    my $sStr=sprintf("%g",$s);
	    my $complete=1;
	    for(my $D=$D0;$D<=$maxD&&$complete;$D=$D+$incrD){
		my $suffix="sL${sStr}_c${c}_N${N}_D${D}";

		my $outfile="${resultsDir}/data_${suffix}";
		my $polfile="${resultsDir}/polX_${suffix}";

		my $JOBNAME="fa.c${c}.N${N}.D${D}.s${sStr}";
		# Check if levels from 0 to nlev-1 exist in data file
		my $levels=`cut -f 1 ${outfile}|grep -v "%"`; # Just numbers, discard headers
		@levels=split('\n',$levels); # list of integers
		for(my $lv=0;$lv<$nlev&&$complete==1;$lv++){
		    # if this level is not there, say not complete and start job
		    my $found=0;
		    #print "Will search for $lv in a list of ",scalar(@levels),"\n";
		    my $ex=0;
		    while($found==0&&$ex<scalar(@levels)){
			if(@levels[$ex]==$lv){
			    $found=1;
			}
			$ex++;
		    }
		    $complete=$found;
		}
		if($complete==0){
		    print "Case not complete! $outfile \n";
    
		    # Decide if I am still (already) running
		    my $alreadyRunning=isRunning($JOBNAME);
		    if(! $alreadyRunning ){
			my $sStrS=sprintf("%g",$s/$N);		
			my $JOBNAMES="fa.c${c}.N${N}.D${D}.s${sStrS}";
			if(!isRunning($JOBNAMES)){
			    my $initmpsDir="${mpsDir}"; 
			    my $args=" -L=${N} -D=${D} -sL=${s} -c=${c} -nlevel=${nlev} -outputfile=${outfile} -polfile=${polfile} -mpsdir=${mpsDir} -initmpsdir=${initmpsDir} -tol=${mytol} -tmpSavingFreq=${saveFreq} -firstN=1 ${reorder} -noise=${noise}  ${refine_} -jobsdir=${jobsdir} -maxD=${maxD} -incrD=${incrD} -knr=${knr} -eigtol=${eigtol} -penaltyZero=$penZero -pbc=${pbc} -maxOverlap=${maxOverlap}";
			    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
			    #print MYFILE "sleep 1\n";
			}
			else{
			    print "$JOBNAMES already running\n";
			}
		    }
		    else{
			print "$JOBNAME already running\n";
		    }
		}
	    } # D loop
	}# s loop
	print MYFILE "cleanjobs \n";
	print MYFILE "sleep 5\n";
	print MYFILE "reLaunchJobs.pl ${jobsdir} 5 30\n";
	print MYFILE "cleanjobs \n";
    }  # End of N loop
    #    print MYFILE "sleep 1800\n";
} # End of c loop
print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 2000 60 &\n";


close(MYFILE);

system("chmod a+x $scriptFile ");


