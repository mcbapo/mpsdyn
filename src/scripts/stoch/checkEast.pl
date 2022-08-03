#!/usr/bin/perl

# Script to prepare a series of jobs using the specSto program. Check whether thy were finish or prepare again

# Name of the script file to be created
my $scriptFile="lanzaEastNeg.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/specSto";

my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $basicDir="/ptmp/mpq/banulsm/stoch/basic";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

my $D0=20; # initial D
my $par="-P 4";
my $incrD=20;
my $maxD=100;
my $app=0;
my $tol=1E-8; # tolerance 
my $eigtol=1E-10; # tolerance for the eigensolver 
my $nlev=2; # how many levels to extract
my $saveFreq=10;
my $refine=""; #" -refine=1 "; # to rerun with increased tolerance
my $knr=1; # fro PRIMME

# lengths of the system
my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(20,40,60,80,100,200);
#my @Ns=(20,40,60,79,80,81,99,100,101,199,200,201,300,400);

# To explore small values of c
my @Cs=(.05); #.1,.05,.02);
#my @Cs=(0.05); #(.1,.5);

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

#print MYFILE "killall reLaunchJobs.pl \n";


foreach $c (@Cs){
    foreach $N (@Ns){
# Directory for results
	my $resultsDir="${basicDir}/c${c}/N${N}/results";
	my $mpsDir="${basicDir}/c${c}/N${N}/MPS";
	if($c==0.02){
	    #$resultsDir="${basicDir}/c${c}/N${N}/results2";
	    $mpsDir="${basicDir}/c${c}/N${N}/MPS2";
	}
	#	print MYFILE "mkdir -p ${resultsDir} \n";
	#	print MYFILE "mkdir -p ${mpsDir} \n";
	## LOOp over s values!
	foreach $s (@sValues){
	    my $complete=1;
	    for(my $D=$D0;$D<=$maxD&&$complete;$D=$D+$incrD){
		#if($s!=0.01){
		my $suffix="s${s}_c${c}_N${N}_D${D}";
		
		my $outfile="${resultsDir}/data_${suffix}";
		my $polfile="${resultsDir}/polX_${suffix}";
		
		my $JOBNAME="E.c${c}.N${N}.D${D}.s${s}";
		
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
		#print "\n look for ${nlev} levels: $complete \n";
		#exit 1;
		if($complete==0){
		    print "Case not complete! $outfile \n";
		    print MYFILE "# TODO: Check for job!!!!! ${JOBNAME} \n";
		    
		    # Decide if I am still (already) running
		    my $alreadyRunning=isRunning($JOBNAME);
		    if(! $alreadyRunning ){
			#my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
			#my $initmpsDir="."; # trick to avoid using D=40!
			my $timeNeeded=""; # -t 5 ";
			#if($N>=300&&$D>=100){$timeNeeded=" -t 60 ";}
			my $initmpsDir="${mpsDir}"; 
			my $args=" -L=${N} -D=${D} -s=${s} -c=${c} -nlevel=${nlev} -outputfile=${outfile} -polfile=${polfile} -mpsdir=${mpsDir} -initmpsdir=${initmpsDir} -tol=${tol} -tmpSavingFreq=${saveFreq} -firstN=1 ${reorder} -noise=${noise} ${reorder} ${refine} -jobsdir=${jobsdir} -maxD=${D} -incrD=${incrD} -knr=${knr} -eigtol=${eigtol}";
			print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
			print MYFILE "sleep 1\n";
			#print MYFILE "${EXEC} ${configFile} ${args} \n";
		    }
		    else{
			print "${JOBNAME} already running!\n";
		    }
		}
	    }
	}
    }
}
#print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 200 1800 &\n";


close(MYFILE);

system("chmod a+x $scriptFile ");


