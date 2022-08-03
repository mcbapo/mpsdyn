#!/usr/bin/perl

# Script to prepare a series of jobs using the specSto program. This time, for s>0

# Name of the script file to be created
my $scriptFile="checkFAPos.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/FASpec";

my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $basicDir="/ptmp/mpq/banulsm/stoch/FAmodel";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

my $D0=20; # initial D
my $par="-P 4";
my $incrD=20;
my $maxD=100;
my $app=0;
my $tol=1E-8; # tolerance 
my $eigtol=1E-12; # tolerance for the eigensolver 
my $nlev=3; # how many levels to extract
my $saveFreq=10;
my $refine=" -refine=1 "; # to rerun with increased tolerance
my $knr=1; # for PRIMME: probably useless
my $penZero=10; # penalty for state will all zeros!

# lengths of the system
#my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(20,40,60,80,100);
#my @Ns=(20,40,60,79,80,81,99,100,101,199,200,201,300,400);

# To explore small values of c
my @Cs=(.5,.2,.1,.05,.02);
my @Cs=(.05);

#my @Cs=(0.5); #(.1,.5);
my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(200,300,400); #(20,40,60,80,100); #300,400); #60,80,100,200);

# Values of sL: FIRST PROBE - at least for N=36, c=0.3 they are in this range
# So I also use this one for c=0.2
my $smin = 1;
my $smax = 100;
my $sstep = 1;
my @sValues2 = map { ($smax-$sstep*$_)/1000 } 0..($smax-$smin)/$sstep;

# For c=0.1 the PT seems around 0.01
$pmin =12;
#$pmax=165;
#$prop=1E-6; # values are $prop times the sequence of powers
$pmax=60;
$prop=1E-3; # values are $prop times the sequence of powers
$base=2;
$dec=10; # powers of base are from $pmin/$dec to $pmax/$dec
my @sValues1 = map $prop*(exp(log($base)*$_/$dec)), $pmin..$pmax;

# For c=0.5 I need larger values (PT around .13?)
my $smin = 60;
my $smax = 180;
my $sstep = 1;
my @sValues5 = map { ($smax-$sstep*$_)/1000 } 0..($smax-$smin)/$sstep;
#$pmin =945;
#$pmax=1045;
#$prop=1E-6; # values are $prop times the sequence of powers
#$base=2;
#$dec=50; # powers of base are from $pmin/$dec to $pmax/$dec
#my @sValues5 = map $prop*(exp(log($base)*$_/$dec)), $pmin..$pmax;

# For c=0.05
$pmin =127;
$pmax=223;
$prop=1E-6; # values are $prop times the sequence of powers
$base=2;
$dec=100./7.; # powers of base are from $pmin/$dec to $pmax/$dec
my @sValues05 = map $prop*(exp(log($base)*$_/$dec)), $pmin..$pmax;

# For c=0.02: smaller
$pmin =97;
$pmax=193;
$prop=1E-6; # values are $prop times the sequence of powers
$base=2;
$dec=100./7.; # powers of base are from $pmin/$dec to $pmax/$dec
my @sValues02 = map $prop*(exp(log($base)*$_/$dec)), $pmin..$pmax;


#print @sValues;exit 1;
my $reorder = " -reorder=1 "; # to force reordering of the levels
#print @sValues ;exit 1;
#my $reorder="";
my $noise= 0.01; # to increase the bond dimension

my $MSUB="msub_slurm";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

#print MYFILE "killall reLaunchJobs.pl \n";

foreach $c (@Cs){
# Really each c needs a range of values:
    my @sValues;
    if($c==0.02){
	@sValues=(0.,@sValues02);
    }
    if($c==0.05){
	@sValues=(0.,@sValues05);
    }
    if($c==0.1){
	@sValues=(0.,@sValues1);
    }
    if($c==0.2){
	@sValues=(0.,@sValues2);
    }
    if($c==0.5){
	@sValues=(0.,@sValues5);
    }
    foreach $N (@Ns){
# Directory for results
	my $resultsDir="${basicDir}/c${c}/N${N}/results";
	my $mpsDir="${basicDir}/c${c}/N${N}/MPS";
	my $mytol=$tol;
	if($N>=200){ $mytol=1E-6;}
#	print MYFILE "mkdir -p ${resultsDir} \n";
#	print MYFILE "mkdir -p ${mpsDir} \n";
## LOOp over s values!
	foreach $s (@sValues){
	    my $complete=1;
	    for(my $D=$D0;$D<=$maxD&&$complete;$D=$D+$incrD){
		#if($s!=0.01){
		my $sStr=sprintf("%.2E",$s);
		#my $suffix="s${s}_c${c}_$N${N}_D${D0}";
		my $suffix="sL${sStr}_c${c}_N${N}_D${D}";
		
		my $outfile="${resultsDir}/data_${suffix}";
		my $polfile="${resultsDir}/polX_${suffix}";

		#my $JOBNAME="sto.c${c}.N${N}.D${D0}.s${s}";
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
		#print "\n look for ${nlev} levels: $complete \n";
		#exit 1;
		if($complete==0){
		    print "Case not complete! $outfile \n";
		    print MYFILE "# TODO: Check for job!!!!! ${JOBNAME} \n";
		    # TODO!!! Decide if I am still (already) running
		    my $initmpsDir="${mpsDir}"; 
		    my $args=" -L=${N} -D=${D} -sL=${s} -c=${c} -nlevel=${nlev} -outputfile=${outfile} -polfile=${polfile} -mpsdir=${mpsDir} -initmpsdir=${initmpsDir} -tol=${mytol} -tmpSavingFreq=${saveFreq} -firstN=1 ${reorder} -noise=${noise}  ${refine} -jobsdir=${jobsdir} -maxD=${maxD} -incrD=${incrD} -knr=${knr} -eigtol=${eigtol} -penaltyZero=$penZero";
		    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
		    print MYFILE "sleep 2\n";
		    #print MYFILE "${EXEC} ${configFile} ${args} \n";
		}
	    }
	}
#	print MYFILE "reLaunchJobs.pl ${jobsdir} 1 1\n";
    }
#    print MYFILE "sleep 300\n";
}
#print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 200 300 &\n";


close(MYFILE);

system("chmod a+x $scriptFile ");


