#!/usr/bin/perl

# Script to prepare a series of jobs using the specSto program. This time, for s>0

# Name of the script file to be created
my $scriptFile="lanzaStochPos.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/specSto0";

my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $basicDir="/ptmp/mpq/banulsm/stoch/basic";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

my $D0=20; # initial D
my $par="-P 10";
my $incrD=20;
my $maxD=100;
my $app=0;
my $tol=1E-8; # tolerance 
my $eigtol=1E-12; # tolerance for the eigensolver 
my $nlev=2; # how many levels to extract
my $saveFreq=10;
my $refine=" -refine=1 "; # to rerun with increased tolerance
my $knr=1; # fro PRIMME

# lengths of the system
#my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(20,40,60,80,100);
#my @Ns=(20,40,60,79,80,81,99,100,101,199,200,201,300,400);

# To explore small values of c
#my @Cs=(.5,.1,.05,.02);
my @Cs=(.2);

#my @Cs=(0.5); #(.1,.5);
my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(200,300,400); #(20,40,60,80,100); #300,400); #60,80,100,200);

# Values of s
my $smin = 1;
my $smax = 100;
my $sstep = 1;
my @sValues0 = map { ($smax-$sstep*$_)/1000 } 0..($smax-$smin)/$sstep;

my $smin = 1;
my $smax = 99;
my $sstep = 1;
my @sValues1 = map { ($smax-$sstep*$_)/100000 } 0..($smax-$smin)/$sstep;

my $smin = 1;
my $smax = 100;
my $sstep = 1;
my @sValues2 = map { ($smax-$sstep*$_)/10000 } 0..($smax-$smin)/$sstep;

my $smin = 1;
my $smax = 100;
my $sstep = 1;
my @sValues3 = map { ($smax-$sstep*$_)/1000000 } 0..($smax-$smin)/$sstep;

# Extra sets for some cases
$pmin =15;
$pmax=65;
$prop=1E-6; # values are $prop times the sequence of powers
$base=2;
$dec=10; # powers of base are from $pmin/$dec to $pmax/$dec
my @sValues6 = map $prop*(exp(log($base)*$_/$dec)), $pmin..$pmax;

$pmin =-15;
$pmax=14;
$prop=1E-6; # values are $prop times the sequence of powers
$base=2;
$dec=10; # powers of base are from $pmin/$dec to $pmax/$dec
my @sValues7 = map $prop*(exp(log($base)*$_/$dec)), $pmin..$pmax;

my $smin = 5;
my $smax = 40;
my $sstep = 5;
my @sValues8 = map { ($smin+$sstep*$_)/10 } 0..($smax-$smin)/$sstep;

# Extra sets for some cases
$pmin =35;
$pmax=125;
$prop=1E-6; # values are $prop times the sequence of powers
$base=2;
$dec=10; # powers of base are from $pmin/$dec to $pmax/$dec
my @sValues9 = map $prop*(exp(log($base)*$_/$dec)), $pmin..$pmax;

my @sValues=(0.,@sValues9,@sValues8); 
#my @sValues=(@sValues8); 

#my @sValues=(@sValues0,@sValues1,@sValues2,@sValues3,@sValues4,@sValues5,@sValues6,@sValues7); #

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
	    my $sStr=sprintf("%.2E",$s);
	    #my $suffix="s${s}_c${c}_$N${N}_D${D0}";
	    my $suffix="s${sStr}_c${c}_$N${N}_D${D0}";

	    my $outfile="${resultsDir}/data_${suffix}";
	    my $polfile="${resultsDir}/polX_${suffix}";

	    #my $JOBNAME="sto.c${c}.N${N}.D${D0}.s${s}";
	    my $JOBNAME="sto.c${c}.N${N}.D${D0}.s${sStr}";
    
	    # TODO!!! Decide if I am still (already) running
	    my $initmpsDir="${mpsDir}"; 
	    my $args=" -L=${N} -D=${D0} -s=${s} -c=${c} -nlevel=${nlev} -outputfile=${outfile} -polfile=${polfile} -mpsdir=${mpsDir} -initmpsdir=${initmpsDir} -tol=${tol} -tmpSavingFreq=${saveFreq} -firstN=1 ${reorder} -noise=${noise}  ${refine} -jobsdir=${jobsdir} -maxD=${maxD} -incrD=${incrD} -knr=${knr} -eigtol=${eigtol}";
	    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
	    print MYFILE "sleep 2\n";
	    #}
	}
	print MYFILE "reLaunchJobs.pl ${jobsdir} 1 1\n";
    }
}
print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 200 300 &\n";


close(MYFILE);

system("chmod a+x $scriptFile ");


