#!/usr/bin/perl

# A few points by hand around some peak, to refine

# Name of the script file to be created
my $scriptFile="lanzaFAPeak.sh";

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
my $nlev=2; # how many levels to extract
my $saveFreq=10;
my $refine=" -refine=0 "; # to rerun with increased tolerance
my $knr=1; # for PRIMME: probably useless
my $penZero=10; # penalty for state will all zeros!

my $pbc=0;

if($pbc==1){
    $basicDir="${basicDir}/PBC";
}

# lengths of the system
#my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(20,40,60,80,100);
#my @Ns=(20,40,60,79,80,81,99,100,101,199,200,201,300,400);

# To explore small values of c
#my @Cs=(.5,.2,.1,.05,.02);
my @Cs=(.5);

#my @Cs=(0.5); #(.1,.5);
#my @Ns=(20,40,60,80,100,200,300,400);

my @Ns=(60);#200,300,400); $tol=1E-6; 
#my @Ns=(200,300,400); #(20,40,60,80,100); #300,400); #60,80,100,200);

# Values of sL: FIRST PROBE - at least for N=36, c=0.3 they are in this range
# So I also use this one for c=0.2
#my @sValues2 = (0.0367,0.0368,0.0369,0.0371,0.0372); # for N=100
my @sValues2 = (0.0366,0.0365); # for N=100
#my @sValues2 = (0.03075,.0307,0.0308,0.0309,0.0311,0.0312); # for N=200
#my @sValues2 = (0.03065,0.0307,0.03075,0.0308);

# For c=0.1 the PT seems around 0.01
$pmin =12;
#$pmax=165;
#$prop=1E-6; # values are $prop times the sequence of powers
#$pmax=60;
#$prop=1E-3; # values are $prop times the sequence of powers
#$base=2;
#$dec=10; # powers of base are from $pmin/$dec to $pmax/$dec
#my @sValues1 = map $prop*(exp(log($base)*$_/$dec)), $pmin..$pmax;
#my @sValues1=(0.01150,0.01170,0.01190,0.01205);
my @sValues1=(0.0135,.0145);

# For c=0.5 I need larger values (PT around .13?)
#my @sValues5 = (0.0945); # N=100
#my @sValues5 = (0.0975); # N=80
my @sValues5 = (0.1035); #(0.1032,0.1034,0.1036,0.1038); # N=60
#my @sValues5 = (0.1142,0.1144,0.1146,0.1148); # N=40

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

foreach $c (@Cs){
# Really each c needs a range of values:
    my @sValues;
    if($c==0.02){
	@sValues=(@sValues02);
    }
    if($c==0.05){
	@sValues=(@sValues05);
    }
    if($c==0.1){
	@sValues=(@sValues1);
    }
    if($c==0.2){
	@sValues=(@sValues2);
    }
    if($c==0.5){
	@sValues=(@sValues5);
    }
    foreach $N (@Ns){
# Directory for results
	my $resultsDir="${basicDir}/c${c}/N${N}/results";
	my $mpsDir="${basicDir}/c${c}/N${N}/MPS";
	print MYFILE "mkdir -p ${resultsDir} \n";
	print MYFILE "mkdir -p ${mpsDir} \n";
	my $mytol=$tol;
	if($N>=200){ $mytol=1E-6; }
## LOOp over s values!
	foreach $s (@sValues){
	    #if($s!=0.01){
	    #my $sStr=sprintf("%.2E",$s);
	    my $sStrG=sprintf("%g",$s);
	    #my $suffix="s${s}_c${c}_$N${N}_D${D0}";
	    #my $suffix="sL${sStr}_c${c}_N${N}_D${D0}";
	    my $suffix="sL${sStrG}_c${c}_N${N}_D${D0}";

	    my $outfile="${resultsDir}/data_${suffix}";
	    my $polfile="${resultsDir}/polX_${suffix}";

	    #my $JOBNAME="sto.c${c}.N${N}.D${D0}.s${s}";
	    my $JOBNAME="fa.c${c}.N${N}.D${D0}.s${sStrG}";
    
	    # TODO!!! Decide if I am still (already) running
	    my $initmpsDir="${mpsDir}"; 
	    my $args=" -L=${N} -D=${D0} -sL=${s} -c=${c} -nlevel=${nlev} -outputfile=${outfile} -polfile=${polfile} -mpsdir=${mpsDir} -initmpsdir=${initmpsDir} -tol=${mytol} -tmpSavingFreq=${saveFreq} -firstN=1 ${reorder} -noise=${noise}  ${refine} -jobsdir=${jobsdir} -maxD=${maxD} -incrD=${incrD} -knr=${knr} -eigtol=${eigtol} -penaltyZero=$penZero -pbc=${pbc} ";
	    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
	    print MYFILE "sleep 1\n";
	    #}
	}
    }
}


close(MYFILE);

system("chmod a+x $scriptFile ");


