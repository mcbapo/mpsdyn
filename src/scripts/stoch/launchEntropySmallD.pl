#!/usr/bin/perl

use File::Basename;

# Script to prepare a series of jobs using the results from the specSto program
# such that entropy is computed for all ground states found



my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

my $D0=2; # initial D
my $par="-P 1";
my $incrD=2;
my $maxD=20;
my $app=0;
my $nlev=4; # how many levels to extract
my $skipDone=1; # skip if the entropy file exists
#$D0=2;$maxD=20;$incrD=2; # for the lim s->-inf in both models

# lengths of the system
my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(300,400);


# To explore small values of c
#my @Cs=(.1,.2,.3,.02,.03,.05,.5);
#my @Cs=(.1,.2,.3,.05,.5);
#my @Cs=(0.1,0.3,0.5,0.02,0.05);
#my @Cs=(0.02); 

# extra points on the right
my @Cs=(.1,.2,.3,.05,.5);
#my @Ns=(20,40,60,80,100);

my $model=0; # East
#my $model=1; # FA
if($model==0){
    $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/eastEntropy";
    $basicDir="/ptmp/mpq/banulsm/stoch/smallD/basic";
    # Name of the script file to be created
    $scriptFile="lanzaEastEntropy.sh";
}
else{
    $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/FAdist";
    $basicDir="/ptmp/mpq/banulsm/stoch/smallD/FAmodel";
    $scriptFile="lanzaFAEntropy.sh";
}
my $MSUB="msub_slurm";

# Pass the list of files and get the sValues inside
sub getSvalues {
    my $N=$_[1];
    my $c=$_[2];
    @sValues=();
    for $f (@{$_[0]}) {
	my $shortN=basename($f);
	$shortN =~ /.*L${N}_s(.+)_c${c}.*/;
	push(@sValues,$1);
    }
    return @sValues;
}

# For FA Pass the list of files and get the sLValues inside
sub getSLvalues {
    my $N=$_[1];
    my $c=$_[2];
    @sValues=();
    for $f (@{$_[0]}) {
	my $shortN=basename($f);
	$shortN =~ /.*L${N}_sL(.+)_eps${c}.*/;
	push(@sValues,$1);
    }
    return @sValues;
}


open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";


foreach $c (@Cs){
    foreach $N (@Ns){
# Directory for results
	my $resultsDir="${basicDir}/c${c}/N${N}/results";
	my $mpsDir="${basicDir}/c${c}/N${N}/MPS";
	if($c==0.02&&$model==0){
	    #	$resultsDir="${basicDir}/c${c}/N${N}/results2";
	    $mpsDir="${basicDir}/c${c}/N${N}/MPS2";
	}
	print MYFILE "mkdir -p ${resultsDir} \n";
	#	    print MYFILE "mkdir -p ${mpsDir} \n"; # Has to exist
	# Which values of s are available??
	opendir(DIR,$mpsDir);
	if($model==0){
	    @files = glob("$mpsDir/MPS_L${N}_s*_c${c}_D${D0}_l0.dat");
	}
	else{
	    @files = glob("$mpsDir/MPS_L${N}_sL*_eps${c}_D${D0}_l0.dat");
	}
	closedir(DIR);
	#print @files ;print "\n";
	if($model==0){
	    @sValues=getSvalues(\@files,$N,$c);
	}
	else{
	    @sValues=getSLvalues(\@files,$N,$c);
	}
	#print @sValues ;
	## LOOp over s values!
	foreach $s (@sValues){
	    
	    if($model==0){
		$suffix="s${s}_c${c}_N${N}";
		$JOBNAME="eastEntropy.c${c}.N${N}.D${D0}.s${s}";
		$outfile="${resultsDir}/entropy_${suffix}";
		$outfileSchm="${resultsDir}/Schmidt_${suffix}";
	    }
	    if($model==1){
		$suffix="sL${s}_eps${c}_N${N}";
		$JOBNAME="faEntropy.c${c}.N${N}.D${D0}.s${s}";
		$outfile="${resultsDir}/distrib_${suffix}";
		$outfileS="${resultsDir}/entropy_${suffix}";
		$outfileSchm="${resultsDir}/Schmidt_${suffix}";
	    }
		
	    # If the file already exists, skip (?)
	    if ($skipDone &&(-f $outfile) ){
		print "Skipping case c=$c$, N=$N s(L)=$s for outputfile exists \n";
	    }
	    else{
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
		my $args=" -s=${s} -L=${N} -D=${D0} -c=${c} -outputfile=${outfile} -mpsdir=${mpsDir} -firstN=1 -maxD=${maxD} -incrD=${incrD} -append=${app} -outputfileSchmidt=${outfileSchm} -append=${app} ";
		if($model==1){
		    $args="${args} -s=${s} ";
		}
		if($model==1){
		    $args="${args} -sL=${s} -outputfileS=${outfileS} "
		}
		print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
		print MYFILE "sleep 1\n";
	    }
	}
	print MYFILE "cleanjobs \n";
	print MYFILE "sleep 600\n";
    }

}

close(MYFILE);

system("chmod a+x $scriptFile ");


