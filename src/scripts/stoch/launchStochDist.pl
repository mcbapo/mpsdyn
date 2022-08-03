#!/usr/bin/perl

use File::Basename;

# Script to prepare a series of jobs using the specSto program

# Name of the script file to be created
my $scriptFile="lanzaStochDist.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/distSto";

my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $basicDir="/ptmp/mpq/banulsm/stoch/basic";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

my $D0=20; # initial D
my $par="-P 4";
my $incrD=20;
my $maxD=100;
my $app=0;
my $nlev=1; # how many levels to extract

# lengths of the system
my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(20,40,60,79,80,81,99,100,101,199,200,201,300,400);

# To explore small values of c
#my @Cs=(.1,.2,.02,.03,.05,.5);
my @Cs=(.1,.2,.05,.5);

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


open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";


    foreach $c (@Cs){
	foreach $N (@Ns){
# Directory for results
	    my $resultsDir="${basicDir}/c${c}/N${N}/results";
	    my $mpsDir="${basicDir}/c${c}/N${N}/MPS";
	    if($c==0.02){
	#	$resultsDir="${basicDir}/c${c}/N${N}/results2";
		$mpsDir="${basicDir}/c${c}/N${N}/MPS2";
	    }
	    print MYFILE "mkdir -p ${resultsDir} \n";
#	    print MYFILE "mkdir -p ${mpsDir} \n"; # Has to exist
	    # Which values of s are available??
	    opendir(DIR,$mpsDir);
	    @files = glob("$mpsDir/MPS_L${N}_s*_c${c}_D${D0}_l0.dat");
	    closedir(DIR);
	    #print @files ;print "\n";
	    @sValues=getSvalues(\@files,$N,$c);
	    #print @sValues ;
## LOOp over s values!
	    foreach $s (@sValues){

		my $suffix="s${s}_c${c}_$N${N}";

		my $outfile="${resultsDir}/distrib_${suffix}";

		my $JOBNAME="dist.c${c}.N${N}.D${D0}.s${s}";
    
		# if the file already exists, skip (?)
		if (-f $outfile ){
		    print "Skipping case c=$c$, N=$N s=$s for outputfile exists \n";
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
		    my $args=" -L=${N} -D=${D0} -s=${s} -c=${c} -nlevel=${nlev} -outputfile=${outfile} -mpsdir=${mpsDir} -firstN=1 -maxD=${maxD} -incrD=${incrD} -append=${app} ";
		    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
		    print MYFILE "sleep 1\n";
		}
	    }
	    print MYFILE "sleep 300\n";
	}

	print MYFILE "reLaunchJobs.pl ${jobsdir} 1 1 \n";
	print MYFILE "sleep 300\n";
	print MYFILE "cleanjobs dist\n";
}

print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 1000 600 & \n";

close(MYFILE);

system("chmod a+x $scriptFile ");


