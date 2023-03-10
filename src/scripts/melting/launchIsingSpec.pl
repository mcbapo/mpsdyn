#!/usr/bin/perl

# Script to prepare a series of jobs using the isingSpec program

# Name of the script file to be created
my $scriptFile="lanzaSpec.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/isingSpec";

my $basicDir="/ptmp/mpq/banulsm/melting/dataIsing";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

#my $D0=20; # initial D
my $D0=200; # initial D
my $par="-P 4";
my $incrD=20;
my $maxD=200;
my $app=0;
my $tol=1E-8; # tolerance 
#my $nlev=2; # how many levels to extract
#my $knr=1; # for PRIMME
#my $maxOverlap=1E-5;

# lengths of the system
my @Ns=(20,24,30);

my $J=-1;
my @Gs=(-0.2,-0.3,-0.4,-0.5);
my $h=-1;

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

#foreach $s (@sValues){$s=-$s;} # for Pos!!

foreach $N (@Ns){
# Directory for results
	my $resultsDir="${basicDir}/results";
	my $mpsDir="${basicDir}/mps";

	print MYFILE "mkdir -p ${resultsDir} \n";
	print MYFILE "mkdir -p ${mpsDir} \n";
## LOOp over s values!
	foreach $g (@Gs){

	    my $suffix="L${N}_J${J}_g${g}_h${h}";

	    my $outfile="${resultsDir}/spec_${suffix}";

	    my $JOBNAME="iS.g${g}.N${N}.D${D0}";
    
	    # Decide if I am still (already) running
	    my $alreadyRunning=isRunning($JOBNAME);
	    if(! $alreadyRunning ){
		my $initmpsDir="${mpsDir}"; 
# ./isingSpec -L=24 -J=-1 -g=-0.2 h=-1 -D=200 -output=/ptmp/mpq/banulsm/melting/dataIsing/spec_L24_Jm1_gm02_hm1.txt -mpsdir=/ptmp/mpq/banulsm/melting/dataIsing/mps


		my $args=" -L=${N} -J=${J} -g=${g} -h=${h} -D=${D0} -output=${outfile} -mpsdir=${mpsDir}";
		print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
		#print MYFILE "sleep 1\n";
	    }
	    else{
#		print "${JOBNAME} already running!\n";
	    }
	}
	print MYFILE "cleanjobs \n";
#	print MYFILE "reLaunchJobs.pl ${jobsdir} 1 1\n";
}
#    print MYFILE "sleep 30\n";
#print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 200 1800 &\n";


close(MYFILE);

system("chmod a+x $scriptFile ");


