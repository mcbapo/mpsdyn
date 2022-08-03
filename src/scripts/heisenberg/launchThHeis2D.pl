#!/usr/bin/perl

# Script to prepare a series of jobs using the isingSpec program

# Name of the script file to be created
my $scriptFile="lanzaHeis2D.sh";

my $EXEC="/u/banulsm/NewSVN/mpsdyn/src/bin/thHeis2";

my $basicDir="/ptmp/mpq/banulsm/heis2D"; #/tests";
my $jobsdir="/u/banulsm/jobsToRun";

#my $D0=20; # initial D
my $D0=40; # initial D
my $par="-P 10";
my $incrD=20;
my $maxD=100;
my $app=0;
my $tol=1E-8; # tolerance 

my $beta=6;
my $delta=0.02;
my $rate=1;

# lengths of the system
#my $Lx=8;
#my $Ly=6;
my $Lx=5;
my $Ly=4;

my $nL=1; # nr of spins per site

my $g=1.;
my @Mus=(0.); #(1.,0.);


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

foreach my $mu (@Mus){
    my $Jx=-2*$g;
    my $Jy=-2*$g;
    my $Jz=0;
    my $h=-$mu;

# Directory for results
    my $resultsDir="${basicDir}/Lx${Lx}_Ly${Ly}/results";
    #my $mpsDir="${basicDir}/mps";

    print MYFILE "mkdir -p ${resultsDir} \n";
    #print MYFILE "mkdir -p ${mpsDir} \n";
    ## LOOp over D values!
    for (my $D=$D0;$D<=$maxD;$D=$D+$incrD){

	my $suffix="Lx${Lx}_Ly${Ly}_g${g}_mu${mu}_D${D}";
	print "mu=${mu}, suffix=${suffix}\n";

	
	my $outfile="${resultsDir}/results_${suffix}";
	
	my $JOBNAME="thH2.g${g}.mu${mu}.${Lx}x${Ly}.D${D}";
	
	# Decide if I am still (already) running
	my $alreadyRunning=isRunning($JOBNAME);
	if(! $alreadyRunning ){
	    #my $initmpsDir="${mpsDir}"; 
	    # ./thHeis2 4 3 1 -.5 -.5 0. 0. 20 6 .02 1 kk

	    my $args=" ${Lx} ${Ly} ${nL} ${Jx} ${Jy} ${Jz} ${h} ${D} ${beta} ${delta} ${rate} ${outfile}";
	    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${args} \n";
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


