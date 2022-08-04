#!/usr/bin/perl -w
use strict;

# Script to prepare a series of jobs using the fhEvol program

# Name of the script file to be created
my $scriptFile="lanzaFHevol.sh";

my $EXEC="./fhEvol";

my $configFile="../config/randomMPO.conf";

my $basicDir="/ptmp/mpq/banulsm/fermihubbard/evol";

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

my $L=20;
my $t=1.;
my @Us=(1.,3.,5.);

#my $allSinglets;
#for my $l (1..$L){
#    $allSinglets=$allsinglets " " $l;
#}
my @configs=("0 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19","0 1 2 3 4 5 6 7 9 10 12 13 14 15 16 17 18 19",
	     "0 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19","0 1 2 3 4 5 6 8 9 11 12 14 15 16 17 18 19",
	     "0 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 19","0 1 2 3 4 6 7 9 10 12 13 15 16 17 18 19");
my @labels=("Ns1","Ns2","Ns2edges","Ns3","Ns3edges","Ns4");

my $app=0;

#$basicDir="${basicDir}/L${L}/";

# Directory for results
my $resultsDir="${basicDir}/results";
my $mpsDir="${basicDir}/MPS";
# Frequency to save intermediate results
my $rateSave=1;

my $delta=0.02;
my $time=10;


my $MSUB="msub_slurm";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

print MYFILE "mkdir -p ${resultsDir} \n";
print MYFILE "mkdir -p ${mpsDir} \n";

my @Ds=(40,60,80);
my @pars=(" -P 4"," -P 6"," -P 10");
my @mems=(" -mem 12G"," -mem 18G"," -mem 30G"); 
#@Ns=(20);

#./fhEvol ../config/randomMPO.conf -L=20 -t=1 -U=3 -D=40 -posSinglet="0 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19" -stepBlock=1 -dt=0.01 -time=10 -output=testFH/testEvol_Ns19.txt -append=0


foreach my $U (@Us){
    my $cnt=0;
    foreach my $paramSet (@configs){
	my $strL=$labels[$cnt];
	
	#print MYFILE "mkdir -p ${resultsDir}/${strL}/results \n";
	
	my $cntD=0;
	foreach my $D (@Ds){
	    my $par=$pars[$cntD];
	    my $mem=$mems[$cntD];
	    
	    my $locpar=$par;
	    my $locmem=$mem;
	    
	    my $suffix="_L${L}t${t}U${U}D${D}_${strL}";
	    # Name of the output file
	    my $outfile="${resultsDir}/data_${suffix}";
	    my $JOBNAME="fhEvol.${strL}.L${L}.U${U}.D${D}";
	    # TODO Check that the job is not yet running!!
	    my $alreadyRunning=isRunning($JOBNAME);
	    if(! $alreadyRunning ){
		my $args=" -L=${L} -t=${t} -U=${U} -D=${D} -dt=${delta} -time=${time} -posSinglet=\\\"${paramSet}\\\" -stepBlock=${rateSave} -output=${outfile} -append=${app}";
		print MYFILE "${MSUB} ${locpar} ${locmem} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
		print MYFILE "sleep 1\n";
	    }
	    else{
		print "Not launching $JOBNAME for it exists\n";
	    }
	    $cntD=$cntD+1;
	}
	$cnt=$cnt+1;
    }
}

close(MYFILE);

system("chmod a+x $scriptFile ");


