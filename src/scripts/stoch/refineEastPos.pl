#!/usr/bin/perl

# Script to prepare a series of jobs using the specSto program. This time, values of s are read from a file
# and relaunches only the positive s values, for a range of s values, to refine previous results

# Name of the script file to be created
my $scriptFile="refineEastPos.sh";

my $EXEC="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/bin/specSto";

my $configFile="/afs/ipp/u/banulsm/NewSVN/mpsdyn/src/config/stochspec.conf";

my $basicDir="/ptmp/mpq/banulsm/stoch/basic";
my $jobsdir="/afs/ipp/u/banulsm/jobsToRun";

my $D0=80; # initial D
my $par="-P 10";
my $incrD=20;
my $maxD=100;
my $app=0;
my $tol=1E-12; # tolerance 
my $eigtol=0;#1E-12; # tolerance for the eigensolver 
my $nlev=3; # how many levels to extract
my $saveFreq=10;
my $refine=" -refine=1 "; # to rerun with increased tolerance
my $knr=1; # for PRIMME

# lengths of the system
#my @Ns=(20,40,60,80,100,200,300,400);
#my @Ns=(20,40,60,80,100);
#my @Ns=(20,40,60,79,80,81,99,100,101,199,200,201,300,400);

# To explore small values of c
#my @Cs=(.5,.2,.1); #
my @Cs=(.1); #,.02);
#my @Cs=(.2);

my @Ns=(20,40,60,80,100,200,300,400);
my @Ns=(200,300,400); #(20,40,60,80,100); #300,400); #60,80,100,200);
my $slimlow=1.625e-6; # only for s between this and slim, read from svals_c0.1
my $slim=9e-6; # only for s between slimlow and this, read from svals_c0.1


my $reorder = " -reorder=1 "; # to force reordering of the levels
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
#	    $resultsDir="${basicDir}/c${c}/N${N}/results2";
	    $mpsDir="${basicDir}/c${c}/N${N}/MPS2";
	}
	print MYFILE "mkdir -p ${resultsDir} \n";
	print MYFILE "mkdir -p ${mpsDir} \n";
## Check all the values of s existing in some file
	my $allVals=`ls ${mpsDir}/MPS_L${N}_s*D${D0}*0.dat |sed  's/^.*MPS_L${N}_s\\(.*\\)_c.*/\\1/g' `;
	#print "ls ${mpsDir}/MPS_L${N}_s*D${D0}*0.dat |sed  's/^.*MPS_L${N}_s\\(.*\\)_c.*/\\1/g'  ";exit 1;
	@sValues=split('\n',$allVals); 

	## LOOp over s values in a file!
	#my $fileSvals="../scripts/stoch/svals_c${c}_small";
	#my $fileSvals="../scripts/stoch/svals_c${c}";
	#$fileSvals =~ s/0\./0o/g ;
	#my $sVals=`cat ${fileSvals}`;

	#@sValues=split('\n',$sVals); 

# print "Values to launch for c=$c are in file $fileSvals: ",$sVals,"\n Split:",@sValues,"\n";
	#print "Values to launch for c=$c are:",@sValues,"\n";
	#exit 1 ;
	foreach $s (@sValues){
	    if($s>$slimlow && $s<=$slim){
		#if($s!=0.01){
		my $sStr=sprintf("%g",$s);
		#my $suffix="s${s}_c${c}_$N${N}_D${D0}";
		my $suffix="s${sStr}_c${c}_N${N}_D${D0}";
		my $suffixLast="s${sStr}_c${c}_N${N}_D${maxD}";
		
		my $outfile="${resultsDir}/data_${suffix}";
		my $polfile="${resultsDir}/polX_${suffix}";
		
		#my $JOBNAME="sto.c${c}.N${N}.D${D0}.s${s}";
		my $JOBNAME="E.c${c}.N${N}.D${D0}.s${sStr}";
		
		# TODO!!! Decide if I am still (already) running
		my $initmpsDir="${mpsDir}"; 

		#my $mpsFiles=`ls -ltr ${mpsDir}/MPS_L${N}_s${s}_c${c}_D${D0}_l1*`;
		my $mpsFiles=`ls -ltr ${mpsDir}/MPS_L${N}_s${s}_c${c}_D*0.dat`;
		my $checkFiles=`cut -f 1,2 ${resultsDir}/data_${suffixLast}`;
		print "======\n${outfile} \n${checkFiles}\n";
		print "$mpsFiles \n";
		#print "rm ${mpsDir}/MPS_L${N}_s${s}_c${c}_D* \n";

		my $args=" -L=${N} -D=${D0} -s=${s} -c=${c} -nlevel=${nlev} -outputfile=${outfile} -polfile=${polfile} -mpsdir=${mpsDir} -initmpsdir=${initmpsDir} -tol=${tol} -tmpSavingFreq=${saveFreq} -firstN=1 ${reorder} -noise=${noise}  ${refine} -jobsdir=${jobsdir} -maxD=${maxD} -incrD=${incrD} -knr=${knr} -eigtol=${eigtol}";
		print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
		print "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
		#print MYFILE "sleep 1\n";
		print "Launching s=$s \n";
	    }
	    else{
		#print "Discarding s=$s \n";
	    }
	}
	print MYFILE "cleanjobs\n";
	print MYFILE "reLaunchJobs.pl ${jobsdir} 1 1\n";
    }
}
print MYFILE "nohup reLaunchJobs.pl ${jobsdir} 200 1800 &\n";


close(MYFILE);

system("chmod a+x $scriptFile ");


