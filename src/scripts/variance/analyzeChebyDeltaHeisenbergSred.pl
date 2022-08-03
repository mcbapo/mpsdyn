#!/usr/bin/perl

# Script to prepare a series of jobs using the minH2 program

# Name of the script file to be created
my $scriptFile="analyzeChebyDeltaHeis.sh";

my $EXEC="./analyzeChebH2";

my $configFile="../config/randomMPO.conf";

##my $basicDir="/ptmp/mpq/banulsm/variance/Ising/ChebyDelta/";
my @models=("Ising","Heisenberg","XY");
#my $model=0; #"Ising";
my $model=1; #"Heisenberg";
#my $model=2; #"XY";

my $modelStr=$models[$model];
my $basicDir="/ptmp/mpq/banulsm/variance/2/${modelStr}/ChebyDelta/Sred";

# Key for initial states
sub initStateName{
    my ($initSt)=@_;
    if(abs($initSt)==1){
	if($initSt>0){	return "Xplus";}
	else{	return "Xminus";}
    }
    else{
	if(abs($initSt)==2){
	    if($initSt>0){return "Yplus";}
	    else{return "Yminus";}
	}
	else{
	    if(abs($initSt)==3){
		if($initSt>0){return "Zplus";}
		else{return "Zminus";}
	    }
	    else{
		if(abs($initSt)==4){
		    if($initSt>0){return "Xst";}
		    else{return "Xstminus";}
		}
		else{
		    if(abs($initSt)==5){
			if($initSt>0){return "Yst";}
			else{return "Ystminus";}
		    }
		    else{
			if(abs($initSt)==6){
			    if($initSt>0){return "Zst";}
			    else{return "Zstminus";}
			}
			else{
			    if($initSt==7){
				return "rndTI";
			    }
			    else{
				if($initSt==8){
				    return "rndTIreal";
				}
				else{
				    if($initSt==9){
					return "YplusRe";
				    }
				    else{
					if(abs($initSt)==10){
					    if($initSt>0){return "Zst2";}
					    else{return "Zst2minus";}
					}
					else{
					    print "Unknown initial state ${initSt}!!\n";exit(1);
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }

};

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

# Check last step reached
sub lastStepComputed{
    my ($outputfile,$model)=@_;
    if($model==2){ # XY --- SOMETHING NOT GOOD! (7+1?)
	$colM=5+1;
    }
    else{
	$colM=5+1; # from 0
    }
    # Read the last line of the file
    my @lines=split "\n",`cut -f $colM $outputfile`;
    chomp(@lines);
    #print "Read @lines; last one is @lines[-1] \n";
    return @lines[-1];
}

# Check last step reached in Sblock!
sub lastStepComputedS{
    my ($outputfileS,$model)=@_;
    $colM=5+1+1; # because of the Dcut
    # Read the last line of the file
    my @lines=split "\n",`cut -f $colM $outputfileS`;
    chomp(@lines);
    #print "Read @lines; last one is @lines[-1] \n";
    return @lines[-1];
}

if($model==0){ # Ising case
    print "Case Ising \n";
    # Hamiltonian parameters
    
    @hamilPars=(" -J=1 -g=-1.05 -h=0.5 ");
    @idStrings=("J1_gm105_h05");

    @Ns=(20,40,60,50,80,100);@Ds=(200,300,400,500,1000);
    @pars=("-P 10","-P 10","-P 20","-P 20","-P 20");
    $initSt=+2; # Yplus
}
else{
    if($model==1){ # Heisenberg
	print "Case Heisenberg \n";
	@hamilPars=(" -Jx=1.1 -Jy=-1 -Jz=.9 -h=0.6 -hx=0. -hy=0. ");#," -Jx=1 -Jy=1 -Jz=1 -h=0.5 -hx=0. -hy=0. ");
	@idStrings=("Jx1o1_Jym1_Jz0o9_h0o6"); #,"Jx1_Jy1_Jz1_h05");
	@labParH=("3");
	#@Ns=(20,40,60,50,80,100);@Ds=(200,300,400,500,1000);
	#@pars=("-P 10","-P 10","-P 20","-P 20","-P 20");
#	@Ns=(20,40,60,50,80,100);@Ds=(400,500,1000);
	@Ns=(20,40,60,80,100);@Ds=(300,400,500,1000);
	#@pars=("-P 20","-P 20","-P 20");
	#@mems=(" -mem 60G "," -mem 60G "," -mem 60G ") ; # -mem ${MEM}G
	@pars=("-P 20","-P 20","-P 20","-P 20");
	@mems=("-mem 60G","-mem 60G","-mem 60G","-mem 60G");
	$initSt=+10; # staggered Z period 2 for HEis
    }
    else{ # XY
	@hamilPars=(" -Jx=1. -Jy=1. -Jz=0. -hx=0.8 -hy=0 -h=0.5 ");
	@idStrings=("Jx1_Jy1_hx0o8_hy0_hz0o5");
	@Ns=(20,40,60,50,80,100);@Ds=(200,300,400,500,1000);
	@pars=("-P 10","-P 10","-P 20","-P 20","-P 20");
	$initSt=+10; # staggered Z period 2 for HEis
    }
}
my $app=0;
my $initStStr=initStateName($initSt);
my $randSeed=117; # to start always the same
my $avrg=1;
my $Lcmax=6;
my $Lb=6;

$basicDir="${basicDir}/${initStStr}/";

# Directory for results
my $resultsDir="${basicDir}/results";
# Frequency to save intermediate results
my $rateSave=100;



my $delta=0.02;


my $MSUB="msub_slurm";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

print MYFILE "mkdir -p ${resultsDir} \n";
#print MYFILE "mkdir -p ${mpsDir} \n";

#my @Ds=(1000);
my $cnt=0;
@Ds=(300,400,500,1000);
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
    print MYFILE "mkdir -p ${resultsDir}/${strL}/results \n";
    # TODO: Set a file for energies (Emin Emax) and check for it
    my $mpsdir="${resultsDir}/${strL}/MPS";
    print MYFILE "mkdir -p ${mpsdir} \n";
    my $cntD=0;
    foreach $D (@Ds){
	my $par=$pars[$cntD];
	my $mem=$mems[$cntD];
	foreach $N (@Ns){
	    my @Ms=(2*$N,int(2*sqrt($N)),5*$N,10*$N,int(5*sqrt($N)),int(10*sqrt($N)),20*$N,int(0.1*($N**2)),int($N*log($N)));
#	    my @Ms=(0);
	    #my @Ms=(int(0.1*($N**2)));
	    #my @Ms=(int($N*log($N)));
	    my $locpar=$par;
	    my $locmem=$mem;
	    foreach $M (@Ms){
#		if($N==100&&$D==1000&&$M>=1000){$locpar=" -P 40 ";$locmem=" -mem 60G ";}
#		else{$locpar=$par;$locmem=$mem;}
		my $suffix="_N${N}M${M}D${D}";
		my $mpsfileroot="mps_${suffix}";
		# Name of the output file
		my $outfileCorr="${resultsDir}/${strL}/results/finalEcorr_${suffix}";
		my $outfileSblock="${resultsDir}/${strL}/results/finalSblock_${suffix}";
		my $JOBNAME="chD.analyze.${modelStr}.N${N}.D${D}.M${M}.par$labParH[${cnt}]";
		# TODO Check that the job is not yet running!!
		my $alreadyRunning=isRunning($JOBNAME);
		if(! $alreadyRunning ){
#		    my $lastM=lastStepComputed($outfileCorr,$model);
		    my $lastM=lastStepComputedS($outfileSblock,$model);
		    if($lastM>$M && $M!=0 ){
			print "# Not launching $JOBNAME for it is done\n";
		    }
		    else{
			my $JOBdata="chD.${modelStr}.N${N}.D${D}.M${M}.par$labParH[${cnt}]";
			my $runningData=isRunning($JOBdata);
			if(! $runningData){
			    my $args=" -ising=0 -L=${N} ${paramSet} -D=${D} -delta=${delta} -initSt=${initSt} -outputSblock=${outfileSblock} -outputCorr=${outfileCorr}  -M=${M} -mpsdir=${mpsdir} -mpsfileroot=${mpsfileroot} -append=${app} ";
			    print MYFILE "${MSUB} ${locpar} ${locmem} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
			    print MYFILE "sleep 1\n";		    
			}
			else{
			    print MYFILE "# Not launching $JOBNAME because the data are still being computed\n";
			}
		    }
		}
		else{
		    print "# Not launching $JOBNAME for it exists\n";
		}
	    }
	}
	$cntD=$cntD+1;
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


