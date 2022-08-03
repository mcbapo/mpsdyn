#!/usr/bin/perl

# Script to prepare a series of jobs using the minH2 program

# Name of the script file to be created
my $scriptFile="lanzaChebyDelta.sh";

my $EXEC="./chebH2";

my $configFile="../config/randomMPO.conf";

##my $basicDir="/ptmp/mpq/banulsm/variance/Ising/ChebyDelta/";
my @models=("Ising","Heisenberg","XY");
my $model=0; #"Ising";
#my $model=1; #"Heisenberg";
#my $model=2; #"XY";

my $modelStr=$models[$model];
my $basicDir="/ptmp/mpq/banulsm/variance/2/${modelStr}/ChebyDelta";
my $extraInitStArgs="";

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
					    if($initSt==11&&$model==0){
						return "stepE";
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
    if($model==2){ # XY
	$colM=7+1;
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

my $initStStr="";
if($model==0){ # Ising case
    print "Case Ising \n";
    # Hamiltonian parameters
    @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 ");
    @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0"); # and an integrable one
    
    @hamilPars=(" -J=1 -g=-1.05 -h=0.5 ");
    @idStrings=("J1_gm105_h05");

    #my @hamilPars=(" 1 0.8 1"," 1 -.5 .5");
    #my @idStrings=("J1_g08_h1","J1_gm05_h05");
    
    #my @hamilPars=(" -J=1 -g=0.905 -h=0.809 ");#
    #my @idStrings=("J1_g0905_h0809");
    #my @hamilPars=(" 1 0.8 0");
    #my @idStrings=("J1_g08_h0"); # an integrable one

    # lengths of the system
    # @Ns=(20,50,100,200);
    # Ds and corresponding parallel options
    # @Ds=(60,100,200,300,400);
    # @pars=("-P 4","-P 10","-P 12","-P 14","-P 20"); # use parallel environment
    
    # short exact ones
    #    @Ns=(12,14,18,20);
    # Ds and corresponding parallel options
    #    @Ds=(500);
    #    @pars=("-P 14");
    #@Ns=(18,20);@Ds=(1024);@pars=("-P 20")
    #@Ns=(22);@Ds=(2048);@pars=("-P 40")

    @Ns=(20,40,60,50,80,100);@Ds=(200,300,400,500,1000);
    @pars=("-P 10","-P 10","-P 20","-P 20","-P 20");
    #@Ns=(14,16,18,22,20,24,26,28,30);


    $initSt=+11;
    $thetaL=0.394;
    $thetaR=0.443;
    $phiL=0;
    $phiR=0;
    $initStStr=initStateName($initSt);
    $initStStr="${initStStr}_thL0o39_thR0o44";
    $extraInitStArgs=" -thetaL=${thetaL} -phiL=${phiL} -thetaR=${thetaR} -phiR=${phiR}";

}
else{
    if($model==1){ # Heisenberg
	print "Case Heisenberg \n";
	#    @hamilPars=(" -Jx=1 -Jy=1 -Jz=1 -h=0. "," -Jx=1 -Jy=1 -Jz=1 -h=0.5 ");
	#    @idStrings=("Jx1_Jy1_Jz1_h0","Jx1_Jy1_Jz1_h05");
	@hamilPars=(" -Jx=1.1 -Jy=-1 -Jz=.9 -h=0. -hx=0. -hy=0. "," -Jx=1 -Jy=1 -Jz=1 -h=0.5 -hx=0. -hy=0. ");
	@idStrings=("Jx1o1_Jym1_Jz0o9_h0","Jx1_Jy1_Jz1_h05");
	# lengths of the system
	#@Ns=(20,50,100,200,300);
	#@Ds=(60,100,200,300,500);
	#@pars=("-P 4","-P 10","-P 12","-P 14","-P 20","-P 20"); # use parallel environment
	
	@Ns=(20,40,60,50,80,100);@Ds=(200,300,400,500,1000);
	@pars=("-P 10","-P 10","-P 20","-P 20","-P 20");
	$initSt=+10; # staggered Z period 2 for HEis
	# sth longer
	#@Ns=(20,30);@Ds=(1024);@pars=("-P 40");
    }
    else{ # XY
	@hamilPars=(" -Jx=1. -Jy=1. -Jz=0. -hx=0.8 -hy=0 -h=0.5 ");
	@idStrings=("Jx1_Jy1_hx0o8_hy0_hz0o5");
	@Ns=(20,40,60,50,80,100);@Ds=(200,300,400,500,1000);
	@pars=("-P 10","-P 10","-P 20","-P 20","-P 20");
	$initSt=+10; # staggered Z period 2 for HEis
    }
    $initStStr=initStateName($initSt);
}
my $app=0;
#my $initStStr=initStateName($initSt);
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


my $cnt=0;
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
    print MYFILE "mkdir -p ${resultsDir}/${strL}/results \n";
    # TODO: Set a file for energies (Emin Emax) and check for it
    my $mpsdir="${resultsDir}/${strL}/MPS";
    print MYFILE "mkdir -p ${mpsdir} \n";
    my $cntD=0;
    foreach $D (@Ds){
	my $par=$pars[$cntD];
	foreach $N (@Ns){
	    #my @Ms=(2*$N,int(2*sqrt($N)),5*$N,10*$N,int(5*sqrt($N)),int(10*sqrt($N)));
	    my @Ms=(2*$N,int(2*sqrt($N)),5*$N,10*$N,int(5*sqrt($N)),int(10*sqrt($N)),20*$N,50*$N,80*$N,100*$N); # just for Ising!!!
	    #my @Ms=(int(.1*($N)**2));
	    #my @Ms=(int($N*log($N)));
	    my $locpar=$par;
	    foreach $M (@Ms){
		if($N==100&&$D==1000&&$M>=1000){$locpar=" -P 40 ";}
		else{$locpar=$par;}
		my $suffix="_N${N}M${M}D${D}";
		my $mpsfileroot="mps_${suffix}";
		# Name of the output file
		my $outfile="${resultsDir}/${strL}/results/data_${suffix}";
		my $outfileS="${resultsDir}/${strL}/results/dataSchmidt_${suffix}";
		my $outfileEd=" -outputEdist=${resultsDir}/${strL}/results/dataE_${suffix} ";
		my $JOBNAME="chD.${modelStr}.N${N}.D${D}.M${M}.par${cnt}";
		if($initSt==11&&$model==0){
		    $JOBNAME="chD.${modelStr}.N${N}.D${D}.M${M}.par${cnt}.st11";
		}
		# TODO Check that the job is not yet running!!
		my $alreadyRunning=isRunning($JOBNAME);
		#if($alreadyRunning && $model==1 && $cnt==0){ # check if the running one is the proper param
		#    my $jobScript=`squeue --noheader -n ${JOBNAME} -o "%o"`;
		#    my $jobdata=`grep srun ${jobScript}`;
		#    $jobdata =~ /.*output=(.*) -outputS.*/ ;
		#    my $usedFile=$1;
		#    my $rightFile= $usedFile eq $outfile ;
		#    print "cadena $jobdata \n match $usedFile outfile $outfile iguales? $rightFile \n";
		#    if(!$rightFile){
		#		print "NOT THE SAME: mark as not running\n";
		#		$alreadyRunning=0;
		#	    }
		#	}
		if(! $alreadyRunning ){
		    my $lastM=lastStepComputed($outfile,$model);
		    if($lastM<($M-1)){
			my $args=" -L=${N} ${paramSet} -D=${D} -delta=${delta} -output=${outfile} -outputS=${outfileS} ${outfileEd} -M=${M} -initSt=${initSt} ${extraInitStArgs} -checkRDM=1 -Lcmax=$Lcmax -average=$avrg -Lb=$Lb -seed=$randSeed -mpsdir=${mpsdir} -mpsfileroot=${mpsfileroot} -savingFreq=${rateSave}";
			if($model!=0){$args="$args -ising=0 ";}
			if(!$alreadyRunning){
			    print MYFILE "${MSUB} ${locpar} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
			    print MYFILE "sleep 1\n";
			}
			else{
			    print MYFILE "CHECK: ${MSUB} ${locpar} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
			}
		    }
		    else{
			print "Not launching $JOBNAME because it already reached M=$lastM \n";
		    }
		}
		else{
		    print "Not launching $JOBNAME for it exists\n";
		}
	    }
	}
	$cntD=$cntD+1;
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


