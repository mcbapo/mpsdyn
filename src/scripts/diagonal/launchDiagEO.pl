#!/usr/bin/perl

# Script to prepare a series of jobs using the minH2 program

# Name of the script file to be created
my $scriptFile="lanzaDiagEO.sh";

my $EXEC="./diagEO";

my $configFile="../config/dissipXY.conf"; #ignored!

my $basicDir="/ptmp/mpq/banulsm/diagonal/diagEO";

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

};


# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=1.05 -h=0.5 "); #," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 ");
my @idStrings=("J1_gm105_h05","J1_g105_h05"); #,"J1_g0905_h0809","J1_g08_h0"); # and an integrable one


#my @hamilPars=(" 1 0.8 0");
#my @idStrings=("J1_g08_h0"); # an integrable one

# lengths of the system
my @Ns=(20,30,40,50);
my @Ds=(80,100,120,200,240);
my @pars=("-P 10","-P 12","-P 14","-P 14","-P 20"); # use parallel environment

#my @Ns=(50);
#my @Ds=(300,400);
#my @pars=("-P 20","-P 20");

my @Ds=(240);
my @pars=("-P 20","-P 20");

my $app=0;
my $initSt=+1; # Xplus
#my $initSt=+2; # Yplus
#my $initSt=+9; # Yplus real part
my $initStStr=initStateName($initSt);
my $randSeed=117; # to start always the same
my $avrg=1;

$basicDir="${basicDir}/${initStStr}/";

# Directory for results
my $resultsDir="${basicDir}"; #/results";
#my $mpsDir="${basicDir}/MPS";

my $rate=50;
my $saveRate=50;
my $delta=0.01;
my $dt=0.01;
my $M=10000;

my $MSUB="msub_slurm";

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

print MYFILE "mkdir -p ${resultsDir} \n";
#print MYFILE "mkdir -p ${mpsDir} \n";


my $cnt=0;
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
    print MYFILE "mkdir -p ${resultsDir}/${strL}/results \n";
    print MYFILE "mkdir -p ${resultsDir}/${strL}/MPS \n";
    my $cntD=0;
    foreach $D (@Ds){
	my $par=$pars[$cntD];
	foreach $N (@Ns){
	    my $suffix="_N${N}D${D}";
	    my $outfile="${resultsDir}/${strL}/results/data_${suffix}";
	    my $mpsdir="${resultsDir}/${strL}/MPS";
	    my $mpsfile="mps_${suffix}";

	    my $JOBNAME="eo.N${N}.D${D}.${cnt}";
	    # Check that the job is not yet running!!
	    #my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
	    my $args=" -L=${N} ${paramSet} -D=${D} -outputfile=${outfile} -mpsdir=${mpsdir} -mpsfile=${mpsfile} -delta=${delta} -dt=${dt} -M=${M} -rate=${rate} -saveFreq=${saveRate} -initSt=${initSt} ";
	    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
	    print MYFILE "sleep 1\n";
	}
	$cntD=$cntD+1;
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


