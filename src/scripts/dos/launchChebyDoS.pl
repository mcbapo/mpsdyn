#!/usr/bin/perl

# Script to prepare a series of jobs using the lowH2 program

# Name of the script file to be created
my $scriptFile="lanzaChebyDoS.sh";

my $EXEC="./dosCh";

my $MSUB="msub_slurm";

my $configFile="../config/randomMPO.conf";

my $mpsdirBas="/ptmp/mpq/banulsm/DoS/Cheby/Ising";

my $delta=0.1;

# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 "," -J=1 -g=1 -h=0 ");
my @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0","J1_g1_h0"); # and an integrable one

my @Ns=(30,50);
my @Ds=(20,40,60,80,100,120);
my @pars=("-P 8","-P 12","-P 12","-P 20","-P 20","-P 20");

#my @Ns=(80,100);
#my @Ns=(50);
#my @Ds=(100,120);
#my @pars=("-P 20","-P 20","-P 20","-P 20");
my $tol=1E-6; # tolerance for convergence


my $Dprobe=20;

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

my $M=200;
my $Mw=200;
my $freqProbe=$M;
my $cntD=0;
foreach $D (@Ds){
    my $par=$pars[$cntD];
    my $cnt=0;
    foreach $paramSet (@hamilPars){
	my $strL=$idStrings[$cnt];
	    foreach $N (@Ns){
		my $JOBNAME="chDoS.${cnt}.N${N}.D${D}";
		my $gendir="${mpsdirBas}/${strL}/N${N}";
		my $mpsdir="${gendir}/dosMPS";
		print MYFILE "mkdir -p ${mpsdir}\n";
		my $mpsprobelist="${gendir}/list_${strL}_N${N}_${Dprobe}"; 
		my $outputfile="${gendir}/test_${strL}_L${N}_D${D}_M${M}_Dpr${Dprobe}";
		my $mpsfile="${mpsdir}/dos_${strL}_L${N}_D${D}_M${M}";
		my $mpsfileW="${mpsdir}/xW_${strL}_L${N}_D${D}_M${Mw}";
		my $args=" -L=${N} ${paramSet} -D=${D} -outputfile=${outputfile} -mpsfile=${mpsfile} -mpsfileW=${mpsfileW} -tol=${tol} -freqprobe=${freqProbe} -M=${M} -Mw=${Mw} -delta=${delta}";
#	    $args="$args -mpsprobelist=${mpsprobelist} ";
		print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} ;\n";
		print MYFILE "sleep 1\n";
	    }
	$cnt=$cnt+1;
    }
    $cntD=$cntD+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


