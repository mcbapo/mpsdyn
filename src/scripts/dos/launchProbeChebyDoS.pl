#!/usr/bin/perl

# Script to prepare a series of jobs using the lowH2 program

# Name of the script file to be created
my $scriptFile="lanzaProbeChebyDoS.sh";

my $EXEC="./dosChprobe";

my $MSUB="msub_slurm";

my $configFile="../config/randomMPO.conf";

my $mpsdirBas="/ptmp/mpq/banulsm/DoS/Cheby/Ising";

my $delta=0.1;

# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 "," -J=1 -g=1 -h=0 ");
my @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0","J1_g1_h0"); # and an integrable one

my @Ds=(20,40,60,80,100,120);
my @pars=("-P 4","-P 10","-P 12","-P 20","-P 20","-P 20");
my $tol=1E-6; # tolerance for convergence

my @Ns=(30,50,80,100);
#my @Ns=(30,50);

my $Dprobe=20;
@Ns=(50);
my @Ds=(100,120);
my @pars=("-P 20","-P 20");

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
	    my $JOBNAME="prDoS.${cnt}.N${N}.D${D}";
	    my $gendir="${mpsdirBas}/${strL}/N${N}";
	    my $listdir="${gendir}/lists";
	    my $mpsdir="${gendir}/dosMPS";
	    my $mpsprobedir="${gendir}/MPS";
	    my $mpsprobelist="${listdir}/list_${strL}_N${N}_${Dprobe}"; 
	    system("mkdir -p ${listdir} ");
	    # have to fill the list file with the files I find in the directory
	    opendir DIR, $mpsprobedir or die "$0: can't opendir $mpsprobedir: $!\n";
	    my $strD="D${Dprobe}";
	    @files=grep { !/^\./ && /$strD/} readdir(DIR);
	    closedir DIR;
	    open LISTFILE, '>', $mpsprobelist or die "Cannot open file $mpsprobelist $!";
	    foreach $f (@files){
# IMPORTANT!!! No space after the name!!!!
		print LISTFILE "${mpsprobedir}/$f\n";
		#my $resLs=`ls -l ${mpsprobedir}/$f` ;  
		#print "\t ls -l $f: $resLs \n" 
	    }
	    close LISTFILE;
	    print "Wrote list to $mpsprobelist \n";
	    my $outputfile="${gendir}/test_${strL}_L${N}_D${D}_M${M}_Dpr${Dprobe}";
	    my $mpsfile="${mpsdir}/dos_${strL}_L${N}_D${D}_M${M}";
	    my $mpsfileW="${mpsdir}/xW_${strL}_L${N}_D${D}_M${Mw}";
	    my $args=" -L=${N} ${paramSet} -D=${D} -outputfile=${outputfile} -mpsfile=${mpsfile} -mpsfileW=${mpsfileW} -mpsprobelist=${mpsprobelist} -tol=${tol} -freqprobe=${freqProbe} -M=${M} -Mw=${Mw} -delta=${delta}";
	    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} ;\n";
	    print MYFILE "sleep 1\n";
	}
	$cnt=$cnt+1;
    }
    $cntD=$cntD+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


