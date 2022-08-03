#!/usr/bin/perl

# Script to prepare a series of jobs using the lowH2 program

# Name of the script file to be created
my $scriptFile="lanzaEvalVars.sh";

my $EXEC="./evalVar";

my $MSUB="msub_slurm";

my $configFile="../config/randomMPO.conf";

my $mpsdirBas="/ptmp/mpq/banulsm/variance/Ising/variational";

my $delta=0.1;

# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 "," -J=1 -g=1 -h=0 ");
my @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0","J1_g1_h0"); # and an integrable one

my $tol=1E-6; # tolerance for convergence

my @Ns=(20,40,50,80,100,200,300);

my $Dmin=20; 
my $Dmax=300;
my $stepD=20; 

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

my $par="-P 10";

my $cnt=0;
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
    foreach $N (@Ns){
	my $JOBNAME="evalVar.${cnt}.N${N}";
	my $gendir="${mpsdirBas}/${strL}/N${N}";
	my $listdir="${gendir}";
	my $mpsdir="${gendir}/MPS";
	my $mpsprobelist="${listdir}/list_${strL}_N${N}"; 
	system("mkdir -p ${listdir} ");
	# have to fill the list file with the files I find in the directory
	opendir DIR, $mpsdir or die "$0: can't opendir $mpsdir: $!\n";
	@files=grep { !/^\./ } readdir(DIR);
	closedir DIR;
	open LISTFILE, '>', $mpsprobelist or die "Cannot open file $mpsprobelist $!";
	foreach $f (@files){
	    # IMPORTANT!!! No space after the name!!!!
	    print LISTFILE "${mpsdir}/$f\n";
	    #my $resLs=`ls -l ${mpsprobedir}/$f` ;  
	    #print "\t ls -l $f: $resLs \n" 
	}
	close LISTFILE;
	print "Wrote list to $mpsprobelist \n";
	my $outputfile="${gendir}/var_${strL}_L${N}";
	my $args=" -L=${N} ${paramSet} -outputfile=${outputfile} -mpslist=${mpsprobelist} -append=0 ";
	print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} ;\n";
	print MYFILE "sleep 1\n";
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


