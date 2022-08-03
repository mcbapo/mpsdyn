#!/usr/bin/perl
use POSIX;
# Script to prepare a series of jobs using the lowH2 program

# Name of the script file to be created
my $scriptFile="lanzaLowVar.sh";

my $EXEC="./lowVar";

my $MSUB="msub_slurm";

my $configFile="../config/randomMPO.conf";

my $mpsdirBas="/ptmp/mpq/banulsm/DoS/Cheby/Ising";


# Hamiltonian parameters
my @hamilPars=(" -J=1 -g=-1.05 -h=0.5 "," -J=1 -g=0.905 -h=0.809 "," -J=1 -g=0.8 -h=0 "," -J=1 -g=1 -h=0 ");
my @idStrings=("J1_gm105_h05","J1_g0905_h0809","J1_g08_h0","J1_g1_h0"); # and an integrable one

my $tol=1E-6; # tolreance for convergence

#my @Ns=(12,20,30,40);
#my @Ns=(30,50);
my @Ns=(80,100);
my $D=30; 
my $Dinit=30; 

open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";

my $cnt=0;
foreach $paramSet (@hamilPars){
    my $strL=$idStrings[$cnt];
    foreach $N (@Ns){
	my $Dn=$D>2**($N/2)?2**($N/2):$D;
	my $gendir="${mpsdirBas}/${strL}/N${N}";
	my $mpsdir="${gendir}/MPS";
	my $outputfile="${gendir}/summary_L${N}_D${Dn}";
	system("mkdir -p ${gendir}");
	#print MYFILE "mkdir -p ${gendir}";
	print MYFILE "mkdir -p ${mpsdir}\n";
	my $listFile="${gendir}/list_${strL}_N${N}_${D}"; 
	open MYLIST, '>', $listFile or die "Cannot open file $listFile $!";
	for (my $E0=-ceil(13*$N/10);$E0<=ceil(13*$N/10);$E0=$E0+1){
	    my $penH=10;
	    if(abs($E0)>1){$penH=$penH*abs($E0);}
	    my $strE0=abs($E0);
	    if($E0<0){
		$strE0="m${strE0}";
	    }
	    my $suffix="N${N}_${strL}_E${strE0}_D${Dn}";
	    my $mpsfile="mps_${suffix}";
	    my $Dninit=$Dinit>2**($N/2)?2**($N/2):$Dinit;
	    my $initsuffix="N${N}_${strL}_E${strE0}_D${Dninit}";
	    my $initmpsfile="mps_${initsuffix}";
	    my $JOBNAME="lowV.${cnt}.N${N}.D${Dn}.E${E0}";
	    my $args=" -L=${N} ${paramSet} -D=${Dn} -initmpsfile=${mpsdir}/${initmpsfile} -mpsfile=${mpsdir}/${mpsfile} -outputfile=${outputfile} -tol=${tol} -E=${E0} -penH=${penH} ";
	    # Check that the job is not yet running!! Add -mem ${MEM}M 
	    #my $MEM=int(((2*$N+1)*$Nb*$D*$D*4+(2*$N+1)*2*$Nb*$Nb*$Nb*$Nb)*16*1E-6)*2+200;
	    # Now write the istruction to launch the job to the script
	    print MYFILE "${MSUB} ${par} -N ${JOBNAME} -- ${EXEC} ${configFile} ${args} \n";
	    print MYFILE "sleep 1\n";
	    print MYLIST "${mpsdir}/${mpsfile}\n";
	}
	close(MYLIST);
    }
    $cnt=$cnt+1;
}

close(MYFILE);

system("chmod a+x $scriptFile ");


