#!/usr/bin/perl


$EXEC="./transLCopSv";



# The ones in the thermalization paper
my $JP=1.;
my $gP=-1.05;
my $hP=0.5;
#my $gP=0.9045;
#my $hP=0.809;
## Params in Frank's paper
#my $gP=1.4;
#my $hP=0.9045;

# Integrable model
#my $JP=1.; my $hP=0.;
#my $gP=.5; # .5; 1.; 1.5
#my $gP=1.;
#my $gP=1.5;


#my $dt=0.01;
my $dt=0.1; # As in Hastings' paper 2014
#my $rate=10;

#my $Lmax=500;

$configfile="/u/banulsm/NewSVN/mpsdyn/src/config/transverseLC.conf";

$basicdir="/ptmp/mpq/banulsm/TransverseLCop/Ising/";#restartEopt/";
$basictmpdir="/ptmp/mpq/banulsm/TransverseLCop/IsingTmp/";#restartEopt/";

#$basicdir="/ptmp/mpq/banulsm/TransverseLCop/Ising_withE/";
#$basictmpdir="/ptmp/mpq/banulsm/TransverseLCop/IsingTmp_withE/";
$computeAll=1;
$computeE=0;


my $Jlab=sprintf("%.2g",$JP); #int(-log($gamma)/log(10));
my $glab=sprintf("%.2g",$gP); #int(-log($gamma)/log(10));
my $hlab=sprintf("%.2g",$hP); #int(-log($gamma)/log(10));
my $dtlab=sprintf("%.2g",$dt); #int(-log($dt)/log(10));

#my @initStNames=("X","Y","Z");

#$initStName=$initStNames[$initSt-1];


$basicdir="${basicdir}/g${glab}/h${hlab}";
$basictmpdir="${basictmpdir}/g${glab}/h${hlab}";


$MSUB="msub_slurm";

my $scriptfile="lanzaTransLCop.sh";
open MYFILE, '>', $scriptfile or die "Cannot open ${scriptfile} $!";
print MYFILE "#!/bin/bash\n\n";
#print MYFILE "mkdir -p ${outputdir} \n";

#my $Ls=2;

#my @Ds=(20,40,60,80); #100,200,300,400,500,1000);
#my @Ds=(80,100,200);
my @Ds=(40,80,100);
my @Ds=(80,100,200);
#my @Ds=(500,1000);
#my @Ds=(1000);
#my @Ds=(40);

my $maxT=20;
my $M=$maxT/$dt;

my $savingfreq=10;

my $append = 0;

my $long = "";

$outputdir="${basicdir}";
$tmpdir="${basictmpdir}";
foreach $Dcut (@Ds){
    if($Dcut>=100){ 
	$parallel=40;
	if($Dcut>=200){
	    $parallel=48;
	#    $long="-l";
	}
#	$savingfreq=50;
    }
    else{    
	if($Dcut>=80){ 
	    $parallel=20;
	    #	    $savingfreq=200;
	}
	else{
	    $parallel=10;
	    #	    $savingfreq=1000;
	}
    }
    if($parallel>1) {
	$parallel="-P $parallel ";
    }
    else{
	$parallel="";
    }
	
    print MYFILE "mkdir -p $outputdir \n";
    print MYFILE "mkdir -p $tmpdir \n";
	
    $NAME="tlcO${jobStr}_g${glab}h${hlab}_dt${dt}_D${Dcut}";
    my $filename="corr_J${Jlab}g${glab}h${hlab}dt${dtlab}_D${Dcut}";
    my $outfile="${outputdir}/${filename}";
    my $args=" ${configfile} -J=${JP} -g=${gP} -h=${hP} -M=${M} -outputfile=${outfile} -delta=${dt} -D=${Dcut} -savingfreq=${savingfreq} -tmpdir=${tmpdir} -append=${append} -computeAll=${computeAll} -computeEnergyCorr=${computeE} ";
    print MYFILE "$MSUB -N $NAME ${long} ${parallel} -- $EXEC ${args} \n";
    print MYFILE "sleep 1\n";
    #    print MYFILE "$EXEC ${args} \n";
}
close(MYFILE);
system("chmod a+x ${scriptfile} ");


