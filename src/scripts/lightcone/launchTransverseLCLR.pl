#!/usr/bin/perl


$EXEC="./transLCLR";


my $initSt=1;

# The ones in the thermalization paper
my $JP=1.;
my $gP=-1.05;
my $hP=0.5;
#my $gP=0.9045;
#my $hP=0.809;

my $JP=1.;
#my $gP=.5; # .5; 1.; 1.5
my $gP=1.;
my $hP=0.;
my $vLR=2*$JP*$gP;
if($gP>1.){
    $vLR=2*$JP;
}

#my $dt=0.01;
my $dt=0.1; # As in Hastings' paper 2014
my $rate=10;

$basicdir="/ptmp/mpq/banulsm/TransverseLCLR/Ising";
$basictmpdir="/ptmp/mpq/banulsm/TransverseLCLR/IsingTmp";
$configfile="/u/banulsm/NewSVN/mpsdyn/src/config/transverseLC.conf";


my $Jlab=sprintf("%.2g",$JP); #int(-log($gamma)/log(10));
my $glab=sprintf("%.2g",$gP); #int(-log($gamma)/log(10));
my $hlab=sprintf("%.2g",$hP); #int(-log($gamma)/log(10));
my $dtlab=sprintf("%.2g",$dt); #int(-log($dt)/log(10));

my @initStNames=("X","Y","Z");

$initStName=$initStNames[$initSt-1];


$basicdir="${basicdir}/g${glab}/h${hlab}/init${initStName}";
$basictmpdir="${basictmpdir}/g${glab}/h${hlab}/init${initStName}";


$MSUB="msub_slurm";

my $scriptfile="lanzaTransLCLR.sh";
open MYFILE, '>', $scriptfile or die "Cannot open ${scriptfile} $!";
print MYFILE "#!/bin/bash\n\n";
#print MYFILE "mkdir -p ${outputdir} \n";

#my $Ls=2;

my @Ds=(20,40,60,80); #100,200,300,400,500,1000);
#my @Ds=(80,100,200);
#my @Ds=(200,500);
#my @Ds=(800);
#my @Ds=(1000);
my @Ds=(100,200);
my @buffersLR=(1,2,3,4);

my $maxT=20;
my $M=$maxT/$dt;

my $savingfreq=100;

my $append = 0;

$outputdir="${basicdir}";
$tmpdir="${basictmpdir}";
foreach $Dcut (@Ds){
    foreach $bufferLR (@buffersLR){
	if($Dcut>=100){ 
	    $parallel=40;
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
	
	$NAME="tlclr${jobStr}${initStName}_g${glab}h${hlab}_dt${dt}_D${Dcut}_${bufferLR}";
	my $filename="init${initSt}_L${L}J${Jlab}g${glab}h${hlab}dt${dtlab}_D${Dcut}_${bufferLR}";
	my $outfile="${outputdir}/${filename}";
	my $args=" ${configfile} -J=${JP} -g=${gP} -h=${hP} -M=${M} -outputfile=${outfile} -delta=${dt} -D=${Dcut} -savingfreq=${savingfreq} -initSt=${initSt} -rate=${rate} -tmpdir=${tmpdir} -append=${append} -bufferLR=${bufferLR} -vLR=${vLR} ";
	print MYFILE "$MSUB -N $NAME ${parallel} -- $EXEC ${args} \n";
	print MYFILE "sleep 1\n";
	#    print MYFILE "$EXEC ${args} \n";
    }
}
close(MYFILE);
system("chmod a+x ${scriptfile} ");


