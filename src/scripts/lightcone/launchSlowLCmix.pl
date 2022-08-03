#!/usr/bin/perl


$EXEC="./slowLCmix";

#my $gP=0.9045;
#my $hP=0.809;

my $initSt=1;

# The ones in the thermalization paper
#my $JP=1.;
#my $gP=-1.05;
#my $hP=0.5;

my $JP=1.;
#my $gP=1.; 1.5
my $gP=0.5;
my $hP=0.;


my $dt=0.01;
my $rate=10;

$basicdir="/ptmp/mpq/banulsm/SlowConeMix/Ising/";
$basictmpdir="/ptmp/mpq/banulsm/SlowConeMix/IsingTmp/";
$configfile="/afs/ipp/home/b/banulsm/SVN/mpsdyn/src/config/lightCone.conf";


my $Jlab=sprintf("%.2g",$JP); #int(-log($gamma)/log(10));
my $glab=sprintf("%.2g",$gP); #int(-log($gamma)/log(10));
my $hlab=sprintf("%.2g",$hP); #int(-log($gamma)/log(10));
my $dtlab=sprintf("%.2g",$dt); #int(-log($dt)/log(10));

my @initStNames=("X","Y","Z");

$initStName=$initStNames[$initSt-1];


$basicdir="${basicdir}/g${glab}/h${hlab}/init${initStName}/${subDir}";
$basictmpdir="${basictmpdir}/g${glab}/h${hlab}/init${initStName}/${subDir}";


$MSUB="msub_slurm";

my $scriptfile="ScriptLanzadorDeLC.sh";
open MYFILE, '>', $scriptfile or die "Cannot open ${scriptfile} $!";
print MYFILE "#!/bin/bash\n\n";
#print MYFILE "mkdir -p ${outputdir} \n";

my $Ls=2;

my @Ds=(100,200);
my @Nxs=(1,2);

my $maxT=20;
my $M=$maxT/$dt;

my $savingfreq=100;


foreach $Nx (@Nxs){
    $outputdir="${basicdir}/Nx${Nx}";
    $tmpdir="${basictmpdir}/Nx${Nx}";
    foreach $Dcut (@Ds){
	if($Dcut>400){ 
	    $parallel=20;
#	$savingfreq=50;
	}
	else{    
	    if($Dcut>=300){ 
		$parallel=12;
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
	
	$NAME="slcM${jobStr}${initStName}_g${glab}h${hlab}_Nx${Nx}D${Dcut}";
	my $filename="init${initSt}_L${L}J${Jlab}g${glab}h${hlab}dt${dtlab}_Ls${Ls}_D${Dcut}";
	my $outfile="${outputdir}/${filename}";
	my $args=" ${configfile} -J=${JP} -g=${gP} -h=${hP} -M=${M} -outputfile=${outfile} -delta=${dt} -Ls=${Ls} -D=${Dcut} -savingfreq=${savingfreq} -initSt=${initSt} -rate=${rate} -Nx=${Nx} -tmpdir=${tmpdir} ";
	print MYFILE "$MSUB -N $NAME ${parallel} -- $EXEC ${args} \n";
	print MYFILE "sleep 1\n";
#    print MYFILE "$EXEC ${args} \n";
    }
}
close(MYFILE);
system("chmod a+x ${scriptfile} ");


