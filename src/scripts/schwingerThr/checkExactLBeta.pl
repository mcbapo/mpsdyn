#!/usr/bin/perl

my $dirData="/ptmp/mpq/banulsm/SchwingerThrm/exactL/";
my $x=$ARGV[0];


my $targetBeta=8;

# I will just check for desired stepBeta
my $stepBeta=0.1;
my $delta=$stepBeta*.25/sqrt($x); 
my @Ks=(1,2,4);
my @Ds=(60,80,100,120,140);
if($x>=256){
    @Ds=(80,100,120,140);
}
if($x>=121&&$x<1000){
    @Ds[scalar(@Ds)]=160;
}
my @Ls=(16,20,24);
my $Lmax=10;
my $tol=1E-8;

foreach $factor (@Ls){
    my $N=int($factor*sqrt($x));
    foreach $k (@Ks){
	my $deltak=$delta/$k;
	foreach $D (@Ds){
	    my $file=sprintf("x%d.N%d.dt%.4E.D%d.Lmax%d.p%d",$x,$N,$deltak,$D,$Lmax,int(-log($tol)/log(10)));
	    my $lastLine=`tail -n 1 ${dirData}/results/x$x/$file `;
	    my @entries=split(/\s+/,$lastLine);
#      print "\t which has entries @entries \n";
	    my $lastBeta=@entries[2];
	    print "x=$x, N=$N, dt=$deltak, D=$D, g*beta=$lastBeta (Lcut=$Lmax) \t";
	    if($lastBeta<$targetBeta){
		# Check if running
		my $NAME=sprintf("x%d.N%d.Lc%d.dt%.4E.D%d",$x,$N,$Lmax,$deltak,$D);
		my $isRunning=`qstat -j $NAME |grep job_number`;
		if($isRunning =~ /(\d+)/){
		    $jobNr=$1;
		    print " running \n";
		}
		else{
		    print " RELAUNCH!!!!!!!********\n ";
		}
	    }
	    else{
		print "Complete\n";
	    }
	}	
    }
}



# my $fileList=`ls ${dirData}/results/x$x/*Lmax${Lcut}*`;
# # Break in lines
# #  print "For x=$x, fileList=$fileList \n";
# my @files=split(/\n/,$fileList);
# # FORMAT: x256.N384.dt7.8125E-04.D80.Lmax15.p8
# foreach my $f (@files){
# #      print " Last line of file $f is: ";
#     my $lastLine=`tail -n 1 $f `;
# # Get the parameters for this file from the file name
#     $f =~ m/.*N(\d+)\.dt(.*)\.D(\d+)\.Lmax(\d+)\.p/; my $N=$1;my $D=$3;my $dt=$2;my $Lmax=$4;
# #      print "$lastLine \n";
#     my @entries=split(/\s+/,$lastLine);
# #      print "\t which has entries @entries \n";
#     my $lastBeta=@entries[2];
#     print "x=$x, N=$N, dt=$dt, D=$D, g*beta=$lastBeta (Lcut=$Lmax) \t";
#     if($lastBeta<$targetBeta){
# 	# Check if running
# 	my $NAME=sprintf("x%d.N%d.Lc%d.dt%.4E.D%d",$x,$N,$Lmax,$dt,$D);
# 	my $isRunning=`qstat -j $NAME |grep job_number`;
# 	if($isRunning =~ /(\d+)/){
# 	    $jobNr=$1;
# 	    print " running \n";
# 	}
# 	else{
# 	    print " RELAUNCH!!!!!!!********\n ";
# 	}
#     }
#     else{
# 	print "Complete\n";
#     }
# }

