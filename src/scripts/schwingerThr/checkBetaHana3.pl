#!/usr/bin/perl

my $dirData="/ptmp/mpq/banulsm/SchwingerThrm/Hana3/";
my @Xs=(55,60,65);


foreach my $x (@Xs){
    my $fileList=`ls ${dirData}/x$x/output*`;
# Break in lines
#  print "For x=$x, fileList=$fileList \n";
  my @files=split(/\n/,$fileList);
  foreach my $f (@files){
#      print " Last line of file $f is: ";
      my $lastLine=`tail -n 1 $f `;
# Get the parameters for this file from the file name
      $f =~ m/.*N(\d+)D(\d+)dlt(.*)_Thrm/; my $N=$1;my $D=$2;my $dt=$3;
#      print "$lastLine \n";
      my @entries=split(/\s+/,$lastLine);
#      print "\t which has entries @entries \n";
      my $lastTime=@entries[0]*2*sqrt($x);
      print "x=$x, N=$N, dt=$dt, D=$D, g*beta=$lastTime \n";
#      print "Last g beta ". 2*sqrt($x)*$lastTime . " for x=$x, N=$N, D=$D, dt=$dt ($f) \n";
  }
}
