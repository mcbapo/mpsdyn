#!/usr/bin/perl

# Launch the calculation of energies and fractal weights from results of searchO

# Name of the script file to be created
my $scriptFile="lanzaSchwFrW.sh";

my $EXEC="/u/banulsm/NewSVN/mpsdyn/src/bin/schFrW";

my $basicDir="/ptmp/mpq/banulsm/schwingerFractal/checks/tol8/";

my $MSUB="msub_slurm";


#masses="0 0.125 0.25 0.5" 
my $x=1;
my @masses=(0,0.05,1.25); #(0.0625, 1.); 
my $x=10;
my @masses=(0,0.05/sqrt(10),1.25/sqrt(10));
my $x=100;
my @masses=(0.,0.005,0.125); # with x=100

# dimensions for which I run the searchO program
#my @Ds=(20, 40, 60, 80, 100, 120, 140, 160);
my @Ds=(20); # Now it continues on its own

my $alpha=0;

my $app=0; # Only for D==20

my $configfile="/u/banulsm/NewSVN/mpsdyn/src/config/phase1Prop.conf";

my @Ls=(12,14,18,20,30,40,50,80,100,150,200,250,306);

my $lev=0; # only GS now

#my $PENALTY=2*sqrt($x); if($PENALTY<100){$PENALTY=100;}


# Check for existing jobs                                                                                                                                                                    
sub isRunning{
    my ($jobname)=@_;
    #print "Checking whether job $jobname exists\n";                                                                                                                                         
    my $result=`squeue --noheader -n $jobname -o "%o" `;
    if (length $result){
        print "Job $jobname exists: $result \n";
	return 1;
    }
    return 0;
};


open MYFILE, '>', $scriptFile or die "Cannot open file $scriptFile $!";

print MYFILE "#!/bin/bash\n\n";


my $app=0; # Only for the first D
foreach my $D (@Ds){
    foreach my $L0 (@Ls){
	#my @Ls_=($L0,$L0-2,$L0-3,$L0-4,$L0-5);
	my @Ls_=($L0,$L0-2,$L0-4); #odd need different mpsdir
	foreach my $L (@Ls_){
	    foreach my $mg (@masses){

		# Require extra memory (more than 3 GB) if 32*L*D*D (storage in Bytes of an MPS)
		# reaches 100 MB.
		# If it goes to 200, move to disk storage, instead, and require long time
		
		my $SIZEMEM=int(16*$L*$D*$D*(10+4*$MINNUMBER+2)*1.1*1E-6)+200;
		
		my $LONGMEM="-mem ${SIZEMEM}M";
		if ($SIZEMEM<3000) {
		    $LONGMEM="";
		}
		
		my $par="";
		if ($D>80){
		    $par="-P 4";
		}
		
#		my $outdir="${basicDir}/M${mg}/x${x}/L${L}/resultsRev"; # projecting beginning of string(mismatch to their conv)
#		my $outdir="${basicDir}/M${mg}/x${x}/L${L}/resultsRev"; # projecting the end of the chain (agrees with conv), but original blocks
#		my $outdir="${basicDir}/M${mg}/x${x}/L${L}/resultsRevNew"; # projecting the end, according to four new blocks 21.10.2021
#		my $outdir="${basicDir}/M${mg}/x${x}/L${L}/resultsNew"; # projecting the beginning, according to four new blocks 24.10.2021
		my $outdir="${basicDir}/M${mg}/x${x}/L${L}/resultsNewDeco"; # projecting the beginning, according to four new blocks 2.11.2021
		
		print MYFILE "mkdir -p $outdir \n";
		
		my $outfile="weightsL${L}x${x}";
		my $outfileH="HeffL${L}x${x}";
		
		my $dirMPSs="${basicDir}/M${mg}/x${x}/L${L}/mps";
		#print MYFILE "mkdir -p $dirMPSs \n";
		
		my $JOBNAME="schW.D${D}.M${mg}.${x}.${L}";
		
		my $args="${configfile} -outfile=${outfile} -outfileHeff=${outfileH} -L=${L} -mg=${mg} -alpha=${alpha} -x=${x} -D=${D} -outputdir=${outdir} -mpsdir=${dirMPSs}  -tol=1E-8 -nLev=${lev} -append=${app} ";

		print MYFILE "${MSUB} -N ${JOBNAME} ${par} ${LONGMEM} -- ${EXEC} ${args}\n"; 
		print MYFILE "sleep 1\n";
	    }
	}
    }
#    $app=1; 
#    print MYFILE "sleep 300\n";
}# end of D loop

