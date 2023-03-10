#!/usr/bin/perl

# Script to prepare a series of jobs using the searchO program, but now for the sizes that are relevant for the recursion

# Name of the script file to be created
my $scriptFile="lanzaSchw.sh";

my $EXEC="/u/banulsm/NewSVN/mpsdyn/src/bin/searchO";

my $basicDir="/ptmp/mpq/banulsm/schwingerFractal/checks/tol8";
my $jobsdir="/u/banulsm/jobsToRun";


my $MSUB="msub_slurm";


#masses="0 0.125 0.25 0.5" 
my $x=1;
my @masses=(0.,0.05,1.25); #(0.0625, 1.); 
#my $x=10;
#my @masses=(0.,0.05/sqrt(10),1.25/sqrt(10)); #(0.0625, 1.); 
#my $x=100;
#my @masses=(0.,0.005,0.125); # with x=100

my $x=10;
my @masses=(-0.5,-0.7,-1.); #(0.0625, 1.); 
my @masses=(-0.7,-1.); 

# dimensions for which I run the searchO program
my @Ds=(20); #, 60, 80, 100, 120, 140, 160);
my $INITD=20;

# min number of states
#MINNUMBER="25"
my $MINNUMBER=2;
#Ds="40 60 80 100 120"
my $ONLYVEC=" -onlyvec=1 ";

my $alpha=0;

my $app=1;

my $configfile="/u/banulsm/NewSVN/mpsdyn/src/config/phase1Prop.conf";

my @Ls=(12,14,18,20,30,40,50,80,100,150,200,250,306);
#my @Ls=(12,14,18,20,30,40,50); my @masses=(1.25);

my $PENALTY=2*sqrt($x); if($PENALTY<100){$PENALTY=100;}




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


foreach my $D (@Ds){

    foreach my $mg (@masses){

	foreach my $L_ (@Ls){

	    #my @Ls_=($L_-2,$L_-3);
	    my @Ls_=($L_,$L_-2,$L_-4,$L_-6,$L_-8);
	    my @Ls_=($L_);

	    foreach my $L (@Ls_){ 
		# Require extra memory (more than 3 GB) if 32*L*D*D (storage in Bytes of an MPS)
		# reaches 100 MB.
		# If it goes to 200, move to disk storage, instead, and require long time
		
		# my $sizeMPS=16*$L*$D*$D*2;
		# my $sizeTmp=16*$N*2*$D*$D*(2*$Lmax+1);
		# my $sizeMPOsL=16*(16*(2*$Lmax+1)*(2*$Lmax+1)*$N);
		# my $sizeMPOsx=16*$N*16*4*4*3;
		# my $sizeAuxMPOs=16*4*(2*$Lmax+1)*(2*$Lmax+1)*$N;
		# my $SIZEMEM=int((3*$sizeMPS+$sizeTmp+2*$sizeMPOsL+3*$sizeMPOsx+$sizeAuxMPOs)*1.2*1E-6)+200; # I will add a 20% AND 200M
		
		my $SIZEMEM=int(16*$L*$D*$D*(10+4*$MINNUMBER+2)*1.1*1E-6)+200;
	
		my $LONGMEM="-mem ${SIZEMEM}M";
		if ($SIZEMEM<3000) {
		    $LONGMEM="";
		}
		
		my $par="";
		if ($D>80){
		    $par="-P 4";
		}
		
		my $outdir="${basicDir}/M${mg}/x${x}/L${L}/results";
		
		print MYFILE "mkdir -p $outdir \n";
		
		my $outfile="L${L}x${x}D${D}";
		
		my $dirMPSs="${basicDir}/M${mg}/x${x}/L${L}/mps";
		print MYFILE "mkdir -p $dirMPSs \n";
		
		my $JOBNAME="D${D}.M${mg}.${x}.${L}";
		
		my $args="${configfile} -outfile=${outfile} -L=${L} -mg=${mg} -x=${x} -D=${D} -outputdir=${outdir} -penalty=${PENALTY} -mpsdir=${dirMPSs} -initmpsdir=${dirMPSs} -initD=${INITD} -tol=1E-8 -minnumber=${MINNUMBER} -jobsdir=${jobsdir} ${ONLYVEC} ";

		print MYFILE "${MSUB} -N ${JOBNAME} ${par} ${LONGMEM} -- ${EXEC} ${args}\n"; 
		print MYFILE "sleep 1\n";
	    }
	}
    }
}

