#!/usr/bin/perl

# Script to convert the binary MPS files from the program to text

# Name of the script file to be created
my $scriptFile="lanzaSchwConv.sh";

my $EXEC="/u/banulsm/NewSVN/mpsdyn/src/bin/mps2txt";

my $basicDir="/ptmp/mpq/banulsm/schwingerFractal/checks/tol8";


my $MSUB="msub_slurm";


#masses="0 0.125 0.25 0.5" 
#my $x=1;
#my @masses=(0.05); #(0.0625, 1.); 
my $x=(1, 100,100);
my @masses=(0.1, 0., 0.125); # with each value of x

# dimensions for which I export
my $D=20; 

my $configfile="/u/banulsm/NewSVN/mpsdyn/src/config/phase1Prop.conf";

my @Ls=(306); #(12, 306);



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


foreach my $L (@Ls){

    foreach my $D (@Ds){

	foreach my $mg (@masses){

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

