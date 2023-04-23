#!/bin/bash

# Set the environment variables required (in the cluster)
# to link to the right libraries

if [[ $HOSTNAME == *tqo* ]]
then
# only in the cluster
    echo "In the Cluster"

# in case it was not done yet
module load mkl
module load intel/19.1.3

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MKL_HOME}/lib/intel64

else
. /opt/intel/bin/compilervars.sh intel64
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/intel/composerxe-2011.2.137/mkl/lib/intel64
fi
