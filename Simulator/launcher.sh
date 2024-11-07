#!/bin/bash

printUsage() {
    echo "---- USAGE OF BASH SCRIPT launcher.sh ----"
    echo "This script must receive the following parameters:"
    echo "1) Configuration file"
    echo "2) Number of replications"
    echo "You can optionally specify the starting seed (default is 1)"
    echo "An example of command is:"
    echo "./launcher.sh cforaging.ini 10"
}

if [ $# -lt 2 ]
  then
    echo "ERROR! Invalid number of parameters!!! Check usage and retry!!!"
    printUsage
    exit -1
fi

CONFIGFILE=$1
NUM_REPLICATIONS=$2
STARTING_SEED=1
if [ $# -ge 3 ]
  then
    STARTING_SEED=$3
fi
GLOBDIR=$(pwd)
if [ $# -ge 4 ]
  then
    GLOBDIR=$4
fi
FILEDIR=$GLOBDIR
if [ $# -eq 5 ]
  then
    FILEDIR=$5
fi

FOLDER=$(pwd)
EXEDIR=$FOLDER
EXEDIR="$(echo $EXEDIR | sed -r 's/cforaging/build/')"
CMD_EXE=$EXEDIR/cforaging

# Run N replications
CREP=1
CURRSEED=$STARTING_SEED
while [ $CREP -le $NUM_REPLICATIONS ]
do
    # Arguments
    ARGS=" $CONFIGFILE $CURRSEED $GLOBDIR $FILEDIR"
    # Run current replication
    $CMD_EXE $ARGS
    # Update current seed
    CURRSEED=$((CURRSEED+1))
    # Update current replication
    CREP=$((CREP+1))
done
