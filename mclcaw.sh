#!/bin/bash
if [ "$#" -ne 5 ]; then
	echo "Usage: <input file> <FS threshold> <CAw alpha> <CAw gamma> <mode:mlr or sr>"

else
    #identify OS
    unameOut="$(uname -s)"
    case "${unameOut}" in
        Linux*)     machine=Linux;;
        Darwin*)    machine=OSX;;
        CYGWIN*)    machine=Cygwin;;
        MINGW*)     machine=MinGw;;
        *)          machine="UNKNOWN:${unameOut}"
    esac

    #assign arguments
    name=$1
    filename=${name##*/}
	threshold=$2
	alpha=$3
	gamma=$4
    mode=$5

    #exe
    WORKDIR=$PWD
    BIN=$WORKDIR/bin/$machine
    FSE=$BIN/FSWeight
    if [ "$mode" == "sr" ]; then MCL=$BIN/srmcl; else MCL=$BIN/mlrmcl; fi
    CAE=$BIN/CAw
    CAEPY=$BIN/CAw.py

    #files
    NETWORKS=$WORKDIR/networks
    NETWORKSFILE="$NETWORKS/${filename}-t-${threshold}"
    NETWORKSFILE4MCL="$NETWORKS/${filename}-mt-${threshold}"
    CLUSTERS=$WORKDIR/clustering
    CLUSTERFILE="$CLUSTERS/${filename}-ct-${mode}-${threshold}"
    COREATTACHMENTS=$WORKDIR/coreattachments
    COREATTCHFILE=$COREATTACHMENTS/${filename}-pt-${mode}-${threshold}-a${alpha}-g${gamma}

    echo "Begin clustering."
    echo "Computing FSWeights at $threshold threshold..."
   
    $FSE "${name}" $NETWORKSFILE $NETWORKSFILE4MCL ${threshold}
        
    echo "Finished computing FSWeights. Now running ${mode}MCL..."

    $MCL $NETWORKSFILE4MCL -o $CLUSTERFILE

    echo "Finished ${mode}MCL. Now running CAw..."

    if [ "$mode" == "sr" ]; then
         python $CAEPY $NETWORKSFILE $CLUSTERFILE $COREATTCHFILE
    else 
        $CAE -a "${alpha}" -g "${gamma}" $NETWORKSFILE $CLUSTERFILE $COREATTCHFILE
    fi

fi
