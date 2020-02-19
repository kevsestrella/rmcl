#!/bin/bash

#identify OS
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=OSX;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="UNKNOWN:${unameOut}"
esac

#exe
WORKDIR=$PWD
BIN=$WORKDIR/bin/$machine
FSE=$BIN/FSWeight
if [ "$mode" == "sr" ]; then MCL=$BIN/srmcl; else MCL=$BIN/mlrmcl; fi
CAE=$BIN/CAw
CAEPY=$BIN/CAw.py

if [ "$#" -eq 0 ]; then
    BUILDDIR='3.package'
    make
    cd $BUILDDIR
    make realclean
    make
    cp srmcl mlrmcl $BIN
else
    if [ "$#" -ne 8 ]; then
        echo "Usage: <input file> <FS threshold> <CAw alpha> <CAw gamma> <mode:mlr or sr> <SRMCL p|MLR c> <SRMCL q|MLR b><SRMCL w|MLR i>"

    else

    #assign arguments
    name=$1
    filename=${name##*/}
    threshold=$2
    alpha=$3
    gamma=$4
    mode=$5
    p=$6
    c=$6
    q=$7
    b=$7
    w=$8
    i=$8

    #files
    NETWORKS=$WORKDIR/networks
    NETWORKSFILE="$NETWORKS/${filename}-t-${threshold}"
    NETWORKSFILE4MCL="$NETWORKS/${filename}-mt-${threshold}"
    CLUSTERS=$WORKDIR/clustering
    CLUSTERFILE="$CLUSTERS/${filename}-ct-${mode}-p$p-q$q-w$w-${threshold}"
    COREATTACHMENTS=$WORKDIR/coreattachments
    COREATTCHFILE=$COREATTACHMENTS/${filename}-pt-${mode}-p$p-q$q-w$w-${threshold}-a${alpha}-g${gamma}

    echo "Begin clustering."
    echo "Computing FSWeights at $threshold threshold..."

    $FSE "${name}" $NETWORKSFILE $NETWORKSFILE4MCL ${threshold}

    echo "Finished computing FSWeights. Now running ${mode}MCL..."

    if [ "$mode" == "sr" ]; then
        $MCL $NETWORKSFILE4MCL -t 30 -p $p -q $q -w $w -o $CLUSTERFILE
    else
        $MCL $NETWORKSFILE4MCL $o -o $CLUSTERFILE
    fi

    echo "Finished ${mode}MCL. Now running CAw..."

    if [ "$mode" == "sr" ]; then
        python $CAEPY -n $NETWORKSFILE -c $CLUSTERFILE -o $COREATTCHFILE -a $alpha -g $gamma
    else 
        $CAE -a "${alpha}" -g "${gamma}" $NETWORKSFILE $CLUSTERFILE $COREATTCHFILE
    fi

    fi
fi
