#!/bin/bash
FSE=FSWeight
BIN=biogrid.in
DIN=dip.in
BMAPPING=BioGRID_AGaser_graph_mapping
DMAPPING=DIP_AGaser_graph_mapping
CAE=CAw
FSCE=FScoring


if [ "$#" -ne 4 ]; then
	echo "Usage: <b or d> <FS threshold> <CAw alpha> <CAw gamma>"

else
	name=$1
	threshold=$2
	alpha=$3
	gamma=$4

	echo "Begin clustering."
	echo "Computing FSWeights at $threshold threshold..."

	if [ "$name" = "b" ]; then
		./$FSE $BIN "${name}t${threshold}" "${name}mt${threshold}" ${threshold}

	else
		./$FSE $DIN "${name}t${threshold}" "${name}mt${threshold}" ${threshold}
	fi
		
	echo "Finished computing FSWeights. Now running MLRMCL..."
	mv "${name}mt${threshold}" mlr_and_mlroutputs
	cd mlr_and_mlroutputs
	./mlrmcl -o "${name}ct${threshold}" "${name}mt${threshold}"
		
	cd ../

	echo "Finished MLRMCL. Now running CAw..."
	./$CAE -a "${alpha}" -g "${gamma}" "${name}t${threshold}" "./mlr_and_mlroutputs/${name}ct${threshold}" "${name}pt${threshold}a${alpha}g${gamma}"
		
	mv "${name}pt${threshold}a${alpha}g${gamma}" ./predicted
	cd ../
	echo "Finished CAw. Clustering complete."

	#clean up!
	#rm "${name}t${threshold}"
	#rm "./mlr_and_mlroutputs/${name}mt${threshold}"
fi
