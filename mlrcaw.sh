#!/bin/bash
FSE=FSWeight
CAE=CAw



if [ "$#" -ne 4 ]; then
	echo "Usage: <input file> <FS threshold> <CAw alpha> <CAw gamma>"

else
	name=$1
	threshold=$2
	alpha=$3
	gamma=$4

	echo "Begin clustering."
	echo "Computing FSWeights at $threshold threshold..."

	mkdir -p clustering
	mkdir -p thresholdfiles

	./$FSE "${name}" "${name}t${threshold}" "./clustering/${name}mt${threshold}" ${threshold}
		
	echo "Finished computing FSWeights. Now running MLRMCL..."

		#mv "${name}mt${threshold}" bin
	cd bin
	#./srmcl -o "../clustering/${name}ct${threshold}" "../clustering/${name}mt${threshold}"
	./srmcl -o "../clustering/${name}ct${threshold}" "../clustering/${name}mt${threshold}"
	cd ../

	echo "Finished MLRMCL. Now running CAw..."
	./$CAE -a "${alpha}" -g "${gamma}" "${name}t${threshold}" "./clustering/${name}ct${threshold}" "./scoring/${name}pt${threshold}a${alpha}g${gamma}"

	echo "Finished CAw. Clustering complete."

	#clean up
	
	mv "${name}t${threshold}" "./thresholdfiles/${name}t${threshold}"

fi
