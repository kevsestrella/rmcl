#!/bin/bash
FSE=FSWeight
BIN=biogrid.in
DIN=dip.in
BMAPPING=BioGRID_AGaser_graph_mapping
DMAPPING=DIP_AGaser_graph_mapping
CLUSTER_MAPPING=clusters_with_proteins_list
BCLUSTER_MAPPING=banno_clusters_with_proteins_list
DCLUSTER_MAPPING=danno_clusters_with_proteins_list
CAE=CAw
FSCE=FScoringParams
FSCEI=FScoringIncParams

array=( 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0 10.5 11.0 11.5 12.0 12.5 13.0 13.5 14.0 14.5 15.0 15.5 16.0 16.5 17.0 17.5 18.0 18.5 19.0 19.5 20.0 )
array2=( 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50 1.55 1.60 1.65 1.70 1.75 1.80 1.85 1.90 1.95 2.00 2.05 2.10 2.15 2.20 2.25 2.30 2.35 2.40 2.45 2.50 2.55 2.60 2.65 2.70 2.75 2.80 2.85 2.90 2.95 3.00 3.05 3.10 3.15 3.20 3.25 3.30 3.35 3.40 3.45 3.50 3.55 3.60 3.65 3.70 3.75 3.80 3.85 3.90 3.95 4.00 4.05 4.10 4.15 4.20 4.25 4.30 4.35 4.40 4.45 4.50 4.55 4.60 4.65 4.70 4.75 4.80 4.85 4.90 4.95 5.00 )
array3=( 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 )

if [ "$#" -ne 4 ]; then
	echo "Usage: <b or d> <a or g> <threshold value> <othervalue>"

else
	name=$1
	parameter=$2
	threshold=$3
	othervalue=$4

	echo "Begin clustering."
	echo "Computing FSWeights for threshold ${threshold}..."

	mkdir -p outputs
	mkdir -p rawoutputs

	if [ "$name" = "b" ]; then
		./$FSE $BIN "${name}t${threshold}" "${name}mt${threshold}" "${threshold}"

	else
		./$FSE $DIN "${name}t${threshold}" "${name}mt${threshold}" "${threshold}"
	fi
			
	echo "Finished computing FSWeights. Now running MLRMCL..."

	./mlrmcl -o "${name}ct${threshold}" "${name}mt${threshold}"


	echo "Finished MLRMCL. Now running CAw..."


	if [ "$parameter" = "a" ]; then
		mkdir "./rawoutputs/${name}varyingalphat${threshold}g${othervalue}"
		for j in "${array[@]}"
			do 
				echo "Running CAw with alpha at ${j}, gamma at ${othervalue}"
				./$CAE -a "$j" -g "${othervalue}" "${name}t${threshold}" "${name}ct${threshold}" "./rawoutputs/${name}varyingalphat${threshold}g${othervalue}/${name}pt${threshold}a${j}g${othervalue}"
			done	

	else
		mkdir "./rawoutputs/${name}varyinggammat${threshold}a${othervalue}"
		for j in "${array2[@]}"
			do
				echo "Running CAw with alpha at ${othervalue}, gamma at ${j}"
				./$CAE -a "${othervalue}" -g "${j}" "${name}t${threshold}" "${name}ct${threshold}" "./rawoutputs/${name}varyinggammat${threshold}a${othervalue}/${name}pt${threshold}a${othervalue}g${j}"
			done
	fi
	

	echo "Finished CAw. Now calculating FScores..."

	if [ "$parameter" = "a" ]; then
		mkdir "./outputs/${name}varyingalphat${threshold}g${othervalue}"
		for k in "${array[@]}"
			do 
				if [ "$name" = "b" ]; then
					./$FSCE "./rawoutputs/${name}varyingalphat${threshold}g${othervalue}/${name}pt${threshold}a${k}g${othervalue}" $BMAPPING $CLUSTER_MAPPING "./outputs/${name}varyingalphat${threshold}g${othervalue}/scenarioD.csv" "${threshold}" "${k}" "${othervalue}"
					./$FSCE "./rawoutputs/${name}varyingalphat${threshold}g${othervalue}/${name}pt${threshold}a${k}g${othervalue}" $BMAPPING $BCLUSTER_MAPPING "./outputs/${name}varyingalphat${threshold}g${othervalue}/scenarioB.csv" "${threshold}" "${k}" "${othervalue}"
					./$FSCEI "./rawoutputs/${name}varyingalphat${threshold}g${othervalue}/${name}pt${threshold}a${k}g${othervalue}" $BMAPPING $CLUSTER_MAPPING "./outputs/${name}varyingalphat${threshold}g${othervalue}/scenarioC.csv" "${threshold}" "${k}" "${othervalue}"
					./$FSCEI "./rawoutputs/${name}varyingalphat${threshold}g${othervalue}/${name}pt${threshold}a${k}g${othervalue}" $BMAPPING $BCLUSTER_MAPPING "./outputs/${name}varyingalphat${threshold}g${othervalue}/scenarioA.csv" "${threshold}" "${k}" "${othervalue}"
					
				else
					./$FSCE "./rawoutputs/${name}varyingalphat${threshold}g${othervalue}/${name}pt${threshold}a${k}g${othervalue}" $DMAPPING $CLUSTER_MAPPING "./outputs/${name}varyingalphat${threshold}g${othervalue}/scenarioD.csv" "${threshold}"  "${k}" "${othervalue}"
					./$FSCE "./rawoutputs/${name}varyingalphat${threshold}g${othervalue}/${name}pt${threshold}a${k}g${othervalue}" $DMAPPING $DCLUSTER_MAPPING "./outputs/${name}varyingalphat${threshold}g${othervalue}/scenarioB.csv" "${threshold}"  "${k}" "${othervalue}"
					./$FSCEI "./rawoutputs/${name}varyingalphat${threshold}g${othervalue}/${name}pt${threshold}a${k}g${othervalue}" $DMAPPING $CLUSTER_MAPPING "./outputs/${name}varyingalphat${threshold}g${othervalue}/scenarioC.csv"  "${threshold}"  "${k}" "${othervalue}"
					./$FSCEI "./rawoutputs/${name}varyingalphat${threshold}g${othervalue}/${name}pt${threshold}a${k}g${othervalue}" $DMAPPING $DCLUSTER_MAPPING "./outputs/${name}varyingalphat${threshold}g${othervalue}/scenarioA.csv" "${threshold}"  "${k}" "${othervalue}"
					
				fi
			done
	else #BRUTE FORCE! HAHA. let's optimize this laterz
		mkdir "./outputs/${name}varyinggammat${threshold}a${othervalue}"
		for k in "${array2[@]}"
			do 
				if [ "$name" = "b" ]; then
					./$FSCE "./rawoutputs/${name}varyinggammat${threshold}a${othervalue}/${name}pt${threshold}a${othervalue}g${k}" $BMAPPING $CLUSTER_MAPPING "./outputs/${name}varyinggammat${threshold}a${othervalue}/scenarioD.csv" "${threshold}" "${othervalue}" "${k}"
					./$FSCE "./rawoutputs/${name}varyinggammat${threshold}a${othervalue}/${name}pt${threshold}a${othervalue}g${k}" $BMAPPING $BCLUSTER_MAPPING "./outputs/${name}varyinggammat${threshold}a${othervalue}/scenarioB.csv" "${threshold}" "${othervalue}" "${k}"
					./$FSCEI "./rawoutputs/${name}varyinggammat${threshold}a${othervalue}/${name}pt${threshold}a${othervalue}g${k}" $BMAPPING $CLUSTER_MAPPING "./outputs/${name}varyinggammat${threshold}a${othervalue}/scenarioC.csv" "${threshold}" "${othervalue}" "${k}"
					./$FSCEI "./rawoutputs/${name}varyinggammat${threshold}a${othervalue}/${name}pt${threshold}a${othervalue}g${k}" $BMAPPING $BCLUSTER_MAPPING "./outputs/${name}varyinggammat${threshold}a${othervalue}/scenarioA.csv" "${threshold}" "${othervalue}" "${k}"
				else
					./$FSCE "./rawoutputs/${name}varyinggammat${threshold}a${othervalue}/${name}pt${threshold}a${othervalue}g${k}" $DMAPPING $CLUSTER_MAPPING "./outputs/${name}varyinggammat${threshold}a${othervalue}/scenarioD.csv" "${threshold}" "${othervalue}" "${k}"
					./$FSCE "./rawoutputs/${name}varyinggammat${threshold}a${othervalue}/${name}pt${threshold}a${othervalue}g${k}" $DMAPPING $DCLUSTER_MAPPING "./outputs/${name}varyinggammat${threshold}a${othervalue}/scenarioB.csv" "${threshold}" "${othervalue}" "${k}"
					./$FSCEI "./rawoutputs/${name}varyinggammat${threshold}a${othervalue}/${name}pt${threshold}a${othervalue}g${k}" $DMAPPING $CLUSTER_MAPPING "./outputs/${name}varyinggammat${threshold}a${othervalue}/scenarioC.csv" "${threshold}" "${othervalue}" "${k}"
					./$FSCEI "./rawoutputs/${name}varyinggammat${threshold}a${othervalue}/${name}pt${threshold}a${othervalue}g${k}" $DMAPPING $DCLUSTER_MAPPING "./outputs/${name}varyinggammat${threshold}a${othervalue}/scenarioA.csv" "${threshold}" "${othervalue}" "${k}"
				
				fi
			done
	fi	
		
	echo "FScores computed. Clustering complete."
fi
