#!/bin/bash
FSE=FSWeight
BIN=biogrid.in
DIN=dip.in
BMAPPING=BioGRID_AGaser_graph_mapping
DMAPPING=DIP_AGaser_graph_mapping
CLUSTER_MAPPING=clusters_with_proteins_list
CAE=CAw
FSCE=FScoring

#echo $1 $2 $3 

if [ "$#" -ne 5 ]; then
	echo "Usage: <b or d> <a or g> <starting value> <end value> <step size>"

else
	name=$1
	parameter=$2
	start=$3
	end=$4
	step=$5	

	mkdir -p threshold_files
	#mkdir -p CAw_output_files
	
	
	echo "Begin clustering."
	echo "Computing FSWeights for default threshold..."

	if [ "$name" = "b" ]; then
		./$FSE $BIN "${name}t0.2" "${name}mt0.2" 0.2

	else
		./$FSE $DIN "${name}t0.2" "${name}mt0.2" 0.2
	fi

	mv "${name}t0.2" threshold_files
		
	echo "Finished computing FSWeights. Now running MLRMCL..."
	mv "${name}mt0.2" mlr_and_mlroutputs
	cd mlr_and_mlroutputs
	./mlrmcl -o "${name}ct0.2" "${name}mt0.2"
	# mv "${name}ct0.2" ../
	cd ../

	echo "Finished MLRMCL. Now running CAw..."

	nums=$(awk -v st="$start" -v en="$end" -v ste="$step" 'BEGIN{for(i=st; i<en; i+=ste)print i}')
	for j in $nums
		do
			if [ "$parameter" = "a" ]; then 
				echo "Running CAw with alpha at ${j}, gamma at 0.75"
				./$CAE -a "$j" "./threshold_files/${name}t0.2" "./mlr_and_mlroutputs/${name}ct0.2" "${name}pt0.2a${j}g0.75"
				mv "${name}pt0.2a${j}g0.75" ./predicted

			else
				echo "Running CAw with alpha at 1.0, gamma at ${j}"
				./$CAE -g "${j}" "./threshold_files/${name}t0.2" "./mlr_and_mlroutputs/${name}ct0.2" "${name}pt0.2a1.0g${j}"
				mv "${name}pt0.2a1.0g${j}" ./predicted
			fi
		done

	echo "Finished CAw. Now calculating FScores..."
	cd ./predicted

	for j in $nums
		do 
		 	if [ "$name" = "b" ]; then

		 		if [ "$parameter" = "a" ]; then
					./$FSCE "${name}pt0.2a${j}g0.75" $BMAPPING $CLUSTER_MAPPING "${name}sct0.2a${j}g0.75"
				else
					./$FSCE "${name}pt0.2a1.0g${j}" $BMAPPING $CLUSTER_MAPPING "${name}sct0.2a1.0g${j}"
				fi
		 	else
		 		if [ "$parameter" = "a" ]; then
					./$FSCE "${name}pt0.2a${j}g0.75" $DMAPPING $CLUSTER_MAPPING "${name}sct0.2a${j}g0.75"
				else
					./$FSCE "${name}pt0.2a1.0g${j}" $DMAPPING $CLUSTER_MAPPING "${name}sct0.2a1.0g${j}"
		 		fi

		 	fi
		done
		cd ../
		echo "FScores computed. Clustering complete."
fi
