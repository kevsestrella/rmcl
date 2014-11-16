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

if [ "$#" -ne 4 ]; then
	echo "Usage: <b or d> <starting threshold> <end threshold> <step size>"

else
	name=$1
	start=$2
	end=$3
	step=$4


	nums=$(awk -v st="$start" -v en="$end" -v ste="$step" 'BEGIN{for(i=st; i<en; i+=ste)print i}')

	mkdir -p threshold_files
	mkdir -p CAw_output_files
	
	for j in $nums
	do
		echo "Begin clustering."
		echo "Computing FSWeights at $j threshold..."

		if [ "$name" = "b" ]; then
			./$FSE $BIN "${name}t${j}" "${name}mt${j}" ${j}

		else
			./$FSE $DIN "${name}t${j}" "${name}mt${j}" ${j}
		fi

		mv "${name}t${j}" threshold_files
		
		echo "Finished computing FSWeights. Now running MLRMCL..."
		mv "${name}mt${j}" mlr_and_mlroutputs
		cd mlr_and_mlroutputs
		./mlrmcl -o "${name}ct${j}" "${name}mt${j}"
		# mv "${name}ct${j}" ../
		cd ../

		echo "Finished MLRMCL. Now running CAw..."
		./$CAE "./threshold_files/${name}t${j}" "./mlr_and_mlroutputs/${name}ct${j}" "${name}pt${j}a1.0g0.75"
		# mv "${name}ct${j}" ./CAw_output_files
		mv "${name}pt${j}a1.0g0.75" ./predicted

		 echo "Finished CAw. Now calculating FScores..."
		 cd ./predicted

		 if [ "$name" = "b" ]; then
			./$FSCE "${name}pt${j}a1.0g0.75" $BMAPPING $CLUSTER_MAPPING "${name}sct${j}a1.0g0.75"
		 else
			./$FSCE "${name}pt${j}a1.0g0.75" $DMAPPING $CLUSTER_MAPPING "${name}sct${j}a1.0g0.75"
		 fi

		cd ../

		echo "FScores computed. Clustering complete."

	done
fi
