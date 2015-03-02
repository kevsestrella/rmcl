#!/bin/bash

BMAPPING=BioGRID_AGaser_graph_mapping
DMAPPING=DIP_AGaser_graph_mapping
CLUSTER_MAPPING=clusters_with_proteins_list
BCLUSTER_MAPPING=banno_clusters_with_proteins_list
DCLUSTER_MAPPING=danno_clusters_with_proteins_list
FSCE=FScoring
FSCIN=FScoringInc

if [ "$#" -ne 2 ]; then
	echo "Usage: <b or d> <input file results to score>"

else
	dataset=$1
	name=$2

	echo "Calculating FScores for all four scenarios..."
	echo "Calculating FScores for Scenario A (Inc NAN, Elim Missing)..."


	if [ "$dataset" = "b" ]; then
		./$FSCIN "${name}" $BMAPPING $BCLUSTER_MAPPING "${name}_scored_scenarioA"
	
	else	
		./$FSCIN "${name}" $DMAPPING $DCLUSTER_MAPPING "${name}_scored_scenarioA"
	fi

	echo "Scenario A computed. Calculating FScores for Scenario B (Ex NAN, Elim Missing)..."

	if [ "$dataset" = "b" ]; then
		./$FSCE "${name}" $BMAPPING $BCLUSTER_MAPPING "${name}_scored_scenarioB"
	
	else	
		./$FSCE "${name}" $DMAPPING $DCLUSTER_MAPPING "${name}_scored_scenarioB"
	fi
	
	echo "Scenario B computed. Calculating FScores for Scenario C (Inc NAN, Use Complete CYC)..."

	if [ "$dataset" = "b" ]; then
		./$FSCIN "${name}" $BMAPPING $CLUSTER_MAPPING "${name}_scored_scenarioC"
	
	else	
		./$FSCIN "${name}" $DMAPPING $CLUSTER_MAPPING "${name}_scored_scenarioC"
	fi

	echo "Scenario C computed. Calculating FScores for Scenario D (Ex NAN, Use Complete CYC)..."

	if [ "$dataset" = "b" ]; then
		./$FSCE "${name}" $BMAPPING $CLUSTER_MAPPING "${name}_scored_scenarioD"
	
	else	
		./$FSCE "${name}" $DMAPPING $CLUSTER_MAPPING "${name}_scored_scenarioD"
	fi

	echo "All scenarios computed. Finished computing FScores."
fi
