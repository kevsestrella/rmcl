<h1 align="center">Hybrid Markov Clustering Algorithms </h1>
<p>
  <img alt="Version" src="https://img.shields.io/badge/version-2.0-blue.svg?cacheSeconds=2592000" />
</p>

> This project consist of 2 Hybrid Clustering Method namely MLR-CAw(1.0) and SR-CAw(2.0). 2.0 was mainly focused on SR-CAw, using SR-MCL as the clustering algorithm to the hybrid FSWeight - SR-MCL - Core Attachment combination, differing from V1 yung MLR-MCL instead.  SR-CAw is not an improvement of MLR-CAw but rather a variant with a different approach. See the accompanying paper for the discussion of the method and the algorithms.

### 🏠 [related page](https://sites.google.com/site/stochasticflowclustering/)

## Install

for python packages

```sh
pip install -r requirement.txt
```

## Usage

basic usage from FSWeighting to Core Attachment

```sh
make 
./mclcaw.sh <unweighted_graph_file> <threshold> <alpha> <gamma> <mode> <redundancy> <beta> <quality>
```
## A typical top-level directory layout

understand you file structure and contents(Team 2.0 learnings #6450)

     .
     └──     3.package/                    #Contains important files such as existing relevant libraries from metis, datasets
     │  └────     Lib/                         #Metis Libraries
     │  └────     gold_standard/               #CYC2008 Dataset
     │  └────     graphs/                      #BioGRID and DIP datasets
     │  └────     realGraphNodeName/           #Protein Mappings
     └──     bin/                          #Contains excecutables
     │  └────     Linux/                       #Executables for Linux
     │  └────     OSX/                         #Executables for macOS
     └──     clustering/                   #output directory of mcl variant
     └──     coreattachments/              #output directory of Core-Attachment
     └──     networks/                     #output directory of FSWeight when mclcaw.sh is used
     └──     scoring/                      #directory of F-Scoring, a separate step, see README 
     └──     src/                          #contains source file except CAw.py(located at bin)
     └──     varyingparameters/            #contains varying parameters files used by V 1.0

## Author

👤 **(1.0) Beltran JC, Montes C, Villar JJ, Valdez AR, (2.0) Cortez JJM, Estrella JKF, Octaviano M, Fabilloren SR, Villar JJ**

## Show your support

Give a ⭐️ if this project helped you!

***
_This README was generated with ❤️ by [readme-md-generator](https://github.com/kefranabg/readme-md-generator)_
