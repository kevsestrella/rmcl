WORKDIR=$PWD
MCLCAW=$WORKDIR/mclcaw.sh
GRAPHS=$WORKDIR/3.package/graphs
BIOGRID=$GRAPHS/BioGRID_unweighted_graph
DIP=$GRAPHS/DIP_unweighted_graph

#default params
threshold=0.2
alpha=1
gamma=0.75
mode="sr"
redundancy=0.6
penalty=1.25
quality=0

for i in $(seq 0 0.1 1)
do
    #=================Biogrid=====================#
    #threshold
    $MCLCAW $BIOGRID $i $alpha $gamma $mode $redundancy $penalty $quality
    #quality
    $MCLCAW $BIOGRID $threshold $alpha $gamma $mode $redundancy $penalty $i
    #redundancy
    $MCLCAW $BIOGRID $threshold $alpha $gamma $mode $i $penalty $quality
    #penalty
    pen=$i*2
    $MCLCAW $BIOGRID $threshold $alpha $gamma $mode $redundancy $pen $quality

    #=================DIP=====================#
    #threshold
    $MCLCAW $DIP $i $alpha $gamma $mode $redundancy $penalty $quality
    #quality
    $MCLCAW $DIP $threshold $alpha $gamma $mode $redundancy $penalty $i
    #redundancy
    $MCLCAW $DIP $threshold $alpha $gamma $mode $i $penalty $quality
    #penalty
    pen=$i*2
    $MCLCAW $DIP $threshold $alpha $gamma $mode $redundancy $pen $quality
done
