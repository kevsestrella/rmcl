import numpy as np
import pandas as pd
from itertools import chain
import sys

clusters = []
with open(sys.argv[2], 'r') as f:
    for row in f:
        clusters.append([int(i) for i in row.split()])

graph = []
with open(sys.argv[1], 'r') as a:
    for row in a:
        graph.append([float(i) for i in row.split()])
V, Vnew, E = [int(i) for i in graph[0]]
graph = graph[1:]

clusters_df = pd.DataFrame(clusters)
clusters_df_filled = clusters_df.fillna(-1).astype(int)

graph_df = pd.DataFrame(graph)
graph_df_filled = graph_df.fillna(-1)
graph_df_filled.loc[:,1::2] = graph_df_filled.loc[:,1::2].astype(int)
graph_protein_index = graph_df_filled.set_index(0)

cluster_mapping = sorted(list(set(chain(*clusters))))
cluster_count = len(cluster_mapping)
cluster_size = []
for i in cluster_mapping:
    cluster_size.append(list(chain(*clusters)).count(i))

cluster_data = pd.DataFrame()
cluster_data['cluster_mapping'] = cluster_mapping
cluster_data['cluster_size'] = cluster_size

connected_proteins = graph_protein_index.iloc[:,0::2]
connection_weights = graph_protein_index.iloc[:,1::2].astype('float64')

p_inconnectivity = pd.DataFrame(np.zeros([Vnew,cluster_count]), index=connected_proteins.index, columns=cluster_mapping, dtype = 'float64')
p_outconnectivity = pd.Series(np.zeros(Vnew),index = connected_proteins.index, dtype = 'float64')

p_inconnectivityC = pd.DataFrame(np.zeros([Vnew,cluster_count]), index=connected_proteins.index, columns=cluster_mapping, dtype = 'float64')
p_outconnectivityC = pd.Series(np.zeros(Vnew),index = connected_proteins.index, dtype = 'float64')

p_interC = pd.DataFrame(np.zeros([Vnew,cluster_count]), index=connected_proteins.index, columns=cluster_mapping, dtype = 'float64')
protein_clusters = clusters_df_filled.set_index(connected_proteins.index)

# Preliminary Cores
for protein_i,rows in connected_proteins.iterrows():
    protein_i_clusters = protein_clusters.loc[protein_i][protein_clusters.loc[protein_i]!=-1]
    protein_i_weights = connection_weights.loc[protein_i][connection_weights.loc[protein_i]!=-1]
    for protein_j_index, protein_j in enumerate(rows[rows!=-1]):
        protein_j_clusters = protein_clusters.loc[protein_j][protein_clusters.loc[protein_j]!=-1]
        cluster_match = protein_i_clusters[protein_i_clusters.isin(protein_j_clusters)]
        if not cluster_match.empty:
            p_inconnectivity.loc[protein_i, cluster_match] += protein_i_weights.iloc[protein_j_index]
        else:
            p_outconnectivity.loc[protein_i] += protein_i_weights.iloc[protein_j_index]

avg_inconnectivity = p_inconnectivity.sum(axis=0)/(2.0*pd.Series(cluster_size, index = p_inconnectivity.columns, dtype = 'float64'))

isCore = pd.DataFrame((p_inconnectivity >= avg_inconnectivity).values & ((p_inconnectivity.T > p_outconnectivity).T).values, index = p_inconnectivity.index, columns= p_inconnectivity.columns)

#Extended Cores
for protein_i,rows in connected_proteins.iterrows():
    protein_i_clusters = protein_clusters.loc[protein_i][protein_clusters.loc[protein_i]!=-1]
    if((isCore.loc[protein_i,protein_i_clusters][isCore.loc[protein_i,protein_i_clusters]==False]).empty):
        continue
    protein_i_clusters = (isCore.loc[protein_i,protein_i_clusters][isCore.loc[protein_i,protein_i_clusters]==False]).index
    protein_i_weights = connection_weights.loc[protein_i][connection_weights.loc[protein_i]!=-1]
    for protein_j_index, protein_j in enumerate(rows[rows!=-1]):
        protein_j_clusters = protein_clusters.loc[protein_j][protein_clusters.loc[protein_j]!=-1]
        cluster_match = protein_i_clusters[protein_i_clusters.isin(protein_j_clusters)]
        if not cluster_match.empty:
            protein_j_clusters = (isCore.loc[protein_j,protein_j_clusters][isCore.loc[protein_j,protein_j_clusters]==True]).index
            cluster_match = protein_i_clusters[protein_i_clusters.isin(protein_j_clusters)]
            if not cluster_match.empty:
                p_inconnectivityC.loc[protein_i, cluster_match] += protein_i_weights.iloc[protein_j_index]
            else:
                p_outconnectivityC.loc[protein_i] += protein_i_weights.iloc[protein_j_index]

avg_inconnectivityC = p_inconnectivityC.sum(axis=0)/(2.0*pd.Series(cluster_size, index = p_inconnectivityC.columns, dtype = 'float64'))

isCore |= (pd.DataFrame((p_inconnectivityC >= avg_inconnectivityC).values & ((p_inconnectivityC.T > p_outconnectivityC).T).values, index = p_inconnectivityC.index, columns= p_inconnectivityC.columns))

#Attachments
cluster_inter = pd.Series(np.zeros(cluster_count), index=cluster_mapping, dtype = 'float64')

for cluster, proteins in isCore.iteritems():
    core_proteins = proteins[proteins == True].index
    for index, core_protein in enumerate(core_proteins):
        for next_core_protein in core_proteins[index+1:]:
            if next_core_protein in connected_proteins.loc[core_protein,:].values:
                cluster_inter[cluster] += connection_weights.loc[core_protein,connected_proteins.loc[core_protein,:][connected_proteins.loc[core_protein,:]==next_core_protein].index[0]+1]

for protein_i,rows in connected_proteins.iterrows():
    protein_i_clusters = protein_clusters.loc[protein_i][protein_clusters.loc[protein_i]!=-1]
    if((isCore.loc[protein_i,protein_i_clusters][isCore.loc[protein_i,protein_i_clusters]==False]).empty):
        continue
    protein_i_clusters = (isCore.loc[protein_i,protein_i_clusters][isCore.loc[protein_i,protein_i_clusters]==False]).index
    protein_i_weights = connection_weights.loc[protein_i][connection_weights.loc[protein_i]!=-1]
    for protein_j_index, protein_j in enumerate(rows[rows!=-1]):
        protein_j_clusters = protein_clusters.loc[protein_j][protein_clusters.loc[protein_j]!=-1]
        cluster_match = protein_i_clusters[protein_i_clusters.isin(protein_j_clusters)]
        if not cluster_match.empty:
            protein_j_clusters = (isCore.loc[protein_j,protein_j_clusters][isCore.loc[protein_j,protein_j_clusters]==True]).index
            cluster_match = protein_i_clusters[protein_i_clusters.isin(protein_j_clusters)]
            if not cluster_match.empty:
                p_interC.loc[protein_i, cluster_match] += protein_i_weights.iloc[protein_j_index]

alpha = 1.0
gamma = 0.75
isAttachment = (p_interC>=alpha*cluster_inter*(np.array(cluster_size, dtype = 'float64')/2.0)**-gamma) #the cluster_size here is based on original code which differs with the definition on the paper

with open('prettycaw.out', 'w+') as prettycaw:
    with open(sys.argv[3], 'w+') as caw:
        cluster_count = 0
        for cluster in cluster_mapping:
            if (isCore.loc[:,cluster].sum() > 0):
                if (isCore.loc[:,cluster].sum() > 0 and (isCore.loc[:,cluster].sum() + isAttachment.loc[p_interC.loc[:,cluster][p_interC.loc[:,cluster]!=0].index,cluster][isAttachment.loc[p_interC.loc[:,cluster][p_interC.loc[:,cluster]!=0].index,cluster]==True].count())>=4 ):
                    prettycaw.write("---------Cluster {}-----------\n".format(cluster_count+1) )
                    prettycaw.write("Core Proteins: {}\n".format(isCore.loc[:,cluster][isCore.loc[:,cluster]==True].index.astype(int).tolist()))
                    caw.write(" ".join(isCore.loc[:,cluster][isCore.loc[:,cluster]==True].index.astype(int).astype(str).tolist()))
                    if not isAttachment.loc[:,cluster][isAttachment.loc[:,cluster]==True].empty:
                        prettycaw.write("Attachment Proteins: {}\n".format(isAttachment.loc[:,cluster][isAttachment.loc[:,cluster]==True].index.tolist()))
                        caw.write(" "+" ".join(isAttachment.loc[p_interC.loc[:,cluster][p_interC.loc[:,cluster]!=0].index,cluster][isAttachment.loc[p_interC.loc[:,cluster][p_interC.loc[:,cluster]!=0].index,cluster]==True].index.astype(int).astype(str).tolist()))
                    caw.write("\n")
                    prettycaw.write("\n\n")
                    cluster_count+=1
