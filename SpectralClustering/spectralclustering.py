import argparse
import numpy as np
from sklearn.cluster import SpectralClustering

ap = argparse.ArgumentParser()
ap.add_argument("--graph", required=True, help="path to the graph file")
ap.add_argument("--output", required=False, help="path to the output file")
args = vars(ap.parse_args())
graph_file = args["graph"]
output_file = args["output"]

contig_cnt = 0
contigs = []
contig_idx = {}
contig_lnk = {}
with open(graph_file) as graph:
    line = graph.readline()
    while line != "":
        if contig_cnt and contig_cnt%1000 == 0:
            print("Read", contig_cnt, "lines ...", flush = True)
        line = line.strip();
        if line[-1] == ':':
            contig = line.split()[0]
            if contig not in contigs:
                contigs.append(contig)
                contig_idx[contig] = contig_cnt
                contig_lnk[contig] = []
                contig_cnt += 1
        elif line[-1] == ';':
            item = line.split()[0]
            if item == "Size":
                line = graph.readline()
                continue
            if item not in contigs:
                contigs.append(item)
                contig_idx[item] = contig_cnt
                contig_lnk[item] = []
                contig_cnt += 1
            contig_lnk[contig].append(item)
            contig_lnk[item].append(contig)
        line = graph.readline()

print("All contigs:", contig_cnt, flush = True)
affinity_matrix = np.zeros((contig_cnt, contig_cnt))
for ctg in contig_lnk:
    for lnk in contig_lnk[ctg]:
        affinity_matrix[contig_idx[ctg]][contig_idx[lnk]] = 1

clustering = SpectralClustering(
    n_clusters=12,
    assign_labels="kmeans", 
    n_init=20, 
    random_state=0, 
    affinity="precomputed",
    n_jobs=250
).fit(affinity_matrix)

cluster_out = open(output_file, "w")
for i in range(clustering.labels_.size):
    cluster_out.write(contigs[i] + "," + str(clustering.labels_[i]) + "\n")
cluster_out.close()