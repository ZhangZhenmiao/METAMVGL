import csv
import math
import argparse
import numpy as np
import scipy as sp
import networkx as nx

def remove_ambiguous_label(G, contigs_bin):
    remove_labels = list()
    for key in contigs_bin:
        closest_neighbours = G.neighbors(key)
        neighbours_have_same_label = True
        for neighbour in closest_neighbours:
            if neighbour in contigs_bin:
                if contigs_bin[key] != contigs_bin[neighbour]:
                    neighbours_have_same_label = False
                    break
        if not neighbours_have_same_label:
            remove_labels.append(key)
    for key in remove_labels:
        contigs_bin.pop(key)

def getClosestLabelledVertices(G, node, contigs_bin):
    queu_l = [list(G.neighbors(node))]
    visited_l = [node]
    labelled = []
    while len(queu_l) > 0:
        active_level = queu_l.pop(0)
        is_finish = False
        visited_l += active_level
        for n in active_level:
            if n in contigs_bin:
                is_finish = True
                labelled.append(n)
        if is_finish:
            return labelled
        else:
            temp = []
            for n in active_level:
                temp += list(G.neighbors(n))
            temp = list(set(temp))
            temp2 = []
            for n in temp:
                if n not in visited_l:
                    temp2.append(n)
            if len(temp2) > 0:
                queu_l.append(temp2)
    return labelled

def remove_ambiguous_label_deeper(G, contigs_bin):
    remove_labels = list()
    for key in contigs_bin:
        closest_neighbours = getClosestLabelledVertices(G, key, contigs_bin)
        if not len(closest_neighbours): continue
        neighbours_have_same_label = True
        for neighbour in closest_neighbours:
            if neighbour in contigs_bin:
                if contigs_bin[key] != contigs_bin[neighbour]:
                    neighbours_have_same_label = False
                    break
        if not neighbours_have_same_label:
            remove_labels.append(key)
    for key in remove_labels:
        contigs_bin.pop(key)

def non_isolated_contigs(G, contigs_bin):
    non_isolated = set()
    for i in G.nodes:
        if i not in non_isolated and i in contigs_bin:
            if not len(list(G.neighbors(i))): continue
            component = set()
            component.add(i)
            length = len(component)
            for neighbor in G.neighbors(i): component.add(neighbor)
            while length!= len(component):
                length = len(component)
                for j in component.copy():
                    for neighbor in G.neighbors(j): component.add(neighbor)
            for j in component: non_isolated.add(j)
    return_list = list()
    for item in non_isolated:
        if item in contigs_bin: return_list.append(item)
    for item in non_isolated:
        if item not in return_list: return_list.append(item)
    return return_list

ap = argparse.ArgumentParser()
ap.add_argument("--assembly_graph", required=True, help="path to the assembly graph file")
ap.add_argument("--PE_graph", required=True, help="path to the PE graph file")
ap.add_argument("--binned", required=True, help="path to the .csv file with the initial binning output from an existing tool")
ap.add_argument("--max_iter", required=True, help="max iteration of AMGL algorithm")
ap.add_argument("--thresh", required=True, help="stop threshold of AMGL algorithm")
ap.add_argument("--remove", required=True, help="rmeove ambiguous contigs or not")

args = vars(ap.parse_args())
assembly_graph_file = args["assembly_graph"]
PE_graph_file = args["PE_graph"]
contig_bins_file = args["binned"]
max_iter = int(args["max_iter"])
thresh = float(args["thresh"])
remove = args["remove"]

all_bins = set() # bin list, start from 0
contigs_bin = dict() # contig to bin, only contain binned contigs

csvfile = open(contig_bins_file, 'r')
readCSV = csv.reader(csvfile, delimiter=',')
for row in readCSV:
    all_bins.add(int(row[1])-1)
    contigs_bin[row[0]] = int(row[1])-1
csvfile.close()
# n_bins = len(all_bins)
n_bins = 0
for i in all_bins:
    if i > n_bins: n_bins = i
n_bins += 1

assembly_graph = nx.Graph()
graph = open(assembly_graph_file, 'r')
line = graph.readline()
while line != "":
    line = line.strip();
    strings = line[:-1].split();
    if line[-1] == ':':
        contig = strings[0]
        assembly_graph.add_node(contig)
    elif line[-1] == ';' and strings[0] != "Size":
        assembly_graph.add_node(strings[0])
        assembly_graph.add_edge(contig, strings[0])
    line = graph.readline()
graph.close()

PE_graph = nx.Graph()
graph = open(PE_graph_file, 'r')
line = graph.readline()
while line != "":
    line = line.strip();
    strings = line[:-1].split();
    if line[-1] == ':':
        contig = strings[0]
        PE_graph.add_node(contig)
        bin_link_left = int(strings[1][3:]) - 1
    elif line[-1] == ';' and strings[0] != "Size":
        assembly_graph.add_node(strings[0])
        bin_link_right = int(strings[1][3:]) - 1
        if bin_link_left == -1 or bin_link_right == -1 or bin_link_left == bin_link_right:
            PE_graph.add_edge(contig, strings[0])
    line = graph.readline()
graph.close()

merged_graph = nx.Graph()
merged_graph.add_nodes_from(assembly_graph.nodes)
merged_graph.add_nodes_from(PE_graph.nodes)
merged_graph.add_edges_from(assembly_graph.edges)
merged_graph.add_edges_from(PE_graph.edges)

print('initial binned contigs:', len(contigs_bin))
print('nodes of assembly graph:', len(assembly_graph.nodes))
print('edges of assembly graph:', len(assembly_graph.edges))
print('nodes of PE graph:', len(PE_graph.nodes))
print('edges of PE graph:', len(PE_graph.edges))
print('nodes of merged graph:', len(merged_graph.nodes))
print('edges of merged graph:', len(merged_graph.edges))

# remove_ambiguous_label(merged_graph, contigs_bin)
# print('binned contigs after remove ambiguous:', len(contigs_bin))
if remove == "True":
    remove_ambiguous_label_deeper(assembly_graph, contigs_bin)
print('binned contigs after remove ambiguous:', len(contigs_bin))

non_isolated = non_isolated_contigs(merged_graph, contigs_bin)
print('non isolated contigs:', len(non_isolated))

binned_cnt = 0
for contig in non_isolated:
    if contig in contigs_bin: binned_cnt += 1
print('non isolated binned contigs:', binned_cnt)

# degree = list()
# for i in range(len(non_isolated)):
#     degree.append(merged_graph.degree[non_isolated[i]])
# merged_graph_degree = np.diag(degree)
# merged_graph_adjacent = nx.adjacency_matrix(merged_graph, nodelist=non_isolated).A
# merged_graph_trans = np.dot(np.linalg.inv(merged_graph_degree), merged_graph_adjacent)
# merged_graph_trans_uu = merged_graph_trans[binned_cnt:, binned_cnt:]
# merged_graph_trans_ul = merged_graph_trans[binned_cnt:, :binned_cnt]
# F = np.zeros([len(non_isolated), n_bins])
# for i in range(len(non_isolated)):
#     if non_isolated[i] in contigs_bin: F[i, contigs_bin[non_isolated[i]]] = 1
# F_l = F[:binned_cnt,]

# F_u = np.dot(np.dot(np.linalg.inv(np.eye(merged_graph_trans_uu.shape[0]) - merged_graph_trans_uu), 
#         merged_graph_trans_ul), F_l)
# F = np.concatenate((F_l, F_u), axis=0)

degree = list()
for i in range(len(non_isolated)):
    degree.append(assembly_graph.degree[non_isolated[i]])
assembly_graph_degree = np.diag(degree)
assembly_graph_adjacent = nx.adjacency_matrix(assembly_graph, nodelist=non_isolated).A

degree = list()
for i in range(len(non_isolated)):
    degree.append(PE_graph.degree[non_isolated[i]])
PE_graph_degree = np.diag(degree)
PE_graph_adjacent = nx.adjacency_matrix(PE_graph, nodelist=non_isolated).A

F = np.zeros([len(non_isolated), n_bins])
for i in range(len(non_isolated)):
    if non_isolated[i] in contigs_bin: F[i, contigs_bin[non_isolated[i]]] = 1
F_l = F[:binned_cnt,]

assembly_graph_L = nx.normalized_laplacian_matrix(assembly_graph, nodelist=non_isolated).A
PE_graph_L = nx.normalized_laplacian_matrix(PE_graph, nodelist=non_isolated).A

Obj_fun = list()
alpha = np.array([0.5, 0.5], dtype=np.float64)
for i in range(max_iter):

    all_degree = alpha[0]*assembly_graph_degree + alpha[1]*PE_graph_degree
    all_adjacant = alpha[0]*assembly_graph_adjacent + alpha[1]*PE_graph_adjacent
    all_trans = np.dot(np.linalg.inv(all_degree), all_adjacant)
    all_trans_uu = all_trans[binned_cnt:, binned_cnt:]
    all_trans_ul = all_trans[binned_cnt:, :binned_cnt]
    F_u = np.dot(np.dot(np.linalg.inv(np.eye(all_trans_uu.shape[0]) - all_trans_uu), all_trans_ul), F_l)
    
    # L = alpha[0]*assembly_graph_L + alpha[1]*PE_graph_L
    # L_uu = L[binned_cnt:, binned_cnt:]
    # L_ul = L[binned_cnt:, :binned_cnt]
    # F_u = -0.5 * np.dot(np.dot(np.linalg.inv(L_uu), L_ul), F_l)
    
    F = np.concatenate((F_l, F_u), axis=0)
    alpha[0] = 0.5/math.sqrt(np.trace(np.dot(np.dot(F.T,assembly_graph_L),F)))
    alpha[1] = 0.5/math.sqrt(np.trace(np.dot(np.dot(F.T,PE_graph_L),F)))
    obj = math.sqrt(np.trace(np.dot(np.dot(F.T,assembly_graph_L),F)))
    obj += math.sqrt(np.trace(np.dot(np.dot(F.T,PE_graph_L),F)))
    Obj_fun.append(obj)
    print("Iteration", i, " Alpla0", alpha[0], " Alpla1", alpha[1], " Obj_value", obj)
    if i >= 1 and (Obj_fun[i-1]-Obj_fun[i])/Obj_fun[i-1] < thresh:
        break

maxCluster = np.argmax(F, axis=1)
maxValue = np.ndarray.max(F, axis=1)
for i in range(len(non_isolated)):
    if maxValue[i] != 0: contigs_bin[non_isolated[i]] = maxCluster[i]
print('final binned contigs:', len(contigs_bin))
if remove == "True":
    remove_ambiguous_label(merged_graph, contigs_bin)
print('final binned contigs after remove ambiguous:', len(contigs_bin))
final_out = open("final.cluster", 'w')
for contig in contigs_bin:
    final_out.write(contig + ',' + str(contigs_bin[contig] + 1) + '\n')
final_out.close()