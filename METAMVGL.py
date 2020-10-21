#!/usr/bin/env python
import csv
import os
import math
import argparse
import numpy as np
import scipy as sp
import networkx as nx
import logging


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
        if not len(closest_neighbours):
            continue
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
            if not len(list(G.neighbors(i))):
                continue
            component = set()
            component.add(i)
            length = len(component)
            for neighbor in G.neighbors(i):
                component.add(neighbor)
            while length != len(component):
                length = len(component)
                for j in component.copy():
                    for neighbor in G.neighbors(j):
                        component.add(neighbor)
            for j in component:
                non_isolated.add(j)
    return_list = list()
    for item in non_isolated:
        if item in contigs_bin:
            return_list.append(item)
    for item in non_isolated:
        if item not in return_list:
            return_list.append(item)
    return return_list


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True, help="path to contigs file")
    ap.add_argument("--assembler", required=True, help="assembler used")
    ap.add_argument("--assembly_graph", required=True, help="path to the .ag file")
    ap.add_argument("--PE_graph", required=True, help="path to the .pe file")
    ap.add_argument(
        "--binned", required=True, help="path to the .csv file as initial binning"
    )
    ap.add_argument("--max_iter", type=int, default=100, help="max iteration (default 100)")
    ap.add_argument(
        "--thresh", default=0.00000001, type=float, help="stop threshold (default 0.00000001)"
    )
    ap.add_argument("--output", required=True, help="output folder")
    ap.add_argument("--quiet", action="store_true", help="supress information", default=False)
    ap.add_argument("--debug", action="store_true", help="Make it more verbose", default=False)

    args = ap.parse_args()

    # define logging
    logLevel = logging.INFO
    if args.quiet:
        logLevel = logging.WARNING
    elif args.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%Y-%m-%dT%H:%M:%S+00:00: ", level=logLevel)

    logging.info("Launching METAMVGL")

    contigs = args.contigs
    assembler = args.assembler
    assembly_graph_file = args.assembly_graph
    PE_graph_file = args.PE_graph
    contig_bins_file = args.binned
    max_iter = args.max_iter
    thresh = args.thresh
    output = args.output

    logging.debug("Loading bin information from csv")
    all_bins = set()
    contigs_bin = dict()
    csvfile = open(contig_bins_file, "r")
    readCSV = csv.reader(csvfile, delimiter=",")
    for row in readCSV:
        all_bins.add(int(row[1]) - 1)
        contigs_bin[row[0]] = int(row[1]) - 1
    csvfile.close()
    n_bins = 0
    for i in all_bins:
        if i > n_bins:
            n_bins = i
    n_bins += 1

    logging.debug("Loading assembly graph")
    assembly_graph = nx.Graph()
    graph = open(assembly_graph_file, "r")
    line = graph.readline()
    while line != "":
        line = line.strip()
        strings = line[:-1].split()
        if line[-1] == ":":
            contig = strings[0]
            assembly_graph.add_node(contig)
        elif line[-1] == ";":
            assembly_graph.add_edge(contig, strings[0])
        line = graph.readline()
    graph.close()

    logging.debug("Loading paired-end graph")
    PE_graph = nx.Graph()
    graph = open(PE_graph_file, "r")
    line = graph.readline()
    while line != "":
        line = line.strip()
        strings = line[:-1].split()
        if line[-1] == ":":
            contig = strings[0]
            PE_graph.add_node(contig)
            if contig in contigs_bin:
                bin_link_left = contigs_bin[contig]
            else:
                bin_link_left = -1
        elif line[-1] == ";":
            PE_graph.add_node(strings[0])
            if strings[0] in contigs_bin:
                bin_link_right = contigs_bin[strings[0]]
            else:
                bin_link_right = -1
            if (
                bin_link_left == -1
                or bin_link_right == -1
                or bin_link_left == bin_link_right
            ):
                PE_graph.add_edge(contig, strings[0])
        line = graph.readline()
    graph.close()

    logging.debug("Merging paired-end with assembly graph")
    merged_graph = nx.Graph()
    logging.debug("Adding assembly graph nodes")
    merged_graph.add_nodes_from(assembly_graph.nodes)
    logging.debug("Adding PE graph nodes")
    merged_graph.add_nodes_from(PE_graph.nodes)
    logging.debug("Adding assembly graph edges")
    merged_graph.add_edges_from(assembly_graph.edges)
    logging.debug("Adding PE graph edges")
    merged_graph.add_edges_from(PE_graph.edges)

    logging.info("Initial binned contigs: {}".format(len(contigs_bin)))
    remove_ambiguous_label_deeper(assembly_graph, contigs_bin)
    logging.info("binned contigs after remove ambiguous: {}".format(len(contigs_bin)))
    non_isolated = non_isolated_contigs(merged_graph, contigs_bin)
    print("non isolated contigs:", len(non_isolated))
    binned_cnt = 0
    for contig in non_isolated:
        if contig in contigs_bin:
            binned_cnt += 1
    print("non isolated binned contigs:", binned_cnt)

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
        if non_isolated[i] in contigs_bin:
            F[i, contigs_bin[non_isolated[i]]] = 1
    F_l = F[
        :binned_cnt,
    ]
    assembly_graph_L = nx.normalized_laplacian_matrix(
        assembly_graph, nodelist=non_isolated
    ).A
    PE_graph_L = nx.normalized_laplacian_matrix(PE_graph, nodelist=non_isolated).A

    Obj_fun = list()
    alpha = np.array([0.5, 0.5], dtype=np.float64)
    for i in range(max_iter):
        all_degree = alpha[0] * assembly_graph_degree + alpha[1] * PE_graph_degree
        all_adjacant = alpha[0] * assembly_graph_adjacent + alpha[1] * PE_graph_adjacent
        all_trans = np.dot(np.linalg.inv(all_degree), all_adjacant)
        all_trans_uu = all_trans[binned_cnt:, binned_cnt:]
        all_trans_ul = all_trans[binned_cnt:, :binned_cnt]
        F_u = np.dot(
            np.dot(
                np.linalg.inv(np.eye(all_trans_uu.shape[0]) - all_trans_uu), all_trans_ul
            ),
            F_l,
        )
        F = np.concatenate((F_l, F_u), axis=0)
        alpha[0] = 0.5 / math.sqrt(np.trace(np.dot(np.dot(F.T, assembly_graph_L), F)))
        alpha[1] = 0.5 / math.sqrt(np.trace(np.dot(np.dot(F.T, PE_graph_L), F)))
        obj = math.sqrt(np.trace(np.dot(np.dot(F.T, assembly_graph_L), F)))
        obj += math.sqrt(np.trace(np.dot(np.dot(F.T, PE_graph_L), F)))
        Obj_fun.append(obj)
        print("Iteration", i, " Alpla0", alpha[0], " Alpla1", alpha[1], " Obj_value", obj)
        if i >= 1 and (Obj_fun[i - 1] - Obj_fun[i]) / Obj_fun[i - 1] < thresh:
            break

    maxCluster = np.argmax(F, axis=1)
    maxValue = np.ndarray.max(F, axis=1)
    for i in range(len(non_isolated)):
        if maxValue[i] != 0:
            contigs_bin[non_isolated[i]] = maxCluster[i]
    print("final binned contigs:", len(contigs_bin))
    remove_ambiguous_label(merged_graph, contigs_bin)
    print("final binned contigs after remove ambiguous:", len(contigs_bin))
    if not os.path.isdir(output):
        os.system("mkdir " + output)
    final_out = open(output + "/binning_result.csv", "w")
    for contig in contigs_bin:
        final_out.write(contig + "," + str(contigs_bin[contig] + 1) + "\n")
    final_out.close()

    contigs_file = open(contigs, "r")
    contigs_map = {}
    header = ""
    content = ""
    for line in contigs_file:
        if line == "":
            continue
        if line[0] == ">":
            if header != "":
                contigs_map[header] = content
            if assembler == "metaSPAdes":
                header = line.split("_")[0][1:] + "_" + line.split("_")[1]
            elif assembler == "MEGAHIT":
                header = line.split()[0][1:]
            content = ""
        else:
            content += line.strip()
    contigs_map[header] = content
    contigs_file.close()

    bin_map = {}
    cluster = output + "/binning_result.csv"
    cluster_file = open(cluster, "r")
    for line in cluster_file:
        if line == "":
            continue
        items = line.strip().split(",")
        if items[1] not in bin_map:
            bin_map[items[1]] = []
        bin_map[items[1]].append(items[0])
    cluster_file.close()

    for file in os.listdir(output):
        if ".fasta" in file:
            os.system("rm " + output + "/" + file)
    for bin in bin_map:
        out = open(output + "/cluster." + bin + ".fasta", "w")
        for header in bin_map[bin]:
            out.write(">" + header + "\n" + contigs_map[header] + "\n")
        out.close()
