import csv
import math
import argparse
import numpy as np
import scipy as sp
import networkx as nx

ap = argparse.ArgumentParser()
ap.add_argument("--assembly_graph", required=True, help="path to the assembly graph file")
ap.add_argument("--PE_graph", required=True, help="path to the PE graph file")

args = vars(ap.parse_args())
assembly_graph_file = args["assembly_graph"]
PE_graph_file = args["PE_graph"]

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
    elif line[-1] == ';' and strings[0] != "Size":
        PE_graph.add_node(strings[0])
        PE_graph.add_edge(contig, strings[0])
    line = graph.readline()
graph.close()

merged_graph = nx.Graph()
merged_graph.add_nodes_from(assembly_graph.nodes)
merged_graph.add_nodes_from(PE_graph.nodes)
merged_graph.add_edges_from(assembly_graph.edges)
merged_graph.add_edges_from(PE_graph.edges)

print("assembly graph")
Max = 0
Sum = 0
compo = nx.connected_components(assembly_graph)
for CP in compo:
    if len(CP) > Max:
        Max = len(CP)
        maxcomp = CP
    Sum += len(CP)
print(Max, Sum, Max/Sum)
print(len(assembly_graph.edges(maxcomp)))
# print([len(c) for c in sorted(nx.connected_components(assembly_graph), key=len, reverse=True)])

print("PE graph")
Max = 0
Sum = 0
compo = nx.connected_components(PE_graph)
for CP in compo:
    if len(CP) > Max:
        Max = len(CP)
        maxcomp = CP
    Sum += len(CP)
print(Max, Sum, Max/Sum)
print(len(PE_graph.edges(maxcomp)))

print("merged graph")
Max = 0
Sum = 0
compo = nx.connected_components(merged_graph)
for CP in compo:
    if len(CP) > Max:
        Max = len(CP)
        maxcomp = CP
    Sum += len(CP)
print(Max, Sum, Max/Sum)
print(len(merged_graph.edges(maxcomp)))