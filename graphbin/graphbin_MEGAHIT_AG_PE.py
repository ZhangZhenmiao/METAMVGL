#!/usr/bin/python3

"""graphbin_MEGAHIT.py: Improved binning of metagenomic contigs using MEGAHIT assembly graphs.

GraphBin is a metagenomic contig binning tool that makes use of the contig 
connectivity information from the assembly graph to bin contigs. It utilizes 
the binning result of an existing binning tool and a label propagation algorithm 
to correct mis-binned contigs and predict the labels of contigs which are 
discarded due to short length.

graphbin_MEGAHIT.py makes use of the assembly graphs produced by MEGAHIT assembler.
"""

import sys
import os
import subprocess
import csv
import operator
import time
import argparse
import re
import logging

from igraph import *
from labelpropagation.labelprop import LabelProp
from bidirectionalmap.bidirectionalmap import BidirectionalMap

# Sample command
# -------------------------------------------------------------------
# python graphbin_MEGAHIT.py    --graph /path/to/graph_file.gfa
#                               --binned /path/to/binning_result.csv
#                               --output /path/to/output_folder
# -------------------------------------------------------------------


# Setup logger
#-----------------------

logger = logging.getLogger('GraphBin 1.0')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
consoleHeader = logging.StreamHandler()
consoleHeader.setFormatter(formatter)
logger.addHandler(consoleHeader)

start_time = time.time()

# Setup argument parser
#-----------------------

ap = argparse.ArgumentParser()

ap.add_argument("--graph1", required=True, help="path to the assembly graph file")
ap.add_argument("--graph2", required=True, help="path to the assembly graph file")
ap.add_argument("--binned", required=True, help="path to the .csv file with the initial binning output from an existing tool")
ap.add_argument("--output", required=False, help="path to the output folder")
ap.add_argument("--prefix", required=False, help="prefix for the output file")
ap.add_argument("--max_iteration", required=False, type=int, help="maximum number of iterations for label propagation algorithm. [default: 100]")
ap.add_argument("--diff_threshold", required=False, type=float, help="difference threshold for label propagation algorithm. [default: 0.1]")

args = vars(ap.parse_args())

assembly_graph_file = args["graph1"]
assembly_graph_file2 = args["graph2"]
contig_bins_file = args["binned"]
output_path = args["output"]
prefix = args["prefix"]
max_iteration = args["max_iteration"]
diff_threshold = args["diff_threshold"]


# Setup output path for log file
#---------------------------------------------------

fileHandler = logging.FileHandler(output_path+"/"+prefix+"graphbin.log")
fileHandler.setLevel(logging.INFO)
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)

logger.info("Welcome to GraphBin: Improved Binning of Metagenomic Contigs using Assembly Graphs.")
logger.info("This version of GraphBin makes use of the assembly graph produced by MEGAHIT which is based on the de Bruijn graph approach.")

logger.info("Assembly graph file: "+assembly_graph_file)
logger.info("Existing binning output file: "+contig_bins_file)
logger.info("Final binning output file: "+output_path)
logger.info("Maximum number of iterations: "+str(max_iteration))
logger.info("Difference threshold: "+str(diff_threshold))

logger.info("GraphBin started")


# Get the number of bins from the initial binning result
#--------------------------------------------------------

n_bins = 0

try:
    all_bins_list = []

    with open(contig_bins_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            all_bins_list.append(row[1])
            
    bins_list = list(set(all_bins_list))
    bins_list.sort()

    # n_bins = len(bins_list)
    for i in bins_list:
        if n_bins < int(i): n_bins = int(i)
    logger.info("Number of bins available in the initial binning result: "+str(n_bins))
except:
    logger.error("Please make sure that the correct path to the initial binning result file is provided and it is having the correct format.")
    logger.info("Exiting GraphBin... Bye...!")
    sys.exit(1)


logger.info("Constructing the assembly graph")


## Construct the assembly graph
#-------------------------------

node_count = 0

links = []

my_map = BidirectionalMap()

try:

    # Get links from .gfa file
    with open(assembly_graph_file) as file:
        line = file.readline()
        while line != "":
            line = line.strip();
            strings = line[:-1].split();
            if line[-1] == ':':
                contig = strings[0]
                if contig not in my_map.values():
                    my_map[node_count] = contig
                    node_count += 1
            elif line[-1] == ';':
                if strings[0] != "Size" and [contig, strings[0]] not in links:
                    links.append([contig, strings[0]])
            line = file.readline()

    with open(assembly_graph_file2) as file:
        line = file.readline()
        while line != "":
            line = line.strip();
            strings = line[:-1].split();
            if line[-1] == ':':
                contig = strings[0]
                bin1 = strings[1]
                if contig not in my_map.values():
                    my_map[node_count] = contig
                    node_count += 1
            elif line[-1] == ';':
                if strings[0] != "Size" and [contig, strings[0]] not in links:
                    bin2 = strings[1]
                    if bin1 == "bin0" or bin2 == "bin0":
                        links.append([contig, strings[0]])
                    elif bin1 == bin2:
                        links.append([contig, strings[0]])
            line = file.readline()

    logger.info("Total number of contigs available: "+str(node_count))

    contigs_map = my_map
    contigs_map_rev = my_map.inverse

    # Create graph
    assembly_graph = Graph()

    # Add vertices
    assembly_graph.add_vertices(node_count)

    # Create list of edges
    edge_list = []

    for i in range(node_count):
        assembly_graph.vs[i]["id"]= i
        assembly_graph.vs[i]["label"]= str(i)
        
    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append((contigs_map_rev[link[0]], contigs_map_rev[link[1]]))

    # Add edges to the graph 
    assembly_graph.add_edges(edge_list)
    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

except:
    logger.error("Please make sure that the correct path to the assembly graph file is provided.")
    logger.info("Exiting GraphBin... Bye...!")
    sys.exit(1)

logger.info("Total number of edges in the assembly graph: "+str(len(edge_list)))

# n_test = assembly_graph.neighbors(contigs_map_rev['k119_2740'], mode='ALL')
# print(n_test)
# n_test = [contigs_map[i] for i in n_test]
# print(n_test)
# n_test = assembly_graph.neighbors(contigs_map_rev['k119_4957'], mode='ALL')
# print(n_test)
# n_test = [contigs_map[i] for i in n_test]
# print(n_test)
# n_test = assembly_graph.neighbors(contigs_map_rev['k119_2644'], mode='ALL')
# print(n_test)
# n_test = [contigs_map[i] for i in n_test]
# print(n_test)
# n_test = assembly_graph.neighbors(contigs_map_rev['k119_326'], mode='ALL')
# print(n_test)
# n_test = [contigs_map[i] for i in n_test]
# print(n_test)
# n_test = assembly_graph.neighbors(contigs_map_rev['k119_7109'], mode='ALL')
# print(n_test)
# n_test = [contigs_map[i] for i in n_test]
# print(n_test)
# print(contigs_map_rev['k119_18477'])
# Get initial binning result
#----------------------------

logger.info("Obtaining the initial binning result")

bins = [[] for x in range(n_bins)]

try:
    with open(contig_bins_file) as contig_bins:
        readCSV = csv.reader(contig_bins, delimiter=',')
        for row in readCSV:
            contig_num = row[0]
            
            bin_num = int(row[1])-1
            bins[bin_num].append(contigs_map_rev[contig_num])
            # print(bin_num, contigs_map[contigs_map_rev[contig_num]])

    for i in range(n_bins):
        bins[i].sort()

except:
    logger.error("Please make sure that you have provided the correct assembler type and the correct path to the binning result file in the correct format.")
    logger.info("Exiting GraphBin... Bye...!")
    sys.exit(1)

binned_contigs = []

for n in range(n_bins):
    binned_contigs = sorted(binned_contigs+bins[n])

logger.info("Number of all contigs: "+str(len(binned_contigs)))

# Remove labels of ambiguous vertices
#-------------------------------------

def getClosestLabelledVertices(graph, node, binned_contigs):
    queu_l = [graph.neighbors(node, mode='ALL')]
    visited_l = [node]
    labelled = []

    while len(queu_l) > 0:
        active_level = queu_l.pop(0)
        is_finish = False
        visited_l += active_level

        for n in active_level:
            if n in binned_contigs:
                is_finish = True
                labelled.append(n)
        if is_finish:
            return labelled
        else:
            temp = []
            for n in active_level:
                temp += graph.neighbors(n, mode='ALL')
            temp = list(set(temp))
            temp2 = []

            for n in temp:
                if n not in visited_l:
                    temp2.append(n)
            if len(temp2) > 0:
                queu_l.append(temp2)
    return labelled


logger.info("Determining ambiguous vertices")

remove_labels = []

neighbours_have_same_label_list = []

for b in range(n_bins):

    for i in bins[b]:

        my_bin = b

        dist = {}

        # Get set of closest labelled vertices with distance = 1
        closest_neighbours = assembly_graph.neighbors(i, mode=ALL)

        # Determine whether all the closest labelled vertices have the same label as its own
        neighbours_have_same_label = True
        
        neighbours_binned = False
        
        for neighbour in closest_neighbours:
            for k in range(n_bins):
                if neighbour in bins[k]:
                    neighbours_binned = True
                    if k != my_bin:
                        neighbours_have_same_label = False
                        break
                        
        if not neighbours_have_same_label:
            remove_labels.append(i)
        elif neighbours_binned:
            neighbours_have_same_label_list.append(i)

for i in remove_labels:
    for n in range(n_bins):
        if i in bins[n]:
            bins[n].remove(i)

# Further remove labels of ambiguous vertices
binned_contigs = []

for n in range(n_bins):
    binned_contigs = sorted(binned_contigs+bins[n])

logger.info("Number of non-ambiguous contigs: "+str(len(binned_contigs)))

for b in range(n_bins):

    for i in bins[b]:
        
        if i not in neighbours_have_same_label_list:

            my_bin = b

            # Get set of closest labelled vertices
            closest_neighbours = getClosestLabelledVertices(assembly_graph, i, binned_contigs)

            if len(closest_neighbours) > 0:

                # Determine whether all the closest labelled vertices have the same label as its own
                neighbours_have_same_label = True

                for neighbour in closest_neighbours:
                    for k in range(n_bins):
                        if neighbour in bins[k]:
                            if k != my_bin:
                                neighbours_have_same_label = False
                                break

                if not neighbours_have_same_label and i not in remove_labels:
                    remove_labels.append(i)

remove_labels.sort()

logger.info("Removing labels of ambiguous vertices")

# Remove labels of ambiguous vertices
for i in remove_labels:
    for n in range(n_bins):
        if i in bins[n]:
            bins[n].remove(i)

binned_contigs = []

for n in range(n_bins):
    binned_contigs = sorted(binned_contigs+bins[n])
logger.info("Number of non-ambiguous contigs: "+str(len(binned_contigs)))

logger.info("Obtaining refined binning result")

bin_out = open("refined_bin", 'w')
for n in range(n_bins):
    for b in bins[n]:
        bin_out.write(contigs_map[b] + ',' + str(n+1) + '\n')
bin_out.close()

# Get vertices which are not isolated and not in components without any labels
#-----------------------------------------------------------------------------

logger.info("Deteremining vertices which are not isolated and not in components without any labels")

non_isolated = []

for i in range(node_count):
    
    if i not in non_isolated and i in binned_contigs:

        component = []
        component.append(i)
        length = len(component)
        neighbours = assembly_graph.neighbors(i, mode=ALL)

        for neighbor in neighbours:
            if neighbor not in component:
                component.append(neighbor)

        component = list(set(component))

        while length!= len(component):

            length = len(component)

            for j in component:

                neighbours = assembly_graph.neighbors(j, mode=ALL)

                for neighbor in neighbours:
                    if neighbor not in component:
                        component.append(neighbor)

        labelled = False
        for j in component:
            if j in binned_contigs:
                labelled = True
                break

        if labelled:
            for j in component:
                if j not in non_isolated:
                    non_isolated.append(j)

logger.info("Number of non-isolated contigs: "+str(len(non_isolated)))

# Run label propagation
#-----------------------

data = []

for contig in range(node_count):
    
    # Consider vertices that are not isolated

    if contig in non_isolated:
        line = []
        line.append(contig)

        assigned = False

        for i in range(n_bins):
            if contig in bins[i]:
                line.append(i+1)
                assigned = True
        
        if not assigned:
            line.append(0)

        neighbours = assembly_graph.neighbors(contig, mode=ALL)

        neighs = []

        for neighbour in neighbours:
            n = []
            n.append(neighbour)
            n.append(1.0)
            neighs.append(n)

        line.append(neighs)

        data.append(line)

lp = LabelProp()

lp.load_data_from_mem(data)

logger.info("Starting label propagation with eps="+str(diff_threshold)+" and max_iteration="+str(max_iteration))

ans = lp.run(diff_threshold, max_iteration, show_log=True, clean_result=False) 
ans.sort()

logger.info("Obtaining Label Propagation result")

for l in ans:
    for i in range(n_bins):
        if l[1]==i+1 and l[0] not in bins[i]:
            bins[i].append(l[0])


# Remove labels of ambiguous vertices
#-------------------------------------

logger.info("Determining ambiguous vertices")

remove_labels = []

for b in range(n_bins):

    for i in bins[b]:

        my_bin = b

        closest_neighbours = assembly_graph.neighbors(i, mode=ALL)

        # Determine whether all the closest labelled vertices have the same label as its own
        neighbours_have_same_label = True
        
        for neighbour in closest_neighbours:
            for k in range(n_bins):
                if neighbour in bins[k]:
                    if k != my_bin:
                        neighbours_have_same_label = False
                        break
                        
        if not neighbours_have_same_label:
            remove_labels.append(i)

remove_labels.sort()

logger.info("Removing labels of ambiguous vertices")

# Remove labels of ambiguous vertices
for i in remove_labels:
    for n in range(n_bins):
        if i in bins[n]:
            bins[n].remove(i)

elapsed_time = time.time() - start_time

logger.info("Obtaining the Final Refined Binning result")

for i in range(n_bins):
    bins[i].sort()

# Print elapsed time for the process
logger.info("Elapsed time: "+str(elapsed_time)+" seconds")


# Write result to output file
#-----------------------------

logger.info("Writing the Final Binning result to file")

output_bins = []

for i in range(node_count):
    for k in range(n_bins):
        if i in bins[k]:
            line = []
            line.append(contigs_map[i])
            line.append(k+1)
            output_bins.append(line)

output_file = output_path + prefix + 'graphbin_output.csv'

with open(output_file, mode='w') as out_file:
    output_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    for row in output_bins:
        output_writer.writerow(row)

logger.info("Final binning results can be found at "+output_file)


unbinned_contigs = []

for i in range(node_count):
    if i in remove_labels or i not in non_isolated:
        line = []
        line.append(contigs_map[i])
        unbinned_contigs.append(line)

if len(unbinned_contigs)!=0:
    unbinned_file = output_path + prefix + 'graphbin_unbinned.csv'

    with open(unbinned_file, mode='w') as out_file:
        output_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        for row in unbinned_contigs:
            output_writer.writerow(row)

    logger.info("Unbinned contigs can be found at "+unbinned_file)


# Exit program
#--------------

logger.info("Thank you for using GraphBin! Bye...!")

logger.removeHandler(fileHandler)
logger.removeHandler(consoleHeader)
