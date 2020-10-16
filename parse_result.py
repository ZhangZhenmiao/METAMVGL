import os
import sys
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("--cluster", required=True, help="path to cluster")
ap.add_argument("--contigs", required=True, help="path to contigs")
ap.add_argument("--output", required=True, help="path to output dir")
ap.add_argument("--assembler", required=True, help="assembler used")

args = vars(ap.parse_args())
cluster = args["cluster"]
contigs = args["contigs"]
output = args["output"]
assembler = args["assembler"]

contigs_file = open(contigs, 'r')
contigs_map = {}
header = ""
content = ""
for line in contigs_file:
	if line == "": continue
	if line[0] == '>':
		if header != "": contigs_map[header] = content
		if assembler == 'spades': header = line.split('_')[0][1:] + '_' + line.split('_')[1]
		elif assembler == 'megahit': header = line.split()[0][1:]
		content = ""
	else: content += line.strip()
contigs_map[header] = content
contigs_file.close()

bin_map = {}
cluster_file = open(cluster, 'r')
for line in cluster_file:
	if line == "": continue
	items = line.strip().split(',')
	if items[1] not in bin_map: bin_map[items[1]] = []
	bin_map[items[1]].append(items[0])
cluster_file.close()

if not os.path.isdir(output): os.system("mkdir " + output)
for file in os.listdir(output):
	if ".fasta" in file: os.system("rm " + output + '/' + file)
for bin in bin_map:
	out = open(output + "/cluster." + bin + ".fasta", 'w')
	for header in bin_map[bin]:
		out.write('>' + header + '\n' + contigs_map[header] + '\n')
	out.close()