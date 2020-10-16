import os
headmap = dict()
mapfile = open('mycc/Alias.txt', 'r')
for line in mapfile:
	headmap[line.strip().split('\t')[1]] = line.strip().split()[0]
for file in os.listdir('mycc'):
	if 'fa' not in file: continue
	contig = open('mycc/' + file, 'r')
	out_contig = open('mycc/modified_' + file, 'w')
	for line in contig:
		if line[0] == '>': 
			out_contig.write(headmap[line.strip()] + '\n')
		else:
			out_contig.write(line)
	contig.close()
	out_contig.close()
os.system('rm mycc/Cluster*.fa')