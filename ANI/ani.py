from sklearn import metrics
input = open("ari.txt", "r")
l1 = []
l2 = []
for line in input:
	line = line.strip().split()
	l1.append(int(line[0]))
	l2.append(line[1])

input.close()
print(metrics.adjusted_rand_score(l1, l2))
