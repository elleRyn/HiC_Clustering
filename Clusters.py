#!/usr/bin/env
#script to separate clusters
#designed by Kohler 2016

#takes output from HiCclust.R script (clusters.txt) as input
#sys.argv[1] i.e. clusters.txt
#sys.argv[2] i.e. longcontigs.fasta
import sys

list_count = 0

#creates an empty list for each cluster
with open (sys.argv[1],'r') as infile:
	for line in infile:
		contig = line.split()[0]
		cluster = line.split()[1]
		if cluster >= list_count:
			list_count = cluster
#outputs how many clusters there are
sys.stdout.write('Number of clusters: ' + list_count + '\n')
lists = [[] for _ in range(int(list_count))]
		

#assigns contig names from clusters.txt 
with open (sys.argv[1],'r') as infile:
	for line in infile:
		contig = line.split()[0]
		cluster = line.split()[1]
		lists[int(cluster)-1].append(contig)
sys.stdout.write('Number of contigs in each cluster: 1-' + list_count + '\n')
for list in lists:
	print len(list)
	

#define this as a function
def contig_grab(name, contig_file, fn):
	print_flag = False
	with open (contig_file,'r') as infile:
		for line in infile:
			if name in line:
				print_flag = True
			elif line[0] == '>':
				print_flag = False
			if print_flag:
				fn.write(line)


#use function on data file to produce one contig file for each genomic cluster		
counter = 0
for list in lists:
	counter += 1
	with open('cluster' + str(counter) + '.fasta', 'a') as fn:
		for name in list:
			contig_grab(name, sys.argv[2], fn)


#not all contigs may have had valid Hi-C reads aligning to them, especially short 
#contigs (~600bp), so they will not show up in any of the individual cluster files as they 
#will not have linked with any of the other contigs