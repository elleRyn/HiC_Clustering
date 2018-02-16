import sys

with open (sys.argv[1],'r') as infile:
	for line in infile:
		column2 = line.split()[2]
		print(column2.split('l')[0])
