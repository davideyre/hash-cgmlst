#!/usr/bin/env python3

import json, glob, sys

fileList = glob.glob('/well/bag/deyre/analysis/spades-flow/replicates_output/*/*.json')
compareOut = '/well/bag/deyre/analysis/spades-flow/replicates_compare.txt'
w = open(compareOut, 'w')
w.write('id1\tid2\tloci_compared\tdifferences\tdist\n')
jsonList = []

sys.stdout.write("Reading in JSON files\n")
for f in fileList:
	with open(f, 'r') as fp:
		j = json.load(fp)
		jsonList.append(j)

for i in range(0, len(jsonList)):
	sys.stdout.write("%s\n"%i)
	a1 = [jsonList[i]['alleles'][k] for k in jsonList[i]['alleles'].keys()]
	for j in range(0, len(jsonList)):
		if i<j:
			a2 = [jsonList[j]['alleles'][k] for k in jsonList[i]['alleles'].keys()]
			compared = 0
			diff = 0
			compared = len([True for l1,l2 in zip(a1,a2) if l1 and l2])
			diff = len([True for l1,l2 in zip(a1,a2) if l1 and l2 and l1!=l2])
			
			if compared>0:
				ratio = "%0.3f"%(diff/compared)
			else:
				ratio = ""			
			w.write('%s\t%s\t%s\t%s\t%s\n'%(jsonList[i]['name'], jsonList[j]['name'], compared, diff, ratio))

w.close()

#total 2270

#ideas
# snps vs cgmlst dist 
# performance in replicates
# performance in reference re-seqeuncing
# ix site of extra / false variation (in all 3 above) - are some genes unreliable

# need to re-assemble with Spades - need to assemble collections of replicates and re-sequenced references