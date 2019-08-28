#!/usr/bin/env python3

import json, glob, sys, time
from optparse import OptionParser


def compareHash(input_folder, compareOut):
	fileList = glob.glob('%s/*/*.json'%input_folder)
	fileList = fileList + glob.glob('%s/*.json'%input_folder)
	w = open(compareOut, 'w')
	w.write('id1\tid2\tloci_compared\tdifferences\tdist\n')
	jsonList = []
	
	sys.stdout.write("Reading in JSON files\n")
	start = time.time()
	for f in fileList:
		with open(f, 'r') as fp:
			j = json.load(fp)
			jsonList.append(j)
	end = time.time()
	sys.stdout.write("Seconds to read in files: %s\n"%(end - start))
	
	sys.stdout.write("Comparing profiles\n")
	start = time.time()
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
	end = time.time()
	sys.stdout.write("Seconds to compare profiles: %s\n"%(end - start))


if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option( '-i', '--input_folder', action = 'store', type='string', dest = 'input_folder', default = '.' )
	parser.add_option( '-o', '--output', action = 'store', type='string', dest = 'output', default = 'output' )
	
	opts, args = parser.parse_args()
	input_folder = opts.input_folder
	compareOut = opts.output
	
	compareHash(input_folder, compareOut)
	
#E.g. compareProfiles.py -i /well/bag/deyre/analysis/spades-flow/replicates_output -o  /well/bag/deyre/analysis/spades-flow/replicates_compare.txt

# compareProfiles.py -i /home/davideyre/hash-cgmlst/comparison_study_data/replicates_output -o  /home/davideyre/hash-cgmlst/comparison_study_data/replicates_compare.txt

