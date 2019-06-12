#!/usr/bin/env python3

# e.g. bin/getReadLengthSummary.py -i comparison_study_data/replicates_output -o  comparison_study_data/replicates_read_length_summary.txt
import json, glob, sys
from optparse import OptionParser

def read_summary(inpath, outfile):
	fileList = glob.glob('%s/*/*_length.txt'%inpath)
	o = open(outfile, 'w')
	o.write('id\tmax_length\n')
	
	for file in fileList:
		id = file.split('/')[-1].split('_')[0][0:8]
		f = open(file, 'r')
		for l in f:
			l = l.strip().split()
		max_length = l[0] #have iterated over all lines the maximum length is the last
		out = '%s\t%s\n'%(id, max_length)
		o.write(out)
	
	o.close()

if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option( '-i', '--inpath', action = 'store', type='string', dest = 'inpath', default = '.' )
	parser.add_option( '-o', '--outfile', action = 'store', type='string', dest = 'outfile', default = 'read_length_summary.txt' )
	
	opts, args = parser.parse_args()
	inpath = opts.inpath
	outfile = opts.outfile
	
	read_summary(inpath, outfile)