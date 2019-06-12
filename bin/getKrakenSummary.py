#!/usr/bin/env python3

# e.g. bin/getKrakenSummary.py -i comparison_study_data/replicates_output -o  comparison_study_data/replicates_kraken_summary.txt
import json, glob, sys
from optparse import OptionParser

def kraken_summary(inpath, outfile):
	fileList = glob.glob('%s/*/*_kraken.txt'%inpath)
	o = open(outfile, 'w')
	o.write('id\ttotal\tunclassified\tcdiff\tother\ttop_species\ttop_species2\ttop_species3\n')
	
	for file in fileList:
		id = file.split('/')[-1].split('_')[0][0:8]
		f = open(file, 'r')
		total = 0
		cdiff = 0
		unclassified = 0
		top_species = ""
		top_species2 = ""
		top_species3 = ""
		for l in f:
			l = l.strip().split('\t')
			taxa = l[5].strip()
			if taxa == "unclassified":
				unclassified = int(l[1])
			if taxa in ("unclassified", "root"):
				total += int(l[1])
			if taxa.startswith("Clostridioides difficile"):
				cdiff += int(l[1])
			if not top_species and l[3]=="S":
				top_species = taxa
			if (not top_species2) and l[3]=="S" and taxa!=top_species:
				top_species2 = taxa
			if (not top_species3) and l[3]=="S" and taxa!=top_species and taxa!=top_species2:
				top_species3 = taxa
		out = '%s\t%s\t%0.3f\t%0.3f\t%0.3f\t%s\t%s\t%s\n'%(id, total, unclassified/total, cdiff/total, 
		(total-unclassified-cdiff)/total, top_species, top_species2, top_species3)
		o.write(out)
	o.close()

if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option( '-i', '--inpath', action = 'store', type='string', dest = 'inpath', default = '.' )
	parser.add_option( '-o', '--outfile', action = 'store', type='string', dest = 'outfile', default = 'kraken_summary.txt' )
	
	opts, args = parser.parse_args()
	inpath = opts.inpath
	outfile = opts.outfile
	
	kraken_summary(inpath, outfile)