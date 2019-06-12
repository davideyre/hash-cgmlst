#!/usr/bin/env python3

# e.g. bin/get_qc_summary.py -i comparison_study_data/replicates_output -o  comparison_study_data/replicates_
import json, glob, sys
from optparse import OptionParser

def read_summary(inpath, outpath):
	fileList = glob.glob('%s/*/*_length.txt'%inpath)
	outfile = '%sread_length_summary.txt'%outpath
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
	
def read_qc(inpath, outpath):
	fileList = glob.glob('%s/*/*_base_qual.txt'%inpath)
	outfile = '%sreads_summary.txt'%outpath
	w = open(outfile, 'w')
	w.write('id\tcoverage\tperc_q30\n')
	
	for file in fileList:
		id = file.split('/')[-1].split('_')[0][0:8]
		f = open(file, "r")
		header = next(f) #Quality	count1	fraction1	count2	fraction2
		total_bases = 0
		high_quality_bases = 0
		for l in f:
			l = l.strip().split()
			total_bases += int(l[1]) + int(l[3])
			if int(l[0])>=30:
				high_quality_bases += int(l[1]) + int(l[3])
		summary = '%s\t%0.1f\t%0.3f\n'%(id, total_bases/4290252, high_quality_bases/total_bases)
		w.write(summary)

def kraken_summary(inpath, outpath):
	fileList = glob.glob('%s/*/*_kraken.txt'%inpath)
	outfile = '%skraken_summary.txt'%outpath
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
	parser.add_option( '-o', '--outpath', action = 'store', type='string', dest = 'outpath', default = '.' )
	
	opts, args = parser.parse_args()
	inpath = opts.inpath
	outpath = opts.outpath
	
	read_summary(inpath, outpath)
	read_qc(inpath, outpath)
	kraken_summary(inpath, outpath)