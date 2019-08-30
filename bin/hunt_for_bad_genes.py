#!/usr/bin/env python3

## find genes differing in replicate pairs
##  - needs a summary of results uploading first: results_filtered.csv


import json
import os
import pandas as pd
from Bio import SeqIO

def geneDiff(id1, id2, sample, path_json1, path_json2, path_fa1, path_fa2):
	o = []
	
	#read in json
	with open(path_json1, 'r') as fp:
		j1 = json.load(fp)
	with open(path_json2, 'r') as fp:
		j2 = json.load(fp)
	
	#read in fasta
		fa1, fa2 = dict(), dict()
	for s in SeqIO.parse(path_fa1, 'fasta'):
		fa1[s.id] = s.seq._data
	for s in SeqIO.parse(path_fa2, 'fasta'):
		fa2[s.id] = s.seq._data
	
	#find non-matching genes
	for k in j1['alleles'].keys():
		if j1['alleles'][k] != j2['alleles'][k] and j1['alleles'][k] and j2['alleles'][k]:
			pos = ", ".join([str(i) for i, (b1, b2) in enumerate(zip(fa1[k], fa2[k])) if b1!=b2])
			if sample:
				o.append("\t".join([id1, id2, sample, k, str(len(fa1[k])), pos]))
			else:
				o.append("\t".join([id1, id2, k, str(len(fa1[k])), pos]))
	return o


def getDiff(infile, search_dir, outfile, assembler):
	df = pd.read_csv(infile)
	w = open(outfile, 'w')
	w.write('id1\tid2\tsample\tgene\tgene_length\tpositions\n')
	for index, row in df.iterrows():
		id1, id2, sample = row['id1'], row['id2'], row['samplename']
		file1_stem = [f.split('_')[0] for f in os.listdir(search_dir+'/'+id1[0:5]) if f[0:8]==id1 and f.split('_')[1]==assembler and f.split('_')[2]=='cgmlst.json'][0]
		file2_stem = [f.split('_')[0] for f in os.listdir(search_dir+'/'+id2[0:5]) if f[0:8]==id2 and f.split('_')[1]==assembler and f.split('_')[2]=='cgmlst.json'][0]
		path_json1 = search_dir+'/'+id1[0:5]+'/'+file1_stem+'_'+assembler+'_cgmlst.json'
		path_json2 = search_dir+'/'+id2[0:5]+'/'+file2_stem+'_'+assembler+'_cgmlst.json'
		path_fa1 = search_dir+'/'+id1[0:5]+'/'+file1_stem+'_'+assembler+'_cgmlst.fa'
		path_fa2 = search_dir+'/'+id2[0:5]+'/'+file2_stem+'_'+assembler+'_cgmlst.fa'
		print (path_json1)
		print (path_json2)
		print (path_fa1)
		print (path_fa2)
		o = geneDiff(id1, id2, sample, path_json1, path_json2, path_fa1, path_fa2)
		if o:
			w.write("\n".join(o)+"\n")
	w.close()

infile = '/home/davideyre/hash-cgmlst/comparison_study_data/replicates_output/results_filtered.csv'
search_dir = '/home/davideyre/hash-cgmlst/comparison_study_data/replicates_output'
outfile = '/home/davideyre/hash-cgmlst/comparison_study_data/replicates_output/genes_with_differences_spades.txt'
getDiff(infile, search_dir, outfile, 'spades')

outfile = '/home/davideyre/hash-cgmlst/comparison_study_data/replicates_output/genes_with_differences_skesa.txt'
getDiff(infile, search_dir, outfile, 'skesa')

# 
# def getDiffSix(infile, search_dir, outfile):
# 	df = pd.read_csv(infile)
# 	w = open(outfile, 'w')
# 	w.write('id1\tid2\tgene\tgene_length\tpositions\n')
# 	for index, row in df.iterrows():
# 		sample = ""
# 		id1, id2 = row['run.1'], row['run.2']
# 		file1_stem = [f.split('_')[0] for f in os.listdir(search_dir+'/'+id1[0:5]) if f[0:8]==id1][0]
# 		file2_stem = [f.split('_')[0] for f in os.listdir(search_dir+'/'+id2[0:5]) if f[0:8]==id2][0]
# 		path_json1 = search_dir+'/'+id1[0:5]+'/'+file1_stem+'_spades_cgmlst.json'
# 		path_json2 = search_dir+'/'+id2[0:5]+'/'+file2_stem+'_spades_cgmlst.json'
# 		path_fa1 = search_dir+'/'+id1[0:5]+'/'+file1_stem+'_spades_contigs.fa'
# 		path_fa2 = search_dir+'/'+id2[0:5]+'/'+file2_stem+'_spades_contigs.fa'
# 		o = geneDiff(id1, id2, sample, path_json1, path_json2, path_fa1, path_fa2)
# 		if o:
# 			w.write("\n".join(o)+"\n")
# 	w.close()
# 
# infile = '/home/davideyre/hash-cgmlst/comparison_study_data/six_hospitals/within_2_snps.csv'
# search_dir = '/home/davideyre/hash-cgmlst/comparison_study_data/six_hospitals'
# outfile = '/home/davideyre/hash-cgmlst/comparison_study_data/six_hospitals/genes_with_differences.txt'
# getDiffSix(infile, search_dir, outfile)

