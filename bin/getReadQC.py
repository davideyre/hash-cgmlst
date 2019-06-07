#!/usr/bin/env python3

import glob, sys

fileList = glob.glob('/well/bag/deyre/analysis/spades-flow/replicates_output/*/*_base_qual.txt')

outfile = '/well/bag/deyre/analysis/spades-flow/replicates_reads_summary.txt'
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
        
    