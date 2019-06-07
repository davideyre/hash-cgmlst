#!/usr/bin/env python3

import json, glob, sys

fileList = glob.glob('/well/bag/deyre/analysis/spades-flow/replicates_output/*/*_length.txt')

outfile = '/well/bag/deyre/analysis/spades-flow/replicates_read_length_summary.txt'
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