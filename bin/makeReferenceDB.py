from Bio import SeqIO
import os

schemePath = '../ridom_scheme/files'
lociList = os.listdir(schemePath)
lociList.sort()

lociSeq = []

for locus in lociList:
	locus_file = "%s/%s"%( schemePath, locus )
	seq = [s for s in SeqIO.parse(locus_file, 'fasta')][0]
	seq.id = locus[:-6] #drop the .fasta from end of locus name
	lociSeq.append(seq)

dbFile = '../ridom_scheme/reference_db.fasta'
SeqIO.write(lociSeq, dbFile, 'fasta')