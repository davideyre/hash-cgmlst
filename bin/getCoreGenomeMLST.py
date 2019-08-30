#!/usr/bin/env python3

#create fasta file containing all cgMLST sites

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import sys, os, hashlib, json
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser

def hashAllele(str):
	m = hashlib.md5()
	m.update(str.encode('utf-8'))
	return (m.hexdigest())

class CgMLST:
	
	"""cgMLST representations for a single genome
	
	Attributes:
		 contigsPath: the file path to the contigs file
		 contigsName: name to use for contigs
		 schemeDB: the multifasta File containing one example of each of the loci to search for
		 schemePath: path to files from scheme, one per locus, with known alleles in
		 verbose: boolean, print more inforamtion to screen
		 matches: dictionary of the alleles for each cgmlst locus
	"""
	
	def __init__(self, contigsPath, contigsName, schemeDB, schemePath, verbose=False):
		self.contigsPath = contigsPath
		self.contigsName = contigsName
		self.schemeDB = schemeDB
		self.schemePath = schemePath
		self.verbose = verbose
		self.checkBlastDB()
		## set up dictionary for blast matches at each locus
		self.matches = { }
		
	def checkBlastDB(self):
		if not os.path.exists( "%s.nhr"%self.contigsPath ):
			#attempt to make blastDB
			cmd = "makeblastdb -dbtype nucl -in %s"%self.contigsPath
			try:
				os.system(cmd)
			except:
				raise RuntimeError('Missing blast database for: %s'%self.contigsPath)
		
	def runBlast(self, blastn_path):
		## Set up the query and execute it - use settings from Mellman paper, but set evalue to 0.01 and 
		# only return one target sequence
		cline = NcbiblastnCommandline( cmd = blastn_path,
									   query=self.schemeDB, db=self.contigsPath, 
									   evalue=0.01, outfmt=5, max_target_seqs=1, perc_identity=90,
									   word_size=11, penalty=-1, reward=1, gapopen=5,gapextend=2)
		if self.verbose:
			print (cline)
		stdout, stderr = cline()
		# Make output string look like file handle to use the parser without writing to file
		fp = StringIO(stdout)
		# Parse the result
		# Looking for a single, exact match for each locus
		blast_records = NCBIXML.parse( fp )
		exact_matches = []
		for record in blast_records:
			locus = record.query[:-2] #get locus name from heading of fasta file
			self.matches[locus] = None #default to empty record
			## should just be a single hit, i.e. one alignemnt, and one hsp within that alignment
			if (len(record.alignments)==1):
				# check if the matched length is >=99% of the query length and
				#   also check that start of query matches start and end matches end to ensure avoid truncation by contig breaks
				if ((record.alignments[ 0 ].hsps[ 0 ].align_length / record.query_length) >=0.99 and
						record.alignments[0].hsps[0].query_start==1 and 
						record.alignments[0].hsps[0].query_end == record.query_length):
					# save the hit - percentage identity is ensured from blast search criteria above as >=90%
					self.matches[locus] = record.alignments[ 0 ].hsps[ 0 ].sbjct
	
	def writeFa(self, outFa):
		seqlist = []
		for k in self.matches.keys():
			seqData = self.matches[k]
			if not seqData:
				seqData = ""
			seq = SeqRecord(seq=Seq(seqData), id=k, description="")
			seqlist.append(seq)
		SeqIO.write(seqlist, outFa, 'fasta')
		
	def writeHash(self, jsonFile):
		hashDict = dict()
		for k in self.matches.keys():
			#print(k)
			seqData = self.matches[k]
			#print(seqData)
			if seqData: #if there is any sequence data...
				#check for any ambiguous characters
				ambiguous = len([True for b in seqData if b not in 'ACGT-'])
				#print("ambiguous: %s"%ambiguous)
				if ambiguous:
					seqData = None
				else:
					#check for premature stop codons
					protein = Seq("".join([b for b in seqData if b in 'ACGT'])).translate()._data
					protein_rev = Seq("".join([b for b in seqData if b in 'ACGT'])).reverse_complement().translate()._data
					stopcodons = min(len([True for aa in protein if aa=="*"]), len([True for aa in protein_rev if aa=="*"]))
					#print("stop codons (strand): %s"%stopcodons)
					if stopcodons>1:
						seqData = None
			#if no data or ambiguous characters or premature stop codons return an empty value
			if not seqData:
				hashDict[k] = ""
			else:
				seqToHash = "".join([b for b in seqData if b in 'ACGT']) #remove any "-" before hashing...
				#    required to avoid false differences for variable placement on "-" in blast alignments
				hashDict[k] = hashAllele(seqToHash)
		cgmlst = {"name": self.contigsName, "file": self.contigsPath, 
					"scheme": self.schemeDB, "alleles": hashDict}
		with open(jsonFile, "w") as f:
			json.dump(cgmlst, fp=f, indent=4)
		missing = len([k for k in hashDict.keys() if not hashDict[k]])
		return(missing)

			
	def matchProfile(self, locusKey):
		locusFa = '%s/%s.fasta'%(self.schemePath, locusKey)
		match = ""
		for s in SeqIO.parse(locusFa, 'fasta'):
			if s.seq._data == self.matches[locusKey]:
				match = s.id
				break
		return(match)
			
	def writeProfile(self, profileFile):
		profiles = dict()
		for k in self.matches.keys():
			if not self.matches[k]:
				profiles[k] = ""
			else:
				profiles[k] = self.matchProfile(k)
		cgmlst = {"name": self.contigsName, "file": self.contigsPath, 
					"scheme": self.schemeDB, "alleles": profiles}
		with open(profileFile, "w") as f:
			json.dump(cgmlst, fp=f, indent=4)


if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option( '-f', '--contig_file', action = 'store', type='string', dest = 'contigsFa', default = '.' )
	parser.add_option( '-n', '--contig_name', action = 'store', type='string', dest = 'contigName', default = '.' )
	parser.add_option( '-s', '--scheme', action = 'store', type='string', dest = 'schemePath', default = '.' )
	parser.add_option( '-d', '--db', action = 'store', type='string', dest = 'schemeDB', default = '.' )
	parser.add_option( '-o', '--output', action = 'store', type='string', dest = 'output', default = '.' )
	parser.add_option( '-b', '--blastn', action = 'store', type='string', dest = 'blastn_path', default = '.' )
	
	opts, args = parser.parse_args()
	contigsFa = opts.contigsFa
	contigName = opts.contigName
	schemePath = opts.schemePath
	schemeDB = opts.schemeDB
	output = opts.output
	blastn_path = opts.blastn_path

	#example
	# cd /users/bag/deyre/analysis/spades-flow/
	# cgmlst/getCoreGenomeMLST.py -f replicates_output/003/003381ae-7ef0-431d-b585-c4c3f02a8399_spades_contigs.fa -n 003381ae -s cgmlst/ridom_scheme/files -d cgmlst/ridom_scheme/ridom_scheme.fasta -o replicates_output/003/003381ae-7ef0-431d-b585-c4c3f02a8399 -b /apps/htseq/ncbi-blast/bin/blastn
	
	print(contigName)
	cgmlstFa = '%s_cgmlst.fa'%output
	jsonFile = '%s_cgmlst.json'%output
	profileFile = '%s_cgmlst.profile'%output
	sample = CgMLST(contigsFa, contigName, schemeDB, schemePath)
	sample.runBlast(blastn_path)
	sample.writeFa(cgmlstFa)
	missing = sample.writeHash(jsonFile)
	sample.writeProfile(profileFile)
	print(missing)
