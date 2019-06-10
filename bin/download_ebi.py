#!/usr/bin/env python3

import pandas
import requests
import sys
from io import StringIO
import ftplib
import hashlib
from optparse import OptionParser

def get_md5(filename, block_size=2**20):
	md5 = hashlib.md5()
	try:
		file = open(filename, 'rb')
		while True:
			data = file.read(block_size)
			if not data:
				break
			md5.update(data)
	except IOError:
		print('File \'' + filename + '\' not found!')
		return None
	except:
		return None
	return md5.hexdigest()


def download_ebi(acc, outdir):
	acc_url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s&result=read_run'%(acc)
	
	#get details of files from EBI
	ena_request = requests.get(acc_url)
	if ena_request.status_code != 200:
		sys.stderr.write("Error: EBI http code %s for accession %s\n"%(ena_request.status_code, accession))
		sys.exit()
	else:
		rio = StringIO(ena_request.text)
		tbl = pandas.read_csv(rio, sep='\t')
	
	#only want one row at a time so error if more than one
	if tbl.shape[0]!=1:
		sys.stderr.write("Error: Wrong number of rows: %s for ENA record for %s\n"%(tbl.shape[0], accession))
		sys.exit()
	
	#get the file names and md5
	for num, row in tbl.iterrows():
		try:
			file_pair = row['fastq_ftp'].split(';')
			md5_pair = row['fastq_md5'].split(';')
		except:
			sys.stderr.write("Error: Wrong number of files ENA record for %s\n"%(accession))
			sys.exit()
	
	#download the file
	with ftplib.FTP('ftp.sra.ebi.ac.uk') as ftp:
		ftp.login()
		for f, md5 in zip(file_pair, md5_pair):
			ftp_url = f.replace('ftp.sra.ebi.ac.uk', '')
			out_file = '%s/%s'%(outdir, f.split('/')[-1]) #keep filename the same
			file_md5 = ""
			i = 0
			while (file_md5 != md5):
				ftp.retrbinary('RETR {0}'.format(ftp_url), open(out_file, 'wb').write)
				file_md5 = get_md5(out_file)
				if i>3:
					sys.stderr.write("Error: Persistent failure to match md5 for %s\n"%(accession))
					sys.exit()
			
if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option( '-a', '--acc', action = 'store', type='string', dest = 'acc', default = '.' )
	parser.add_option( '-o', '--outdir', action = 'store', type='string', dest = 'outdir', default = 'output' )
	
	opts, args = parser.parse_args()
	acc = opts.acc
	outdir = opts.outdir
	
	download_ebi(acc, outdir)
	
	#E.g. ./download_ebi.py -a ERS199799 -o ~/Desktop
