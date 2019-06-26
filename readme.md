# hash-cgMLST
A tool for core-genome MLST typing for bacterial data. This has been initially developed for Clostridium difficile, but could be adpated to other bacteria.

## Dependencies
* Singularity - instructions can be found here https://github.com/sylabs/singularity/blob/master/INSTALL.md
* Java version 8 or later (required for nextflow)
* Nextflow - https://www.nextflow.io

## Installation

### Singularity

Build image
```
cd singularity
sudo singularity build hash-cgmlst.img Singularity
```

If you do not have sudo access to your compute machine, the image can be built on another installation and copied across.

The image can be mounted and tested if desired:
```
singularity shell hash-cgmlst.img
```

### miniKraken2
Fetch the miniKraken2 database, from the root of the repository
```
mkdir minikraken2
cd minikraken2
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v1_8GB_201904_UPDATE.tgz
tar -xzf minikraken2_v1_8GB_201904_UPDATE.tgz
rm minikraken2_v1_8GB_201904_UPDATE.tgz
cd minikraken2_v1_8GB
mv * ../
cd ..
rmdir minikraken2_v1_8GB
```

### cgmlst.org scheme
A version of this is included - it can be updated if desired, but this is not required for hash-cgMLST.

To update this, visit cgmlst.org and download the allele files for each gene to `ridom_scheme/files` and then run:
```
cd bin
python makeReferenceDB.py
```


## Running hash-cgMLST
At present the scripts are designed to be run on two specific server set-ups, each with its own profile. Updates can be made to the nextflow.config file for other set-ups. 

The two profiles are:
* cluster - an example implementation for a SGE based cluster
* ophelia - an example for a local bare-metal server (would adapt to a stand-alone machine if the amount of memory available per core is changed)

To run the cgMLST analysis
```
nextflow hash-cgMLST.nf --seqlist example_data/example_input.csv --outputPath comparison_study_data/example_output -resume -profile ophelia
```

This will download 2 example pairs of fastq files from EBI and process these through the pipeline. Please see `example_data/example_input.csv` for an example of an input file. The file type column should be set to ebi to download from EBI and the file_name column should be a sample or run identifier. 

Using the `example_six_hospitals.csv` file as an input would allow you to run hash-cgMLST on the 973 samples used in the study describing hash-cgMLST.

## Outputs
Outputs provided include:
 - QC data: `*_base_qual.txt`, `*_kraken.txt`, `*_length.txt`, `*.raw_fastqc.html`, `*.clean_fastqc.html`
 - Spades contigs and assembly stats: `*_spades_contigs.fa`, `*_cgmlst.stats`
 - cgMLST genes as a multifasta file: `*_cgmlst.fa`
 - standard cgMLST calls as a json file: `*_cgmlst.profile` (no tracking of novel alleles is done, just recorded missing for now)
 - hash-cgMLST calls as a json file: `*_cgmlst.json`
 - standard MLST calls: `*_mlst.txt`


## Comparison scripts
To run a comparison for hash-cgMLST profiles after the nextflow pipeline above is complete use the `bin/compareProfiles.py` script:
```
bin/compareProfiles.py -i comparison_study_data/example_output -o  comparison_study_data/example_compare.txt
```

The `bin/compareProfilesExclude.py` script ignores the 26 genes likely prone to mis-assembly.


## Limitations
### File size
Each downloaded set of gzipped fastq files will be stored in the working directory, as will a pair of gzipped cleaned fastq files. You will need to periodically delete your working directory to prevent it from growing too large for very large projects.