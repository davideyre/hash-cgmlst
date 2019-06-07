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
cd minikraken2
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v1_8GB_201904_UPDATE.tgz
tar -xzf minikraken2_v1_8GB_201904_UPDATE.tgz
cd minikraken2_v1_8GB
cp * ../
rmdir minikraken2_v1_8GB
```

### cgmlst.org scheme
A version of this is included - it can be updated if desired. [Instructions to follow.]


## Running hash-cgMLST
At present the scripts are designed to be run on a specific server set-up. More generalised versions of the scripts will follow.

To run the cgMLST analysis