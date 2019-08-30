#!/usr/bin/env nextflow

// parameters 
params.seqlist = "example_input.csv"
params.outputPath = "example_output"
params.krakendb = "minikraken2" //location of minikrakenDB
params.runSkesa = 0
params.aoSpades = 1

def firstFive( str ) {
    return str.substring(0,5)
}



// initial logging
log.info "\n" 
log.info "hash-cgMLST -- version 0.1"
log.info "Input sequence list    :  ${params.seqlist}"
log.info "Output path            :  ${params.outputPath}"
log.info "Container engine       :  ${workflow.containerEngine}"
log.info "\n"

// rename input parameters
outputPath = file(params.outputPath)
krakendb = file(params.krakendb)
runSkesa = params.runSkesa
aoSpades = params.aoSpades

spadesReadCorrection = ""
if (aoSpades==1) {
	spadesReadCorrection = "--only-assembler"
}

//location for bbduk adapter sequences
bbduk_adapaters = "/opt/conda/opt/bbmap-38.22-1/resources/adapters.fa" //path within docker/singularity image


// set up initial channel based on CSV file
Channel
    .fromPath(params.seqlist)
    .splitCsv(header:true)
    .map{ row-> tuple(row.file_type, row.file_name, row.fq1, row.fq2) }
    .set { samples_ch }


process fetchReads {
	
	input:
		set file_type, file_name, fq1, fq2 from samples_ch
	output:
		set file_type, file_name, file("*") into reads_ch

	
	executor = 'local'
	tag "$file_name"

	
	script:
		if (file_type=="bam") { 
			"""
			scp -P 8081 ana:/mnt/microbio/ndm-hicf/ogre/pipeline_output/${file_name}/MAPPING/103e39d6-096c-46da-994d-91c5acbda565_R00000003/STD/${file_name}_v3.bam in.bam || scp -P 8081 ana:/mnt/microbio/ndm-hicf/ogre/pipeline_output/${file_name}/MAPPING/103e39d6-096c-46da-994d-91c5acbda565_R00000003/STD/${file_name}_v2.bam in.bam
			"""
		}
		else if (file_type=="ebi") {
			"""
			${baseDir}/bin/download_ebi.py -a ${file_name} -o .
			mv *_1.fastq.gz in.1.fq.gz
			mv *_2.fastq.gz in.2.fq.gz
			"""
		}
		else if (file_type=="local") {
			f1 = file(fq1)
			f2 = file(fq2)
			
			"""
			ln -s $f1 in.1.fq.gz
			ln -s $f2 in.2.fq.gz
			"""
			
		}
		
}


process makeFastQ {
	input:
		set file_type, file_name, file("*") from reads_ch
	output:
		set file_name, file("in.1.fq.gz"), file("in.2.fq.gz") into fq_ch
	
	tag "$file_name"
	
	script:
		if (file_type=="bam") {
		"""
		samtools sort -@${task.cpus} -n -o sorted.bam in.bam
		bedtools bamtofastq -i sorted.bam \
						  -fq in.1.fq \
						  -fq2 in.2.fq
		rm sorted.bam
		gzip in.1.fq
		gzip in.2.fq
		"""
		}
		else {
		"""
		
		"""
		}
}

//split raw reads into 3 channels - for QC and assembly and kraken2
fq_ch.into { fq_ch1; fq_ch2; fq_ch3}

process kraken2 {

	input:
    	set file_name, file("in.1.fq.gz"), file("in.2.fq.gz") from fq_ch3
	
	output:
		file "*"
	
	tag "$file_name"
    
    publishDir "${outputPath}/${firstFive(file_name)}", mode: 'copy', pattern: "${file_name}*"

	"""
	kraken2 --report ${file_name}_kraken.txt --db ${krakendb} --paired in.1.fq.gz in.2.fq.gz > /dev/null
	"""

}


process rawFastQC {
	
	input:
    	set file_name, file("in.1.fq.gz"), file("in.2.fq.gz") from fq_ch1
	
	output:
		file "*"
	
	tag "$file_name"
    
	publishDir "${outputPath}/${firstFive(file_name)}", mode: 'copy', pattern: "${file_name}*"
	
	"""
	cat in.1.fq.gz in.2.fq.gz > ${file_name}.raw.fq.gz
	fastqc --threads ${task.cpus} ${file_name}.raw.fq.gz > ${file_name}_raw_fastqc_log.txt
	rm ${file_name}.raw.fq.gz
	"""

}


//adapter trimming with bbDuk
process bbDuk {
	
	input:
		set file_name, file("in.1.fq.gz"), file("in.2.fq.gz") from fq_ch2
	
	output:
		set file_name, file("clean.1.fq.gz"), file("clean.2.fq.gz") into bbduk_out_ch
		file("${file_name}_base_qual.txt")
		file("${file_name}_length.txt")
	
	tag "$file_name"
	memory '8 GB' //not able to use 16GB standard
    
    publishDir "${outputPath}/${firstFive(file_name)}", mode: 'copy', pattern: "${file_name}*"
	
	"""
	bbduk.sh in1=in.1.fq.gz in2=in.2.fq.gz out1=clean.1.fq out2=clean.2.fq \
				ref=$bbduk_adapaters ktrim=r k=23 mink=11 hdist=1 tpe tbo \
				qtrim=rl trimq=30 \
				qchist=${file_name}_base_qual.txt \
				lhist=${file_name}_length.txt \
				-Xmx${task.memory.toGiga()}g threads=${task.cpus}
	gzip clean.1.fq
	gzip clean.2.fq
	"""
}

//split cleaned reads into 2 channels - for QC and assembly
bbduk_out_ch.into { bbduk_out_ch1; bbduk_out_ch2; bbduk_out_ch3;}


process cleanFastQC {
	
	input:
    	set file_name, file("clean.1.fq.gz"), file("clean.2.fq.gz") from bbduk_out_ch1
	
	output:
		file "*"
	
	tag "$file_name"
    
	publishDir "${outputPath}/${firstFive(file_name)}", mode: 'copy', pattern: "${file_name}*"
	
	"""
	cat clean.1.fq.gz clean.2.fq.gz > ${file_name}.clean.fq.gz
	fastqc --threads ${task.cpus} ${file_name}.clean.fq.gz > ${file_name}_clean_fastqc_log.txt
	rm ${file_name}.clean.fq.gz
	"""

}


process spades {

	input:
		set file_name, file("clean.1.fq.gz"), file("clean.2.fq.gz") from bbduk_out_ch2
	
	output:
		set file_name, file("${file_name}_spades_contigs.fa") into spades_out
			
	tag "$file_name"
	
	publishDir "${outputPath}/${firstFive(file_name)}", mode: 'copy', pattern: "${file_name}_*"
    
	"""
	spades.py --careful ${spadesReadCorrection} -o spades -1 clean.1.fq.gz -2 clean.2.fq.gz \
		-t ${task.cpus} -m ${task.memory.toGiga()}
	cp spades/contigs.fasta ${file_name}_spades_contigs.fa
	cp spades/assembly_graph.fastg ${file_name}_spades_assembly_graph.fastg
	cp spades/assembly_graph_with_scaffolds.gfa ${file_name}_spades_assembly_graph_with_scaffolds.gfa
	cp spades/spades.log ${file_name}_spades.log
	#remove spades directory to save space
	rm -rf spades*
	"""
}

//split assembly into 2 channels - for cgmlst and mslt
spades_out.into { spades_out_ch1; spades_out_ch2; }

process cgmlst {
	input:
	 	set file_name, file("${file_name}_spades_contigs.fa") from spades_out_ch1
	output:
		file "${file_name}_spades_cgmlst.*"
	
	tag "$file_name"
	publishDir "${outputPath}/${firstFive(file_name)}", mode: 'copy', pattern: "${file_name}_spades_cgmlst.*"
	
	"""
	#get stats
	/opt/conda/opt/bbmap-38.22-1/stats.sh in=${file_name}_spades_contigs.fa > ${file_name}_spades_cgmlst.stats
	/opt/conda/opt/bbmap-38.22-1/statswrapper.sh in=${file_name}_spades_contigs.fa > ${file_name}_spades_cgmlst.statlog
	#run hash cgmlst 
	${baseDir}/bin/getCoreGenomeMLST.py -f ${file_name}_spades_contigs.fa \
		-n ${file_name}_spades \
		-s ${baseDir}/ridom_scheme/files \
		-d ${baseDir}/ridom_scheme/ridom_scheme.fasta \
		-o ${file_name}_spades \
		-b blastn 
	"""	
}

process mlst {
	input:
	 	set file_name, file("${file_name}_spades_contigs.fa") from spades_out_ch2
	output:
		file "${file_name}_mlst.txt"
	
	tag "$file_name"
	publishDir "${outputPath}/${firstFive(file_name)}", mode: 'copy', pattern: "${file_name}_mlst.txt"
	
	"""
	mlst --scheme cdifficile --legacy  ${file_name}_spades_contigs.fa > ${file_name}_mlst.txt
	"""	
}

if (runSkesa==1) {

	//run skesa and generate cgmlst from that too
	process skesa {
	input:
		set file_name, file("clean.1.fq.gz"), file("clean.2.fq.gz") from bbduk_out_ch3
	
	output:
		set file_name, file("${file_name}_skesa_contigs.fa") into skesa_out
			
	tag "$file_name"
	
	publishDir "${outputPath}/${firstFive(file_name)}", mode: 'copy', pattern: "${file_name}_*"
	
	"""
	skesa --fastq clean.1.fq.gz,clean.2.fq.gz \
	      --cores ${task.cpus} --memory ${task.memory.toGiga()} > ${file_name}_skesa_contigs.fa
	"""
	
	}
	
	process cgmlst_skesa {
	input:
		set file_name, file("${file_name}_skesa_contigs.fa") from skesa_out
	output:
		file "${file_name}_skesa_cgmlst.*"
	
	tag "$file_name"
	publishDir "${outputPath}/${firstFive(file_name)}", mode: 'copy', pattern: "${file_name}_skesa_cgmlst.*"
	
	"""
	#get stats
	/opt/conda/opt/bbmap-38.22-1/stats.sh in=${file_name}_skesa_contigs.fa > ${file_name}_skesa_cgmlst.stats
	/opt/conda/opt/bbmap-38.22-1/statswrapper.sh in=${file_name}_skesa_contigs.fa > ${file_name}_skesa_cgmlst.statlog
	#run hash cgmlst 
	${baseDir}/bin/getCoreGenomeMLST.py -f ${file_name}_skesa_contigs.fa \
		-n ${file_name}_skesa \
		-s ${baseDir}/ridom_scheme/files \
		-d ${baseDir}/ridom_scheme/ridom_scheme.fasta \
		-o ${file_name}_skesa \
		-b blastn 
	"""	
	}
}