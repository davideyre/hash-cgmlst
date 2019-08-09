#!/usr/bin/env nextflow

// parameters 
params.seqlist = "example_input.csv"
params.outputPath = "example_output"
params.krakendb = "minikraken2" //location of minikrakenDB

def firstFive( str ) {
    return str.substring(0,5)
}



// initial logging
log.info "\n" 
log.info "Spadeflow -- version 0.1"
log.info "Input sequence list    :  ${params.seqlist}"
log.info "Output path            :  ${params.outputPath}"
log.info "Container engine       :  ${workflow.containerEngine}"
log.info "\n"

// rename input parameters
outputPath = file(params.outputPath)
krakendb = file(params.krakendb)

//location for bbduk adapter sequences
bbduk_adapaters = "/opt/conda/opt/bbmap-38.22-1/resources/adapters.fa" //path within docker/singularity image


// set up initial channel based on CSV file
Channel
    .fromPath(params.seqlist)
    .splitCsv(header:true)
    .map{ row-> tuple(row.file_type, row.file_name) }
    .set { samples_ch }


process fetchReads {
	
	input:
		set file_type, file_name from samples_ch
	output:
		set file_name, file_type, file("*") into reads_ch
	tag "$file_name"
	
	executor = 'local'
	
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
		


}


process makeFastQ {
	input:
		set file_name, file_type, file("*") from reads_ch
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
bbduk_out_ch.into { bbduk_out_ch1; bbduk_out_ch2; }


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
	spades.py --careful --only-assembler -o spades -1 clean.1.fq.gz -2 clean.2.fq.gz \
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
		file "${file_name}_cgmlst.*"
	
	tag "$file_name"
	publishDir "${outputPath}/${firstFive(file_name)}", mode: 'copy', pattern: "${file_name}_cgmlst.*"
	
	"""
	#get stats
	/opt/conda/opt/bbmap-38.22-1/stats.sh in=${file_name}_spades_contigs.fa > ${file_name}_cgmlst.stats
	/opt/conda/opt/bbmap-38.22-1/statswrapper.sh in=${file_name}_spades_contigs.fa > ${file_name}_cgmlst.statlog
	#run hash cgmlst
	${baseDir}/bin/getCoreGenomeMLST.py -f ${file_name}_spades_contigs.fa \
		-n ${file_name} \
		-s ${baseDir}/ridom_scheme/files \
		-d ${baseDir}/ridom_scheme/ridom_scheme.fasta \
		-o ${file_name} \
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
