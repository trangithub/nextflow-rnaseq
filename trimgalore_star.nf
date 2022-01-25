nextflow.enable.dsl=2

baseDir = "/home/tran/rnaseq2"
params.reads = "$baseDir/rawdata/*{1,2}.fastq.gz"
params.trimgaloreOut = "$baseDir/output/trimgalore/run4"
params.index = "$baseDir/ref/index_10/"
params.starOut = "$baseDir/output/star/run3"

println "reads: $params.reads"
println "index: $params.index"

process trimgalore {
    tag "Trimming from $pair_id"
    publishDir "$params.trimgaloreOut", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*"), emit: trimmed_reads 
        
    script:
    """
    trim_galore --paired --gzip ${pair_id}1.fastq.gz ${pair_id}2.fastq.gz
    """
}

process star {
    tag "STAR mapping for $pair_id"
    publishDir "$params.starOut", mode: 'copy'

    input:
    tuple val(pair_id), path(trimmed_reads)
    path index
       
    output:
    tuple val(pair_id), path("${pair_id}_Aligned.sortedByCoord.out.bam"), path("${pair_id}_Aligned.sortedByCoord.out.bam.bai") 
    
    script:
    """
    STAR --readFilesIn ${pair_id}1_val_1.fq.gz ${pair_id}2_val_2.fq.gz \
    --genomeDir ${params.index} --outFileNamePrefix ${pair_id}_ \
    --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat

    samtools index ${pair_id}_Aligned.sortedByCoord.out.bam
    """
}

workflow {
    Channel
    .fromFilePairs( params.reads, checkIfExists: true )
    .set { read_pairs_ch } 

    trimgalore_ch = trimgalore(read_pairs_ch)
    star_ch = star(trimgalore.out, params.index)	
}