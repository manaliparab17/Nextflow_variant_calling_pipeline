#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reads       = "Sequences/*_{1,2}.fastq.gz"
params.reference   = "Reference/chr11.fa"
params.snpeff_db   = "GRCm39.115"
params.outdir      = "results"
params.threads     = 4
params.qual        = 30
params.depth       = 10

workflow {

    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    reference_ch  = Channel.fromPath(params.reference, checkIfExists: true)

    INDEX_REFERENCE(reference_ch)

    FASTQC(read_pairs_ch)

    ALIGN(read_pairs_ch, INDEX_REFERENCE.out.fasta, INDEX_REFERENCE.out.bwa_index)

    SORT_BAM(ALIGN.out.sam)

    INDEX_BAM(SORT_BAM.out.sorted_bam)

    FLAGSTAT(INDEX_BAM.out.bam_bai)

    CALL_VARIANTS(INDEX_BAM.out.bam_bai, INDEX_REFERENCE.out.fasta, INDEX_REFERENCE.out.fai)

    FILTER_VARIANTS(CALL_VARIANTS.out.vcf)

    FILTER_SNPS(FILTER_VARIANTS.out.filtered_vcf)
    FILTER_INDELS(FILTER_VARIANTS.out.filtered_vcf)

    ANNOTATE(FILTER_VARIANTS.out.filtered_vcf)

    EXTRACT_HIGH_IMPACT(ANNOTATE.out.annotated_vcf)
}

process INDEX_REFERENCE {
    tag "bwa index + faidx"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path ref

    output:
    path ref,                           emit: fasta
    path "${ref}.fai",                  emit: fai
    path "${ref}.{amb,ann,bwt,pac,sa}", emit: bwa_index

    script:
    """
    bwa index ${ref}
    samtools faidx ${ref}
    """
}

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/QC", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}"

    script:
    """
    fastqc ${reads[0]} ${reads[1]} -o .
    """
}

process ALIGN {
    tag "$sample_id"
    cpus params.threads
    maxForks 1

    input:
    tuple val(sample_id), path(reads)
    path ref
    path bwa_index

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam

    script:
    """
    bwa mem -t ${task.cpus} ${ref} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}

process SORT_BAM {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: sorted_bam

    script:
    """
    samtools view -Sb ${sam} > ${sample_id}.bam
    samtools sort ${sample_id}.bam -o ${sample_id}.sorted.bam
    """
}

process INDEX_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/Alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path(sorted_bam), path("${sorted_bam}.bai"), emit: bam_bai

    script:
    """
    samtools index ${sorted_bam}
    """
}

process FLAGSTAT {
    tag "$sample_id"
    publishDir "${params.outdir}/Alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(bai)

    output:
    path "${sample_id}.alignment_stats.txt"

    script:
    """
    samtools flagstat ${sorted_bam} > ${sample_id}.alignment_stats.txt
    """
}

process CALL_VARIANTS {
    tag "$sample_id"
    publishDir "${params.outdir}/Variants", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(bai)
    path ref
    path fai

    output:
    tuple val(sample_id), path("${sample_id}.variants.vcf.gz"), emit: vcf

    script:
    """
    bcftools mpileup -f ${ref} ${sorted_bam} -Ou -o ${sample_id}.output.bcf
    bcftools call -mv -Oz -o ${sample_id}.variants.vcf.gz ${sample_id}.output.bcf
    """
}

process FILTER_VARIANTS {
    tag "$sample_id"
    publishDir "${params.outdir}/Variants", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"), emit: filtered_vcf

    script:
    """
    bcftools filter -i 'QUAL>${params.qual} && DP>${params.depth}' ${vcf} -Oz -o ${sample_id}.filtered.vcf.gz
    """
}

process FILTER_SNPS {
    tag "$sample_id"
    publishDir "${params.outdir}/Variants", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf)

    output:
    path "${sample_id}.snps.vcf.gz"

    script:
    """
    bcftools view -v snps ${filtered_vcf} -Oz -o ${sample_id}.snps.vcf.gz
    """
}

process FILTER_INDELS {
    tag "$sample_id"
    publishDir "${params.outdir}/Variants", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf)

    output:
    path "${sample_id}.indels.vcf.gz"

    script:
    """
    bcftools view -v indels ${filtered_vcf} -Oz -o ${sample_id}.indels.vcf.gz
    """
}

process ANNOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/Variants", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf)

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf"), emit: annotated_vcf
    path "${sample_id}.snpEff_summary.html"
    path "${sample_id}.snpEff_summary.genes.txt"

    script:
    """
    snpEff -Xmx6g -stats ${sample_id}.snpEff_summary.html ${params.snpeff_db} ${filtered_vcf} > ${sample_id}.annotated.vcf
    """
}

process EXTRACT_HIGH_IMPACT {
    tag "$sample_id"
    publishDir "${params.outdir}/Variants", mode: 'copy'

    input:
    tuple val(sample_id), path(annotated_vcf)

    output:
    path "${sample_id}.high_impact.vcf"

    script:
    """
    grep "HIGH" ${annotated_vcf} > ${sample_id}.high_impact.vcf || touch ${sample_id}.high_impact.vcf
    """
}
