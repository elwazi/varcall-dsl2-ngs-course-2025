#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Get memory safely
def getMem(task_mem) {
    // Calculate memory safely - ensure at least 1GB for Xmx
    def total_mem = task_mem.toGiga().doubleValue()
    def reserve_mem = Math.min(4.0, Math.max(1.0, total_mem / 4.0))
    def mem = Math.max(1.0, total_mem - reserve_mem)
    return mem
}

process run_combine_gvcfs {
    tag { "${params.project_name}.${chr}.rCG" }
    label 'gatk'
    label 'extra_large_memory'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.g.vcf*"

    input:
    path(gvcf_list)
    each chr

    output:
    tuple path("${cohort}.${chr}.g.vcf.gz"), path("${cohort}.${chr}.g.vcf.gz.tbi"), emit: cohort_chr_calls

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def mem = getMem(task.memory)

    cohort = params.cohort_id ?: params.project_name
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms1g -Xmx${mem.intValue()}g" CombineGVCFs \
        -R ${params.reference_files.ref} \
        --intervals ${chr} \
        --variant ${gvcf_list} \
        -O ${cohort}.${chr}.g.vcf.gz
    """
}

process run_concat_gvcfs {
    label 'gatk'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.{g.vcf.gz,g.vcf.gz.tbi}"

    input:
    tuple path(gvcfs), path(indices)
    
    output:
    path "cohort.g.vcf.gz", emit: combined_gvcf
    path "cohort.g.vcf.gz.tbi"
    val true, emit: concat_done

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def mem = task.memory.toGiga() < 1 ? 1 : task.memory.toGiga().intValue()
    def vcf_list_file = 'gvcf_files.list'

    """
    # Create list of VCF files to combine
    for gvcf in ${gvcfs}; do
        echo \$gvcf >> ${vcf_list_file}
    done

    # Combine the VCFs 
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms1g -Xmx${mem}g" \\
        MergeVcfs \\
        --INPUT ${vcf_list_file} \\
        --OUTPUT cohort.g.vcf.gz
    """
}
