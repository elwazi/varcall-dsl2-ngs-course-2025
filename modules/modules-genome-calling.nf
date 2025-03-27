#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
// ref                       = file(params.ref, type: 'file')
// known_indels_1            = file(params.known_indels_1, type: 'file')
// known_indels_2            = file(params.known_indels_2, type: 'file')
// dbsnp                     = file(params.dbsnp, type: 'file')
// outdir                    = file(params.outdir, type: 'dir')

// build                     = params.build
// sample_coverage           = params.sample_coverage
// project_name              = params.project_name
// cohort_id                 = params.cohort_id

// db                        = file(params.db_path, type: 'dir')

// hapmap                    = file(params.hapmap, type: 'file')
// omni                      = file(params.omni, type: 'file')
// phase1_snps               = file(params.phase1_snps, type: 'file')
// golden_indels             = file(params.golden_indels, type: 'file')

// ts_filter_level_snps      = params.ts_filter_level_snps
// ts_filter_level_indels    = params.ts_filter_level_indels
// max_gaussians_snps        = params.max_gaussians_snps
// max_gaussians_indels      = params.max_gaussians_indels

// outdir.mkdir()

// CALL_CONF
// if ( sample_coverage == "high" ) {
//     call_conf = 30
// } else if ( sample_coverage == "low" ) {
//     call_conf = 10
// } else {
//     call_conf = 30
// }

process run_genotype_gvcf_on_genome_db {
    tag { "${params.project_name}.${interval}.rGGoG" }
    label 'gatk'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/vcfs/chroms", mode: 'copy', overwrite: true

    input:
    tuple val(interval), val(db_path)

    output:
    tuple val(interval), path("${params.cohort_id}.${interval}.vcf.gz"), path("${params.cohort_id}.${interval}.vcf.gz.tbi"), emit: db_gg_vcf_set

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx, but not more than 75% of available memory
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    def min_mem = Math.min(1, max_mem)
    
    """
    echo "Using reference: ${params.reference_files.ref}"
    ls -la ${params.reference_files.ref}
    
    # Check if the database directory exists
    if [ ! -d "${db_path}" ]; then
        echo "Error: GenomicsDB directory not found: ${db_path}"
        exit 1
    fi
    
    echo "Running genotyping on interval: ${interval}"
    echo "Using database: ${db_path}"
    
    gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \\
        --reference ${params.reference_files.ref} \\
        --intervals ${interval} \\
        --variant gendb://${db_path} \\
        -stand-call-conf 30 \\
        --annotation Coverage \\
        --annotation FisherStrand \\
        --annotation StrandOddsRatio \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation RMSMappingQuality \\
        --annotation ReadPosRankSumTest \\
        --allow-old-rms-mapping-quality-annotation-data \\
        --output "${params.cohort_id}.${interval}.vcf.gz"
    """
}

// process run_genotype_gvcf_on_genome_gvcf {
//     tag { "${project_name}.${cohort_id}.${chr}.rGGoG" }
//     label 'gatk'
//     memory { 16.GB * task.attempt }
//     publishDir "${outdir}/${params.workflow}/${project_name}", mode: 'copy', overwrite: true
    
//     input:
//     tuple path(gvcf), path(index)
//     each interval
  
//     output:
//     tuple path("${project_name}.${chr}.vcf.gz"), path("${project_name}.${chr}.vcf.gz.tbi"), emit: gg_vcf_set
  
//     script:
//     mem = task.memory.toGiga() - 4
//     if (db_import == "no") {
//         variants = "-V ${gvcf}"
//     } else if (db_import == "yes") {
//         variants = "-V gendb://${db}/${chr}.gdb"
//     } else {
//         exit 1
//     }

//     """
//     gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" GenotypeGVCFs \
//         --reference ${ref} \
//         --intervals ${chr} \
//         ${variants} \
//         -stand-call-conf ${call_conf} \
//         --annotation Coverage \
//         --annotation FisherStrand \
//         --annotation StrandOddsRatio \
//         --annotation MappingQualityRankSumTest \
//         --annotation QualByDepth \
//         --annotation RMSMappingQuality \
//         --annotation ReadPosRankSumTest \
//         --allow-old-rms-mapping-quality-annotation-data \
//         -output "${project_name}.${chr}.vcf.gz"
//     """
// }

process run_concat_vcf {
    tag { "${params.project_name}.rGCC" }
    label 'gatk'
    label 'extra_large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/vcfs", mode: 'copy',
        pattern: "*.{vcf.gz,vcf.gz.tbi}"

    input:
    tuple val(interval), path(vcfs), path(indices)
    
    output:
    tuple path("${params.project_name}.vcf.gz"), path("${params.project_name}.vcf.gz.tbi"), emit: combined_vcf
    
    script:
    // Calculate memory safely - ensure at least 1GB for Xmx, but not more than 75% of available memory
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    
    """
    # Display the VCF files to be concatenated
    echo "VCF files to concatenate:"
    ls -la *vcf.gz
    
    # Create list of VCF files to combine - handle single or multiple files
    if [ \$(ls -1 *vcf.gz | wc -l) -eq 1 ]; then
        echo "Only one VCF file found, copying instead of concatenating"
        cp \$(ls *vcf.gz) ${params.project_name}.vcf.gz
        cp \$(ls *vcf.gz.tbi) ${params.project_name}.vcf.gz.tbi
    else
        echo "Multiple VCF files found, concatenating"
        ls -1 *vcf.gz > vcf.list
        
        # Combine the VCFs
        gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms1g -Xmx${max_mem}g" \\
            GatherVcfs \\
            --INPUT=vcf.list \\
            --OUTPUT=${params.project_name}.vcf.gz
            
        # Index the combined VCF
        gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${max_mem}g" \\
            IndexFeatureFile \\
            --input ${params.project_name}.vcf.gz
    fi
    """
}

process run_vqsr_on_snps {
    tag { "${params.project_name}.${params.cohort_id}.rVoS" }
    label 'gatk'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/vcfs/vqsr", mode: 'copy', overwrite: true
    
    input:
    tuple path(vcf), path(vcf_index)
 
    output:
    tuple path(vcf), path(vcf_index),
        path("${params.project_name}.vcf.recal-SNP.recal"),
        path("${params.project_name}.vcf.recal-SNP.recal.idx"),
        path("${params.project_name}.vcf.recal-SNP.tranches"), emit: snps_vqsr_recal

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx, but not more than 75% of available memory
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    def min_mem = Math.min(1, max_mem)
    def max_gaussians = params.max_gaussians_snps ?: "4"

    """
    echo "Running VQSR for SNPs using reference files:"
    echo "Reference: ${params.reference_files.ref}"
    echo "HapMap: ${params.reference_files.hapmap}"
    echo "Omni: ${params.reference_files.omni}"
    echo "1000G SNPs: ${params.reference_files.phase1_snps}"
    echo "dbSNP: ${params.reference_files.dbsnp}"
    
    # Run VQSR for SNPs - exclude potentially problematic annotations for small datasets
    gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g" VariantRecalibrator \\
        --reference ${params.reference_files.ref} \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.reference_files.hapmap} \\
        --resource:omni,known=false,training=true,truth=true,prior=12.0 ${params.reference_files.omni} \\
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.reference_files.phase1_snps} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.reference_files.dbsnp} \\
        --use-annotation DP \\
        --use-annotation FS \\
        --use-annotation SOR \\
        --use-annotation MQ \\
        --use-annotation QD \\
        --ignore-filter FAIL_variant_filter \\
        --mode SNP \\
        --max-gaussians ${max_gaussians} \\
        --variant ${vcf} \\
        --output "${params.project_name}.vcf.recal-SNP.recal" \\
        --tranches-file "${params.project_name}.vcf.recal-SNP.tranches"
    """
}
 
process apply_vqsr_on_snps {
    tag { "${params.project_name}.${params.cohort_id}.aVoS" }
    label 'gatk'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/vcfs", mode: 'copy', overwrite: true
    
    input:
    tuple path(vcf), path(vcf_index),
        path(snp_recal), path(snp_recal_index), path(snp_tranches)
    
    output:
    tuple path("${params.project_name}.recal-SNP.vcf.gz"),
        path("${params.project_name}.recal-SNP.vcf.gz.tbi"), emit: snps_recalibrated

    script:
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    def min_mem = Math.min(1, max_mem)
    def filter_level = params.snp_filter_level ?: "99.0"

    """
    echo "Applying VQSR for SNPs using filter level: ${filter_level}"
    echo "Input VCF: ${vcf}"
    echo "SNP Recal: ${snp_recal}"
    echo "SNP Tranches: ${snp_tranches}"
    
    gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g" ApplyVQSR \\
        --reference ${params.reference_files.ref} \\
        --recal-file ${snp_recal} \\
        --tranches-file ${snp_tranches} \\
        --truth-sensitivity-filter-level ${filter_level} \\
        --create-output-variant-index true \\
        --mode SNP \\
        --variant ${vcf} \\
        --output ${params.project_name}.recal-SNP.vcf.gz
    """
}

// Add a new process for hard filtering of INDELs
process hard_filter_indels {
    tag { "${params.project_name}.${params.cohort_id}.hFI" }
    label 'gatk'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/vcfs", mode: 'copy', overwrite: true
    
    input:
    tuple path(vcf), path(vcf_index)
    
    output:
    tuple path("${params.project_name}.hard-filtered-indels.vcf.gz"),
        path("${params.project_name}.hard-filtered-indels.vcf.gz.tbi"), emit: hard_filtered_indels

    script:
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    def min_mem = Math.min(1, max_mem)
    
    """
    echo "Applying hard filtering for INDELs - small dataset detected"
    
    # Extract INDELs
    gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g" SelectVariants \\
        -R ${params.reference_files.ref} \\
        -V ${vcf} \\
        --select-type-to-include INDEL \\
        -O ${params.project_name}.indels.vcf.gz
    
    # Apply hard filters to INDELs
    gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g" VariantFiltration \\
        -R ${params.reference_files.ref} \\
        -V ${params.project_name}.indels.vcf.gz \\
        --filter-name "QD_filter" --filter-expression "QD < 2.0" \\
        --filter-name "FS_filter" --filter-expression "FS > 200.0" \\
        --filter-name "SOR_filter" --filter-expression "SOR > 10.0" \\
        -O ${params.project_name}.hard-filtered-indels.vcf.gz
    """
}

// Modify the run_vqsr_on_indels process to check indel count
process run_vqsr_on_indels {
    tag { "${params.project_name}.${params.cohort_id}.rVoI" }
    label 'large_mem'
    label 'gatk'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/vcfs/vqsr", mode: 'copy', overwrite: true
    errorStrategy { task.exitStatus == 2 ? 'ignore' : 'retry' }
    
    input:
    tuple path(vcf), path(vcf_index)
    
    output:
    tuple path(vcf), path(vcf_index),
        path("${params.project_name}.vcf.recal-INDEL.recal"),
        path("${params.project_name}.vcf.recal-INDEL.recal.idx"),
        path("${params.project_name}.vcf.recal-INDEL.tranches"),
        emit: indels_vqsr_recal, optional: true
    tuple path(vcf), path(vcf_index), emit: vcf_for_hard_filtering, optional: true
    
    script:
    // Calculate memory safely - ensure at least 1GB for Xmx, but not more than 75% of available memory
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    def min_mem = Math.min(1, max_mem)
    def max_gaussians = params.max_gaussians_indels ?: "4"

    """
    echo "Running VQSR for INDELs using reference files:"
    echo "Reference: ${params.reference_files.ref}"
    echo "Mills: ${params.reference_files.mills}"
    echo "dbSNP: ${params.reference_files.dbsnp}"
    
    # Run VQSR for INDELs - exclude potentially problematic annotations for small datasets
    gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator \\
        --reference ${params.reference_files.ref} \\
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ${params.reference_files.mills} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.reference_files.dbsnp} \\
        --use-annotation DP \\
        --use-annotation FS \\
        --use-annotation SOR \\
        --use-annotation QD \\
        --ignore-filter FAIL_variant_filter \\
        --mode INDEL \\
        --max-gaussians ${max_gaussians} \\
        --variant ${vcf} \\
        --output "${params.project_name}.vcf.recal-INDEL.recal" \\
        --tranches-file "${params.project_name}.vcf.recal-INDEL.tranches"
    """
    
    // Return tuple for hard filtering if VQSR fails
    exec:
    if (task.exitStatus == 2) {
        log.warn "VQSR failed for INDELs. Switching to hard filtering."
    }
}

// Now modify the apply_vqsr_on_indels process to handle either VQSR or hard filtering
process apply_vqsr_on_indels {
    tag { "${params.project_name}.${params.cohort_id}.aVoI" }
    label 'gatk'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/vcfs/vqsr", mode: 'copy', overwrite: true
    
    input:
    tuple path(vcf), path(vcf_index),
        path(indel_recal), path(indel_recal_index), path(indel_tranches)
    
    output:
    tuple path("${params.project_name}.recal-INDEL.vcf.gz"),
        path("${params.project_name}.recal-INDEL.vcf.gz.tbi"), emit: indels_recalibrated

    script:
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    def min_mem = Math.min(1, max_mem)
    def filter_level = params.indel_filter_level ?: "99.0"

    """
    echo "Applying VQSR for INDELs using filter level: ${filter_level}"
    echo "Input VCF: ${vcf}"
    echo "INDEL Recal: ${indel_recal}"
    echo "INDEL Tranches: ${indel_tranches}"
    
    gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g" ApplyVQSR \\
        --reference ${params.reference_files.ref} \\
        --recal-file ${indel_recal} \\
        --tranches-file ${indel_tranches} \\
        --truth-sensitivity-filter-level ${filter_level} \\
        --create-output-variant-index true \\
        --mode INDEL \\
        --variant ${vcf} \\
        --output ${params.project_name}.recal-INDEL.vcf.gz
    """
}

process run_genotype_gvcfs {
    tag { "${params.project_name}.${interval}.rGGVCF" }
    label 'gatk'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/vcfs/genotyped", mode: 'copy',
        pattern: "*.{vcf.gz,vcf.gz.tbi}"
        
    input:
    tuple val(interval), val(database_path)
    
    output:
    tuple val(interval), path("${params.project_name}.${interval}.vcf.gz"), path("${params.project_name}.${interval}.vcf.gz.tbi"), emit: genotyped_vcf
    
    script:
    // Calculate memory safely - ensure at least 1GB for Xmx, but not more than 75% of available memory
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    def min_mem = Math.min(1, max_mem)
    
    """
    echo "Database path: ${database_path}"
    echo "Reference: ${params.reference_files.ref}"
    
    # Genotype GVCFs
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms${min_mem}g -Xmx${max_mem}g" \\
        GenotypeGVCFs \\
        --reference ${params.reference_files.ref} \\
        --variant gendb://${database_path} \\
        --output ${params.project_name}.${interval}.vcf.gz \\
        --dbsnp ${params.reference_files.dbsnp} \\
        --include-non-variant-sites
    """
}

