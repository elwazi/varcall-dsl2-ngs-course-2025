#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Function to set build-specific parameters
def setBuildSpecificParams(build) {
    def params = [:]
    
    if (build == "b37") {
        params.x = "X"
        params.y = "Y"
        params.x_par1 = "60001-2699520"
        params.x_par1_target_default = "60001-60001"
        params.x_par2 = "154931044-155260560"
        params.x_par2_target_default = "154931044-154931044"
        params.y_par1 = "10001-2649520"
        params.y_par1_target_default = "10001-10001"
        params.y_par2 = "59034050-59363566"
        params.y_par2_target_default = "59034050-59034050"
        params.mt = "MT"
    } else if (build == "b38") {
        params.x = "chrX"
        params.y = "chrY"
        params.x_par1 = "10001-2781479"
        params.x_par1_target_default = "10001-10001"
        params.x_par2 = "155701383-156030895"
        params.x_par2_target_default = "155701383-155701383"
        params.y_par1 = "10001-2781479"
        params.y_par1_target_default = "10001-10001"
        params.y_par2 = "56887903-57217415"
        params.y_par2_target_default = "56887903-56887903"
        params.mt = "chrM"
    } else {
        error "Please specify a genome build (b37 or b38)!"
    }
    return params
}

// Set call confidence based on coverage
def getCallConf(coverage) {
    return coverage == "low" ? 10 : 30
}

// Get memory safely
def getMem(task_mem) {
    // Calculate memory safely - ensure at least 1GB for Xmx
    def total_mem = task_mem.toGiga().doubleValue()
    def reserve_mem = Math.min(4.0, Math.max(1.0, total_mem / 4.0))
    def mem = Math.max(1.0, total_mem - reserve_mem)
    return mem
}

// Move validateChromosomes outside the workflow
process validateChromosomes {
    label 'small_mem'
    
    input:
    val autosomes
    
    output:
    path "valid_chromosomes.txt", emit: valid_chrs
    
    script:
    """
    # Create a list of valid chromosomes
    touch valid_chromosomes.txt
    
    for chr in ${autosomes.join(' ')}; do
        if samtools faidx ${params.reference_files.ref} \$chr >/dev/null 2>&1; then
            echo "\$chr" >> valid_chromosomes.txt
        else
            echo "WARNING: Chromosome \$chr not found in reference genome" >&2
        fi
    done
    """
}

// RUN HAPLOTYPE CALLER ON AUTOSOMES
process run_haplotype_caller_auto {
    tag { "${sample_id}.${autosome}.rHCoA" }
    label 'gatk'
    label 'gatk_haplotypecaller'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.{g.vcf.gz,g.vcf.gz.tbi}"

    input:
    tuple val(sample_id), val(gender), path(bam)
    each autosome
    
    output:
    tuple val(sample_id), path("${sample_id}.${autosome}.g.vcf.gz"), path("${sample_id}.${autosome}.g.vcf.gz.tbi"), emit: auto_calls

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def mem = getMem(task.memory)
    
    // Get call confidence based on coverage
    def call_conf = getCallConf(params.sample_coverage)

    def intervals = params.type == "wes" ? 
        "--intervals ${params.target_regions} --intervals ${autosome} --interval-set-rule INTERSECTION" :
        "--intervals ${autosome}"
    
    """
    echo "Debug information:"
    echo "Build: ${params.build}"
    echo "Autosome: ${autosome}"
    echo "Reference: ${params.reference_files.ref}"
    echo "Intervals: ${intervals}"
    
    # Directly check if chromosome exists in reference
    if [[ "${params.reference_files.ref}" == *".${autosome}."* ]]; then
        echo "Using chromosome-specific reference file for ${autosome}"
    elif ! samtools faidx ${params.reference_files.ref} ${autosome} >/dev/null 2>&1; then
        echo "ERROR: Chromosome ${autosome} not found in reference genome"
        exit 1
    fi
    
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms1g -Xmx${mem.intValue()}g" \\
        HaplotypeCaller \\
        --reference ${params.reference_files.ref} \\
        --input ${bam} \\
        --emit-ref-confidence GVCF \\
        --dbsnp ${params.reference_files.dbsnp} \\
        ${intervals} \\
        -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \\
        -stand-call-conf ${call_conf} \\
        --sample-ploidy 2 \\
        -O ${sample_id}.${autosome}.g.vcf.gz
    """
}

// RUN HAPLOTYPE CALLER ON MALES
process run_haplotype_caller_males {
    tag { "${sample_id}:${nonautosome}.rHCoNAmales" }
    label 'gatk'
    label 'gatk_haplotypecaller'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.{g.vcf.gz,g.vcf.gz.tbi}"

    input:
    tuple val(sample_id), val(gender), path(bam)
    each nonautosome
    
    output:
    tuple val(sample_id), path("${sample_id}.*.g.vcf.gz"), path("${sample_id}.*.g.vcf.gz.tbi"), emit: male_calls

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def mem = getMem(task.memory)
    
    // Get call confidence based on coverage
    def call_conf = getCallConf(params.sample_coverage)
    
    def base
    def ploidy
    def intervals
    
    if (nonautosome == 'x_par1_male') {
        base = "X_PAR1"
        ploidy = "2"
        intervals = params.type == "wes" ? 
            "--intervals ${params.target_regions} --intervals ${params.x}:${params.x_par1} --interval-set-rule INTERSECTION" :
            "--intervals ${params.x}:${params.x_par1}"
    } else if (nonautosome == 'x_par2_male') {
        base = "X_PAR2"
        ploidy = "2"
        intervals = params.type == "wes" ? 
            "--intervals ${params.target_regions} --intervals ${params.x}:${params.x_par2} --interval-set-rule INTERSECTION" :
            "--intervals ${params.x}:${params.x_par2}"
    } else if (nonautosome == 'x_nonpar_male') {
        base = "X_nonPAR"
        ploidy = "1"
        intervals = params.type == "wes" ? 
            "--intervals ${params.target_regions} --intervals ${params.x} --exclude-intervals ${params.x}:${params.x_par1} --exclude-intervals ${params.x}:${params.x_par2} --interval-set-rule INTERSECTION" :
            "--intervals ${params.x} --exclude-intervals ${params.x}:${params.x_par1} --exclude-intervals ${params.x}:${params.x_par2}"
    } else if (nonautosome == 'y_par1_male') {
        base = "Y_PAR1"
        ploidy = "2"
        intervals = params.type == "wes" ? 
            "--intervals ${params.target_regions} --intervals ${params.y}:${params.y_par1} --interval-set-rule INTERSECTION" :
            "--intervals ${params.y}:${params.y_par1}"
    } else if (nonautosome == 'y_par2_male') {
        base = "Y_PAR2"
        ploidy = "2"
        intervals = params.type == "wes" ? 
            "--intervals ${params.target_regions} --intervals ${params.y}:${params.y_par2} --interval-set-rule INTERSECTION" :
            "--intervals ${params.y}:${params.y_par2}"
    } else if (nonautosome == 'y_nonpar_male') {
        base = "Y_nonPAR"
        ploidy = "1"
        intervals = params.type == "wes" ? 
            "--intervals ${params.target_regions} --intervals ${params.y} --exclude-intervals ${params.y}:${params.y_par1} --exclude-intervals ${params.y}:${params.y_par2} --interval-set-rule INTERSECTION" :
            "--intervals ${params.y} --exclude-intervals ${params.y}:${params.y_par1} --exclude-intervals ${params.y}:${params.y_par2}"
    } else {
        error "Unknown nonautosome value: ${nonautosome}"
    }

    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms1g -Xmx${mem.intValue()}g" \\
        HaplotypeCaller \\
        --reference ${params.reference_files.ref} \\
        --input ${bam} \\
        --emit-ref-confidence GVCF \\
        --dbsnp ${params.reference_files.dbsnp} \\
        ${intervals} \\
        -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \\
        -stand-call-conf ${call_conf} \\
        --sample-ploidy ${ploidy} \\
        -O ${sample_id}.${base}.g.vcf.gz
    """
}

// RUN HAPLOTYPE CALLER ON FEMALES
process run_haplotype_caller_females {
    tag { "${sample_id}:${params.x}.rHCoNAfemales" }
    label 'gatk'
    label 'gatk_haplotypecaller'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.{g.vcf.gz,g.vcf.gz.tbi}"

    input:
    tuple val(sample_id), val(gender), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.${params.x}.g.vcf.gz"), path("${sample_id}.${params.x}.g.vcf.gz.tbi"), emit: female_calls

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def mem = getMem(task.memory)
    
    // Get call confidence based on coverage
    def call_conf = getCallConf(params.sample_coverage)
    
    def intervals
    if (params.type == "wes") {
        intervals = "--intervals ${params.target_regions} --intervals ${params.x} --interval-set-rule INTERSECTION"
    } else {
        intervals = "--intervals ${params.x}"
    }
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms1g -Xmx${mem.intValue()}g" \\
        HaplotypeCaller \\
        --reference ${params.reference_files.ref} \\
        --input ${bam} \\
        --emit-ref-confidence GVCF \\
        --dbsnp ${params.reference_files.dbsnp} \\
        ${intervals} \\
        -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \\
        -stand-call-conf ${call_conf} \\
        --sample-ploidy 2 \\
        -O ${sample_id}.${params.x}.g.vcf.gz
    """
}

// RUN HAPLOTYPE CALLER ON MT
process run_haplotype_caller_mt {
    tag { "${sample_id}:${params.mt}.rHCoNAmt" }
    label 'gatk'
    label 'gatk_haplotypecaller'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.{g.vcf.gz,g.vcf.gz.tbi}"

    input:
    tuple val(sample_id), val(gender), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.${params.mt}.g.vcf.gz"), path("${sample_id}.${params.mt}.g.vcf.gz.tbi"), emit: mt_calls

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def mem = getMem(task.memory)
    
    // Get call confidence based on coverage
    def call_conf = getCallConf(params.sample_coverage)
    
    def intervals
    if (params.type == "wes") {
        intervals = "--intervals ${params.mt}"
    } else {
        intervals = "--intervals ${params.mt}"
    }
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms1g -Xmx${mem.intValue()}g" \\
        HaplotypeCaller \\
        --reference ${params.reference_files.ref} \\
        --input ${bam} \\
        --emit-ref-confidence GVCF \\
        --dbsnp ${params.reference_files.dbsnp} \\
        ${intervals} \\
        -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \\
        -stand-call-conf ${call_conf} \\
        --sample-ploidy 2 \\
        -O ${sample_id}.${params.mt}.g.vcf.gz
    """
}

process run_sort_male_gvcfs {
    tag { "${sample_id}.sMgVCF" }
    label 'gatk'
    label 'medium_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.{g.vcf.gz,g.vcf.gz.tbi}"

    input:
    tuple val(sample_id), path(vcf), path(index) 
    
    output:
    tuple val(sample_id), path("${sample_id}.X.g.vcf.gz"), path("${sample_id}.Y.g.vcf.gz"), path("${sample_id}.X.g.vcf.gz.tbi"), path("${sample_id}.Y.g.vcf.gz.tbi"), emit: male_calls_combined
    
    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def mem = getMem(task.memory)
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem.intValue()}g" SortVcf \
        --INPUT ${sample_id}.X_PAR1.g.vcf.gz \
        --INPUT ${sample_id}.X_PAR2.g.vcf.gz \
        --INPUT ${sample_id}.X_nonPAR.g.vcf.gz \
        --OUTPUT ${sample_id}.X.g.vcf.gz
  
    gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem.intValue()}g" SortVcf \
        --INPUT ${sample_id}.Y_PAR1.g.vcf.gz \
        --INPUT ${sample_id}.Y_PAR2.g.vcf.gz \
        --INPUT ${sample_id}.Y_nonPAR.g.vcf.gz \
        --OUTPUT ${sample_id}.Y.g.vcf.gz
    """
}

process run_combine_sample_gvcfs {
    tag { "${sample_id}.cCgVCF" }
    label 'gatk'
    label 'extra_large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.{g.vcf.gz,g.vcf.gz.tbi}"
    
    input:
    tuple val(sample_id), path(gvcfs), path(indexes)
  
    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), emit: combined_gvcf
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), emit: gvcf_out_loc
    
    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def mem = getMem(task.memory)

    def gvcf_list = "${sample_id}.gvcf.list"
    """
    # Create list of GVCFs to combine
    for gvcf in ${gvcfs}; do
        echo "\$gvcf" >> ${gvcf_list}
    done
    
    gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem.intValue()}g" GatherVcfs \\
        --INPUT ${gvcf_list} \\
        --OUTPUT ${sample_id}.g.vcf.gz
    """
}

process run_combine_sample_gvcfs_index {
    tag { "${sample_id}.indexGVCF" }
    label 'bcftools'
    label 'small_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.{g.vcf.gz,g.vcf.gz.tbi}"

    input:
    tuple val(sample_id), path(gvcf), path(index)
  
    output:
    tuple val(sample_id), path(gvcf), path("${gvcf}.tbi"), emit: indexed_gvcf
    tuple val(sample_id), path(gvcf), emit: gvcf_out_loc

    script:
    """
    # Index the GVCF file
    tabix -p vcf ${gvcf}
    """
}

process run_index_gvcf {
    tag { "${sample_id}.index" }
    label 'bcftools'
    label 'small_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.{g.vcf.gz,g.vcf.gz.tbi}"

    input:
    tuple val(sample_id), path(gvcf)

    output:
    tuple val(sample_id), path(gvcf), path("${gvcf}.tbi"), emit: indexed_gvcf

    script:
    """
    tabix -p vcf ${gvcf}
    """
}

process run_create_gvcf_md5sum {
    tag { "${sample_id}.cGMD5" }
    label 'small_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs/md5", mode: 'copy',
        pattern: "*.md5"

    input:
    tuple val(sample_id), path(gvcf), path(index)

    output:
    tuple val(sample_id), path("${gvcf}.md5"), path("${index}.md5"), emit: gvcf_md5sum

    script:
    """
    md5sum ${gvcf} > ${gvcf}.md5
    md5sum ${index} > ${index}.md5
    """
}

process run_validate_gvcf {
    tag { "${sample_id}" }
    label 'gatk'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs", mode: 'copy',
        pattern: "*.validation"

    input:
    tuple val(sample_id), path(gvcf), path(index)
    each autosome

    output:
    tuple val(sample_id), path("${gvcf}.validation"), path("${index}.validation"), emit: gvcf_validation

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def mem = getMem(task.memory)

    def intervals = params.type == "wes" ? 
        "--intervals ${params.target_regions} --intervals ${autosome} --interval-set-rule INTERSECTION" :
        "--intervals ${autosome}"

    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms1g -Xmx${mem.intValue()}g" ValidateVariants \
        --reference ${params.reference_files.ref} \
        --input ${gvcf} \
        ${intervals} \
        --dbsnp ${params.reference_files.dbsnp} \
        --output ${gvcf}.validation
    """
}

// Add this workflow at the bottom of your file
workflow GENERATE_GVCFS {
    take:
    bam_files // Channel of sample_id, gender, bam tuples
    
    main:
    // Initialize build-specific parameters within the workflow
    def buildParams = setBuildSpecificParams(params.build)
    
    // First, create a channel to check which chromosomes exist in the reference
    Channel
        .of(1..22)
        .map { it.toString() }
        .map { autosome -> 
            def chr_prefix = params.build == "b38" ? "chr" : ""
            return "${chr_prefix}${autosome}"
        }
        .set { all_autosomes }
    
    // Run chromosome validation
    validateChromosomes(all_autosomes.collect())
    
    // Create a channel of valid autosomes
    validateChromosomes.out.valid_chrs
        .splitText()
        .map { it.trim() }
        .filter { it != "" }
        .set { valid_autosomes }
    
    // Run HaplotypeCaller on valid autosomes
    run_haplotype_caller_auto(bam_files, valid_autosomes)
    
    // Sex chromosome handling depends on gender
    bam_files
        .branch {
            male: it[1].toLowerCase() == "male" || it[1].toLowerCase() == "m"
            female: it[1].toLowerCase() == "female" || it[1].toLowerCase() == "f"
        }
        .set { sex_specific_bams }
    
    // Male-specific processing
    if (sex_specific_bams.male) {
        Channel
            .of('x_par1_male', 'x_par2_male', 'x_nonpar_male', 'y_par1_male', 'y_par2_male', 'y_nonpar_male')
            .set { male_regions }
            
        run_haplotype_caller_males(sex_specific_bams.male, male_regions)
        run_sort_male_gvcfs(run_haplotype_caller_males.out.male_calls.groupTuple())
    }
    
    // Female-specific processing
    if (sex_specific_bams.female) {
        run_haplotype_caller_females(sex_specific_bams.female)
    }
    
    // MT processing for all samples
    run_haplotype_caller_mt(bam_files)
    
    // Combine all gVCFs
    autosome_calls = run_haplotype_caller_auto.out.auto_calls.groupTuple()
    
    // Combine different types of calls based on gender
    bam_files
        .map { it[0] } // Get just the sample IDs
        .combine(autosome_calls, by: 0) // Join with autosome calls
        .map { sample_id, gvcfs, indexes -> 
            def gender = bam_files.find { it[0] == sample_id }[1]
            return [sample_id, gender, gvcfs.flatten(), indexes.flatten()]
        }
        .branch {
            male: it[1].toLowerCase() == "male" || it[1].toLowerCase() == "m"
            female: it[1].toLowerCase() == "female" || it[1].toLowerCase() == "f"
        }
        .set { combined_by_gender }
    
    // For males, add X and Y
    if (combined_by_gender.male) {
        combined_by_gender.male
            .join(run_sort_male_gvcfs.out.male_calls_combined, by: 0)
            .map { sample_id, gender, auto_gvcfs, auto_indexes, x_gvcf, y_gvcf, x_index, y_index ->
                [sample_id, auto_gvcfs + [x_gvcf, y_gvcf], auto_indexes + [x_index, y_index]]
            }
            .set { male_combined }
    }
    
    // For females, add X
    if (combined_by_gender.female) {
        combined_by_gender.female
            .join(run_haplotype_caller_females.out.female_calls, by: 0)
            .map { sample_id, gender, auto_gvcfs, auto_indexes, x_gvcf, x_index ->
                [sample_id, auto_gvcfs + [x_gvcf], auto_indexes + [x_index]]
            }
            .set { female_combined }
    }
    
    // Add MT for all samples
    def all_combined_ch = channel.empty()
    if (combined_by_gender.male) {
        all_combined_ch = all_combined_ch.mix(male_combined)
    }
    if (combined_by_gender.female) {
        all_combined_ch = all_combined_ch.mix(female_combined)
    }
    
    all_combined_ch
        .join(run_haplotype_caller_mt.out.mt_calls, by: 0)
        .map { sample_id, gvcfs, indexes, mt_gvcf, mt_index ->
            [sample_id, gvcfs + [mt_gvcf], indexes + [mt_index]]
        }
        .set { final_combined }
    
    // Run final process to combine all gVCFs - split into two steps
    run_combine_sample_gvcfs(final_combined)
    run_combine_sample_gvcfs_index(run_combine_sample_gvcfs.out.combined_gvcf)
    
    // Use the indexed output for downstream processes
    run_create_gvcf_md5sum(run_combine_sample_gvcfs_index.out.indexed_gvcf)
    run_validate_gvcf(run_combine_sample_gvcfs_index.out.indexed_gvcf)

    emit:
    gvcfs = run_combine_sample_gvcfs_index.out.indexed_gvcf
    md5sums = run_create_gvcf_md5sum.out.gvcf_md5sum
    validation = run_validate_gvcf.out.gvcf_validation
}
