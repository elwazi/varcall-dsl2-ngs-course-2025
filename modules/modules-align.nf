#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// process validate_fastq {
//     tag "validating ${sample_id}"
//     publishDir "${params.outdir}/${params.project_name}/${params.workflow}/validation", mode: 'copy', pattern: '*.validation.txt'
    
//     label 'small_mem'
    
//     input:
//     tuple val(sample_id), path(fastq_r1), path(fastq_r2), val(flowcell), val(lane)
    
//     output:
//     path "${sample_id}.validation.txt", emit: validation_report
    
//     script:
//     def force_pass = params.containsKey('force_validation_pass') && params.force_validation_pass ? 'true' : 'false'
//     """
//     # Function to write to both file and stdout with clear formatting
//     write_report() {
//         echo -e "\$1" | tee -a ${sample_id}.validation.txt >&2
//     }

//     # Clear the file if it exists
//     > ${sample_id}.validation.txt

//     # Print header with clear separation
//     write_report "
// ============================================================"
//     write_report "                Validation Report for ${sample_id}"
//     write_report "============================================================
// "

//     # Check if files exist and are readable
//     write_report "1. File Access Check:"
//     if [ -r "${fastq_r1}" ] && [ -r "${fastq_r2}" ]; then
//         write_report "   ✓ Both FASTQ files are readable"
//     else
//         write_report "   ✗ One or both FASTQ files are not readable"
//     fi
//     write_report ""

//     # Check if files are gzipped
//     write_report "2. Compression Check:"
//     if [[ "${fastq_r1}" == *.gz ]] && [[ "${fastq_r2}" == *.gz ]]; then
//         write_report "   ✓ Files are gzipped"
        
//         # Test gzip integrity
//         write_report "3. Gzip Integrity Check:"
//         if gzip -t "${fastq_r1}" 2>/dev/null && gzip -t "${fastq_r2}" 2>/dev/null || [ "${force_pass}" = "true" ]; then
//             write_report "   ✓ Gzip integrity check passed ${force_pass == 'true' ? '(forced)' : ''}"
//         else
//             write_report "   ✗ One or both files failed gzip integrity check"
//         fi
        
//         # Count reads
//         READ_COUNT_1=\$(zcat "${fastq_r1}" 2>/dev/null | awk 'NR%4==1' | wc -l)
//         READ_COUNT_2=\$(zcat "${fastq_r2}" 2>/dev/null | awk 'NR%4==1' | wc -l)
//     else
//         write_report "   ! Files are not gzipped"
//         write_report ""
        
//         # Count reads
//         READ_COUNT_1=\$(awk 'NR%4==1' "${fastq_r1}" 2>/dev/null | wc -l)
//         READ_COUNT_2=\$(awk 'NR%4==1' "${fastq_r2}" 2>/dev/null | wc -l)
//     fi

//     # Report read counts
//     write_report "4. Read Count Check:"
//     write_report "   R1 reads: \$READ_COUNT_1"
//     write_report "   R2 reads: \$READ_COUNT_2"
//     if [ "\$READ_COUNT_1" = "\$READ_COUNT_2" ] || [ "${force_pass}" = "true" ]; then
//         write_report "   ✓ Read counts match between pairs ${force_pass == 'true' && "\$READ_COUNT_1" != "\$READ_COUNT_2" ? '(forced)' : ''}"
//     else
//         write_report "   ✗ Read counts do not match between pairs"
//     fi
//     write_report ""

//     # Add file sizes
//     write_report "5. File Size Information:"
//     write_report "   R1 size: \$(ls -lh ${fastq_r1} | awk '{print \$5}')"
//     write_report "   R2 size: \$(ls -lh ${fastq_r2} | awk '{print \$5}')"
//     write_report ""

//     # Add file paths
//     write_report "6. File Paths:"
//     write_report "   R1: ${fastq_r1}"
//     write_report "   R2: ${fastq_r2}"
//     write_report "
// ============================================================
// "

//     # Print a copy to stdout explicitly
//     cat ${sample_id}.validation.txt >&2
//     """
// }

process run_bwa {
    tag { sample_id }
    label 'bwa'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/bams", mode: 'copy',
        pattern: "*.sam"
    
    input:
    tuple val(sample_id), path(fastq_r1), path(fastq_r2), val(flowcell), val(lane)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam"), val(flowcell), val(lane), emit: bam_file
    
    script:
    """
    bwa mem \
        -R "@RG\\tID:${flowcell}.${lane}\\tPL:ILLUMINA\\tPU:${flowcell}.${lane}.${sample_id}\\tLB:${sample_id}\\tSM:${sample_id}" \
        -t ${task.cpus} \
        ${params.reference_files.ref} \
        ${fastq_r1} ${fastq_r2} \
        > ${sample_id}.sam
    """
}

process run_bam_sort {
    tag { sample_id }
    label 'samtools'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/bams", mode: 'copy',
        pattern: "*.bam"
    
    input:
    tuple val(sample_id), path(sam_file), val(flowcell), val(lane)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: raw_bam
    
    script:
    """
    # Convert SAM to BAM and sort
    samtools sort \
        -@ ${task.cpus-1} \
        -o ${sample_id}.bam \
        ${sam_file}
    
    # Index the BAM file
    samtools index -@ ${task.cpus-1} ${sample_id}.bam
    
    # Remove the SAM file to save space
    rm ${sam_file}
    """
}

process run_mark_duplicates {
    tag { "${sample_id}" }
    label 'gatk'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/bams", mode: 'copy',
        pattern: "*.md.*"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.md.bam"), path("${sample_id}.md.bai"), emit: md_bam

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def total_mem = task.memory.toGiga().doubleValue()
    def reserve_mem = Math.min(4.0, Math.max(1.0, total_mem / 4.0))
    def mem = Math.max(1.0, total_mem - reserve_mem)
    
    """
    # Run MarkDuplicates
    gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem.intValue()}g" \
        MarkDuplicates \
        --INPUT ${bam} \
        --OUTPUT ${sample_id}.md.bam \
        --METRICS_FILE ${sample_id}.md.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --CLEAR_DT false \
        --ADD_PG_TAG_TO_READS false
    """
}

process run_create_recalibration_table {
    tag { "${sample_id}" }
    label 'gatk'
    label 'large_mem'

    input:
    tuple val(sample_id), path(bam), path(index)
    
    output:
    tuple val(sample_id), path(bam), path(index), path("${sample_id}.recal.table"), emit: recal_table

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def total_mem = task.memory.toGiga().doubleValue()
    def reserve_mem = Math.min(4.0, Math.max(1.0, total_mem / 4.0))
    def mem = Math.max(1.0, total_mem - reserve_mem)

    """
    # Run BaseRecalibrator  
    echo "Running BaseRecalibrator..."
    gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem.intValue()}g" \
        BaseRecalibrator \
        --input ${bam} \
        --output ${sample_id}.recal.table \
        --tmp-dir . \
        -R ${params.reference_files.ref} \
        --known-sites ${params.reference_files.dbsnp} \
        --known-sites ${params.reference_files.known_indels_1} \
        --known-sites ${params.reference_files.known_indels_2} \
        --use-original-qualities
    """
}

process run_recalibrate_bam {
    tag { "${sample_id}" }
    label 'gatk'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/bams", mode: 'copy',
        pattern: "*.md.recal.*"

    input:
    tuple val(sample_id), path(bam), path(index), file(recal_table)

    output:
    tuple val(sample_id), path("${sample_id}.md.recal.bam"), path("${sample_id}.md.recal.bai"), emit: recal_bam

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def total_mem = task.memory.toGiga().doubleValue()
    def reserve_mem = Math.min(4.0, Math.max(1.0, total_mem / 4.0))
    def mem = Math.max(1.0, total_mem - reserve_mem)

    """
    echo "Running ApplyBQSR..."
    gatk --java-options "-XX:+UseSerialGC -Xms1g -Xmx${mem.intValue()}g" \
        ApplyBQSR \
        --input ${bam} \
        --output ${sample_id}.md.recal.bam \
        --tmp-dir . \
        -R ${params.reference_files.ref} \
        --create-output-bam-index true \
        --bqsr-recal-file ${recal_table} \
        --use-original-qualities \
        --static-quantized-quals 10 \
        --static-quantized-quals 20 \
        --static-quantized-quals 30
    """
}

process bam_to_cram {
    tag { "${sample_id}" }
    label 'samtools'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/crams", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(index)

    output:
    tuple val(sample_id), path(cram_out), path(cram_index), emit: cram_file
    tuple val(sample_id), val("${params.outdir}/${params.project_name}/${params.workflow}/crams/${cram_out}"), emit: cram_out_loc

    script:
    cram_out = "${bam.baseName}.cram"
    cram_index = "${bam.baseName}.cram.crai"
    """
    samtools view \
        --threads ${task.cpus-1} \
        --reference ${params.reference_files.ref} \
        --output-fmt cram,version=3.0 \
        -o ${cram_out} ${bam}
        
    samtools index \
        -@ ${task.cpus-1} \
        ${cram_out} \
        ${cram_index}
    """
}

process run_cram_flagstat {
    tag { "${sample_id}" }
    label 'samtools'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/crams", mode: 'copy',
        pattern: "*.{stats,md5}"

    input:
    tuple val(sample_id), path(cram), path(index)

    output:
    tuple val(sample_id), path(cram), path(index), path("${cram}.flagstat"), emit: cram_stats

    script:
    """
    samtools flagstat \
        --threads ${task.cpus-1} \
        ${cram} > ${cram}.flagstat
    """
}

process create_cram_md5sum {
    tag { "${sample_id}" }
    label 'samtools'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/crams", mode: 'copy',
        pattern: "*.{stats,md5}"

    input:
    tuple val(sample_id), path(cram), path(index), path(flagstat)

    output:
    tuple val(sample_id), path("${cram}.md5"), path("${index}.md5"), path("${flagstat}.md5"), emit: cram_all_md5sum

    script:
    """
    md5sum ${cram} > ${cram}.md5
    md5sum ${index} > ${index}.md5
    md5sum ${flagstat} > ${flagstat}.md5
    """
}

// process combine_validation_reports {
//     publishDir "${params.outdir}/${params.project_name}/${params.workflow}/validation", mode: 'copy'
    
//     label 'small_mem'
    
//     input:
//     tuple val(output_name), path(reports)
    
//     output:
//     path "${output_name}", emit: combined_report
//     tuple path("success_count.txt"), path("failed_count.txt"), path("validation_result.txt"), emit: validation_status
    
//     script:
//     """
//     # Create combined report header
//     echo "================================================" > ${output_name}
//     echo "          COMBINED FASTQ VALIDATION REPORT      " >> ${output_name}
//     echo "================================================" >> ${output_name}
//     echo "" >> ${output_name}
//     echo "Report generated on \$(date)" >> ${output_name}
//     echo "" >> ${output_name}
    
//     # Add individual reports
//     for report in ${reports}; do
//         echo "--------------------------------" >> ${output_name}
//         cat \$report >> ${output_name}
//         echo "" >> ${output_name}
//     done
    
//     # Add summary
//     echo "================================================" >> ${output_name}
//     echo "                   SUMMARY                      " >> ${output_name}
//     echo "================================================" >> ${output_name}
    
//     # Count successful files (no "✗" marks)
//     TOTAL=\$(ls -1 ${reports} | wc -l)
//     SUCCESS=\$(grep -L "✗" ${reports} | wc -l)
//     FAILED=\$((TOTAL - SUCCESS))
    
//     echo "Total files validated: \$TOTAL" >> ${output_name}
//     echo "Files with no issues: \$SUCCESS" >> ${output_name}
//     echo "Files with issues: \$FAILED" >> ${output_name}
//     echo "" >> ${output_name}
    
//     # Write summary values to separate files for use in the workflow
//     echo "\$SUCCESS" > success_count.txt
//     echo "\$FAILED" > failed_count.txt
    
//     if [ \$FAILED -eq 0 ]; then
//         echo "✅ All FASTQ files passed validation" >> ${output_name}
//         echo "PASS" > validation_result.txt
//     else
//         echo "❌ Some FASTQ files failed validation" >> ${output_name}
//         echo "   Please check individual reports for details" >> ${output_name}
//         echo "FAIL" > validation_result.txt
//     fi
    
//     echo "" >> ${output_name}
//     echo "================================================" >> ${output_name}
    
//     # Output summary to stdout
//     cat ${output_name}
//     """
// }

// workflow VALIDATE_FASTQ {
//     take:
//     samplesheet

//     main:
//     log.info """
//     ===============================================================================================
//     🧬 FASTQ Validation Workflow 🧬
//     ===============================================================================================
//     This workflow performs comprehensive validation on FASTQ files:
//     • Checking FASTQ file format and compression
//     • Validating read pairs and read counts
//     • Verifying sequence quality encoding
//     • Creating a combined validation summary report
//     ===============================================================================================
//     """
    
//     // Create input channel from samplesheet
//     Channel
//         .fromPath(samplesheet)
//         .splitCsv(header:true, sep:'\t')
//         .map { row -> 
//             def fastq_r1 = file(row.FastqR1)
//             def fastq_r2 = file(row.FastqR2)
//             def flowcell = row.containsKey('Flowcell') ? row.Flowcell : 'unknown'
//             def lane = row.containsKey('Lane') ? row.Lane : 'unknown'
            
//             if (!fastq_r1.exists()) {
//                 error """
//                 ERROR: FastqR1 file not found: ${fastq_r1}
//                 """
//                 System.exit(1)
//             }
//             if (!fastq_r2.exists()) {
//                 error """
//                 ERROR: FastqR2 file not found: ${fastq_r2}
//                 """
//                 System.exit(1)
//             }
            
//             return [row.SampleID, fastq_r1, fastq_r2, flowcell, lane]
//         }
//         .set { fastq_ch }
    
//     // Run validation on each sample
//     validate_fastq(fastq_ch)
    
//     // Combine validation reports
//     validate_fastq.out.validation_report
//         .collect()
//         .map { reports -> 
//             ["combined_validation_report.txt", reports]
//         }
//         .set { all_reports }
    
//     // COMBINE_VALIDATION_REPORTS process
//     combine_validation_reports(all_reports)
    
//     // Log completion and results message
//     combine_validation_reports.out.validation_status.last().subscribe { successFile, failedFile, resultFile ->
//         def successCount = file(successFile).text.trim()
//         def failedCount = file(failedFile).text.trim()
//         def result = file(resultFile).text.trim()
        
//         log.info """
//     ===============================================================================================
//     🧬 FASTQ Validation Workflow Completed 🧬
//     ===============================================================================================
//     • ${successCount} files validated with NO issues ${result == "PASS" ? "✅" : ""}
//     • ${failedCount} files validated with issues ${result == "FAIL" ? "❌" : ""}
//     • Validation reports are available in: ${params.outdir}/${params.project_name}/${params.workflow}/validation/
//     • Combined validation report: ${params.outdir}/${params.project_name}/${params.workflow}/validation/combined_validation_report.txt
//     ===============================================================================================
//     """
//     }

//     emit:
//     validation_reports = validate_fastq.out.validation_report
//     combined_report = combine_validation_reports.out.combined_report
// }

// workflow ALIGN {
//     take:
//         sample_sheet_ch

//     main:
//         // First validate FASTQ files
//         VALIDATE_FASTQ(sample_sheet_ch)

//         // Run BWA alignment with validated FASTQs
//         run_bwa(VALIDATE_FASTQ.out.validated_fastqs)
        
//         // Sort BAM
//         run_bam_sort(run_bwa.out.bam_file)
        
//         // Mark duplicates
//         run_mark_duplicates(run_bam_sort.out.raw_bam)
        
//         // Create recalibration table
//         run_create_recalibration_table(run_mark_duplicates.out.md_bam)
        
//         // Apply BQSR
//         run_recalibrate_bam(run_create_recalibration_table.out.recal_table)
        
//         // Convert to CRAM
//         bam_to_cram(run_recalibrate_bam.out.recal_bam)
        
//         // Generate CRAM statistics
//         run_cram_flagstat(bam_to_cram.out.cram_file)
        
//         // Create MD5 checksums
//         create_cram_md5sum(run_cram_flagstat.out.cram_stats)

//     emit:
//         cram_file = bam_to_cram.out.cram_file
//         cram_stats = run_cram_flagstat.out.cram_stats
//         cram_md5 = create_cram_md5sum.out.cram_all_md5sum
//         cram_out_loc = bam_to_cram.out.cram_out_loc
// }

// workflow {
//     // Check if sample sheet is provided
//     if (!params.sample_sheet) {
//         error "Please provide a sample sheet with --sample_sheet"
//     }

//     // Parse sample sheet and create channel
//     Channel
//         .fromPath(params.sample_sheet)
//         .splitCsv(header:true, sep:'\t')
//         .map { row -> 
//             // Validate required fields
//             if (!row.SampleID || !row.FastqR1 || !row.FastqR2) {
//                 error "Missing required fields in sample sheet. Each row must have SampleID, FastqR1, and FastqR2"
//             }
            
//             // Create tuple with default values for NA fields
//             tuple(
//                 row.SampleID,
//                 file(row.FastqR1),
//                 file(row.FastqR2),
//                 row.Flowcell == 'NA' ? 'default' : row.Flowcell,
//                 row.Lane == 'NA' ? '1' : row.Lane
//             )
//         }
//         .set { sample_sheet_ch }

//     // Choose workflow based on params
//     if (params.workflow == "validate_fastq") {
//         VALIDATE_FASTQ(sample_sheet_ch)
//     } else if (params.workflow == "align") {
//         ALIGN(sample_sheet_ch)
//     } else {
//         error "Please specify --workflow as either 'validate_fastq' or 'align'"
//     }
// }
