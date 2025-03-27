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

process validate_gvcf {
    tag { "${sample_id}.${gvcf.name}" }
    label 'gatk'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/validation", mode: 'copy',
        pattern: "*.{validation,summary.txt}"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(gvcf), path(index)

    output:
    tuple val(sample_id), path("${gvcf.simpleName}.validation"), path("${gvcf.simpleName}.summary.txt"), emit: validation_results

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx
    def mem = getMem(task.memory)
    def file_basename = gvcf.simpleName

    """
    # Run GATK ValidateVariants with safer options
    set +e  # Don't exit on error
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms1g -Xmx${mem.intValue()}g" ValidateVariants \\
        --reference ${params.reference_files.ref} \\
        --variant ${gvcf} \\
        --dbsnp ${params.reference_files.dbsnp} \\
        > ${file_basename}.validation 2>&1
    validate_exit=\$?
    set -e
    
    # Generate summary statistics
    echo "Validation summary for ${sample_id} (${file_basename}):" > ${file_basename}.summary.txt
    echo "---------------------------------------------------" >> ${file_basename}.summary.txt
    echo "Date: \$(date)" >> ${file_basename}.summary.txt
    echo "" >> ${file_basename}.summary.txt
    
    # Extract variant counts using bcftools (if available) or grep
    echo "Variant statistics:" >> ${file_basename}.summary.txt
    if command -v bcftools >/dev/null 2>&1; then
        bcftools stats ${gvcf} 2>/dev/null | grep "number of" >> ${file_basename}.summary.txt || echo "Could not extract variant statistics" >> ${file_basename}.summary.txt
    else
        echo "bcftools not available - skipping variant count statistics" >> ${file_basename}.summary.txt
    fi
    
    # Check for validation status
    if [ \${validate_exit} -ne 0 ] || grep -q "ERROR" ${file_basename}.validation; then
        echo "" >> ${file_basename}.summary.txt
        echo "VALIDATION ERRORS FOUND:" >> ${file_basename}.summary.txt
        grep -E "ERROR|Exception|fail" ${file_basename}.validation >> ${file_basename}.summary.txt 2>/dev/null || echo "No specific error message found" >> ${file_basename}.summary.txt
        echo "" >> ${file_basename}.summary.txt
        echo "VALIDATION FAILED" >> ${file_basename}.summary.txt
    else
        echo "" >> ${file_basename}.summary.txt
        echo "No validation errors found." >> ${file_basename}.summary.txt
        echo "VALIDATION PASSED" >> ${file_basename}.summary.txt
    fi
    
    # Ensure the process succeeds even if validation fails
    exit 0
    """
}

process generate_validation_report {
    tag { "validation_report" }
    label 'small_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/validation", mode: 'copy',
        pattern: "*.{html,txt}"
    errorStrategy 'ignore'

    input:
    path summaries

    output:
    path "validation_report.html", emit: html_report
    path "validation_summary.txt", emit: text_summary

    script:
    """
    # Debug info
    echo "Current directory: \$(pwd)"
    echo "Listing summary files:"
    ls -la 
    
    # Create a summary text file
    echo "gVCF Validation Summary Report" > validation_summary.txt
    echo "============================" >> validation_summary.txt
    echo "Generated: \$(date)" >> validation_summary.txt
    echo "" >> validation_summary.txt
    
    # Process all summary files to count passes and failures first
    passed=0
    failed=0
    for summary in *.summary.txt; do
        if [ -f "\$summary" ]; then
            # Count passes and failures inline
            if grep -q "VALIDATION PASSED" "\$summary"; then
                passed=\$((passed + 1))
            elif grep -q "VALIDATION FAILED" "\$summary"; then
                failed=\$((failed + 1))
            fi
        fi
    done
    
    # Calculate total
    total=\$((passed + failed))
    
    # Add summary at the top
    echo "SUMMARY STATISTICS" >> validation_summary.txt
    echo "===================" >> validation_summary.txt
    echo "Total files validated: \$total" >> validation_summary.txt
    echo "Passed: \$passed" >> validation_summary.txt
    echo "Failed: \$failed" >> validation_summary.txt
    echo "" >> validation_summary.txt
    echo "" >> validation_summary.txt
    
    # Process all summary files with simpler approach
    echo "DETAILED VALIDATION RESULTS" >> validation_summary.txt
    echo "=========================" >> validation_summary.txt
    echo "" >> validation_summary.txt
    for summary in *.summary.txt; do
        if [ -f "\$summary" ]; then
            echo "Adding \$summary to report" >> validation_summary.txt
            echo "File: \$summary" >> validation_summary.txt
            cat "\$summary" >> validation_summary.txt
            
            echo "" >> validation_summary.txt
            echo "------------------------------------" >> validation_summary.txt
            echo "" >> validation_summary.txt
        fi
    done
    
    # Generate HTML report with simpler approach
    cat <<EOF > validation_report.html
    <!DOCTYPE html>
    <html>
    <head>
        <title>gVCF Validation Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            h1 { color: #2c3e50; }
            h2 { color: #3498db; }
            .summary { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; }
            .passed { color: green; }
            .failed { color: red; }
            .warning { color: orange; }
            table { border-collapse: collapse; width: 100%; margin-top: 20px; }
            th, td { text-align: left; padding: 8px; border-bottom: 1px solid #ddd; }
            th { background-color: #3498db; color: white; }
            tr:nth-child(even) { background-color: #f2f2f2; }
            .details { white-space: pre-wrap; font-family: monospace; background-color: #f0f0f0; padding: 10px; }
        </style>
    </head>
    <body>
        <h1>gVCF Validation Report</h1>
        <p>Generated: \$(date)</p>
        
        <div class="summary">
            <h2>SUMMARY STATISTICS</h2>
            <p>Total files validated: <strong>\$total</strong></p>
            <p>Passed: <strong class="passed">\$passed</strong></p>
            <p>Failed: <strong class="failed">\$failed</strong></p>
        </div>
        
        <h2>Validation Details</h2>
        <table>
            <tr>
                <th>Sample</th>
                <th>Status</th>
                <th>Details</th>
            </tr>
EOF
    
    # Add each file's details to the HTML with simpler approach
    filesProcessed=0
    for summary in *.summary.txt; do
        if [ -f "\$summary" ]; then
            filesProcessed=\$((filesProcessed + 1))
            # Extract sample name from filename
            sample=\$(basename "\$summary" .summary.txt)
            
            if grep -q "VALIDATION PASSED" "\$summary"; then
                status="PASSED"
                status_class="passed"
            elif grep -q "VALIDATION FAILED" "\$summary"; then
                status="FAILED"
                status_class="failed"
            else
                status="UNKNOWN"
                status_class="warning"
            fi
            
            # Get content safely
            content=\$(cat "\$summary" | sed 's/&/\\&amp;/g' | sed 's/</\\&lt;/g' | sed 's/>/\\&gt;/g')
            
            cat <<EOF >> validation_report.html
        <tr>
            <td>\$sample</td>
            <td class="\$status_class">\$status</td>
            <td><details>
                <summary>Click to view details</summary>
                <div class="details">\$content</div>
            </details></td>
        </tr>
EOF
        fi
    done
    
    # If no files were processed, add a notice row
    if [ "\$filesProcessed" -eq 0 ]; then
        cat <<EOF >> validation_report.html
        <tr>
            <td colspan="3" class="warning" style="text-align: center;">No validation files were found or processed</td>
        </tr>
EOF
    fi
    
    # Close the HTML file
    cat <<EOF >> validation_report.html
        </table>
    </body>
    </html>
EOF
    """
}

workflow VALIDATE_GVCF {
    take:
    gvcf_files // Channel of sample_id, gvcf, index tuples
    
    main:
    log.info """
    ===============================================================================================
    🧬 gVCF Validation Workflow 🧬
    ===============================================================================================
    This workflow validates gVCF files and generates validation reports:
    • Running GATK ValidateVariants on each gVCF
    • Generating summary statistics for each file
    • Creating a consolidated validation report
    • Inputs:
        • gVCF files
        • Reference genome: ${params.reference_files.ref}
        • dbSNP: ${params.reference_files.dbsnp}
        • Validation Report folder: ${params.outdir}/${params.project_name}/${params.workflow}/validation
        • Validation Report: ${params.outdir}/${params.project_name}/${params.workflow}/validation/validation_report.html
        • Validation Summary: ${params.outdir}/${params.project_name}/${params.workflow}/validation/validation_summary.txt
    ===============================================================================================
    """
    
    // Validate each gVCF file
    validate_gvcf(gvcf_files)
    
    // Collect all validation results
    validate_gvcf.out.validation_results
        .map { it[2] } // Extract just the summary files
        .collect()
        .set { all_summaries }
    
    // Generate validation report
    generate_validation_report(all_summaries)
    
    emit:
    validation_results = validate_gvcf.out.validation_results
    html_report = generate_validation_report.out.html_report
    text_summary = generate_validation_report.out.text_summary
} 