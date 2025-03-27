#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process print_sample_info {
    tag { sample_id }
    debug true

    input:
    tuple val(sample_id), val(gender), path(fastq_r1), path(fastq_r2), val(flowcell), val(lane), path(bam), path(gvcf)

    script:
    """
    printf "[Sample Info]"
    printf "Sample   : ${sample_id}"
    printf "Fastq R1 : ${fastq_r1}"
    printf "Fastq R2 : ${fastq_r2}"
    printf "Flowcell : ${flowcell}"
    printf "Lane     : ${lane}\n"
    printf "BAM      : ${lane}\n"
    printf "GVCF     : ${lane}\n"
    printf "--------------------------------------------------
    """
}

// Split the combined process into separate processes for each tool
process log_tool_version_samtools {
    tag "tool_ver"
    label 'bwa_samtools'
    label 'small_mem'
    
    output:
    path("tool.samtools.version"), emit: samtools_version
    
    script:
    """
    samtools --version > tool.samtools.version
    """
}

process log_tool_version_bwa {
    tag "tool_ver"
    label 'bwa_samtools'
    label 'small_mem'
    
    output:
    path("tool.bwa.version"), emit: bwa_version
    
    script:
    """
    bwa 2>&1 | head -n 5 > tool.bwa.version
    """
}

// Keep the original process for backward compatibility but make it call the new processes
process log_tool_version_samtools_bwa {
    tag "tool_ver"
    label 'small_mem'
    
    output:
    path("tool.samtools.bwa.version"), emit: version
    
    script:
    """
    echo "WARNING: This process is deprecated. Please use log_tool_version_samtools and log_tool_version_bwa instead." > tool.samtools.bwa.version
    echo "============================================================" >> tool.samtools.bwa.version
    samtools --version >> tool.samtools.bwa.version
    echo "============================================================" >> tool.samtools.bwa.version
    bwa 2>&1 | head -n 5 >> tool.samtools.bwa.version
    """
}

// Log tool versions
process log_tool_version_gatk {
    tag { "tool_ver" }
    label 'gatk'
    debug true
    publishDir "${params.outdir}/${params.workflow}", mode: 'copy', overwrite: true
    
    output:
    path("tool.gatk.version"), emit: tool_version_gatk
    
    script:
    mem = task.memory.toGiga() - 1
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms500m -Xmx${mem}g" --version > tool.gatk.version 2>&1
    """
}

// CHECK IF THE FILES EXIST
def check_files(file_list) {
    file_list.each { myfile ->
        if (!myfile) {
            exit 1, "File path is null or empty"
        }
        
        def f = file(myfile)
        if (!f.exists()) {
            exit 1, "File not found: ${myfile}"
        }
        if (!f.isFile()) {
            exit 1, "Path exists but is not a file: ${myfile}" 
        }
    }
}


// Helper function to generate updated samplesheet name
def get_updated_samplesheet_name(original_name, workflow_step) {
    def base_name = original_name.replaceAll(/\.[^.]+$/, '')
    return "${base_name}_${workflow_step}.tsv"
}

// Update the samplesheet
process update_samplesheet {
    tag "update_samplesheet_${workflow_type}"
    label 'process_low'
    // Publish to a dedicated directory in the output folder
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/samplesheets", mode: 'copy', overwrite: true
    // Also publish to the original sample sheet directory for backward compatibility
    publishDir { "${file(params.sample_sheet).parent}" }, mode: 'copy', overwrite: true
    cache 'deep'
    
    input:
    path samplesheet
    val workflow_type
    val suffix
    path paths_file
    
    output:
    path "${samplesheet.baseName}_${suffix}.csv", emit: updated_samplesheet
    
    script:
    """
    # Log input files for debugging
    echo "Content of paths_file (${paths_file}):"
    cat ${paths_file} | head -n 10
    echo ""
    
    # Run the update script with explicit python call
    python3 ${workflow.projectDir}/templates/update_samplesheet.py "${samplesheet}" "${workflow_type}" "${suffix}" "${paths_file}" "${samplesheet.baseName}_${suffix}.csv"
    
    # Check the created file
    echo "Content of generated file (${samplesheet.baseName}_${suffix}.csv):"
    cat "${samplesheet.baseName}_${suffix}.csv" | head -n 10
    echo ""
    """
}


process check_resources {
    label 'process_low'
    
    output:
    path 'resource_check.txt'
    
    exec:
    def maxCpus = Runtime.runtime.availableProcessors()
    def maxMemory = Runtime.runtime.maxMemory()
    def maxMemoryGB = maxMemory / (1024 * 1024 * 1024)
    
    // Get process resource configs from workflow object
    def warnings = []
    def critical = []
    
    // Examine process configurations from the workflow context
    workflow.session.config.process.each { name, value ->
        if (name == 'cpus' && value > maxCpus) {
            critical << "Default process CPU setting (${value}) exceeds available CPUs (${maxCpus})"
        }
        
        if (value instanceof Map) {
            if (value.cpus && value.cpus instanceof Integer && value.cpus > maxCpus) {
                critical << "Process '${name}' CPU setting (${value.cpus}) exceeds available CPUs (${maxCpus})"
            }
            
            if (value.memory) {
                def memStr = value.memory.toString()
                if (memStr.endsWith('GB')) {
                    def mem = memStr.replace('GB', '').trim().toFloat()
                    if (mem > maxMemoryGB) {
                        critical << "Process '${name}' memory setting (${mem} GB) exceeds available memory (${String.format("%.2f", maxMemoryGB)} GB)"
                    }
                }
            }
        }
    }
    
    // Check withLabel process configurations
    workflow.session.config.process.each { name, value ->
        if (name == 'withLabel' && value instanceof Map) {
            value.each { label, settings ->
                if (settings.cpus && settings.cpus instanceof Integer && settings.cpus > maxCpus) {
                    critical << "Process with label '${label}' CPU setting (${settings.cpus}) exceeds available CPUs (${maxCpus})"
                }
                
                if (settings.memory) {
                    def memStr = settings.memory.toString()
                    if (memStr.endsWith('GB')) {
                        def mem = memStr.replace('GB', '').trim().toFloat()
                        if (mem > maxMemoryGB) {
                            critical << "Process with label '${label}' memory setting (${mem} GB) exceeds available memory (${String.format("%.2f", maxMemoryGB)} GB)"
                        }
                    }
                }
            }
        }
    }
    
    // Check withName process configurations
    workflow.session.config.process.each { name, value ->
        if (name == 'withName' && value instanceof Map) {
            value.each { procName, settings ->
                if (settings.cpus && settings.cpus instanceof Integer && settings.cpus > maxCpus) {
                    critical << "Process '${procName}' CPU setting (${settings.cpus}) exceeds available CPUs (${maxCpus})"
                }
                
                if (settings.memory) {
                    def memStr = settings.memory.toString()
                    if (memStr.endsWith('GB')) {
                        def mem = memStr.replace('GB', '').trim().toFloat()
                        if (mem > maxMemoryGB) {
                            critical << "Process '${procName}' memory setting (${mem} GB) exceeds available memory (${String.format("%.2f", maxMemoryGB)} GB)"
                        }
                    }
                }
            }
        }
    }
    
    // Display resource information only if explicitly requested or if there are warnings/errors
    // if (params.display_resource_check || !warnings.isEmpty() || !critical.isEmpty()) {
    //     log.info """
    // ===============================================================================================
    // System Resource Check
    // ===============================================================================================
    // Available CPUs: ${maxCpus}
    // Available Memory: ${String.format("%.2f", maxMemoryGB)} GB
    // ===============================================================================================
    // """
    // }
    
    // Display warnings and critical issues
    if (!warnings.isEmpty()) {
        log.warn """
    ===============================================================================================
    Resource Configuration Warnings
    ===============================================================================================
    Available CPUs: ${maxCpus}
    Available Memory: ${String.format("%.2f", maxMemoryGB)} GB
    -----------------------------------------------------------------------------------------------
    ${warnings.collect { "• ${it}" }.join('\n    ')}
    ===============================================================================================
    """
    }
    
    if (!critical.isEmpty()) {
        log.error """
    ===============================================================================================
    CRITICAL Resource Configuration Issues
    ===============================================================================================
    Available CPUs: ${maxCpus}
    Available Memory: ${String.format("%.2f", maxMemoryGB)} GB
    -----------------------------------------------------------------------------------------------
    ${critical.collect { "• ${it}" }.join('\n    ')}
    -----------------------------------------------------------------------------------------------
    Recommendations:
    1. Reduce process requirements in nextflow.config or your custom config
    2. Run on a machine with more resources
    3. Use a cluster/cloud executor instead of 'local'
    ===============================================================================================
    """
        if (params.enforce_resource_check) {
            error "Aborting due to insufficient resources. Set params.enforce_resource_check = false to continue anyway."
        }
    }
    
    def resourceFile = task.workDir.resolve('resource_check.txt')
    resourceFile.text = """
    CPUs: ${maxCpus}
    Memory: ${String.format("%.2f", maxMemoryGB)} GB
    """
}