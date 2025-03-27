#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_genomics_db_import_new {
    tag { "${interval}.rGDIN" }
    label 'gatk'
    label 'gatk_genomicsdb'
    publishDir "${params.db_path}", mode: 'copy', overwrite: true
    
    input:
    path(gvcf_list)
    each interval

    output:
    path("${interval}.gdb"), emit: interval_db

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx, but not more than 75% of available memory
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    def min_mem = Math.min(1, max_mem)
    
    """
    set -e  # Exit on error
    
    gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g" \
        GenomicsDBImport \
        --reference ${params.reference_files.ref} \
        --intervals ${interval} \
        --sample-name-map ${gvcf_list} \
        --batch-size 50 \
        --bypass-feature-reader \
        --consolidate \
        --genomicsdb-shared-posixfs-optimizations \
        --genomicsdb-workspace-path ${interval}.gdb
    """
}

process run_genomics_db_import_update {
    tag { "${interval}.rGDIU" }
    label 'gatk'
    label 'gatk_genomicsdb'
    publishDir "${params.db_path}", mode: 'copy', overwrite: true
    
    input:
    path(gvcf_list)
    each interval
    val(backup_status)

    output:
    path("${interval}.gdb"), emit: interval_db

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx, but not more than 75% of available memory
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    def min_mem = Math.min(1, max_mem)
    
    """
    set -e  # Exit on error
    
    # Create and setup temporary directories
    mkdir -p "\${PWD}/tmp"
    chmod 777 "\${PWD}/tmp"
    
    # Set environment variables for Java and native libraries
    export TMPDIR="\${PWD}/tmp"
    export _JAVA_OPTIONS="-Djava.io.tmpdir=\${PWD}/tmp"
    export LD_LIBRARY_PATH=/gatk/libs:\$LD_LIBRARY_PATH
    
    # Check if the workspace already exists
    if [ -d "${params.db_path}/${interval}.gdb" ] && [ -f "${params.db_path}/${interval}.gdb/callset.json" ]; then
        echo "Existing genomics database found. Performing an update operation."
        
        # Copy the existing workspace from the database path
        cp -r ${params.db_path}/${interval}.gdb .
        
        # Extract existing sample names from the callset.json file
        grep -o '"name":"[^"]*"' ${interval}.gdb/callset.json | cut -d':' -f2 | tr -d '"' > existing_samples.txt
        echo "Existing samples:"
        cat existing_samples.txt
        
        # Create a filtered sample map with only new samples
        awk 'NR==FNR{existing[\$1]=1; next} !existing[\$1]{print \$0}' existing_samples.txt ${gvcf_list} > new_samples.txt
        
        # Check if we have any new samples to add
        if [ ! -s new_samples.txt ]; then
            echo "No new samples to add. Workspace is already up to date."
            exit 0
        fi
        
        echo "New samples to add:"
        cat new_samples.txt
        
        # Create properly formatted sample map file for new samples
        awk '{printf "%s\\t%s\\n", \$1, \$2}' new_samples.txt > sample_map.txt
        echo "Sample map content:"
        cat sample_map.txt
        
        # Run GenomicsDBImport in update mode
        gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \\
            GenomicsDBImport \\
            --reference ${params.reference_files.ref} \\
            --intervals ${interval} \\
            --sample-name-map sample_map.txt \\
            --batch-size 50 \\
            --consolidate \\
            --reader-threads ${task.cpus} \\
            --genomicsdb-update-workspace-path ${interval}.gdb \\
            --genomicsdb-shared-posixfs-optimizations true \\
            --tmp-dir "\${PWD}/tmp"
            
    else
        echo "No existing genomics database found. Creating a new workspace."
        
        # Create properly formatted sample map file for all samples
        cp ${gvcf_list} sample_map.txt
        echo "Sample map content:"
        cat sample_map.txt
        
        # Run GenomicsDBImport in create mode
        gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \\
            GenomicsDBImport \\
            --reference ${params.reference_files.ref} \\
            --intervals ${interval} \\
            --sample-name-map sample_map.txt \\
            --batch-size 50 \\
            --bypass-feature-reader \\
            --consolidate \\
            --reader-threads ${task.cpus} \\
            --genomicsdb-workspace-path ${interval}.gdb \\
            --genomicsdb-shared-posixfs-optimizations true \\
            --tmp-dir "\${PWD}/tmp"
    fi
    """
}

process run_backup_genomic_db {
    tag { "rBGDB" }
    label 'small_mem'

    output:
    val(true), emit: backup_status

    script:
    """
    set -e  # Exit on error
    
    # Create backup directory with timestamp
    backup_dir="genomicsdb_backup/${params.db_path}_backup_\$(date +%Y%m%d_%H%M%S)"
    mkdir -p "\${backup_dir}"
    
    # Check if source directory exists
    if [ -d "${params.db_path}" ]; then
        # Check if directory has content
        if [ "\$(ls -A ${params.db_path} 2>/dev/null)" ]; then
            echo "Backing up existing data from ${params.db_path}"
            # Copy all files from genomics DB to backup
            cp -r ${params.db_path}/* "\${backup_dir}/"
            
            # Verify backup was successful
            if [ \$? -ne 0 ]; then
                echo "Error: Backup failed"
                exit 1
            fi
            
            echo "Backup completed successfully to \${backup_dir}"
        else
            echo "Source directory ${params.db_path} exists but is empty. No backup needed."
        fi
    else
        echo "Source directory ${params.db_path} does not exist yet. Creating it."
        mkdir -p "${params.db_path}"
    fi
    
    # Always succeed if we get to this point
    touch "\${backup_dir}/backup_manifest.txt"
    echo "Backup process completed at \$(date)" > "\${backup_dir}/backup_manifest.txt"
    """
}

process run_genomics_db_import {
    tag { "import-${interval}" }
    label 'gatk'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/genomicsdb", mode: 'copy',
        pattern: "*.tar"

    input:
    path(gvcf_list)
    each interval

    output:
    path("${interval}.gdb"), emit: interval_db

    script:
    // Calculate memory safely - ensure at least 1GB for Xmx, but not more than 75% of available memory
    def max_mem = Math.max(1, Math.min(task.memory.toGiga()-1, (task.memory.toGiga() * 0.75).intValue()))
    def min_mem = Math.min(1, max_mem)
    
    """
    set -e  # Exit on error
    
    gatk --java-options "-XX:+UseSerialGC -Xms${min_mem}g -Xmx${max_mem}g" \
        GenomicsDBImport \
        --reference ${params.reference_files.ref} \
        --intervals ${interval} \
        --sample-name-map ${gvcf_list} \
        --batch-size 50 \
        --bypass-feature-reader \
        --consolidate \
        --genomicsdb-shared-posixfs-optimizations \
        --genomicsdb-workspace-path ${interval}.gdb
    """
}
