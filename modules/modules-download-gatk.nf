#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process download_reference {
    tag { "dDL_full" }
    label 'small_mem'
    errorStrategy { task.exitStatus == 1 ? 'retry' : 'terminate' }
    maxRetries 3
    cache 'deep'
    storeDir "${params.ref_dir}"

    output:
    path("*.fasta.gz"), optional: true, emit: compressed_ref
    path("*.fasta"), emit: full_ref

    script:
    def resources = params.build == 'b37' ? params.b37_resources : params.b38_resources
    def ref_name = params.build == 'b37' ? 'human_g1k_v37_decoy' : 'Homo_sapiens_assembly38'
    """
    # Function to handle download with retries
    function download_with_retry() {
        local url=\$1
        local output=\$2
        local max_attempts=5
        local attempt=1
        local wait_time=10
        
        # Validate URL
        if [[ -z "\$url" || "\$url" == "null" ]]; then
            echo "Error: Empty or null URL provided"
            return 1
        fi
        
        echo "Downloading reference from: \$url"
        while [ \$attempt -le \$max_attempts ]; do
            echo "Download attempt \$attempt of \$max_attempts"
            if wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 --tries=5 "\$url" -O "\$output"; then
                if [ ! -s "\$output" ]; then
                    echo "Error: Downloaded file is empty"
                    rm -f "\$output"
                    return 1
                fi
                if grep -q '<?xml' "\$output" || grep -q '<html' "\$output"; then
                    echo "Error: Received HTML/XML instead of FASTA file"
                    rm -f "\$output"
                    return 1
                fi
                return 0
            fi
            attempt=\$((attempt + 1))
            sleep \$wait_time
            wait_time=\$((wait_time * 2))
        done
        return 1
    }

    # Download full reference genome
    if [[ "${resources.ref}" == *.gz ]]; then
        echo "Downloading compressed reference genome"
        if ! download_with_retry "${resources.ref}" "${ref_name}.fasta.gz"; then
            echo "Failed to download reference after multiple attempts"
            exit 1
        fi
        
        # Validate the downloaded file
        if ! gunzip -t "${ref_name}.fasta.gz"; then
            echo "Error: Downloaded file is not a valid gzip file"
            exit 1
        fi
        
        echo "Decompressing reference genome"
        if ! gunzip -c "${ref_name}.fasta.gz" > "${ref_name}.fasta"; then
            echo "Error: Failed to decompress reference genome"
            exit 1
        fi
        
        # Verify the uncompressed file exists and has content
        if [ ! -s "${ref_name}.fasta" ]; then
            echo "Error: Uncompressed reference file is empty"
            exit 1
        fi
        
        echo "Successfully downloaded and processed reference genome"
    else
        echo "Downloading uncompressed reference genome"
        if ! download_with_retry "${resources.ref}" "${ref_name}.fasta"; then
            echo "Failed to download reference after multiple attempts"
            exit 1
        fi
        
        # Verify the file exists and has content
        if [ ! -s "${ref_name}.fasta" ]; then
            echo "Error: Downloaded reference file is empty"
            exit 1
        fi
        
        echo "Successfully downloaded reference genome"
        
        # Don't compress the reference - per user request
        # bgzip -f ${ref_name}.fasta
    fi
    """
}

process samtools_index {
    tag { "dIDX_full" }
    label 'small_mem'
    label 'samtools'
    cache 'deep'
    storeDir "${params.ref_dir}"

    input:
    path(fasta)

    output:
    tuple path(fasta), path("*.fai"), emit: indexed_ref

    script:
    """
    samtools faidx ${fasta}
    """
}

process bwa_index {
    tag { "dBWA_full" }
    label 'small_mem'
    label 'bwa'
    cache 'deep'
    storeDir "${params.ref_dir}"

    input:
    tuple path(fasta), path(fai)

    output:
    path("*.{amb,ann,bwt,pac,sa}"), emit: bwa_indices

    script:
    """
    # Create BWA index
    bwa index ${fasta}
    """
}

process samtools_dict {
    tag { "dSD_full" }
    label 'small_mem'
    label 'samtools'
    cache 'deep'
    storeDir "${params.ref_dir}"

    input:
    path(fasta)

    output:
    path("*.samtools.dict"), emit: dict_file

    script:
    def ref_name = params.build == 'b37' ? 'human_g1k_v37_decoy' : 'Homo_sapiens_assembly38'
    """
    samtools dict ${fasta} > ${ref_name}.samtools.dict
    """
}

process gatk_dict {
    tag { "dSD_full" }
    label 'small_mem'
    label 'gatk'
    cache 'deep'
    storeDir "${params.ref_dir}"

    input:
    path(fasta)

    output:
    path("*.dict"), emit: dict_file

    script:
    def ref_name = params.build == 'b37' ? 'human_g1k_v37_decoy' : 'Homo_sapiens_assembly38'
    """
    gatk CreateSequenceDictionary -R ${fasta} -O ${ref_name}.dict
    ln -sf ${ref_name}.dict ${fasta}.dict
    """
}

process download_full_vcf {
    tag { "dl_vcf" }
    label 'small_mem'
    publishDir "${params.ref_dir}", mode: 'copy', overwrite: true, failOnError: false
    cache 'deep'
    storeDir "${params.ref_dir}"
    errorStrategy { task.exitStatus == 1 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(res_name), val(res_url)
    
    output:
    tuple val(res_name), path("${res_name}__*.vcf.gz"), emit: vcf_file
    
    script:
    if (res_url == null || res_url.toString().trim().isEmpty()) {
        error "Error: Empty or null URL provided for ${res_name}"
    }
    
    """
    #!/bin/bash

    # Store input variables
    RESOURCE_NAME="${res_name}"
    RESOURCE_URL="${res_url}"
    
    # Extract original filename from URL
    ORIGINAL_FILENAME=\$(basename "\$RESOURCE_URL")
    
    # Remove .gz extension if present
    if [[ "\$ORIGINAL_FILENAME" == *.gz ]]; then
        BASE_NAME=\${ORIGINAL_FILENAME%.gz}
    else
        BASE_NAME="\$ORIGINAL_FILENAME"
    fi
    
    # Remove .vcf extension if present
    if [[ "\$BASE_NAME" == *.vcf ]]; then
        BASE_NAME=\${BASE_NAME%.vcf}
    fi
    
    # Create output filename
    OUTPUT_FILE="\${RESOURCE_NAME}__\${BASE_NAME}.vcf.gz"
    
    echo "Processing resource: \$RESOURCE_NAME"
    echo "URL: \$RESOURCE_URL"
    echo "Output file: \$OUTPUT_FILE"
    
    # Function to handle download with retries
    function download_with_retry() {
        local url=\$1
        local output=\$2
        local max_attempts=5
        local attempt=1
        local wait_time=10
        
        # Validate URL before attempting download
        if [[ -z "\$url" || "\$url" == "null" ]]; then
            echo "Error: Empty or null URL provided"
            return 1
        fi
        
        while [ \$attempt -le \$max_attempts ]; do
            echo "Download attempt \$attempt of \$max_attempts for \$url to \$output"
            if wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 --tries=5 "\$url" -O "\$output"; then
                # Verify file was downloaded properly
                if [ ! -s "\$output" ]; then
                    echo "Error: Downloaded file is empty"
                    rm -f "\$output"
                    return 1
                fi
                # Check if the file is an HTML/XML error page
                if grep -q '<?xml' "\$output" || grep -q '<html' "\$output"; then
                    echo "Error: Received HTML/XML instead of VCF file"
                    rm -f "\$output"
                    return 1
                fi
                return 0
            fi
            attempt=\$((attempt + 1))
            sleep \$wait_time
            wait_time=\$((wait_time * 2))
        done
        return 1
    }
    
    # Download VCF file
    if ! download_with_retry "\$RESOURCE_URL" "temp_download.vcf"; then
        echo "Failed to download VCF file after multiple attempts: \$RESOURCE_URL"
        exit 1
    fi
    
    # Ensure proper BGZF compression for all files
    if [[ "\$RESOURCE_URL" == *.gz ]]; then
        # For gzipped files, decompress and recompress with bgzip
        gunzip -c "temp_download.vcf" | bgzip > "\$OUTPUT_FILE"
        rm "temp_download.vcf"
    else
        # For uncompressed files, compress with bgzip
        bgzip -c "temp_download.vcf" > "\$OUTPUT_FILE"
        rm "temp_download.vcf"
    fi
    """
}

process index_vcf {
    tag { "idx_vcf" }
    label 'small_mem'
    publishDir "${params.ref_dir}", mode: 'copy', overwrite: true, failOnError: false
    cache 'deep'
    storeDir "${params.ref_dir}"
    
    input:
    tuple val(res_name), path(vcf_file)
    
    output:
    tuple val(res_name), path(vcf_file), path("${vcf_file}.tbi"), emit: indexed_vcf
    
    script:
    """
    #!/bin/bash
    
    # Store input variables
    RESOURCE_NAME="${res_name}"
    VCF_FILE="${vcf_file}"
    
    echo "Indexing VCF file: \$VCF_FILE for resource: \$RESOURCE_NAME"
    
    # Index VCF file
    if [[ "\$RESOURCE_NAME" == "dbsnp" ]]; then
        # For dbsnp, ensure it's properly compressed first
        echo "Special handling for dbSNP resource"
        mv "\$VCF_FILE" temp.vcf.gz
        gunzip -c temp.vcf.gz | bgzip > "\$VCF_FILE"
        rm temp.vcf.gz
    fi
    
    if ! tabix -p vcf "\$VCF_FILE"; then
        echo "Failed to index VCF file"
        exit 1
    fi
    
    echo "Successfully indexed VCF file: \$VCF_FILE"
    """
}

workflow download_and_index_reference_workflow {
    main:
    // Download reference
    download_reference()

    // Create samtools index
    samtools_index(download_reference.out.full_ref)

    // Create BWA index
    bwa_index(samtools_index.out.indexed_ref)

    // Create samtools dictionary
    samtools_dict(download_reference.out.full_ref)

    // Create GATK dictionary
    gatk_dict(download_reference.out.full_ref)

    emit:
    done = true
}

workflow download_and_index_vcfs {
    main:
    def resources = params.build == 'b37' ? params.b37_resources : params.b38_resources
    
    // Create list of resources to download (only include those with non-null URLs)
    def vcf_configs = [
        [name: 'dbsnp', url: resources.dbsnp],
        [name: 'hapmap', url: resources.hapmap],
        [name: 'omni', url: resources.omni],
        [name: 'phase1_snps', url: resources.phase1_snps],
        [name: 'golden_indels', url: resources.golden_indels],
        [name: 'known_indels', url: resources.known_indels]
    ]
    
    // For mills resource, only add if it exists in the configuration
    if (resources.containsKey('mills')) {
        vcf_configs.add([name: 'mills', url: resources.mills])
    } else {
        log.warn "Note: 'mills' resource URL not found in configuration. Skipping download."
    }
    
    // Filter out any null URLs and create the download channel
    def vcf_list = vcf_configs.findAll { it.url != null }.collect { [it.name, it.url] }
    
    if (vcf_list.size() == 0) {
        log.warn "Warning: No valid VCF URLs found in configuration. Skipping VCF downloads."
        return
    }
    
    // Log resources being downloaded
    log.info "Downloading VCF resources: ${vcf_list.collect { it[0] }.join(', ')}"
    
    // Create the channel for VCF downloads
    def vcf_ch = Channel.fromList(vcf_list)

    // Download VCFs
    download_full_vcf(vcf_ch)
    
    // Index VCFs
    index_vcf(download_full_vcf.out.vcf_file)
    vcf_files = index_vcf.out.indexed_vcf.collect()

    emit:
    vcf_files
}

workflow {
    if (params.workflow == 'download_gatk') {
        // Download full reference and VCFs
        download_and_index_reference_workflow()
        download_and_index_vcfs()
    }
} 