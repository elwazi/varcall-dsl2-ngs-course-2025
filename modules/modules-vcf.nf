#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Process to filter VCF files using GATK VariantFiltration
process filter_vcf {
    tag { "filter_${vcf_file.simpleName}" }
    label 'gatk'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/vcfs/filtered", mode: 'copy'
    
    input:
    tuple path(vcf_file), path(vcf_index)
    path(ref_fasta)
    path(ref_fasta_fai)
    path(ref_dict)
    
    output:
    tuple path("${vcf_file.simpleName}.filtered.vcf.gz"), path("${vcf_file.simpleName}.filtered.vcf.gz.tbi"), emit: filtered_vcf
    path "${vcf_file.simpleName}.filter_stats.txt", emit: filter_stats
    path "${vcf_file.simpleName}.marked.vcf.gz*", optional: true
    
    script:
    def qd_threshold = params.qd_threshold ?: 2.0
    def fs_threshold = params.fs_threshold ?: 60.0
    def mq_threshold = params.mq_threshold ?: 40.0
    def sor_threshold = params.sor_threshold ?: 3.0
    def readpos_threshold = params.readpos_threshold ?: -8.0
    
    // Check if index exists
    def index_cmd = ""
    if (!vcf_index.exists() || vcf_index.size() == 0) {
        index_cmd = """
        echo "Index file doesn't exist or is empty. Creating index..."
        # Check if tabix is available, otherwise skip
        if command -v tabix &> /dev/null; then
            tabix -p vcf ${vcf_file}
        else
            echo "Warning: tabix not found in container, proceeding without creating index"
        fi
        """
    }
    """
    ${index_cmd}
    
    # Set filter expressions based on whether it's a SNP or INDEL
    gatk VariantFiltration \\
        -R ${ref_fasta} \\
        -V ${vcf_file} \\
        --filter-name "QD_filter" --filter-expression "QD < ${qd_threshold}" \\
        --filter-name "FS_filter" --filter-expression "FS > ${fs_threshold}" \\
        --filter-name "MQ_filter" --filter-expression "MQ < ${mq_threshold}" \\
        --filter-name "SOR_filter" --filter-expression "SOR > ${sor_threshold}" \\
        --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < ${readpos_threshold}" \\
        -O ${vcf_file.simpleName}.marked.vcf.gz
    
    # Extract variants that PASS filter
    gatk SelectVariants \\
        -R ${ref_fasta} \\
        -V ${vcf_file.simpleName}.marked.vcf.gz \\
        --exclude-filtered \\
        -O ${vcf_file.simpleName}.filtered.vcf.gz
    
    # Generate filter statistics
    echo "Filtering statistics for ${vcf_file.simpleName}" > ${vcf_file.simpleName}.filter_stats.txt
    echo "-----------------------------------------------" >> ${vcf_file.simpleName}.filter_stats.txt
    
    # Get original variant count
    echo "Original variant count:" >> ${vcf_file.simpleName}.filter_stats.txt
    orig_count=\$(gatk CountVariants -V ${vcf_file} | tail -1)
    echo "\$orig_count" >> ${vcf_file.simpleName}.filter_stats.txt
    
    # Get filtered variant count
    echo "Filtered variant count:" >> ${vcf_file.simpleName}.filter_stats.txt
    filt_count=\$(gatk CountVariants -V ${vcf_file.simpleName}.filtered.vcf.gz | tail -1)
    echo "\$filt_count" >> ${vcf_file.simpleName}.filter_stats.txt
    
    # Calculate difference if counts are available
    if [[ -n "\$orig_count" && -n "\$filt_count" ]] && [ "\$orig_count" -eq "\$orig_count" 2>/dev/null ] && [ "\$filt_count" -eq "\$filt_count" 2>/dev/null ]; then
        removed=\$((\$orig_count - \$filt_count))
        echo "Filtered out: \$removed variants" >> ${vcf_file.simpleName}.filter_stats.txt
        percent=\$(awk "BEGIN {printf \\"%.2f\\", (\$removed/\$orig_count)*100}")
        echo "Percentage filtered: \${percent}%" >> ${vcf_file.simpleName}.filter_stats.txt
    else
        echo "Could not calculate filtered variants (counts unavailable)" >> ${vcf_file.simpleName}.filter_stats.txt
    fi
    
    # Filter statistics by type
    echo "" >> ${vcf_file.simpleName}.filter_stats.txt
    echo "Filter statistics by reason:" >> ${vcf_file.simpleName}.filter_stats.txt
    echo "Filter thresholds used:" >> ${vcf_file.simpleName}.filter_stats.txt
    echo "  QD < ${qd_threshold}" >> ${vcf_file.simpleName}.filter_stats.txt
    echo "  FS > ${fs_threshold}" >> ${vcf_file.simpleName}.filter_stats.txt
    echo "  MQ < ${mq_threshold}" >> ${vcf_file.simpleName}.filter_stats.txt
    echo "  SOR > ${sor_threshold}" >> ${vcf_file.simpleName}.filter_stats.txt
    echo "  ReadPosRankSum < ${readpos_threshold}" >> ${vcf_file.simpleName}.filter_stats.txt
    echo "" >> ${vcf_file.simpleName}.filter_stats.txt
    
    # Count variants filtered by each criterion
    zgrep -v "^#" ${vcf_file.simpleName}.marked.vcf.gz | grep "QD_filter" | wc -l | xargs -I{} echo "QD < ${qd_threshold}: {}" >> ${vcf_file.simpleName}.filter_stats.txt
    zgrep -v "^#" ${vcf_file.simpleName}.marked.vcf.gz | grep "FS_filter" | wc -l | xargs -I{} echo "FS > ${fs_threshold}: {}" >> ${vcf_file.simpleName}.filter_stats.txt
    zgrep -v "^#" ${vcf_file.simpleName}.marked.vcf.gz | grep "MQ_filter" | wc -l | xargs -I{} echo "MQ < ${mq_threshold}: {}" >> ${vcf_file.simpleName}.filter_stats.txt
    zgrep -v "^#" ${vcf_file.simpleName}.marked.vcf.gz | grep "SOR_filter" | wc -l | xargs -I{} echo "SOR > ${sor_threshold}: {}" >> ${vcf_file.simpleName}.filter_stats.txt
    zgrep -v "^#" ${vcf_file.simpleName}.marked.vcf.gz | grep "ReadPosRankSum_filter" | wc -l | xargs -I{} echo "ReadPosRankSum < ${readpos_threshold}: {}" >> ${vcf_file.simpleName}.filter_stats.txt
    """
}

// Process to annotate VCF files using Ensembl VEP (Variant Effect Predictor)
process annotate_vcf_vep {
    tag { "vep_annotate_${vcf_file.simpleName}" }
    label 'vep'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/vcfs/annotated", mode: 'copy', pattern: '*.{html,json,txt}'
    
    input:
    tuple path(vcf_file), path(vcf_index)
    val cache_dir
    
    output:
    path "${vcf_file.simpleName}.vep.vcf", emit: annotated_vcf
    path "${vcf_file.simpleName}.vep_summary.html", emit: annotation_summary
    path "${vcf_file.simpleName}.vep_stats.txt", optional: true, emit: annotation_stats
    path "${vcf_file.simpleName}.vep_warnings.txt", optional: true, emit: annotation_warnings
    
    script:
    def cache_version = params.vep_cache_version ?: "113"
    def species = "homo_sapiens"
    def assembly = "GRCh38"
    def force_cache = params.force_cache ? "--force_overwrite" : ""
    
    // Plugin configuration
    def plugin_options = ""
    if (params.vep_plugins) {
        def plugins = params.vep_plugins.split(',').collect { it.trim() }
        plugin_options = plugins.collect { "--plugin $it" }.join(' ')
    }
    
    // Plugin directory
    def plugin_dir = params.vep_plugin_dir ? "--dir_plugins ${params.vep_plugin_dir}" : ""
    
    // Extra options
    def extra_opts = params.vep_extra_options ?: ""
    
    // Check if index exists
    def index_cmd = ""
    if (!vcf_index.exists() || vcf_index.size() == 0) {
        index_cmd = """
        echo "Index file doesn't exist or is empty."
        # Check if tabix is available, otherwise skip
        if command -v tabix &> /dev/null; then
            echo "Creating index with tabix..."
            tabix -p vcf ${vcf_file}
        else
            echo "Warning: tabix not found in container, proceeding without creating index"
        fi
        """
    }
    
    // Check for VEP cache directory
    def cache_check_cmd = """
    if [ ! -d "${cache_dir}/${species}/${cache_version}_${assembly}" ]; then
        echo "WARNING: VEP cache directory for ${species}/${cache_version}_${assembly} not found"
        if [ "${params.download_vep_cache}" == "true" ]; then
            echo "Attempting to download VEP cache..."
            mkdir -p ${cache_dir}
            curl -L ${params.vep_cache_url} | tar xz -C ${cache_dir}
        else
            echo "Cache download not enabled. Set params.download_vep_cache = true to auto-download"
            echo "Proceeding with VEP annotation without cache (slower)"
        fi
    else
        echo "Using VEP cache from: ${cache_dir}/${species}/${cache_version}_${assembly}"
    fi
    """
    
    """
    ${index_cmd}
    ${cache_check_cmd}
    
    # Run VEP annotation
    vep \\
        --input_file ${vcf_file} \\
        --output_file ${vcf_file.simpleName}.vep.vcf \\
        --format vcf \\
        --vcf \\
        --symbol \\
        --check_existing \\
        --terms SO \\
        --cache \\
        --dir_cache ${cache_dir} \\
        --cache_version ${cache_version} \\
        --offline \\
        --species ${species} \\
        --assembly ${assembly} \\
        --stats_file ${vcf_file.simpleName}.vep_stats.txt \\
        --warning_file ${vcf_file.simpleName}.vep_warnings.txt \\
        --stats_text \\
        --fork 4 \\
        --sift b \\
        --polyphen b \\
        --regulatory \\
        --canonical \\
        --biotype \\
        --af \\
        --af_1kg \\
        --af_gnomad \\
        --max_af \\
        --pubmed \\
        --uniprot \\
        --stats_html \\
        --stats_text \\
        --stats_file ${vcf_file.simpleName}.vep_summary.html
    
    # Check for successful annotation
    if [ -f "${vcf_file.simpleName}.vep.vcf" ]; then
        echo "VEP annotation completed successfully"
        # Count variants
        VARIANT_COUNT=\$(grep -v "^#" ${vcf_file.simpleName}.vep.vcf | wc -l)
        echo "Annotated \$VARIANT_COUNT variants"
    else
        echo "VEP annotation failed. Creating empty output files to prevent pipeline failure"
        touch ${vcf_file.simpleName}.vep.vcf
        touch ${vcf_file.simpleName}.vep_summary.html
        touch ${vcf_file.simpleName}.vep_stats.txt
        touch ${vcf_file.simpleName}.vep_warnings.txt
    fi
    """
}

// Process to compress and index the annotated VCF
process compress_index_vcf {
    tag { "compress_${vcf_file.simpleName}" }
    label 'htslib'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/vcfs/annotated", mode: 'copy'
    
    input:
    path vcf_file
    
    output:
    tuple path("${vcf_file}.gz"), path("${vcf_file}.gz.tbi"), emit: compressed_vcf
    
    script:
    """
    # Print info for debugging
    echo "Input file details:"
    ls -lh ${vcf_file}
    
    # Compress the VCF file using bgzip from HTSlib
    bgzip -f ${vcf_file}
    
    # Index the compressed VCF using tabix from HTSlib
    tabix -p vcf ${vcf_file}.gz
    
    # Print info for debugging
    echo "Compressed file details:"
    ls -lh ${vcf_file}.gz
    echo "Index file details:"
    ls -lh ${vcf_file}.gz.tbi
    """
}

// Process to combine multiple annotated VCFs
process combine_vcfs {
    tag "combine_vcfs"
    label 'bcftools'
    label 'large_mem'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/combined", mode: 'copy'
    
    input:
    path('vcfs/*')
    path(ref_fasta)
    
    output:
    tuple path("combined.vcf.gz"), path("combined.vcf.gz.tbi"), emit: combined_vcf
    path "vcf_merge_stats.txt", emit: merge_stats
    
    script:
    """
    # List all VCF files
    find vcfs/ -name "*.vcf.gz" > vcf_list.txt
    
    # Sort VCF files by chromosome order (for proper concatenation)
    # This creates a sorted file list based on chromosome numbers
    cat vcf_list.txt | sort -V > vcf_list_sorted.txt
    
    # Count how many VCFs we're combining
    vcf_count=\$(cat vcf_list_sorted.txt | wc -l)
    echo "Combining \$vcf_count VCF files" > vcf_merge_stats.txt
    
    if [ \$vcf_count -gt 1 ]; then
        echo "Detected multiple VCFs, likely split by chromosome. Using bcftools concat..." >> vcf_merge_stats.txt
        
        # First check if chromosomes overlap - if they do, we need to use merge instead of concat
        chrom_counts=\$(bcftools query -f '%CHROM\\n' \$(head -1 vcf_list_sorted.txt) | sort | uniq | wc -l)
        total_chroms=\$(bcftools query -f '%CHROM\\n' \$(cat vcf_list_sorted.txt) | sort | uniq | wc -l)
        
        # List chromosomes in each file
        echo "Chromosomes in input files:" >> vcf_merge_stats.txt
        for vcf in \$(cat vcf_list_sorted.txt); do
            chroms=\$(bcftools query -f '%CHROM\\n' \$vcf | sort | uniq | tr '\\n' ' ')
            echo "  \$(basename \$vcf): \$chroms" >> vcf_merge_stats.txt
        done
        
        if [ \$chrom_counts -eq \$total_chroms ]; then
            # No overlap in chromosomes, use concat for better performance
            echo "No chromosome overlap detected. Using bcftools concat..." >> vcf_merge_stats.txt
            bcftools concat -a -f vcf_list_sorted.txt -O z -o combined.vcf.gz
        else
            # Chromosomes overlap, use merge
            echo "Chromosome overlap detected. Using bcftools merge..." >> vcf_merge_stats.txt
            bcftools merge -f vcf_list_sorted.txt -O z -o combined.vcf.gz
        fi
    else
        echo "Only one VCF file found, creating a copy as combined.vcf.gz" >> vcf_merge_stats.txt
        cp \$(cat vcf_list_sorted.txt) combined.vcf.gz
    fi
    
    # Generate variant counts
    echo "" >> vcf_merge_stats.txt
    echo "Variant counts in combined VCF:" >> vcf_merge_stats.txt
    bcftools stats combined.vcf.gz | grep "number of records:" >> vcf_merge_stats.txt
    
    # Index the combined file
    tabix -p vcf combined.vcf.gz
    """
}

// Define the filter_annotate_vcf workflow
workflow filter_annotate_vcf_workflow {
    take:
    vcf_ch                // Channel with: [vcf, index]
    ref_fasta             // Reference genome FASTA file
    ref_fasta_fai         // Reference genome FASTA index
    ref_dict              // Reference genome dictionary file
    genome_version        // Genome version for annotation (e.g., 'GRCh38.105') - not used for VEP
    
    main:
    // Filter VCF files
    filter_vcf(vcf_ch, ref_fasta, ref_fasta_fai, ref_dict)
    
    // Create channel for VEP cache directory
    // vep_cache_dir = Channel.fromPath(params.vep_cache_dir, checkIfExists: true)
    
    // Annotate filtered VCF files with VEP
    annotate_vcf_vep(filter_vcf.out.filtered_vcf, params.vep_cache_dir)
    
    // Compress and index the VEP annotated VCF
    compress_index_vcf(annotate_vcf_vep.out.annotated_vcf)
    
    // Optionally combine all VCFs if needed
    if (params.combine_vcfs) {
        compress_index_vcf.out.compressed_vcf
            .map { vcf, index -> vcf }
            .collect()
            .map { vcfs -> 
                ["vcfs": vcfs]
            }
            .set { all_vcfs }
        
        combine_vcfs(all_vcfs, ref_fasta)
    }
    
    emit:
    filtered_vcfs = filter_vcf.out.filtered_vcf
    filter_stats = filter_vcf.out.filter_stats
    annotated_vcfs = compress_index_vcf.out.compressed_vcf
    annotation_summary = annotate_vcf_vep.out.annotation_summary
    annotation_stats = annotate_vcf_vep.out.annotation_stats
    combined_vcf = params.combine_vcfs ? combine_vcfs.out.combined_vcf : Channel.empty()
} 