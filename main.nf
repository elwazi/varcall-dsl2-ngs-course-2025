#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { print_sample_info; log_tool_version_samtools; log_tool_version_bwa; log_tool_version_gatk; check_files; update_samplesheet; check_resources } from './modules/modules-general.nf'
include { validate_fastq; VALIDATE_FASTQ; fastq_qc_workflow; bam_qc_workflow; vcf_qc_workflow } from './modules/modules-qc.nf'
include { run_bwa; run_bam_sort; run_mark_duplicates; run_create_recalibration_table; run_recalibrate_bam; bam_to_cram; run_cram_flagstat; create_cram_md5sum; } from './modules/modules-align.nf'
include { run_haplotype_caller_auto as run_haplotype_caller_auto_nosex; run_haplotype_caller_auto as run_haplotype_caller_auto_males; run_haplotype_caller_auto as run_haplotype_caller_auto_females; run_haplotype_caller_males; run_haplotype_caller_females; run_haplotype_caller_mt as run_haplotype_caller_mt_nosex; run_haplotype_caller_mt as run_haplotype_caller_mt_males; run_haplotype_caller_mt as run_haplotype_caller_mt_females; run_sort_male_gvcfs; run_combine_sample_gvcfs as run_combine_sample_gvcfs_nosex; run_combine_sample_gvcfs as run_combine_sample_gvcfs_males; run_combine_sample_gvcfs as run_combine_sample_gvcfs_females; run_create_gvcf_md5sum as run_create_gvcf_md5sum_nosex; run_create_gvcf_md5sum as run_create_gvcf_md5sum_males; run_create_gvcf_md5sum as run_create_gvcf_md5sum_females; run_index_gvcf as run_index_gvcf_nosex; run_index_gvcf as run_index_gvcf_males; run_index_gvcf as run_index_gvcf_females; run_validate_gvcf } from './modules/modules-generate-gvcf.nf'
include { run_combine_gvcfs; run_concat_gvcfs } from './modules/modules-combine-gvcfs.nf'
include { run_genomics_db_import_new; run_backup_genomic_db; run_genomics_db_import_update } from './modules/modules-genomics-db-import.nf'
include { run_genotype_gvcf_on_genome_db; run_concat_vcf; run_vqsr_on_snps; apply_vqsr_on_snps; hard_filter_indels } from './modules/modules-genome-calling.nf'
include { download_and_index_reference_workflow; download_and_index_vcfs } from './modules/modules-download-gatk.nf'
include { VALIDATE_GVCF } from './modules/modules-validate-gvcf.nf'
include { filter_annotate_vcf_workflow } from './modules/modules-vcf'


// Helper functions
def initializeWorkflowParams() {
    // Set chromosomes from config
    def chroms_par = params.chroms_par
    def chroms_auto = params.build == "b37" ? params.b37_chroms_auto : params.b38_chroms_auto
    def chroms_all = params.build == "b37" ? params.b37_chroms_all : params.b38_chroms_all

    if(params.test) {
        chroms_auto = chroms_all
        chroms_par = []
    }

    return [chroms_par, chroms_auto, chroms_all]
}

def validateWorkflow(workflow_name) {
    def available_workflows = [
        'download_gatk': 'Download and prepare GATK resources',
        'validate_fastq': 'Validate FASTQ files (check format, compression, and read counts)',
        'fastq_qc': 'Run QC on FASTQ files (FastQC, Trimmomatic, MultiQC)',
        'bam_qc': 'Run QC on BAM files (samtools stats, flagstat, idxstats, MultiQC)',
        'validate-gvcf': 'Validate gVCF files and generate validation reports',
        'vcf_qc': 'Run QC on VCF files (bcftools stats, MultiQC)',
        'align': 'Align FASTQ files to reference genome (BWA-MEM, GATK)',
        'generate-gvcfs': 'Generate GVCFs from BAM files (GATK HaplotypeCaller)',
        'combine-gvcfs': 'Combine GVCFs from multiple samples',
        'genomics-db-import': 'Import GVCFs into GenomicsDB',
        'genome-calling': 'Joint genotyping from GenomicsDB',
        'filter-vcf': 'Filter and annotate VCF files'
    ]

    if (!workflow_name) {
        log.error """
============================================================
ERROR: No workflow specified!
============================================================
Please specify a workflow using --workflow WORKFLOW_NAME

Available workflows:
${available_workflows.collect { name, desc -> "  - ${name}: ${desc}" }.join('\n')}

Example usage:
    nextflow run main.nf --workflow fastq_qc [other parameters]
============================================================
        """.stripIndent()
        System.exit(1)
    }

    if (!available_workflows.containsKey(workflow_name)) {
        log.error """
        ============================================================
        ERROR: Unknown workflow '${workflow_name}'
        ============================================================
        Available workflows:
        ${available_workflows.collect { name, desc -> "  - ${name}: ${desc}" }.join('\n')}
        
        Example usage:
          nextflow run main.nf --workflow fastq_qc [other parameters]
        ============================================================
        """.stripIndent()
        System.exit(1)
    }

    return workflow_name
}

def checkRequiredFile(file_path, file_type, workflow_name) {
    if (!file_path) {
        log.error """
        ============================================================
        ERROR: ${file_type} is required for ${workflow_name}
        ============================================================
        Please provide the required file path.
        """.stripIndent()
        System.exit(1)
    }
    if (!file(file_path).exists()) {
        log.error """
        ============================================================
        ERROR: ${file_type} not found for ${workflow_name}
        ============================================================
        Could not find file: ${file_path}
        Please check that the file exists and the path is correct.
        """.stripIndent()
        System.exit(1)
    }
}


// Check if the reference file exists
def check_reference_file(reference_file) {
    def reference_files = params.reference_files
    if (!reference_files.containsKey(reference_file)) {
        log.error """
        ============================================================
        ERROR: Reference file not found: ${reference_file}
        ============================================================
        Please check that the file exists and the path is correct.
        """
        System.exit(1)
    }
}

// Add this function after the imports at the top of the file
def checkReferenceFilesExist() {
    // Check reference file and index files
    def ref_path = params.reference_files.ref
    log.info "Checking reference files for: ${ref_path}"
    
    // Check main reference file
    def ref_file = file(ref_path)
    if (!ref_file.exists()) {
        log.error """
        ============================================================
        ERROR: Reference file not found: ${ref_path}
        ============================================================
        Please check that the file exists and the path is correct.
        """
        System.exit(1)
    }
    
    // Check if BWA index files exist
    def index_extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
    def missing_indexes = []
    
    index_extensions.each { ext ->
        def idx_file = file("${ref_path}.${ext}")
        if (!idx_file.exists()) {
            missing_indexes << "${ref_path}.${ext}"
        }
    }
    
    if (missing_indexes.size() > 0) {
        log.warn """
        ============================================================
        WARNING: Missing BWA index files for reference: ${ref_path}
        Missing index files:
        ${missing_indexes.collect { "  - ${it}" }.join('\n')}
        
        BWA requires these index files in the same location as the reference.
        You can create them with: bwa index ${ref_path}
        ============================================================
        """
        System.exit(1)
    }
    
    // Check for FASTA index (.fai) file
    def fai_file = file("${ref_path}.fai")
    if (!fai_file.exists()) {
        log.error """
        ============================================================
        ERROR: FASTA index file (.fai) not found: ${ref_path}.fai
        This file is required for the alignment workflow.
        You can create it with: samtools faidx ${ref_path}
        ============================================================
        """
        System.exit(1)
    }
    
    // Check for FASTA dictionary (.dict) file
    def dict_path = ref_path.replaceAll(/\.fa(sta)?(\.gz)?$/, '.dict')
    if (!file(dict_path).exists()) {
        log.warn """
        ============================================================
        WARNING: FASTA dictionary file (.dict) not found: ${dict_path}
        This file is required for some GATK operations.
        You can create it with: gatk CreateSequenceDictionary -R ${ref_path}
        ============================================================
        """
    }
    
    // Check other reference files
    if (params.workflow == "align" || params.workflow == "generate-gvcfs") {
        def required_files = [
            'dbsnp': params.reference_files.dbsnp,
            'known_indels_1': params.reference_files.known_indels_1,
            'known_indels_2': params.reference_files.known_indels_2
        ]
        
        def missing_files = [:]
        required_files.each { name, path ->
            def file_obj = file(path)
            if (!file_obj.exists()) {
                missing_files[name] = path
            }
        }
        
        if (missing_files.size() > 0) {
            log.warn """
            ============================================================
            WARNING: Missing reference files:
            ${missing_files.collect { name, path -> "  - ${name}: ${path}" }.join('\n')}
            
            Required reference files are missing. Please ensure all reference files exist.
            ============================================================
            """
            System.exit(1)
        }
    }
    
    log.info "All required reference files exist and are properly indexed."
}

// Define workflow functions
workflow FASTQ_QC {
    take:
    samplesheet

    main:
    log.info """
    ===============================================================================================
    🧬 FASTQ Quality Control Workflow 🧬
    ===============================================================================================
    This workflow performs comprehensive quality control on FASTQ files:
    • Running FastQC to assess read quality metrics
    • Generating quality reports for each sample
    • ${params.trim_reads ? 'Trimming low-quality bases and adapter sequences with Trimmomatic' : 'Skipping trimming (--trim_reads is disabled)'}
    • Creating MultiQC report to summarize all QC metrics
    • Inputs:
        • Sample Sheet: ${samplesheet}
        • Output Directory: ${params.outdir}
        • Trim Reads: ${params.trim_reads}
        • Adapter File: ${params.adapter_file}
    ===============================================================================================
    """
    
    // Create a channel from the samplesheet
    Channel
        .fromPath(samplesheet)
        .map { sheet ->
            // First check if file exists and has content
            if (!file(sheet).exists()) {
                error """
                ERROR: Sample sheet file not found: ${sheet}
                """.stripIndent()
                System.exit(1)
            }
            if (file(sheet).isEmpty()) {
                error """
                ERROR: Sample sheet file is empty: ${sheet}
                """.stripIndent()
                System.exit(1)
            }
            
            // Validate header line
            def header = file(sheet).readLines()[0]
            def expected_columns = ['SampleID', 'Gender', 'FastqR1', 'FastqR2', 'Flowcell', 'Lane', 'BAM', 'gVCF']
            def header_columns = header.split('\t')
            
            // Check for required columns
            def missing_columns = []
            expected_columns.each { required_col -> 
                if (!header_columns.contains(required_col)) {
                    missing_columns << required_col
                }
            }
            
            if (missing_columns.size() > 0) {
                log.error """
                ============================================================
                ERROR: Missing required columns in sample sheet: ${sheet}
                ============================================================
                Missing columns: ${missing_columns.join(', ')}
                
                Expected format:
                ${expected_columns.join('\t')}
                ============================================================
                """
                exit 1
            }
            
            // Continue with the file if validation passes
            return sheet
        }
        .splitCsv(header:true, sep: '\t')
        .map { row -> 
            // Check if all expected columns are present in this row
            if (row.size() < 8) {
                log.error """
                ============================================================
                ERROR: Incomplete row in sample sheet for sample ${row.SampleID}
                ============================================================
                Row has ${row.size()} columns but 8 are expected.
                Please check the sample sheet format and ensure all columns are present.
                
                Row content: ${row}
                ============================================================
                """
                return null
            }
            
            // Check if FastqR1 and FastqR2 columns exist and contain valid file paths
            def fastq_r1 = row.FastqR1 ? file(row.FastqR1, checkIfExists: true) : null
            def fastq_r2 = row.FastqR2 ? file(row.FastqR2, checkIfExists: true) : null
            
            if (!fastq_r1) {
                log.error 'FastqR1 file not found for sample ' + row.SampleID + ': ' + row.FastqR1
                return null
            }
            
            if (!fastq_r2 && params.paired_end) {
                log.error 'FastqR2 file not found for sample ' + row.SampleID + ': ' + row.FastqR2
                return null
            }
            
            def flowcell = row.containsKey('Flowcell') ? row.Flowcell : 'unknown'
            def lane = row.containsKey('Lane') ? row.Lane : 'unknown'
            
            return [row.SampleID, fastq_r1, fastq_r2, flowcell, lane]
        }
        .filter { it != null }
        .tap { reads_ch }
    
    // Log the number of samples
    reads_ch.count().subscribe { count ->
        log.info "Processing ${count} samples for FASTQ QC"
    }
    
    // Run the FASTQ QC workflow
    fastq_qc_workflow(reads_ch)

    // Log completion message with results
    if (params.trim_reads) {
        fastq_qc_workflow.out.updated_samplesheet.map { it.toString() }.subscribe { updated_sheet ->
            log.info """
        ===============================================================================================
        🧬 FASTQ Quality Control Workflow Completed 🧬
        ===============================================================================================
        • Outputs:
            • FASTQ QC reports: ${params.outdir}/${params.project_name}/${params.workflow}/fastqc/
            • MultiQC reports: ${params.outdir}/${params.project_name}/${params.workflow}/multiqc/multiqc_report.html
            • Trimmed reads available in: ${params.outdir}/${params.project_name}/${params.workflow}/trimmed/
            • Updated samplesheet: 
              - ${updated_sheet}
              - ${params.outdir}/${params.project_name}/${params.workflow}/samplesheets/${file(updated_sheet).name}
        ===============================================================================================
        """
        }
    } else {
        fastq_qc_workflow.out.multiqc_report.map { report ->
            log.info """
        ===============================================================================================
        🧬 FASTQ Quality Control Workflow Completed 🧬
        ===============================================================================================
        • Outputs:
            • FASTQ QC reports: ${params.outdir}/${params.project_name}/${params.workflow}/fastqc/
            • MultiQC reports: ${params.outdir}/${params.project_name}/${params.workflow}/multiqc/multiqc_report.html
            • Trimming was disabled
        ===============================================================================================
        """
        }
    }

    emit:
    qc_reports = fastq_qc_workflow.out.qc_reports
    qc_summary = fastq_qc_workflow.out.qc_summary
    reads = fastq_qc_workflow.out.reads
    trim_reports = fastq_qc_workflow.out.trim_reports
    multiqc_report = fastq_qc_workflow.out.multiqc_report
    multiqc_data = fastq_qc_workflow.out.multiqc_data
    updated_samplesheet = fastq_qc_workflow.out.updated_samplesheet
}

workflow BAM_QC {
    take:
    samplesheet

    main:
    log.info """
    ===============================================================================================
    🧬 BAM Quality Control Workflow 🧬
    ===============================================================================================
    This workflow performs comprehensive quality control on BAM files:
    • Validating BAM file format and integrity
    • Calculating alignment statistics and coverage metrics
    • Generating per-sample quality reports
    • Creating MultiQC report to summarize all alignment metrics
    • Inputs:
        • Sample Sheet: ${samplesheet}
        • Output Directory: ${params.outdir}
    ===============================================================================================
    """
    
    // Create input channel from samplesheet
    Channel.fromPath(file(samplesheet))
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def sample_id = row['SampleID']
            def bam = file(row['BAM'])
            def bai = bam.name.endsWith('.cram') ? file("${bam}.crai") : file("${bam}.bai")

            // Check if the BAM files exist
            if (!bam.exists() || !bai.exists()) {
                throw new RuntimeException("BAM/BAI files for sample ${sample_id} do not exist: ${bam}, ${bai}")
            }

            return [sample_id, bam, bai]
        }
        .set { bam_ch }

    // Run BAM QC workflow
    bam_qc_workflow(bam_ch)

    // Log completion message once after all samples are processed
    bam_qc_workflow.out.bam_qc_results.collect().subscribe { bam_qc_results ->
        log.info """
    ===============================================================================================
    🧬 BAM Quality Control Workflow Completed 🧬
    ===============================================================================================
    • BAM QC results are available in: ${params.outdir}/${params.project_name}/${params.workflow}/bams/
    • MultiQC report: ${params.outdir}/${params.project_name}/${params.workflow}/bams/multiqc_report.html
    ===============================================================================================
    """
    }

    emit:
    bam_qc_results = bam_qc_workflow.out.bam_qc_results
    multiqc_report = bam_qc_workflow.out.multiqc_report
    multiqc_data = bam_qc_workflow.out.multiqc_data
}

workflow VCF_QC {
    main:
    log.info """
    ===============================================================================================
    🧬 VCF Quality Control Workflow 🧬
    ===============================================================================================
    This workflow performs comprehensive quality control on VCF files:
    • Validating VCF file format and integrity
    • Calculating variant statistics and quality metrics
    • Generating per-sample variant reports
    • Creating MultiQC report to summarize all variant metrics
    • Inputs:
        • VCF Input: ${params.vcf_input ?: 'Not provided'}
        • VCF Directory: ${params.vcf_dir ?: 'Not provided'}
        • Output Directory: ${params.outdir}
    ===============================================================================================
    """
    
    // Create input channel based on provided parameters
    def vcf_ch = Channel.empty()
    
    if (params.vcf_input) {
        // Single VCF file provided
        def vcf_file = file(params.vcf_input)
        
        if (!vcf_file.exists()) {
            error """
            ============================================================
            ERROR: VCF file not found: ${params.vcf_input}
            ============================================================
            Please check that the file exists and the path is correct.
            """
        }
        
        def index = file("${params.vcf_input}.tbi")
        if (!index.exists()) {
            log.warn """
            ============================================================
            WARNING: VCF index (.tbi) not found for ${params.vcf_input}
            ============================================================
            Will attempt to create index during processing.
            """
        }
        
        vcf_ch = Channel.of([vcf_file, index])
    } 
    else if (params.vcf_dir) {
        // Directory of VCFs provided
        def vcf_dir = file(params.vcf_dir)
        
        if (!vcf_dir.exists() || !vcf_dir.isDirectory()) {
            error """
            ============================================================
            ERROR: VCF directory not found or is not a directory: ${params.vcf_dir}
            ============================================================
            Please check that the directory exists and the path is correct.
            """
        }
        
        // Find all VCF files in the directory
        vcf_ch = Channel.fromPath("${vcf_dir}/*.vcf.gz")
            .map { vcf_file -> 
                def vcf_index = file("${vcf_file}.tbi")
                if (!vcf_index.exists()) {
                    log.warn """
                    ============================================================
                    WARNING: VCF index (.tbi) not found for ${vcf_file}
                    ============================================================
                    Will attempt to create index during processing.
                    """
                }
                return [vcf_file, vcf_index]
            }
    }
    else if (params.sample_sheet) {
        // Legacy mode: create input channel from samplesheet
        log.warn "Using sample sheet mode for backward compatibility. Consider using --vcf_input or --vcf_dir instead."
        
        Channel.fromPath(file(params.sample_sheet))
            .splitCsv(header: true, sep: '\t')
            .map { row -> 
                def sample_id = row['SampleID']
                def vcf_path = row['VCF']
                def vcf = file(vcf_path)
                
                if (!vcf.exists()) {
                    log.error """
                    ============================================================
                    ERROR: VCF file not found for sample ${sample_id}
                    ============================================================
                    Path: ${vcf_path}
                    Please check that the file exists and the path is correct.
                    """
                    return null
                }
                
                def index = file("${vcf}.tbi")
                if (!index.exists()) {
                    log.warn """
                    ============================================================
                    WARNING: VCF index (.tbi) not found for sample ${sample_id}
                    ============================================================
                    Path: ${index}
                    Will attempt to create index during processing.
                    """
                }

                return [vcf, index]
            }
            .filter { it != null }
            .set { vcf_ch }
    }
    else {
        error """
        ===============================================================================================
        ERROR: No VCF input provided!
        ===============================================================================================
        Please provide one of the following:
          --vcf_input    : Path to a single VCF file
          --vcf_dir      : Path to a directory containing VCF files
          --sample_sheet : Legacy mode - path to a sample sheet with VCF file paths
        ===============================================================================================
        """
    }
    
    // Count and log the number of VCF files
    vcf_ch.count().subscribe { count ->
        log.info "Processing ${count} VCF files for QC"
    }

    // Run VCF QC workflow
    vcf_qc_workflow(vcf_ch)

    // Log completion message when the workflow completes
    vcf_qc_workflow.out.stats.collectFile(name: 'vcf_qc_complete.txt').subscribe { 
        log.info """
    ===============================================================================================
    🧬 VCF Quality Control Workflow Completed 🧬
    ===============================================================================================
    • VCF QC results are available in: ${params.outdir}/${params.project_name}/${params.workflow}/
    • MultiQC report: ${params.outdir}/${params.project_name}/${params.workflow}/multiqc_report.html
    ===============================================================================================
    """
    }

    emit:
    vcf_results = vcf_qc_workflow.out.stats
    multiqc_report = vcf_qc_workflow.out.multiqc_report
}

workflow ALIGN {
    take:
    samplesheet

    main:
    log.info """
    ===============================================================================================
    🧬 Sequence Alignment Workflow 🧬
    ===============================================================================================
    This workflow aligns FASTQ reads to the reference genome:
    • Aligning reads using BWA-MEM algorithm
    • Sorting and marking duplicates in BAM files
    • Adding read groups and other metadata
    • Generating alignment statistics and quality metrics
    • Inputs:
        • Sample Sheet: ${samplesheet}
        • Trim Reads: ${params.trim_reads ? "${params.trim_reads}" : ""}
        • Input type: ${params.trim_reads ? 'FastqTrimmedR1/R2' : 'FastqR1/R2'}
        • Output Directory: ${params.outdir}
    ===============================================================================================
    """
    
    // Log tool versions - use separate processes
    log_tool_version_samtools()
    log_tool_version_bwa()
    log_tool_version_gatk()

    // Check if reference files exist
    checkReferenceFilesExist()

    // Check if sample sheet is provided
    if (!params.sample_sheet)  {
        error "Please provide a sample sheet with --sample_sheet"
    }

    // Process samplesheet
    Channel.fromPath(file(samplesheet))
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_id = row['SampleID']
            
            // Validate required columns exist based on trim_reads parameter
            if (params.trim_reads) {
                if (!row.containsKey('FastqTrimmedR1') || !row.containsKey('FastqTrimmedR2')) {
                    error """
                    ERROR: Missing required columns for trimmed FASTQ files.
                    Sample ${sample_id} is missing FastqTrimmedR1 and/or FastqTrimmedR2 columns.
                    These columns are required when --trim_reads is enabled.
                    Please ensure the sample sheet contains these columns or disable --trim_reads.
                    """
                }
            }
            
            // Get FASTQ files based on trim_reads parameter
            def fastq_r1_file = params.trim_reads ? 
                (row['FastqTrimmedR1'] ? file(row['FastqTrimmedR1']) : null) : 
                (row['FastqR1'] ? file(row['FastqR1']) : null)
            def fastq_r2_file = params.trim_reads ? 
                (row['FastqTrimmedR2'] ? file(row['FastqTrimmedR2']) : null) : 
                (row['FastqR2'] ? file(row['FastqR2']) : null)
            def flowcell = row.containsKey('Flowcell') ? row['Flowcell'] : 'unknown'
            def lane = row.containsKey('Lane') ? row['Lane'] : 'unknown'

            // Check if the FASTQ files exist
            if (!fastq_r1_file || !fastq_r2_file) {
                error """
                ERROR: FASTQ files not found for sample ${sample_id}
                ${params.trim_reads ? 'Trimmed' : 'Input'} FASTQ files:
                R1: ${fastq_r1_file ?: 'missing'}
                R2: ${fastq_r2_file ?: 'missing'}
                Please check that the files exist and the paths are correct.
                """
            }

            if (!fastq_r1_file.exists() || !fastq_r2_file.exists()) {
                error """
                ERROR: FASTQ files do not exist for sample ${sample_id}
                ${params.trim_reads ? 'Trimmed' : 'Input'} FASTQ files:
                R1: ${fastq_r1_file}
                R2: ${fastq_r2_file}
                Please check that the files exist and the paths are correct.
                """
            }

            return [sample_id, fastq_r1_file, fastq_r2_file, flowcell, lane]
        }
        .set { samples_fastq }

    // Run alignment workflow with split processes
    run_bwa(samples_fastq)
    run_bam_sort(run_bwa.out.bam_file)
    run_mark_duplicates(run_bam_sort.out.raw_bam)
    
    run_create_recalibration_table(run_mark_duplicates.out.md_bam)
    run_recalibrate_bam(run_create_recalibration_table.out.recal_table)
    bam_to_cram(run_recalibrate_bam.out.recal_bam)
    run_cram_flagstat(bam_to_cram.out.cram_file)
    create_cram_md5sum(run_cram_flagstat.out.cram_stats)

    // Collect CRAM paths and update samplesheet
    bam_to_cram.out.cram_file
        .map { sample_id, cram, crai -> 
            def published_path = "${params.outdir}/${params.project_name}/${params.workflow}/crams/${cram.name}"
            // Check if the CRAM file exists and has content
            if (!cram.exists() || cram.size() == 0) {
                log.error "Error: CRAM file is empty or not created properly for sample ${sample_id}"
                exit 1
            }
            "${sample_id}\t${published_path}"
        }
        .collectFile(
            name: 'cram_paths.tsv',
            newLine: true,
            seed: "SampleID\tFilePath\n"
        )
        .set { cram_paths }

    // Update samplesheet with CRAM paths
    if (params.update_samplesheet) {
    //     log.info """
    // Updating samplesheet with CRAM paths...
    //     """
        update_samplesheet(
            file(params.sample_sheet),
            "align",
            "crams",
            cram_paths
        )
    }

    // Log completion message
    update_samplesheet.out.updated_samplesheet.map { it.toString() }.subscribe { updated_sheet ->
        def input_dir = file(samplesheet).parent
        def updated_filename = file(updated_sheet).name
        def published_samplesheet = "${input_dir}/${updated_filename}"
        log.info """
    ===============================================================================================
    🧬 Sequence Alignment Workflow Completed 🧬
    ===============================================================================================
    • Outputs:
        • CRAM files are available in: ${params.outdir}/${params.project_name}/${params.workflow}/crams/
        • MultiQC report: ${params.outdir}/${params.project_name}/${params.workflow}/multiqc/multiqc_report.html
        • Updated samplesheet: ${published_samplesheet}
    ===============================================================================================
    """
    }

    emit:
    cram_paths
}

workflow GENERATE_GVCFS {
    main:
    // Define local chroms variables internally rather than taking them as params
    def (chroms_par, chroms_auto, chroms_all) = initializeWorkflowParams()
    
    log.info """
    ===============================================================================================
    🧬 GVCF Generation Workflow 🧬
    ===============================================================================================
    This workflow generates Genomic VCF files from aligned reads:
    • Performing per-sample variant calling with GATK HaplotypeCaller
    • Processing ${chroms_auto.size()} autosomal chromosomes and ${chroms_par.size()} PAR regions
    • Creating GVCF files for downstream joint genotyping
    • Validating output files for data integrity
    • Inputs:
        • Sample Sheet: ${params.sample_sheet}
        • Chromosomes: ${chroms_auto}
        • PAR Chromosomes: ${chroms_par}
        • All Chromosomes: ${chroms_all}
        • MT Chromosome: ${params.mt}
        • Output Directory: ${params.outdir}
    ===============================================================================================
    """
    
    log_tool_version_gatk()

    // Process samples from the samplesheet - using params directly
    def all_samples = Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_id = row['SampleID']
            def gender = row['Gender']
            def bam = file(row['BAM'])

            // Check if the BAM file exists
            check_files([bam])
            return [sample_id, gender, bam]
        }

    // Filter samples by gender
    def samples_male = all_samples.filter { it[1] == 'M' }
    def samples_female = all_samples.filter { it[1] == 'F' }
    def samples_nosex = all_samples.filter { !(it[1] =~ 'F|M') }
    
    // Placeholder for collecting all GVCF paths
    def all_gvcf_paths = Channel.empty()

    // Placeholder for collecting all GVCF files
    def all_gvcfs = Channel.empty()
    
    // Process male samples
    // if (params.containsKey('generate_gvcfs') && params.generate_gvcfs) {
        // Run haplotype caller on each autosome
        run_haplotype_caller_auto_males(samples_male, chroms_auto)
        
        // Group results by sample ID and prepare for combining
        def male_calls = run_haplotype_caller_auto_males.out.auto_calls
            .groupTuple()
            .map { sample_id, gvcf_files, idx_files -> 
                [sample_id, gvcf_files.flatten(), idx_files.flatten()]
            }
        
        // Combine GVCFs for males
        run_combine_sample_gvcfs_males(male_calls)
        run_index_gvcf_males(run_combine_sample_gvcfs_males.out.combined_gvcf)
        run_create_gvcf_md5sum_males(run_index_gvcf_males.out.indexed_gvcf)
        
        // Collect all GVCF files
        all_gvcfs = all_gvcfs.mix(
            run_index_gvcf_males.out.indexed_gvcf
        )

        // Collect paths for output
        all_gvcf_paths = all_gvcf_paths.mix(
            run_combine_sample_gvcfs_males.out.gvcf_out_loc
                .map { sample_id, gvcf -> 
                    def published_path = "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs/${file(gvcf).name}"
                    return "${sample_id}\t${published_path}"
                }
        )
    // }
    
    // Process female samples
    // if (params.containsKey('generate_gvcfs') && params.generate_gvcfs) {
        // Run haplotype caller on each autosome
        run_haplotype_caller_auto_females(samples_female, chroms_auto)
        
        // Group results by sample ID and prepare for combining
        def female_calls = run_haplotype_caller_auto_females.out.auto_calls
            .groupTuple()
            .map { sample_id, gvcf_files, idx_files -> 
                [sample_id, gvcf_files.flatten(), idx_files.flatten()]
            }
        
        // Combine GVCFs for females
        run_combine_sample_gvcfs_females(female_calls)
        run_index_gvcf_females(run_combine_sample_gvcfs_females.out.combined_gvcf)
        run_create_gvcf_md5sum_females(run_index_gvcf_females.out.indexed_gvcf)
        
        // Collect all GVCF files
        all_gvcfs = all_gvcfs.mix(
            run_index_gvcf_females.out.indexed_gvcf
        )

        // Collect paths for output
        all_gvcf_paths = all_gvcf_paths.mix(
            run_combine_sample_gvcfs_females.out.gvcf_out_loc
                .map { sample_id, gvcf -> 
                    def published_path = "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs/${file(gvcf).name}"
                    return "${sample_id}\t${published_path}"
                }
        )
    // }
    
    // Process samples without gender distinction
    // if (params.containsKey('generate_gvcfs') && params.generate_gvcfs) {
        // Run haplotype caller on each autosome
        run_haplotype_caller_auto_nosex(samples_nosex, chroms_auto)
        
        // Group results by sample ID and prepare for combining
        def nosex_calls = run_haplotype_caller_auto_nosex.out.auto_calls
            .groupTuple()
            .map { sample_id, gvcf_files, idx_files -> 
                [sample_id, gvcf_files.flatten(), idx_files.flatten()]
            }
        
        // Combine GVCFs for nosex
        run_combine_sample_gvcfs_nosex(nosex_calls)
        run_index_gvcf_nosex(run_combine_sample_gvcfs_nosex.out.combined_gvcf)
        run_create_gvcf_md5sum_nosex(run_index_gvcf_nosex.out.indexed_gvcf)
        
        // Collect all GVCF files
        all_gvcfs = all_gvcfs.mix(
            run_index_gvcf_nosex.out.indexed_gvcf
        )

        // Collect paths for output
        all_gvcf_paths = all_gvcf_paths.mix(
            run_combine_sample_gvcfs_nosex.out.gvcf_out_loc
                .map { sample_id, gvcf -> 
                    def published_path = "${params.outdir}/${params.project_name}/${params.workflow}/gvcfs/${file(gvcf).name}"
                    return "${sample_id}\t${published_path}"
                }
        )
    // }

    // Validate GVCF files
    // run_validate_gvcf(all_gvcfs, chroms_auto)
    
    // Create output file with all paths
    def gvcf_paths_file = all_gvcf_paths
        .collectFile(
            name: 'gvcf_paths.tsv',
            newLine: true,
            seed: "SampleID\tFilePath\n"  // Add header line
        )
    
    // Update samplesheet if requested
    def updated_samplesheet_ch = Channel.value(file(params.sample_sheet))
    if (params.update_samplesheet) {
        updated_samplesheet_ch = update_samplesheet(
            file(params.sample_sheet),
            "generate-gvcfs",
            "gvcfs",
            gvcf_paths_file
        )
    }

    // Log completion message
    if (params.update_samplesheet) {
        updated_samplesheet_ch.map { it.toString() }.subscribe { updated_sheet ->
            def input_dir = file(params.sample_sheet).parent
            def updated_filename = file(updated_sheet).name
            def published_path = "${input_dir}/${updated_filename}"
            log.info """
        ===============================================================================================
        🧬 GVCF Generation Workflow Completed 🧬
        ===============================================================================================
        • GVCF files are available in: ${params.outdir}/${params.project_name}/${params.workflow}/gvcfs/
        • MultiQC report: ${params.outdir}/${params.project_name}/${params.workflow}/multiqc/multiqc_report.html
        • Updated samplesheet: ${published_path}
        ===============================================================================================
        """
        }
    }

    emit:
    gvcf_paths = all_gvcf_paths
}

workflow COMBINE_GVCFS {
    take:
    samplesheet
    chroms_all
    
    main:
    log.info """
    ===============================================================================================
    🧬 GVCF Combination Workflow 🧬
    ===============================================================================================
    This workflow combines sample GVCFs for joint genotyping:
    • Consolidating GVCFs from multiple samples
    • Processing ${chroms_all.size()} chromosomal regions
    • Preparing data for GenomicsDB import
    • Optimizing for efficient downstream analysis
    • Inputs:
        • Sample Sheet: ${samplesheet}
        • Chromosomes: ${chroms_all}
        • Output Directory: ${params.outdir}
    ===============================================================================================
    """
    
    log_tool_version_gatk()

    // Create input channel from samplesheet with validation
    Channel.fromPath(samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_id = row['SampleID']
            // Check both gVCF and GVCF column names
            def gvcf_path = row.containsKey('gVCF') ? row['gVCF'] : 
                           row.containsKey('GVCF') ? row['GVCF'] : null
            
            if (!gvcf_path) {
                error """
                ERROR: Missing GVCF path for sample ${sample_id}
                The sample sheet must contain either a 'gVCF' or 'GVCF' column.
                Available columns: ${row.keySet().join(', ')}
                """
            }
            
            def gvcf = file(gvcf_path)
            if (!gvcf.exists()) {
                error """
                ERROR: GVCF file not found for sample ${sample_id}
                Path: ${gvcf_path}
                Please check that the file exists and the path is correct.
                """
            }
            return [sample_id, gvcf]
        }
        .set { samples_gvcfs }

    // Create GVCF list and run combine process
    samples_gvcfs
        .collectFile() { item -> 
            ["gvcf.list", "${item[1]}\n"]
        }
        .set { gvcf_list }

    run_combine_gvcfs(gvcf_list, chroms_all)
    
    // Collect and process the combined GVCFs
    run_combine_gvcfs.out.cohort_chr_calls
        .flatten()
        .collect()
        .map { files -> 
            def vcfs = files.findAll { it.toString().endsWith('.g.vcf.gz') }
            def indices = files.findAll { it.toString().endsWith('.g.vcf.gz.tbi') }
            if (vcfs.isEmpty()) {
                error "No GVCF files found in output"
            }
            [vcfs, indices]
        }
        .set { cohort_calls }

    cohort_calls.view()
    run_concat_gvcfs(cohort_calls)

    // Create completion channel
    run_concat_gvcfs.out.concat_done
        .last()
        .subscribe { 
            log.info """
    ===============================================================================================
    🧬 GVCF Combination Workflow Completed 🧬
    ===============================================================================================
    • Combined GVCF files are available in: ${params.outdir}/${params.project_name}/${params.workflow}/gvcfs/
    • MultiQC report: ${params.outdir}/${params.project_name}/${params.workflow}/multiqc/multiqc_report.html
    ===============================================================================================
    """
        }

    emit:
    combined_gvcfs = run_concat_gvcfs.out.concat_done
}

workflow GENOMICS_DB_IMPORT {
    take:
    samplesheet
    chroms_all
    
    main:
    log.info """
    ===============================================================================================
    🧬 GenomicsDB Import Workflow 🧬
    ===============================================================================================
    This workflow consolidates sample gVCFs into a GenomicsDB:
    • Sample Sheet: ${samplesheet}
    • Chromosomes: ${chroms_all.size()} regions
    • Database Path: ${params.db_path ?: "${params.outdir}/genomicsdb"}
    • Database Update: ${params.db_update ? 'Yes' : 'No'}
    • Output Directory: ${params.outdir}
    ===============================================================================================
    """
    
    // Set default database path if not provided
    if (!params.db_path) {
        params.db_path = "${params.outdir}/genomicsdb"
    }
    
    // Create the directory if it doesn't exist
    new File(params.db_path).mkdirs()
    
    // Check if DB update is needed
    def db_update = params.db_update ?: false
    
    // Create input channel from samplesheet with validation
    gvcf_list_ch = Channel.fromPath(samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_id = row['SampleID']
            // Check both gVCF and GVCF column names
            def gvcf_path = row.containsKey('gVCF') ? row['gVCF'] : 
                           row.containsKey('GVCF') ? row['GVCF'] : null
            
            if (!gvcf_path || gvcf_path == 'NA') {
                log.warn """
                WARNING: Missing GVCF path for sample ${sample_id}
                The sample sheet must contain either a 'gVCF' or 'GVCF' column with valid paths.
                Available columns: ${row.keySet().join(', ')}
                Skipping this sample.
                """
                return null
            }
            
            def gvcf_file = new File(gvcf_path)
            if (!gvcf_file.exists()) {
                log.warn """
                WARNING: GVCF file not found for sample ${sample_id}
                Path: ${gvcf_path}
                Skipping this sample.
                """
                return null
            }
            return [sample_id, gvcf_path]
        }
        .filter { it != null }
        .collectFile(name: "gvcf.list", newLine: true) { item -> 
            "${item[0]}\t${item[1]}"
        }
    
    // Create channels for DB path and chromosomes
    db_channel = Channel.empty()
    
    if (db_update) {
        run_backup_genomic_db()
        run_genomics_db_import_update(gvcf_list_ch, chroms_all, run_backup_genomic_db.out.backup_status)
        run_genomics_db_import_update.out.interval_db.set { db_channel }
    } else {
        run_genomics_db_import_new(gvcf_list_ch, chroms_all)
        run_genomics_db_import_new.out.interval_db.set { db_channel }
    }
    
    db_channel
        .ifEmpty { error "No GenomicsDB content generated" }
        .subscribe {
            log.info """
    ===============================================================================================
    🧬 GenomicsDB Import Workflow Completed 🧬
    ===============================================================================================
    • GenomicsDB is available in: ${params.db_path}
    ===============================================================================================
    """
        }

    emit:
        genomics_db = db_channel
}

workflow GENOME_CALLING {
    take:
    chroms_all
    
    main:
    log.info """
    ===============================================================================================
    🧬 Genome Calling Workflow 🧬
    ===============================================================================================
    This workflow performs joint genotyping across samples:
    • Calling variants from the GenomicsDB database
    • Processing ${chroms_all.size()} chromosomal regions
    • Applying VQSR for SNPs and hard filtering for INDELs
    • Generating population-level VCF files
    • Inputs:
        • Database Path: ${params.db_path}
        • Chromosomes: ${chroms_all}
        • Output Directory: ${params.outdir}
    ===============================================================================================
    """
    
    log_tool_version_gatk()
    
    // Ensure database path is set
    if (!params.db_path) {
        error "Database path not provided. Please specify --db_path parameter."
    }
    
    // Check if database exists
    def db_dir = file(params.db_path)
    if (!db_dir.exists()) {
        error "Database directory ${params.db_path} does not exist. Please run the genomics-db-import workflow first."
    }
    
    // Create channel with database path and intervals
    Channel.fromList(chroms_all)
        .map { interval -> 
            def db_path = "${params.db_path}/${interval}.gdb"
            if (!file(db_path).exists()) {
                log.warn "Database for interval ${interval} not found at ${db_path}"
                return null
            }
            return tuple(interval, db_path)
        }
        .filter { it != null }
        .set { interval_db_ch }
    
    // Run genotyping on each interval
    run_genotype_gvcf_on_genome_db(interval_db_ch)
    
    // Group by interval and prepare for concatenation
    run_genotype_gvcf_on_genome_db.out.db_gg_vcf_set
        .groupTuple()
        .set { concat_ready }
        
    // Concatenate VCFs
    run_concat_vcf(concat_ready)
    
    // Apply VQSR on SNPs
    run_vqsr_on_snps(run_concat_vcf.out.combined_vcf)
    apply_vqsr_on_snps(run_vqsr_on_snps.out.snps_vqsr_recal)
    
    // For small datasets, always use hard filtering for INDELs
    hard_filter_indels(apply_vqsr_on_snps.out.snps_recalibrated)
    
    // Log completion message
    hard_filter_indels.out.hard_filtered_indels
        .subscribe {
            log.info """
    ===============================================================================================
    🧬 Genome Calling Workflow Completed 🧬
    ===============================================================================================
    • Final VCF is available in: ${params.outdir}/${params.project_name}/${params.workflow}/vcf/hard_filter/
    • Combined before VQSR and Hard filtering: ${params.outdir}/${params.project_name}/${params.workflow}/vcfs/${params.project_name}.vcf.gz
    • Final VCF SNPs: ${params.outdir}/${params.project_name}/${params.workflow}/vcfs/${params.project_name}.recal-SNP.vcf.gz
    • Final VCF INDELs: ${params.outdir}/${params.project_name}/${params.workflow}/vcfs/${params.project_name}.recal-INDEL.vcf.gz
    ===============================================================================================
    """
        }
        
    emit:
    final_vcf = hard_filter_indels.out.hard_filtered_indels
}

// Add new workflow definition for resource checking
workflow CHECK_RESOURCES {
    main:
    // Silent operation - no initial log message
    
    // Check resources
    def maxCpus = Runtime.runtime.availableProcessors()
    def maxMemory = Runtime.runtime.maxMemory()
    def maxMemoryGB = maxMemory / (1024 * 1024 * 1024)
    def critical = []
    
    // Check CPU and memory requirements
    workflow.session.config.process.each { name, value ->
        // Check default process settings
        if (name == 'cpus' && value > maxCpus) {
            critical << "Default process CPU setting (${value}) exceeds available CPUs (${maxCpus})"
        }
        
        // Check specific process settings
        if (value instanceof Map) {
            if (value.cpus && value.cpus instanceof Integer && value.cpus > maxCpus) {
                critical << "Process '${name}' CPU setting (${value.cpus}) exceeds available CPUs (${maxCpus})"
            }
            
            // Check memory settings
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
    
    // Check labeled processes
    workflow.session.config.process.each { name, value ->
        if (name == 'withLabel' && value instanceof Map) {
            value.each { label, settings ->
                if (settings.cpus && settings.cpus instanceof Integer && settings.cpus > maxCpus) {
                    critical << "Process 'withLabel:${label}' CPU setting (${settings.cpus}) exceeds available CPUs (${maxCpus})"
                }
                
                if (settings.memory) {
                    def memStr = settings.memory.toString()
                    if (memStr.endsWith('GB')) {
                        def mem = memStr.replace('GB', '').trim().toFloat()
                        if (mem > maxMemoryGB) {
                            critical << "Process 'withLabel:${label}' memory setting (${mem} GB) exceeds available memory (${String.format("%.2f", maxMemoryGB)} GB)"
                        }
                    }
                }
            }
        }
    }
    
    // Check named processes
    workflow.session.config.process.each { name, value ->
        if (name == 'withName' && value instanceof Map) {
            value.each { procName, settings ->
                if (settings.cpus && settings.cpus instanceof Integer && settings.cpus > maxCpus) {
                    critical << "Process 'withName:${procName}' CPU setting (${settings.cpus}) exceeds available CPUs (${maxCpus})"
                }
                
                if (settings.memory) {
                    def memStr = settings.memory.toString()
                    if (memStr.endsWith('GB')) {
                        def mem = memStr.replace('GB', '').trim().toFloat()
                        if (mem > maxMemoryGB) {
                            critical << "Process 'withName:${procName}' memory setting (${mem} GB) exceeds available memory (${String.format("%.2f", maxMemoryGB)} GB)"
                        }
                    }
                }
            }
        }
    }
    
    // Only display critical issues
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
        
        // Exit in a controlled way based on enforce_resource_check parameter
        if (params.enforce_resource_check) {
            System.exit(1)  // This exits without the full stack trace
        }
    }

    emit:
    checked = true
}

// Add the rest of the workflows
workflow DOWNLOAD_GATK {
    main:
    log.info """
    ===============================================================================================
    🧬 Reference Data Download Workflow 🧬
    ===============================================================================================
    This workflow downloads and prepares reference data:
    • Retrieving reference genome and index files
    • Downloading variant databases for GATK
    • Processing chromosome-specific reference files
    • Validating downloaded files for integrity
    • Specifying chromosomes: ${params.chromosomes}
    • Work directory: ${workflow.workDir}
    • Saving reference data to: ${params.ref_dir}
    ===============================================================================================
    """
    // Download reference and VCFs using the new workflow
    download_and_index_reference_workflow()
    download_and_index_vcfs()

    emit:
    done = true
}

// Add this at the appropriate place where other workflows are defined
workflow WORKFLOW_VALIDATE_GVCF {
    // Get gVCF files from sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def gvcf = row.gVCF ? file(row.gVCF) : null
            def sample_id = row.SampleID
            
            if (!gvcf || !gvcf.exists()) {
                log.error "gVCF file not found for sample ${sample_id}: ${gvcf}"
                return null
            }
            
            // Assume index exists in same location with .tbi extension
            def index = file("${gvcf}.tbi")
            if (!index.exists()) {
                log.error "gVCF index file not found for sample ${sample_id}: ${index}"
                return null
            }
            
            return [sample_id, gvcf, index]
        }
        .filter { it != null }
        .set { gvcf_files }
    
    // Run the gVCF validation workflow
    VALIDATE_GVCF(gvcf_files)
}

workflow FILTER_VCF {
    take:
    samplesheet

    main:
    log.info """
    ===============================================================================================
    🧬 VCF Filtering and Annotation Workflow 🧬
    ===============================================================================================
    This workflow filters and annotates VCF files:
    • Applying quality filters using GATK VariantFiltration
    • Annotating variants with Ensembl VEP (v${params.vep_cache_version})
    • Generating summary statistics and reports
    • Optionally combining chromosome-split VCFs into a single file
    • Inputs:
        • VCF Input: ${params.vcf_input ?: 'Not provided'}
        • VCF Directory: ${params.vcf_dir ?: 'Not provided'}
        • Output Directory: ${params.outdir}
        • Genome Version: ${params.genome_version ?: 'GRCh38.105'}
        • SnpEff Data Dir: ${params.snpeff_data_dir ?: 'Auto-generated'}
        • Combine VCFs: ${params.combine_vcfs ?: 'false'}
    ===============================================================================================
    """
    
    // Check if reference files exist
    checkReferenceFilesExist()
    
    // Get reference files
    def ref_fasta = file(params.reference_files.ref)
    def ref_fasta_fai = file("${params.reference_files.ref}.fai")
    def ref_dict = file(params.reference_files.ref.replaceFirst(/\.fa(sta)?$/, '.dict'))
    
    // Set genome version for annotation
    def genome_version = params.genome_version ?: 'GRCh38.105'
    
    // Create input channel based on provided parameters
    def vcf_ch = Channel.empty()
    
    if (params.vcf_input) {
        // Single VCF file provided
        def vcf_file = file(params.vcf_input)
        def vcf_index = file("${params.vcf_input}.tbi")
        
        if (!vcf_file.exists()) {
            error "VCF file not found: ${vcf_file}"
        }
        
        if (!vcf_index.exists()) {
            log.warn "VCF index not found: ${vcf_index}. Will attempt to create it."
            // Index will be created in the process
        }
        
        vcf_ch = Channel.of([vcf_file, vcf_index])
    } 
    else if (params.vcf_dir) {
        // Directory of split VCFs provided
        def vcf_dir = file(params.vcf_dir)
        
        if (!vcf_dir.exists() || !vcf_dir.isDirectory()) {
            error "VCF directory not found or is not a directory: ${vcf_dir}"
        }
        
        // Find all VCF files in the directory
        vcf_ch = Channel.fromPath("${vcf_dir}/*.vcf.gz")
            .map { vcf_file -> 
                def vcf_index = file("${vcf_file}.tbi")
                if (!vcf_index.exists()) {
                    log.warn "VCF index not found for ${vcf_file}. Will attempt to create it."
                    // Index will be created in the process
                }
                return [vcf_file, vcf_index]
            }
    }
    else if (samplesheet) {
        // Legacy mode: create input channel from samplesheet
        log.warn "Using sample sheet mode for compatibility. Consider using --vcf_input or --vcf_dir instead."
        
        Channel.fromPath(file(samplesheet))
            .splitCsv(header: true, sep: '\t')
            .map { row -> 
                def vcf = file(row['VCF'])
                def index = file("${vcf}.tbi")

                // Check if the VCF files exist
                if (!vcf.exists()) {
                    throw new RuntimeException("VCF file does not exist: ${vcf}")
                }
                
                if (!index.exists()) {
                    log.warn "VCF index not found for ${vcf}. Will attempt to create it."
                }

                return [vcf, index]
            }
            .set { vcf_ch }
    }
    else {
        error """
        ===============================================================================================
        ERROR: No VCF input provided!
        ===============================================================================================
        Please provide one of the following:
          --vcf_input     : Path to a single joint-called VCF file
          --vcf_dir       : Path to a directory containing chromosome-split VCF files
          --sample_sheet  : Legacy mode - path to a sample sheet with VCF file paths
        ===============================================================================================
        """
    }

    // Run VCF filtering and annotation workflow
    filter_annotate_vcf_workflow(vcf_ch, ref_fasta, ref_fasta_fai, ref_dict, genome_version)
    
    // Log completion message
    filter_annotate_vcf_workflow.out.annotated_vcfs.collect().subscribe() { filtered_vcf ->
    log.info """
    ===============================================================================================
    🧬 VCF Filtering and Annotation Workflow Completed 🧬
    ===============================================================================================
    • Filtered VCFs are available in: ${params.outdir}/${params.project_name}/${params.workflow}/filtered/
    • VEP Annotated VCFs are available in: ${params.outdir}/${params.project_name}/${params.workflow}/annotated/
    """ + (params.combine_vcfs ? "• Combined VCF: ${params.outdir}/${params.project_name}/${params.workflow}/combined/combined.vcf.gz" : "") + """
    ===============================================================================================
    """
    }

    emit:
    filtered_vcfs = filter_annotate_vcf_workflow.out.filtered_vcfs
    annotated_vcfs = filter_annotate_vcf_workflow.out.annotated_vcfs
    filter_stats = filter_annotate_vcf_workflow.out.filter_stats
    annotation_summary = filter_annotate_vcf_workflow.out.annotation_summary
    annotation_stats = filter_annotate_vcf_workflow.out.annotation_stats
    combined_vcf = filter_annotate_vcf_workflow.out.combined_vcf
}

// Main workflow
workflow {
    // Initialize workflow parameters
    def (chroms_par, chroms_auto, chroms_all) = initializeWorkflowParams()
    
    // Validate workflow
    validateWorkflow(params.workflow)
    
    // Check system resources only if requested or errors exist
    if (params.display_resource_check == null || params.display_resource_check) {
        CHECK_RESOURCES()
    }
    // Initialize reference directory for download_gatk workflow
    if (params.workflow == 'download_gatk') {
        if (!params.ref_dir) {
            params.ref_dir = "${params.outdir}/${params.project_name}/${params.build}"
        }
        file(params.ref_dir).mkdirs()
    }

    // Ensure that Singularity cache is set and exists
    if (!params.singularity_cache) {
        params.singularity_cache = "${params.outdir}/singularity_cache"
    }
    file(params.singularity_cache).mkdirs()
    
    // Ensure that Singularity temp directory exists
    if (!params.singularity_tempdir) {
        params.singularity_tempdir = "${params.outdir}/singularity_tmp"
    }
    file(params.singularity_tempdir).mkdirs()
    
    // Set Singularity environment variables if not already set
    if (!workflow.containerEngine && workflow.profile.contains('singularity')) {
        if (System.getenv('SINGULARITY_CACHEDIR') == null) {
            System.setProperty('SINGULARITY_CACHEDIR', params.singularity_cache)
        }
        if (System.getenv('SINGULARITY_TMPDIR') == null) {
            System.setProperty('SINGULARITY_TMPDIR', params.singularity_tempdir)
        }
    }

    // List of workflows that require a sample sheet
    def samplesheet_required = ['validate_fastq', 'fastq_qc', 'bam_qc', 'align', 'generate-gvcfs', 
                               'combine-gvcfs', 'genomics-db-import', 'validate-gvcf']
    
    // Check if sample sheet is required
    if (samplesheet_required.contains(params.workflow)) {
        if (!params.sample_sheet) {
            error """
            Please provide a sample sheet with --sample_sheet for the ${params.workflow} workflow
            """
        }
    }

    // List of workflows that require reference files
    def reference_required = ['align', 'generate-gvcfs', 'combine-gvcfs', 'genomics-db-import', 'genome-calling', 'validate-gvcf']
    
    // Check reference files exist for workflows that need them
    if (reference_required.contains(params.workflow)) {
        checkReferenceFilesExist()
    }

    // Execute workflow based on selection
    if (params.workflow == 'download_gatk') {
        checkRequiredFile(params.ref_dir, "Reference directory", "Reference data download workflow")
        DOWNLOAD_GATK()
    }
    else if (params.workflow == 'validate_fastq') {
        checkRequiredFile(params.sample_sheet, "Sample sheet", "FASTQ validation workflow")
        VALIDATE_FASTQ(params.sample_sheet)
    }
    else if (params.workflow == 'fastq_qc') {
        checkRequiredFile(params.sample_sheet, "Sample sheet", "FASTQ QC workflow")
        FASTQ_QC(params.sample_sheet)
    }
    else if (params.workflow == 'bam_qc') {
        checkRequiredFile(params.sample_sheet, "Sample sheet", "BAM QC workflow")
        BAM_QC(params.sample_sheet)
    }
    else if (params.workflow == 'validate-gvcf') {
        checkRequiredFile(params.sample_sheet, "Sample sheet", "gVCF validation workflow")
        WORKFLOW_VALIDATE_GVCF()
    }
    else if (params.workflow == 'vcf_qc') {
        // No need to check for sample sheet anymore
        VCF_QC()
    }
    else if (params.workflow == 'align') {
        checkRequiredFile(params.sample_sheet, "Sample sheet", "Alignment workflow")
            ALIGN(params.sample_sheet)
    }
    else if (params.workflow == 'generate-gvcfs') {
        checkRequiredFile(params.sample_sheet, "Sample sheet", "GVCF generation workflow")
        GENERATE_GVCFS()
    }
    else if (params.workflow == 'combine-gvcfs') {
        checkRequiredFile(params.sample_sheet, "Sample sheet", "GVCF combination workflow")
        COMBINE_GVCFS(params.sample_sheet, chroms_all)
    }
    else if (params.workflow == 'genomics-db-import') {
        checkRequiredFile(params.sample_sheet, "Sample sheet", "GenomicsDB import workflow")
        GENOMICS_DB_IMPORT(params.sample_sheet, chroms_all)
    }
    else if (params.workflow == 'genome-calling') {
        // Ensure database path is set
        if (!params.db_path) {
            error "Database path not provided. Please specify --db_path parameter."
        }
        
        // Check if database exists
        def db_dir = file(params.db_path)
        if (!db_dir.exists()) {
            error "Database directory ${params.db_path} does not exist. Please run the genomics-db-import workflow first."
        }
        
        GENOME_CALLING(chroms_all)
    }
    else if (params.workflow == 'filter-vcf') {
        if (!params.vcf_input && !params.vcf_dir && !params.sample_sheet) {
            error """
            ERROR: No VCF input provided for filter-vcf workflow.
            Please specify one of:
              --vcf_input     : Path to a single joint-called VCF file
              --vcf_dir       : Path to a directory containing chromosome-split VCF files
              --sample_sheet  : Legacy mode - path to a sample sheet with VCF file paths
            """
        }
        
        // Run the workflow with whatever input was provided
        FILTER_VCF(params.sample_sheet)
    }
    
    
}
