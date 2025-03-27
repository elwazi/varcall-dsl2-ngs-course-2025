#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// FASTQ QC Process for Raw Reads
process raw_fastqc {
    tag { "fQCraw_${sample_id}" }
    label 'big_mem'
    label 'fastqc'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/fastqc/raw", mode: 'copy',
        pattern: "*_fastqc.{html,zip}"
    cache 'deep'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(fastq_r1), path(fastq_r2), val(flowcell), val(lane)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html"), emit: fastqc_reports
    tuple val(sample_id), path("${sample_id}_fastqc_summary.txt"), emit: fastqc_summary
    path "*_fastqc.zip", emit: zip_files

    script:
    """
    # Run FastQC on both read files
    echo "Running FastQC on raw reads..."
    fastqc -t 2 --noextract ${fastq_r1} ${fastq_r2}
    
    # Create summary report
    echo "Sample: ${sample_id}" > ${sample_id}_fastqc_summary.txt
    echo "Flowcell: ${flowcell}" >> ${sample_id}_fastqc_summary.txt
    echo "Lane: ${lane}" >> ${sample_id}_fastqc_summary.txt
    echo "" >> ${sample_id}_fastqc_summary.txt
    
    # Extract basic statistics and per base sequence quality from FastQC data
    for zip in *_fastqc.zip; do
        unzip -p \$zip */fastqc_data.txt | awk '/>>Basic Statistics/,/>>END_MODULE/' >> ${sample_id}_fastqc_summary.txt
        echo "" >> ${sample_id}_fastqc_summary.txt
        unzip -p \$zip */fastqc_data.txt | awk '/>>Per base sequence quality/,/>>END_MODULE/' >> ${sample_id}_fastqc_summary.txt
        echo "" >> ${sample_id}_fastqc_summary.txt
    done
    """
}

// FASTQ QC Process for Trimmed Reads
process trimmed_fastqc {
    tag { "fQCtrim_${sample_id}" }
    label 'big_mem'
    label 'fastqc'
    publishDir "${params.outdir}/${params.project_name}/${params.workflow}/fastqc/trimmed", mode: 'copy',
        pattern: "*_fastqc.{html,zip}"
    cache 'deep'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(fastq_r1), path(fastq_r2), val(flowcell), val(lane)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html"), emit: fastqc_reports
    tuple val(sample_id), path("${sample_id}_fastqc_summary.txt"), emit: fastqc_summary
    path "*_fastqc.zip", emit: zip_files

    script:
    """
    # Run FastQC on both read files
    echo "Running FastQC on trimmed reads..."
    fastqc -t 2 --noextract ${fastq_r1} ${fastq_r2}
    
    # Create summary report
    echo "Sample: ${sample_id}" > ${sample_id}_fastqc_summary.txt
    echo "Flowcell: ${flowcell}" >> ${sample_id}_fastqc_summary.txt
    echo "Lane: ${lane}" >> ${sample_id}_fastqc_summary.txt
    echo "" >> ${sample_id}_fastqc_summary.txt
    
    # Extract basic statistics and per base sequence quality from FastQC data
    for zip in *_fastqc.zip; do
        unzip -p \$zip */fastqc_data.txt | awk '/>>Basic Statistics/,/>>END_MODULE/' >> ${sample_id}_fastqc_summary.txt
        echo "" >> ${sample_id}_fastqc_summary.txt
        unzip -p \$zip */fastqc_data.txt | awk '/>>Per base sequence quality/,/>>END_MODULE/' >> ${sample_id}_fastqc_summary.txt
        echo "" >> ${sample_id}_fastqc_summary.txt
    done
    """
} 