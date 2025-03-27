# VarCall DSL2 eLwazi NGS Variant Calling Course 2025

VarCall DSL2 is a modular Nextflow pipeline for variant calling following GATK Best Practices. It supports a variety of workflows from initial FASTQ quality control to joint genotyping of VCF files.

**Key Features:**
- Modular design using Nextflow DSL2
- Support for both GRCh37 (b37) and GRCh38 (b38) reference genomes
- Sex-aware variant calling for autosomal, X, Y, and MT chromosomes
- Container support via Singularity and Docker
- Comprehensive QC reporting with MultiQC
- Cluster execution support (SLURM, PBS)

## Pipeline Overview

The pipeline consists of multiple workflows that can be run independently or sequentially:

1. **fastq_qc**: Quality control and optional trimming of raw reads
2. **align**: BWA-MEM alignment, duplicate marking, and BQSR
3. **generate-gvcfs**: Per-sample variant calling with GATK HaplotypeCaller
4. **genomics-db-import**: Consolidate gVCFs into GenomicsDB
5. **genome-calling**: Joint variant calling using GenomicsDB
6. **filter-vcf**: Variant filtration using GATK and annotation using VEP
9. **download_gatk**: Download and prepare GATK resources



## Installation

### Prerequisites

- Nextflow (version 21.10.3 or later)
- Java 8 or later
- Singularity (recommended) or Docker
- Sufficient disk space for reference files and outputs

### Setup Steps

1. Clone the repository:

```bash
git clone https://github.com/elwazi/varcall-dsl2-ngs-course-2025
cd varcall-dsl2
```

2. Install Nextflow if not already installed:

```bash
curl -s https://get.nextflow.io | bash
```

3. Test the installation:

```bash
./nextflow run main.nf --help
```

## Configuration

The pipeline uses three configuration files:

1. **nextflow.config**: Default configuration with standard resource settings
2. **user.config**: User-specific configuration overrides
3. **profiles**: Predefined execution environments (standard, singularity, docker, ilifu, wits, cbio)

### Basic Configuration

Create a copy of the user.config file and modify it for your environment:

```bash
cp user.config my_user.config
```

Key parameters to configure:

```nextflow
params {
    // Project parameters
    project_dir = "/path/to/your/project"  // Base directory for your project
    outdir = "/path/to/output/directory"   // Directory where results will be stored
    sample_sheet = "/path/to/sample_sheet.csv"  // Path to your sample information file
    project_name = "your_project_name"     // Used for naming output directories and reports
    
    // Reference settings
    build = "b38"    // Reference genome build: "b37" or "b38"
    resources_dir = "/path/to/reference/data"  // Directory containing reference files
    
    // Analysis settings
    type = "wgs"     // Analysis type: "wgs" (whole genome) or "wes" (whole exome)
    trim_reads = true  // Whether to trim adapters and low-quality bases from reads
    
    // Resource checking options
    system_resource_check = true  // Verify system resources meet pipeline requirements
    enforce_resource_check = true  // Stop execution if resource requirements not met
    
    // FASTQ validation options
    force_validation_pass = true  // Continue despite FASTQ validation failures
    
    // System resources
    max_cpus = 4     // Maximum number of CPUs to allocate to any process
    max_memory = 14.GB  // Maximum memory to allocate to any process
}
```