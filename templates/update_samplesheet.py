#!/usr/bin/env python3
import csv
import os
import sys

# Debug information
print("Current working directory:", os.getcwd())
print("Samplesheet path:", "${samplesheet}")
print("Paths file:", "${paths_file}")
print("Workflow type:", "${workflow_type}")
print("Suffix:", "${suffix}")

# Generate output filename using baseName and suffix
output_file = "${samplesheet.baseName}_${suffix}.csv"
print("Output file:", output_file)

# Check if files exist
if not os.path.exists("${samplesheet}"):
    print("ERROR: Samplesheet file not found: ${samplesheet}")
    sys.exit(1)

# Detect if samplesheet is tab-delimited or comma-delimited
delimiter = '\\t'  # Default to tab
with open("${samplesheet}", 'r') as f:
    first_line = f.readline().strip()
    if ',' in first_line and '\\t' not in first_line:
        delimiter = ','
        print("Assuming comma-delimited samplesheet")
    else:
        print("Detected tab-delimited samplesheet")

# Read the samplesheet
samplesheet_data = []
headers = []
try:
    with open("${samplesheet}", 'r') as f:
        reader = csv.reader(f, delimiter=delimiter)
        headers = next(reader)
        samplesheet_data = list(reader)
    
    print("Original samplesheet columns:", headers)
    print("Original samplesheet rows:", len(samplesheet_data))
except Exception as e:
    print("ERROR: Failed to read samplesheet:", e)
    sys.exit(1)

# Check if paths file exists and has content
if not os.path.exists("${paths_file}"):
    print("WARNING: Paths file not found: ${paths_file}")
    print("Creating an empty paths file with headers")
    with open("${paths_file}", 'w') as f:
        if "${workflow_type}" == "fastq_qc":
            f.write("SampleID\\tFilePathR1\\tFilePathR2\\n")
        else:
            f.write("SampleID\\tFilePath\\n")

# Check if paths file is empty (only header)
file_size = os.path.getsize("${paths_file}")
if file_size == 0:
    print("WARNING: Paths file is empty: ${paths_file}")
    with open("${paths_file}", 'w') as f:
        if "${workflow_type}" == "fastq_qc":
            f.write("SampleID\\tFilePathR1\\tFilePathR2\\n")
        else:
            f.write("SampleID\\tFilePath\\n")

# Read the paths file
paths_data = []
paths_headers = []
try:
    with open("${paths_file}", 'r') as f:
        reader = csv.reader(f, delimiter='\\t')
        paths_headers = next(reader)
        paths_data = list(reader)
    
    print("Paths file columns:", paths_headers)
    print("Paths file rows:", len(paths_data))
except Exception as e:
    print("WARNING: Failed to read paths file:", e)
    if "${workflow_type}" == "fastq_qc":
        paths_headers = ['SampleID', 'FilePathR1', 'FilePathR2']
    else:
        paths_headers = ['SampleID', 'FilePath']
    paths_data = []

# Ensure SampleID column exists in the samplesheet
if 'SampleID' not in headers:
    print("ERROR: SampleID column not found in samplesheet")
    print("The samplesheet must contain a 'SampleID' column")
    sys.exit(1)

sampleid_index = headers.index('SampleID')

# Update the samplesheet based on workflow type
if "${workflow_type}" == "fastq_qc":
    # Check if paths has separate R1 and R2 columns and has data
    if 'FilePathR1' in paths_headers and 'FilePathR2' in paths_headers and len(paths_data) > 0:
        print("Found separate R1 and R2 columns in paths file")
        
        # Create FastqTrimmedR1 and FastqTrimmedR2 columns if they don't exist
        if 'FastqTrimmedR1' not in headers:
            print("Creating FastqTrimmedR1 column")
            headers.append('FastqTrimmedR1')
            for row in samplesheet_data:
                row.append("")
        
        if 'FastqTrimmedR2' not in headers:
            print("Creating FastqTrimmedR2 column")
            headers.append('FastqTrimmedR2')
            for row in samplesheet_data:
                row.append("")
        
        r1_index = headers.index('FastqTrimmedR1')
        r2_index = headers.index('FastqTrimmedR2')
        
        # Get indices for paths file
        paths_sampleid_idx = paths_headers.index('SampleID')
        paths_r1_idx = paths_headers.index('FilePathR1')
        paths_r2_idx = paths_headers.index('FilePathR2')
        
        # Update the columns with the trimmed FASTQ paths
        for path_row in paths_data:
            sample_id = path_row[paths_sampleid_idx] if len(path_row) > paths_sampleid_idx else ""
            if not sample_id:
                print("WARNING: Row has no valid SampleID, skipping")
                continue
                
            r1_path = path_row[paths_r1_idx] if len(path_row) > paths_r1_idx else ""
            r2_path = path_row[paths_r2_idx] if len(path_row) > paths_r2_idx else ""
            
            # Find the sample in the samplesheet
            updated = False
            for i, row in enumerate(samplesheet_data):
                if row[sampleid_index] == sample_id:
                    row[r1_index] = r1_path
                    row[r2_index] = r2_path
                    print(f"Updated sample {sample_id} with trimmed FastQ paths")
                    updated = True
            
            if not updated:
                print(f"WARNING: Sample {sample_id} not found in samplesheet")
    
    # For backward compatibility - handle single FilePath column
    elif 'FilePath' in paths_headers and len(paths_data) > 0:
        print("Found single FilePath column in paths file")
        
        # Create FastqTrimmed column if it doesn't exist
        if 'FastqTrimmed' not in headers:
            print("Creating FastqTrimmed column")
            headers.append('FastqTrimmed')
            for row in samplesheet_data:
                row.append("")
        
        trimmed_index = headers.index('FastqTrimmed')
        
        # Get indices for paths file
        paths_sampleid_idx = paths_headers.index('SampleID')
        paths_file_idx = paths_headers.index('FilePath')
        
        # Update the column with the trimmed FASTQ paths
        for path_row in paths_data:
            sample_id = path_row[paths_sampleid_idx] if len(path_row) > paths_sampleid_idx else ""
            if not sample_id:
                print("WARNING: Row has no valid SampleID, skipping")
                continue
                
            file_path = path_row[paths_file_idx] if len(path_row) > paths_file_idx else ""
            
            # Find the sample in the samplesheet
            updated = False
            for i, row in enumerate(samplesheet_data):
                if row[sampleid_index] == sample_id:
                    row[trimmed_index] = file_path
                    print(f"Updated sample {sample_id} with path: {file_path}")
                    updated = True
            
            if not updated:
                print(f"WARNING: Sample {sample_id} not found in samplesheet")
    
    else:
        print("WARNING: Paths file does not have expected columns or is empty")
        print("Creating FastqTrimmedR1 and FastqTrimmedR2 columns without updating values")
        if 'FastqTrimmedR1' not in headers:
            headers.append('FastqTrimmedR1')
            for row in samplesheet_data:
                row.append("")
        if 'FastqTrimmedR2' not in headers:
            headers.append('FastqTrimmedR2')
            for row in samplesheet_data:
                row.append("")

elif "${workflow_type}" == "align":
    # Create BAM column if it doesn't exist
    if 'BAM' not in headers:
        print("Creating BAM column")
        headers.append('BAM')
        for row in samplesheet_data:
            row.append("")
    
    # Update the column with the BAM paths if paths file has data
    if 'FilePath' in paths_headers and len(paths_data) > 0:
        bam_index = headers.index('BAM')
        
        # Get indices for paths file
        paths_sampleid_idx = paths_headers.index('SampleID')
        paths_file_idx = paths_headers.index('FilePath')
        
        for path_row in paths_data:
            sample_id = path_row[paths_sampleid_idx] if len(path_row) > paths_sampleid_idx else ""
            if not sample_id:
                print("WARNING: Row has no valid SampleID, skipping")
                continue
                
            file_path = path_row[paths_file_idx] if len(path_row) > paths_file_idx else ""
            
            # Find the sample in the samplesheet
            updated = False
            for i, row in enumerate(samplesheet_data):
                if row[sampleid_index] == sample_id:
                    row[bam_index] = file_path
                    print(f"Updated sample {sample_id} with BAM path: {file_path}")
                    updated = True
            
            if not updated:
                print(f"WARNING: Sample {sample_id} not found in samplesheet")
    else:
        print("WARNING: Paths file does not have expected columns or is empty")

elif "${workflow_type}" == "generate-gvcfs":
    # Create GVCF column if it doesn't exist
    if 'gVCF' not in headers:
        print("Creating GVCF column")
        headers.append('gVCF')
        for row in samplesheet_data:
            row.append("")
    
    # Update the column with the GVCF paths if paths file has data
    if 'FilePath' in paths_headers and len(paths_data) > 0:
        gvcf_index = headers.index('gVCF')
        
        # Get indices for paths file
        paths_sampleid_idx = paths_headers.index('SampleID')
        paths_file_idx = paths_headers.index('FilePath')
        
        for path_row in paths_data:
            sample_id = path_row[paths_sampleid_idx] if len(path_row) > paths_sampleid_idx else ""
            if not sample_id:
                print("WARNING: Row has no valid SampleID, skipping")
                continue
                
            file_path = path_row[paths_file_idx] if len(path_row) > paths_file_idx else ""
            
            # Find the sample in the samplesheet
            updated = False
            for i, row in enumerate(samplesheet_data):
                if row[sampleid_index] == sample_id:
                    row[gvcf_index] = file_path
                    print(f"Updated sample {sample_id} with GVCF path: {file_path}")
                    updated = True
            
            if not updated:
                print(f"WARNING: Sample {sample_id} not found in samplesheet")
    else:
        print("WARNING: Paths file does not have expected columns or is empty")

elif "${workflow_type}" == "vcf":
    # Create VCF column if it doesn't exist
    if 'VCF' not in headers:
        print("Creating VCF column")
        headers.append('VCF')
        for row in samplesheet_data:
            row.append("")
    
    # Update the column with the VCF paths if paths file has data
    if 'FilePath' in paths_headers and len(paths_data) > 0:
        vcf_index = headers.index('VCF')
        
        # Get indices for paths file
        paths_sampleid_idx = paths_headers.index('SampleID')
        paths_file_idx = paths_headers.index('FilePath')
        
        for path_row in paths_data:
            sample_id = path_row[paths_sampleid_idx] if len(path_row) > paths_sampleid_idx else ""
            if not sample_id:
                print("WARNING: Row has no valid SampleID, skipping")
                continue
                
            file_path = path_row[paths_file_idx] if len(path_row) > paths_file_idx else ""
            
            # Find the sample in the samplesheet
            updated = False
            for i, row in enumerate(samplesheet_data):
                if row[sampleid_index] == sample_id:
                    row[vcf_index] = file_path
                    print(f"Updated sample {sample_id} with VCF path: {file_path}")
                    updated = True
            
            if not updated:
                print(f"WARNING: Sample {sample_id} not found in samplesheet")
    else:
        print("WARNING: Paths file does not have expected columns or is empty")

# Save the updated samplesheet
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\\t')
    writer.writerow(headers)
    writer.writerows(samplesheet_data)

print(f"Updated samplesheet saved to {output_file}")
print("Updated samplesheet columns:", headers)
print("Updated samplesheet rows:", len(samplesheet_data))