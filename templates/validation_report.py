#!/usr/bin/env python3

import sys
from datetime import datetime
from pathlib import Path

def write_section(f, title, char='='):
    separator = char * 47
    f.write(f"{separator}\n")
    if title:
        f.write(f"{title.center(47)}\n")
        f.write(f"{separator}\n")

def generate_validation_summary(report_files, samplesheet, output_file):
    # Initialize counters and lists
    total_samples = 0
    samples_with_errors = []
    validation_data = []

    # Process each validation file
    for report_file in report_files:
        if Path(report_file).exists():
            total_samples += 1
            with open(report_file) as f:
                content = f.read()
                validation_data.append(content)
                
                # Check for errors
                if '✗' in content:
                    # Extract sample ID from the validation report header
                    for line in content.split('\n'):
                        if 'Validation Report for' in line:
                            sample_id = line.split()[-1]
                            samples_with_errors.append(sample_id)
                            break

    # Write directly to the output file
    with open(output_file, 'w') as f:
        # Write header
        write_section(f, "Validation Summary")
        f.write(f"Sample Sheet: {samplesheet}\n")
        f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        write_section(f, "", char='-')
        f.write('\n')
        
        # Write statistics
        f.write(f"Total Samples: {total_samples}\n")
        f.write(f"Samples with Errors: {len(samples_with_errors)}\n")
        f.write(f"Samples Passed: {total_samples - len(samples_with_errors)}\n")
        write_section(f, "", char='-')
        f.write('\n')
        
        # List samples with errors
        if samples_with_errors:
            f.write("Samples with validation errors:\n")
            for sample in sorted(samples_with_errors):
                f.write(f"- {sample}\n")
            write_section(f, "", char='-')
            f.write('\n')
        
        # Write individual reports
        f.write("Individual Validation Reports:\n")
        write_section(f, "", char='-')
        f.write('\n')
        
        for report in validation_data:
            f.write(report)
            f.write('\n')
            write_section(f, "", char='-')
            f.write('\n')

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: validation_report.py <samplesheet> <output_file> <report_files...>")
        sys.exit(1)
    
    samplesheet = sys.argv[1]
    output_file = sys.argv[2]
    report_files = sys.argv[3:]
    
    try:
        generate_validation_summary(report_files, samplesheet, output_file)
    except Exception as e:
        print(f"Error generating validation summary: {str(e)}", file=sys.stderr)
        sys.exit(1) 