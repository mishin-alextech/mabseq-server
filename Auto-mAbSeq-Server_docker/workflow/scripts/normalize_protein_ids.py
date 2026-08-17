#!/usr/bin/env python3
"""
Step 7: Normalize protein sequence IDs from ANARCI results

"""

import yaml
import argparse
import re
from pathlib import Path
from Bio import SeqIO
import os


def get_data_dirs():
    """
    Получить директории из переменных окружения или использовать defaults.
    Это позволяет изолировать каждый job в отдельную директорию.
    """
    data_dir = Path(os.environ.get('DATA_DIR', 'data'))
    raw_data_dir = Path(os.environ.get('RAW_DATA_DIR', data_dir / 'raw' / 'sanger'))
    intermediate_dir = Path(os.environ.get('INTERMEDIATE_DIR', data_dir / 'intermediate'))
    results_dir = Path(os.environ.get('RESULTS_DIR', data_dir / 'results'))
    logs_dir = Path(os.environ.get('LOGS_DIR', data_dir / 'logs'))
    
    return {
        'results_dir': results_dir,   
        'logs_dir': logs_dir,  
    }


def load_config(config_file):
    """Load YAML configuration"""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def normalize_protein_id(original_id):
    """
    Remove _fwd_frame* or _rev_frame* suffix from protein ID
    
    Example:
        VHH_cDNA_Cm_TN95089_fwd_frame1 → VHH_cDNA_Cm_TN95089
        VHH_cDNA_Cm_TN95090_rev_frame2 → VHH_cDNA_Cm_TN95090
    """
    # Remove _fwd_frame[123] or _rev_frame[123] suffix
    normalized = re.sub(r'_(fwd|rev)_frame[123]$', '', original_id)
    return normalized

def main():
    parser = argparse.ArgumentParser(
        description='Normalize protein sequence IDs from ANARCI results'
    )
    parser.add_argument('--config', required=True, help='Path to config.yaml')
    parser.add_argument('--dry-run', action='store_true', help='Perform dry run without writing')
    
    args = parser.parse_args()
    
    config = load_config(args.config)
    dirs = get_data_dirs()  
    
    # Используем только env vars 
    input_fasta = dirs['results_dir'] / 'dna_files' / 'protein_valid_anarci.fasta'
    output_fasta = dirs['results_dir'] / 'dna_files' / 'protein_valid_normalized.fasta'
    log_dir = dirs['logs_dir']  
    
    # Create output directories if needed
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Read input FASTA
    if not input_fasta.exists():
        print(f" Input file not found: {input_fasta}")
        return False
    
    # Process sequences
    output_records = []
    normalization_log = []
    normalized_count = 0
    
    print(f"\n Reading sequences from: {input_fasta}")
    
    for record in SeqIO.parse(str(input_fasta), 'fasta'):
        original_id = record.id
        normalized_id = normalize_protein_id(original_id)
        
        # Log: original_id -> normalized_id
        normalization_log.append(f"{original_id} -> {normalized_id}")
        
        # Create normalized record
        record.id = normalized_id
        record.description = ''  # Clear description to avoid duplication
        output_records.append(record)
        normalized_count += 1
    
    # Write output
    if not args.dry_run:
        SeqIO.write(output_records, str(output_fasta), 'fasta')
        print(f" Written {len(output_records)} sequences to: {output_fasta}")
    else:
        print(f" DRY RUN: Would write {len(output_records)} sequences to: {output_fasta}")
    
    # Write simple log (only ID mappings)
    log_file = log_dir / 'protein_normalization.log'
    with open(log_file, 'w') as f:
        for log_entry in normalization_log:
            f.write(log_entry + '\n')
    
    print(f" Log file: {log_file}")
    print(f"\n Summary:")
    print(f"   Normalized: {normalized_count}")
    print(f"   Total:      {len(output_records)}")
    print(f"\n Step 7 completed!")
    
    return True

if __name__ == '__main__':
    main()

