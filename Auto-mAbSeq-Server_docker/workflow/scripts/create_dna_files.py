#!/usr/bin/env python3

"""
Create individual .dna files from verified nucleotide sequences.
Simplifies sequence IDs by removing strand and frame information.

On input:
- verify_frames.fasta — verified nucleotide sequences with IDs like:
  VHH_raw_cDNA_Cm_TN95089_fwd_frame1

On output:
- results/dna_files/*.dna — individual .dna files with simplified names:
  VHH_cDNA_Cm_TN95089.dna
"""

import yaml
import re
import os
from pathlib import Path
from Bio import SeqIO

def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def get_data_dirs():
    """
    Получить директории из переменных окружения или использовать defaults.
    Это позволяет изолировать каждый job в отдельную директорию.
    """
    data_dir = Path(os.environ.get('DATA_DIR', 'data'))
    intermediate_dir = Path(os.environ.get('INTERMEDIATE_DIR', data_dir / 'intermediate'))
    results_dir = Path(os.environ.get('RESULTS_DIR', data_dir / 'results'))
    
    
    return {
        'data_dir': data_dir,
        'intermediate_dir': intermediate_dir,
        'results_dir': results_dir
    }

def simplify_sequence_id(complex_id):
    """
    Simplify sequence ID by removing strand and frame information.
    Examples:
    - VHH_raw_cDNA_Cm_TN95089_fwd_frame1 → VHH_cDNA_Cm_TN95089
    - VHeavy_raw_cDNA_Cm_TN95077_rev_frame2 → VHeavy_cDNA_Cm_TN95077
    - VKappa_raw_cDNA_Ms_lib#1_1_fwd_frame3 → VKappa_cDNA_Ms_lib#1_1
    """
    # Remove _raw from the middle (raw_cDNA → cDNA)
    simplified = complex_id.replace("_raw_", "_")
    
    # Remove _fwd_frame*, _rev_frame* from the end
    simplified = re.sub(r'_(fwd|rev)_frame\d+$', '', simplified)
    
    return simplified

def create_dna_files(nucleotide_fasta, output_dir):
    """
    Create individual .dna files for each verified nucleotide sequence.
    
    Args:
    - nucleotide_fasta: Path to FASTA file with nucleotide sequences
    - output_dir: Directory for output .dna files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    count = 0
    id_mapping = []
    
    print(f"\n{'='*80}")
    print(f"Creating Individual .dna Files")
    print(f"{'='*80}\n")
    print(f"{'#':<4} {'Original ID':<45} {'→':<2} {'File Name':<35}")
    print(f"{'-'*86}")
    
    for idx, record in enumerate(SeqIO.parse(nucleotide_fasta, "fasta"), 1):
        original_id = record.id
        simplified_id = simplify_sequence_id(original_id)
        filename = f"{simplified_id}.dna"
        filepath = output_dir / filename
        
        # Write to .dna file
        with open(filepath, 'w') as f:
            f.write(f">{simplified_id}\n")
            f.write(f"{str(record.seq)}\n")
        
        print(f"{idx:<4} {original_id:<45} → {filename:<35}")
        id_mapping.append({
            "original": original_id,
            "simplified": simplified_id,
            "filename": filename
        })
        count += 1
    
    print(f"\n✓ Created {count} .dna files in {output_dir}")
    return count, id_mapping

def save_id_mapping_log(id_mapping, log_file):
    """Save ID mapping for reference"""
    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(log_path, 'w') as f:
        f.write("DNA File ID Simplification Log\n")
        f.write("=" * 90 + "\n\n")
        f.write(f"{'Original ID':<45} → {'Simplified ID':<40}\n")
        f.write("-" * 87 + "\n")
        for entry in id_mapping:
            f.write(f"{entry['original']:<45} → {entry['simplified']:<40}\n")

def main(config_file="config/config.yaml"):
    config = load_config(args.config)
    dirs = get_data_dirs()
    
    # Используем переменные окружения для путей
    nucleotide_fasta = dirs['intermediate_dir'] / "verify_frames.fasta"
    output_dir = dirs['results_dir'] / "dna_files"
    mapping_log = dirs['data_dir'] / "logs" / "dna_files_id_mapping.txt"
    
    # Check input file exists
    if not Path(nucleotide_fasta).exists():
        print(f"✗ Nucleotide FASTA not found: {nucleotide_fasta}")
        print(f"✗ Run extract_verified_nucleotides.py first")
        return False
    
    # Create .dna files
    count, id_mapping = create_dna_files(nucleotide_fasta, output_dir)
    
    if count == 0:
        print("✗ No .dna files created")
        return False
    
    # Save mapping log
    save_id_mapping_log(id_mapping, mapping_log)
    print(f"✓ ID mapping log: {mapping_log}")
    
    print(f"\n{'='*80}")
    print(f"✓ COMPLETE")
    print(f"✓ Created {count} .dna files in {output_dir}")
    print(f"{'='*80}\n")
    
    return True

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Create .dna files from verified nucleotide sequences"
    )
    parser.add_argument("--config", default="config/config.yaml",
                        help="Path to config file")
    args = parser.parse_args()
    
    success = main(args.config)
    exit(0 if success else 1)
