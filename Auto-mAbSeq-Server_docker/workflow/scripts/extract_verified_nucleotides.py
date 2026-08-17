#!/usr/bin/env python3

"""
На входе:
- vhh_valid_anarci.fasta — верифицированные БЕЛКИ (используем только ID)
- raw_seqs.fasta — ИСХОДНЫЕ НУКЛЕОТИДНЫЕ последовательности (без фреймов)

На выходе:
- verify_frames.fasta — НУКЛЕОТИДНЫЕ последовательности c нормализованными по ID названиями
"""

import yaml
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

def extract_verified_ids(protein_fasta):
    """
    Извлекает ID верифицированных аминокислотных последовательностей.
    ID содержат информацию об исходном файле и фрейме, например:
    - sample1_fwd_frame1
    
    Нужно получить исходный ID файла (sample1), чтобы найти его в raw_seqs.fasta
    """
    verified_ids = {}
    
    for record in SeqIO.parse(protein_fasta, "fasta"):
        parts = record.id.split("_")
        
        if "fwd" in parts or "rev" in parts:
            if "fwd" in parts:
                idx = parts.index("fwd")
            else:
                idx = parts.index("rev")
            original_id = "_".join(parts[:idx])
        else:
            # Если нет fwd/rev, используем весь ID
            original_id = record.id
        
        verified_ids[original_id] = record.id
    
    print(f"Found {len(verified_ids)} verified sequence IDs")
    return verified_ids

def extract_nucleotide_sequences(raw_seqs_fasta, verified_ids, output_fasta):
    """
    Извлекает нуклеотидные последовательности.
    
    Args:
    - raw_seqs_fasta: путь к файлу с исходными нуклеотидными последовательностями
    - verified_ids: dict {original_id: protein_id} верифицированных последовательностей
    - output_fasta: путь к выходному файлу
    """
    count = 0
    not_found = []
    
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(raw_seqs_fasta, "fasta"):
            if record.id in verified_ids:
                verified_protein_id = verified_ids[record.id]
                record.id = verified_protein_id
                record.description = ""
                
                SeqIO.write(record, out_f, "fasta")
                count += 1
                print(f"  ✓ {record.id}: {len(record.seq)} nt")
            else:
                if record.id not in verified_ids:
                    not_found.append(record.id)
    
    if not_found:
        print(f"\n✗ Not verified ({len(not_found)}): {', '.join(not_found[:5])}...")
    
    print(f"✓ Wrote {count} nucleotide sequences to {output_fasta}")
    
    return count

def main(config_file="config/config.yaml"):
    
    config = load_config(config_file)
    dirs = get_data_dirs()
    
    protein_fasta = dirs['results_dir'] / "dna_files" / "protein_valid_anarci.fasta"
    raw_seqs_fasta = dirs['intermediate_dir'] / "raw_seqs.fasta"
    output_fasta = dirs['intermediate_dir'] / "verify_frames.fasta"
    
    # Проверяем входные файлы
    if not Path(protein_fasta).exists():
        print(f"✗ Protein FASTA not found: {protein_fasta}")
        return False
    
    if not Path(raw_seqs_fasta).exists():
        print(f"✗ Raw sequences FASTA not found: {raw_seqs_fasta}")
        print(f"✗ Run convert_ab1_multiframe.py first")
        return False
    
    print(f"\n{'='*60}")
    print(f"Step 1: Extracting verified sequence IDs from proteins")
    print(f" Source: {Path(protein_fasta).name}")
    print(f"{'='*60}")
    
    verified_ids = extract_verified_ids(protein_fasta)
    
    print(f"\n{'='*60}")
    print(f"Step 2: Extracting matching nucleotide sequences")
    print(f" Source: {Path(raw_seqs_fasta).name}")
    print(f" Output: {Path(output_fasta).name}")
    print(f"{'='*60}")
    
    count = extract_nucleotide_sequences(raw_seqs_fasta, verified_ids, output_fasta)
    
    if count == 0:
        print("✗ No matching sequences found")
        return False
    
    print(f"\n{'='*60}")
    print(f"✓ COMPLETE: {count} verified nucleotide sequences extracted")
    print(f" Output: {output_fasta}")
    print(f"{'='*60}\n")
    
    return True

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Extract verified nucleotide sequences from raw_seqs matching ANARCI results"
    )
    parser.add_argument("--config", default="config/config.yaml",
                        help="Path to config file")
    args = parser.parse_args()
    
    success = main(args.config)
    exit(0 if success else 1)
