#!/usr/bin/env python3

"""
Convert .ab1 files to FASTA with 3-frame translation.
Assumes IDs are already normalized by normalize_ids.py
Uses normalized file names as sequence IDs.

Input:
- .ab1 files from data/raw/sanger with normalized names

Output:
- data/intermediate/translated_frames.fasta — protein sequences, 3 frames
- data/intermediate/raw_seqs.fasta — original nucleotide sequences
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
    raw_data_dir = Path(os.environ.get('RAW_DATA_DIR', data_dir / 'raw' / 'sanger'))
    intermediate_dir = Path(os.environ.get('INTERMEDIATE_DIR', data_dir / 'intermediate'))
    
    return {
        'data_dir': data_dir,
        'raw_data_dir': raw_data_dir,
        'intermediate_dir': intermediate_dir
    }

def translate_three_frames(normalized_id, seq, forward=True):
    """
    Translate nucleotide sequence into 3 reading frames.
    
    Args:
    - normalized_id: The normalized sequence ID from file name
    - seq: BioPython Seq object with nucleotide sequence
    - forward: If True, use forward strand; if False, use reverse complement
    
    Returns:
    - list of tuples (frameid, aasequence)
    """
    if not forward:
        seq = seq.reverse_complement()
    
    frames = []
    
    for frame in range(3):
        subseq = seq[frame:]
        aa_seq = subseq.translate(to_stop=False)
        
        frame_id = frame + 1
        strand_label = 'fwd' if forward else 'rev'
        
        new_id = f"{normalized_id}_{strand_label}_frame{frame_id}"
        frames.append((new_id, aa_seq))
    
    return frames

def main(config_file="config/config.yaml"):
    dirs = get_data_dirs()
    config = load_config(config_file)
    
    input_pattern = config["paths"]["ab1_pattern"]
    
    # Используем переменные окружения для путей
    output_fasta = dirs['intermediate_dir'] / "translated_frames.fasta"
    raw_seqs_fasta = dirs['intermediate_dir'] / "raw_seqs.fasta"
    forward_line = config['ab1_conversion']['forward_line']
    
    # Создать директории если их нет
    dirs['intermediate_dir'].mkdir(parents=True, exist_ok=True)
    
    # Найти входные .ab1 файлы
    input_files = sorted(Path(dirs['raw_data_dir']).glob('*.ab1'))
    
    if not input_files:
        raise SystemExit(f"No input .ab1 files found in {dirs['raw_data_dir']}")
    
    strand_label = 'forward' if forward_line else 'reverse complement'
    
    print("="*80)
    print(f"AB1 to FASTA Conversion and 3-Frame Translation")
    print("="*80)
    print(f"Found {len(input_files)} .ab1 files")
    print(f"Processing strand: {strand_label}")
    print(f"{'#':<4} {'File Name':<40} {'Normalized ID':<30} {'Seq Length':<15}")
    print("-"*89)
    
    frames_records = []
    raw_records = []
    
    for idx, ab1_file in enumerate(input_files, 1):
        # Get normalized ID from FILE NAME (not from internal ID)
        normalized_id = ab1_file.stem  # Remove .ab1 extension
        
        # Read .ab1 file to get sequence
        record = SeqIO.read(ab1_file, 'abi')
        seq = record.seq
        
        print(f"{idx:<4} {ab1_file.name:<40} {normalized_id:<30} {len(seq):<15} nt")
        
        # Extract nucleotide sequence
        orig_seq = seq
        if not forward_line:
            orig_seq = orig_seq.reverse_complement()
        
        # Create raw nucleotide record with NORMALIZED ID
        raw_record = SeqIO.SeqRecord(orig_seq, id=normalized_id, description='')
        raw_records.append(raw_record)
        
        # Generate 3-frame translations using NORMALIZED ID
        frames = translate_three_frames(normalized_id, seq, forward=forward_line)
        for frame_id, aa_seq in frames:
            rec = SeqIO.SeqRecord(aa_seq, id=frame_id, description='')
            frames_records.append(rec)
    
    print("="*80)
    print("Writing output files...")
    print("="*80)
    
    # Write raw nucleotide sequences
    count_raw = SeqIO.write(raw_records, raw_seqs_fasta, 'fasta')
    print(f"✓ Raw nucleotide sequences: {count_raw} sequences")
    print(f"  File: {raw_seqs_fasta}")
    
    # Write translated protein sequences
    count_frames = SeqIO.write(frames_records, output_fasta, 'fasta')
    expected = len(raw_records) * 3
    
    print(f"✓ Translated protein sequences: {count_frames} sequences")
    print(f"  File: {output_fasta}")
    print(f"  Expected: {len(raw_records)} files × 3 frames = {expected}")
    
    if count_frames == expected:
        print(f"✓ All frames successfully translated!")
    else:
        print(f"✗ Warning: expected {expected}, got {count_frames}")
    
    print("="*80)
    print("✓ COMPLETE")
    print(f"  Raw sequences: {raw_seqs_fasta}")
    print(f"  Protein frames: {output_fasta}")
    print("="*80)
    
    return True

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Convert .ab1 files to multi-frame translated FASTA"
    )
    parser.add_argument("--config", default="config/config.yaml",
                        help="Path to config file")
    args = parser.parse_args()
    
    try:
        main(args.config)
    except Exception as e:
        print(f"✗ Error: {e}")
        exit(1)
