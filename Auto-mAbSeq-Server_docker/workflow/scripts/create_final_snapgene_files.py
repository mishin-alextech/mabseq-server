#!/usr/bin/env python3
"""
Step 8: Create Final SnapGene Files
Merges Protein Sequence, Nucleotide Sequence, and Trace Data into a single .dna file.
Uses autosnapgene library and staden io_lib for complete integration.
"""

import yaml
import argparse
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import autosnapgene as snap
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
        'data_dir': data_dir,
        'raw_data_dir': raw_data_dir,
        'intermediate_dir': intermediate_dir,
        'results_dir': results_dir,
        'logs_dir': logs_dir
    }


def load_config(config_file):
    
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def find_cds_match(dna_seq_str, protein_seq_str):
    """
    Translates DNA in 3 forward frames to find the Protein sequence.
    Returns: (dna_substring, frame_number) or (None, None)
    
    dna_substring includes one nucleotide prefix for proper positioning.
    
    Example:
        find_cds_match("ATGGCTAGC...", "MASR...") → ("NATGGCTAGC...", 1)
    """
    dna_seq = Seq(dna_seq_str)
    
    # Clean protein sequence (remove stop * if present at end)
    clean_protein = protein_seq_str.rstrip('*')
    
    # Try each of 3 forward reading frames
    for frame in range(3):
        translated = dna_seq[frame:].translate(table=1)
        translated_str = str(translated)
        
        # Find protein in translation
        idx = translated_str.find(clean_protein)
        
        if idx != -1:
            # Calculate start and end in DNA coordinates
            start = frame + (idx * 3)
            end = start + (len(clean_protein) * 3)
            
            # Add one nucleotide prefix for proper positioning in add_feature
            start_with_prefix = max(0, start - 1)
            dna_substring = dna_seq_str[start_with_prefix:end]
            return dna_substring, frame + 1
            
    return None, None


def parse_feature_name_from_id(seq_id):
    """
    Parse feature name from sequence ID.
    
    Example:
        VHH_cDNA_Cm_TN95089 → VHH_TN95089
        VHeavy_cDNA_Ov_lib#1_ABC123 → VHeavy_lib#1_ABC123
    
    Logic: Remove "_cDNA" and the animal code (single component after _cDNA)
    """
    parts = seq_id.split('_')
    
    # Format: [chain, 'cDNA', animal, ...rest]
   
    
    if len(parts) >= 4 and parts[1] == 'cDNA':
        chain = parts[0]
        sample_parts = parts[3:]
        sample = '_'.join(sample_parts)
        feature_name = f"{chain}_{sample}"
        return feature_name
    else:
        return seq_id

def main():
    parser = argparse.ArgumentParser(
        description='Create Final SnapGene Files with integrated traces and CDS annotations'
    )
    parser.add_argument('--config', required=True, help='Path to config.yaml')
    parser.add_argument('--dry-run', action='store_true', help='Perform dry run without writing')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    dirs = get_data_dirs()

    protein_fasta = dirs['results_dir'] / 'dna_files' / 'protein_valid_normalized.fasta'
    dna_input_dir = dirs['results_dir'] / 'dna_files'
    ab1_input_dir = dirs['raw_data_dir']
    output_dir = dirs['results_dir'] / 'final_snapgene'
    log_dir = dirs['logs_dir']
    
    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    
    print(f"\n Starting Final SnapGene Assembly (Step 8)...")
    print(f"   Proteins source: {protein_fasta}")
    print(f"   DNA source:      {dna_input_dir}")
    print(f"   Traces source:   {ab1_input_dir}")
    print(f"   Output dir:      {output_dir}")

    # Validate inputs
    if not protein_fasta.exists():
        print(f" Error: Protein file not found: {protein_fasta}")
        return False

    # Load Proteins
    proteins = {rec.id: str(rec.seq) for rec in SeqIO.parse(protein_fasta, 'fasta')}
    print(f"   Loaded {len(proteins)} protein sequences.\n")
    
    success_count = 0
    failed_count = 0
    assembly_log = []
    
    # Process each protein
    for seq_id, protein_seq in proteins.items():
        # Find matching DNA file
        dna_file = dna_input_dir / f"{seq_id}.dna"
        ab1_file = ab1_input_dir / f"{seq_id}.ab1"
        
        if not dna_file.exists():
            print(f"  [{seq_id}] Missing DNA file, skipping.")
            print(f"  [{seq_id}] Missing DNA file: {dna_file}")
            assembly_log.append(f"{seq_id} -> FAILED: Missing DNA file")
            failed_count += 1
            continue
            
        # Read raw DNA sequence
        try:
            with open(dna_file, 'r') as f:
                content = f.read().strip()
                if content.startswith('>'):
                    # Simple FASTA parsing
                    dna_seq = ''.join(content.split('\n')[1:])
                else:
                    dna_seq = content
            # Remove whitespace
            dna_seq = ''.join(dna_seq.split())
        except Exception as e:
            print(f" [{seq_id}] Error reading DNA: {e}")
            assembly_log.append(f"{seq_id} -> FAILED: {str(e)}")
            failed_count += 1
            continue

        # Find CDS match (translate in 3 frames)
        dna_substring, frame = find_cds_match(dna_seq, protein_seq)
        
        if not dna_substring:
            print(f" [{seq_id}] CDS match failed in all frames!")
            assembly_log.append(f"{seq_id} -> FAILED: CDS not found")
            failed_count += 1
            continue
            
        # Create SnapGene object
        try:
            snap_dna = snap.SnapGene()
            snap_dna.sequence = dna_seq
            snap_dna.topology = 'linear'
            snap_dna.strandedness = 'double'
        except Exception as e:
            print(f" [{seq_id}] Error creating SnapGene object: {e}")
            assembly_log.append(f"{seq_id} -> FAILED: {str(e)}")
            failed_count += 1
            continue

        # Parse feature name from ID
        feature_name = parse_feature_name_from_id(seq_id)

        # Create CDS Feature WITHOUT segment (will be set by add_feature)
        try:
            cds_feat = snap.Feature()
            cds_feat.name = feature_name
            cds_feat.type = 'CDS'
            cds_feat.qualifiers = {
                'translation': protein_seq,
                'note': f"Predicted VHH coding sequence (Frame {frame})"
            }
            
            # Add feature to sequence with explicit dna_substring
            # add_feature will find dna_substring and set coordinates automatically
            snap_dna.add_feature(cds_feat, dna_substring)
            
        except Exception as e:
            print(f" [{seq_id}] Error adding CDS feature: {e}")
            assembly_log.append(f"{seq_id} -> FAILED: {str(e)}")
            failed_count += 1
            continue

        # Embed Trace if available
        trace_status = "no trace"
        if ab1_file.exists():
            try:
                snap_dna.add_trace(str(ab1_file))
                trace_status = "trace added"
            except Exception as e:
                print(f"    [{seq_id}] Warning: Could not embed trace: {e}")
                trace_status = "trace failed"
        
        # Write Final File
        if not args.dry_run:
            try:
                out_path = output_dir / f"{seq_id}.dna"
                snap_dna.write(str(out_path))
                success_count += 1
                assembly_log.append(f"{seq_id} -> SUCCESS (Frame {frame}, {trace_status})")
                print(f" [{seq_id}] Assembled (Frame {frame}, {trace_status})")
            except Exception as e:
                print(f" [{seq_id}] Error writing file: {e}")
                assembly_log.append(f"{seq_id} -> FAILED: {str(e)}")
                failed_count += 1
        else:
            success_count += 1
            assembly_log.append(f"{seq_id} -> DRY RUN (Frame {frame}, {trace_status})")
            print(f" [{seq_id}] DRY RUN (Frame {frame}, {trace_status})")

    # Write assembly log
    log_file = log_dir / 'snapgene_assembly.log'
    with open(log_file, 'w') as f:
        for log_entry in assembly_log:
            f.write(log_entry + '\n')

    print(f"\n📊 Summary:")
    print(f"   Success:  {success_count}")
    print(f"   Failed:   {failed_count}")
    print(f"   Total:    {len(proteins)}")
    print(f"📝 Log file: {log_file}")
    print(f"\n Step 8 completed!")
    
    return success_count > 0

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
