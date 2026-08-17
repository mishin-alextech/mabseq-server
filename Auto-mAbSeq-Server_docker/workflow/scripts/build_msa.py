#!/usr/bin/env python3
"""
Build multiple sequence alignment (MSA) for verified VHH sequences.
Groups sequences by chain type and creates separate alignments.

On input:
  - vhh_valid_anarci.fasta — verified protein sequences
    (may contain VHH, VHeavy, VKappa, VLambda)

On output:
  - results/alignment/[CHAIN]_aligned.fasta — aligned sequences per chain
  - results/alignment/[CHAIN]_alignment.clustal — Clustal format per chain
  - results/alignment/alignment_summary.txt — summary of all alignments
"""

import yaml
import re
import subprocess
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict
import os


def get_data_dirs():
    """
    Получить директории из переменных окружения или использовать defaults.
    Это позволяет изолировать каждый job в отдельную директорию.
    """
    data_dir = Path(os.environ.get('DATA_DIR', 'data'))
    results_dir = Path(os.environ.get('RESULTS_DIR', data_dir / 'results'))
    logs_dir = Path(os.environ.get('LOGS_DIR', data_dir / 'logs'))
    
    return {
        'results_dir': results_dir,
        'logs_dir': logs_dir,
    }


def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def extract_chain_type(sequence_id):
    """
    Extract chain type from sequence ID.
    
    Examples:
      VHH_raw_cDNA_Cm_TN95089 → VHH
      VHeavy_raw_cDNA_Cm_TN95077 → VHeavy
      VKappa_raw_cDNA_Ms_TN2201 → VKappa
      VLambda_raw_cDNA_Rb_TN3101 → VLambda
   
    """
    # Match chain types at the beginning
    match = re.match(r'^(VHH|VHeavy|VKappa|VLambda)', sequence_id)
    if match:
        return match.group(1)
    return None

def group_sequences_by_chain(protein_fasta):
    """
    Group sequences by chain type.
    
    Returns:
        dict {chain_type: [SeqRecord, ...]}
    """
    chain_groups = defaultdict(list)
    unclassified = []
    
    for record in SeqIO.parse(protein_fasta, "fasta"):
        chain_type = extract_chain_type(record.id)
        
        if chain_type:
            chain_groups[chain_type].append(record)
        else:
            unclassified.append(record.id)
    
    if unclassified:
        print(f"\n Warning: {len(unclassified)} sequences could not be classified:")
        for seq_id in unclassified[:5]:
            print(f"    - {seq_id}")
    
    return dict(chain_groups)

def run_clustalo(protein_fasta, output_fasta, output_clustal, num_threads=4):
    """Run Clustal Omega for alignment"""
    cmd = [
        "clustalo",
        "-i", protein_fasta,
        "-o", output_fasta,
        "--outfmt=fasta",
        "--verbose",
        "--force",
        "--threads", str(num_threads)
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"    ✓ FASTA alignment: {Path(output_fasta).name}")
        
        # Also create Clustal format
        cmd_clustal = [
            "clustalo",
            "-i", protein_fasta,
            "-o", output_clustal,
            "--outfmt=clustal",
            "--force",
            "--threads", str(num_threads)
        ]
        
        result_clustal = subprocess.run(cmd_clustal, capture_output=True, text=True, check=True)
        print(f"     Clustal format: {Path(output_clustal).name}")
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"     Clustal Omega failed: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"     clustalo command not found")
        return False

def get_alignment_stats(aligned_fasta):
    """Get alignment statistics"""
    records = list(SeqIO.parse(aligned_fasta, "fasta"))
    
    if not records:
        return None
    
    num_seqs = len(records)
    alignment_length = len(records[0].seq)
    
    gaps_per_seq = {}
    for record in records:
        gap_count = str(record.seq).count("-")
        gaps_per_seq[record.id] = gap_count
    
    avg_gaps = sum(gaps_per_seq.values()) / len(gaps_per_seq) if gaps_per_seq else 0
    
    return {
        "num_sequences": num_seqs,
        "alignment_length": alignment_length,
        "avg_gaps": avg_gaps,
        "min_gaps": min(gaps_per_seq.values()) if gaps_per_seq else 0,
        "max_gaps": max(gaps_per_seq.values()) if gaps_per_seq else 0
    }

def main(config_file="config/config.yaml"):
    config = load_config(args.config)
    dirs = get_data_dirs()
    
    protein_fasta = dirs['results_dir'] / 'dna_files' / 'protein_valid_normalized.fasta'
    output_base_dir = dirs['results_dir'] / 'alignment' 
    log_dir = dirs['logs_dir']
    
    alignment_config = config.get("alignment", {})
    tool = alignment_config.get("tool", "clustalo")
    num_threads = alignment_config.get("num_threads", 4)
    group_by_chain = alignment_config.get("group_by_chain", True)
    
    # Check input file
    if not Path(protein_fasta).exists():
        print(f"✗ Input FASTA not found: {protein_fasta}")
        return False
    
    # Read and group sequences
    print(f"\n{'='*80}")
    print(f"Building Multiple Sequence Alignments (MSA)")
    print(f"{'='*80}")
    print(f"Tool: {tool}")
    print(f"Group by chain: {group_by_chain}")
    print(f"Threads: {num_threads}")
    print(f"Input: {Path(protein_fasta).name}\n")
    
    chain_groups = group_sequences_by_chain(protein_fasta)
    
    if not chain_groups:
        print(f" No sequences found in input file")
        return False
    
    print(f"Found {sum(len(v) for v in chain_groups.values())} sequences")
    print(f"Chain types: {', '.join(sorted(chain_groups.keys()))}\n")
    
    # Create output directory
    output_base_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each chain type
    alignment_results = {}
    
    print(f"{'='*80}")
    print(f"Running alignments...")
    print(f"{'='*80}\n")
    
    for chain_type in sorted(chain_groups.keys()):
        records = chain_groups[chain_type]
        num_seqs = len(records)
        
        print(f"{chain_type}:")
        print(f"  Sequences: {num_seqs}")
        
        # Need at least 2 sequences for alignment
        if num_seqs < 2:
            print(f"   Skipped: need at least 2 sequences")
            continue
        
        # Create temporary FASTA for this chain
        temp_fasta = output_base_dir / f"temp_{chain_type}.fasta"
        SeqIO.write(records, str(temp_fasta), "fasta")
        
        # Create output file names
        output_fasta = output_base_dir / f"{chain_type}_aligned.fasta"
        output_clustal = output_base_dir / f"{chain_type}_alignment.clustal"
        
        # Run alignment
        if tool == "clustalo":
            success = run_clustalo(str(temp_fasta), str(output_fasta), str(output_clustal), num_threads)
        else:
            print(f"   Unknown tool: {tool}")
            success = False
        
        # Get statistics
        if success:
            stats = get_alignment_stats(str(output_fasta))
            if stats:
                print(f"  Alignment length: {stats['alignment_length']} aa")
                print(f"  Avg gaps: {stats['avg_gaps']:.1f}")
                alignment_results[chain_type] = {
                    "fasta": str(output_fasta),
                    "clustal": str(output_clustal),
                    "stats": stats,
                    "success": True
                }
        else:
            alignment_results[chain_type] = {"success": False}
        
        # Clean up temporary file
        temp_fasta.unlink(missing_ok=True)
        print()
    
    # Write summary report
    summary_file = output_base_dir / "alignment_summary.txt"
    write_summary_report(alignment_results, summary_file)
    
    # Final output
    print(f"{'='*80}")
    print(f"✓ COMPLETE")
    print(f"{'='*80}")
    
    successful = sum(1 for r in alignment_results.values() if r.get("success", False))
    print(f"\nAlignments created: {successful}/{len(alignment_results)}")
    
    for chain_type in sorted(alignment_results.keys()):
        if alignment_results[chain_type].get("success", False):
            stats = alignment_results[chain_type]["stats"]
            print(f"  {chain_type:<10} {stats['num_sequences']:>2} seqs, {stats['alignment_length']:>3} aa")
    
    print(f"\nOutput directory: {output_base_dir}")
    print(f"Summary report: {summary_file}\n")
    
    return True

def write_summary_report(alignment_results, summary_file):
    """Write summary report of all alignments"""
    with open(summary_file, 'w') as f:
        f.write("VHH Pipeline - Multiple Sequence Alignment Summary\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("Alignments by Chain Type:\n")
        f.write("-" * 80 + "\n\n")
        
        for chain_type in sorted(alignment_results.keys()):
            result = alignment_results[chain_type]
            
            f.write(f"Chain Type: {chain_type}\n")
            
            if result.get("success", False):
                stats = result["stats"]
                f.write(f"  Status: ✓ Completed\n")
                f.write(f"  Number of sequences: {stats['num_sequences']}\n")
                f.write(f"  Alignment length: {stats['alignment_length']} aa\n")
                f.write(f"  Average gaps per sequence: {stats['avg_gaps']:.1f}\n")
                f.write(f"  Min gaps: {stats['min_gaps']}\n")
                f.write(f"  Max gaps: {stats['max_gaps']}\n")
                f.write(f"  Output files:\n")
                f.write(f"    - {Path(result['fasta']).name}\n")
                f.write(f"    - {Path(result['clustal']).name}\n")
            else:
                f.write(f"  Status: ✗ Failed or skipped\n")
            
            f.write("\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Build per-chain multiple sequence alignments"
    )
    parser.add_argument("--config", default="config/config.yaml",
                        help="Path to config file")
    
    args = parser.parse_args()
    success = main(args.config)
    exit(0 if success else 1)

