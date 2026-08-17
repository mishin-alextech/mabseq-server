#!/usr/bin/env python3
import yaml
import subprocess
import csv
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
        
        'intermediate_dir': intermediate_dir,
        'results_dir': results_dir,
     
    }


def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def run_anarci(input_fasta, output_prefix, scheme="imgt"):
   
    cmd = [
        "ANARCI",
        "-i", str(input_fasta),
        "--scheme", scheme,
        "--csv",
        "--outfile", str(output_prefix)
    ]
    
    print(f"Running command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f" ANARCI completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f" ANARCI failed with error:")
        print(f"  stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f" ANARCI command not found in PATH")
        print(f"  ANARCI_DIR should be: $HOME/ANARCI-master")
        print(f"  Current PATH: {os.environ.get('PATH', 'NOT SET')}")
        return False

def parse_anarci_csv(anarci_output_prefix):
    
    sequences = {}
    
    # Ищем все CSV файлы с этим префиксом
    prefix_str = str(anarci_output_prefix)
    csv_files = sorted(Path(anarci_output_prefix.parent).glob(f"{anarci_output_prefix.name}*.csv"))
    

    
    for csv_file in csv_files:
        print(f"Parsing: {csv_file.name}")
        
        try:
            with open(csv_file, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                header = next(reader)  
                
                meta_cols = 13  
                
                
                for row in reader:
                    if len(row) < meta_cols:
                        continue
                    
                    seq_id = row[0].strip()
                    
                    aa_sequence = ""
                    for i in range(meta_cols, len(row)):
                        aa = row[i].strip()
                        
                        
                        if aa and aa != "-" and aa != ".":
                            aa_sequence += aa
                    
                    if aa_sequence:
                        sequences[seq_id] = aa_sequence
                        print(f"   {seq_id}: {len(aa_sequence)} aa")
                    else:
                        print(f"   {seq_id}: no sequence extracted")
        
        except Exception as e:
            print(f"   Error parsing {csv_file}: {e}")
    
    print(f" Parsed {len(sequences)} sequences from ANARCI output")
    return sequences

def write_output_fasta(sequences, output_fasta):
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    """
    Записывает последовательности в FASTA без пробелов и пропусков.
    """
    with open(output_fasta, 'w') as f:
        for seq_id, sequence in sequences.items():
            # Очищаем от всех пробелов, дефисов и точек (на случай если что-то проскочило)
            clean_seq = sequence.replace(" ", "").replace("-", "").replace(".", "")
            f.write(f">{seq_id}\n{clean_seq}\n")
    
    print(f" Wrote {len(sequences)} sequences to {output_fasta}")

def main(config_file="config/config.yaml"):
    config = load_config(config_file)
    dirs = get_data_dirs()
    
    input_fasta = dirs['intermediate_dir'] / "translated_frames.fasta"
    output_fasta = dirs['results_dir'] / "dna_files" / "protein_valid_anarci.fasta"
    anarci_temp_prefix = dirs['intermediate_dir'] / "anarci_temp"
    
    anarci_config = config["anarci"]
    scheme = anarci_config["scheme"]
    min_aa = anarci_config["min_aa"]
    max_aa = anarci_config["max_aa"]
    
    # Проверяем входной файл
    if not Path(input_fasta).exists():
        print(f" Input file not found: {input_fasta}")
        return False
    
    print(f"\n{'='*60}")
    print(f"Step 1: Running ANARCI on {input_fasta}")
    print(f"{'='*60}")
    
    if not run_anarci(input_fasta, anarci_temp_prefix, scheme):
        return False
    
    print(f"\n{'='*60}")
    print(f"Step 2: Parsing ANARCI output")
    print(f"{'='*60}")
    
    sequences = parse_anarci_csv(anarci_temp_prefix)
    
    if not sequences:
        print(" No sequences parsed from ANARCI output")
        return False
    
    print(f"\n{'='*60}")
    print(f"Step 3: Filtering by length ({min_aa}-{max_aa} aa)")
    print(f"{'='*60}")
    
    filtered_sequences = {}
    for seq_id, seq in sequences.items():
        length = len(seq)
        
        if length < min_aa or length > max_aa:
            print(f"   {seq_id}: length {length} aa (outside range, skipped)")
            continue
        
        print(f"   {seq_id}: {length} aa")
        filtered_sequences[seq_id] = seq
    
    if not filtered_sequences:
        print(" No sequences passed length filter")
        return False
    
    print(f" Filtered to {len(filtered_sequences)} valid sequences")
    
    print(f"\n{'='*60}")
    print(f"Step 4: Writing output FASTA")
    print(f"{'='*60}")
    
    write_output_fasta(filtered_sequences, output_fasta)
    
    # Очистка временных файлов
    temp_dir = anarci_temp_prefix.parent
    temp_pattern = f"{anarci_temp_prefix.name}*"
    for temp_file in temp_dir.glob(temp_pattern):
        temp_file.unlink(missing_ok=True)
    
    print(f"\n{'='*60}")
    print(f"✓ COMPLETE: {len(filtered_sequences)} sequences written to {output_fasta}")
    print(f"{'='*60}\n")
    
    return True

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Run ANARCI and parse output to FASTA"
    )
    parser.add_argument("--config", default="config/config.yaml",
                        help="Path to config file (default: config/config.yaml)")
    
    args = parser.parse_args()
    
    success = main(args.config)
    
    exit(0 if success else 1)

