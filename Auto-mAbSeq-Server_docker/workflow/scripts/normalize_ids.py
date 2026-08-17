#!/usr/bin/env python3
"""
Normalize .ab1 file IDs and rename files accordingly.
Creates standardized file names for downstream processing.

On input:
  - .ab1 files in data/raw/sanger/ with various naming formats
  
On output:
  - Renamed .ab1 files with standardized IDs
  - ID mapping log for reference
  - Detailed processing logs
"""

import yaml
import re
import shutil
from pathlib import Path
from datetime import datetime
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
        'raw_data_dir': raw_data_dir,
        'logs_dir': logs_dir
    }


def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def setup_logging(log_dir):

    log_path = Path(log_dir)
    log_path.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    return {
        "mapping": log_path / "id_mapping.txt",
        "detailed": log_path / f"id_processing_detailed_{timestamp}.log",
        "summary": log_path / f"id_processing_summary_{timestamp}.log",
        "errors": log_path / f"id_processing_errors_{timestamp}.log"
    }

class IDProcessor:
    def __init__(self, config):
        self.config = config
        self.id_config = config.get("id_processing", {})
        self.results = []
        self.stats = {
            "total": 0,
            "success": 0,
            "failed": 0,
            "patterns_used": {},
            "renamed": 0,
            "skipped": 0
        }
    
    def parse_id_with_regex(self, original_id):
        """Parse ID using regex patterns"""
        regex_patterns = self.id_config.get("regex_patterns", [])
        
        for pattern_config in regex_patterns:
            pattern_name = pattern_config["name"]
            pattern = pattern_config["pattern"]
            groups = pattern_config["groups"]
            
            match = re.search(pattern, original_id)
            if match:
                result = {"pattern_used": pattern_name}
                for key, group_num in groups.items():
                    result[key] = match.group(group_num)
                
                return result, pattern_name, True
        
        return None, None, False
    
    def validate_parsed_id(self, parsed):
        """Validate parsed ID against config"""
        if not parsed:
            return False, "Failed to parse ID"
        
        chain = parsed.get("chain")
        animal = parsed.get("animal")
        sample = parsed.get("sample")
        
        valid_chains = self.id_config.get("valid_chains", [])
        valid_animals = self.id_config.get("valid_animals", [])
        
        if chain not in valid_chains:
            return False, f"Invalid chain '{chain}'. Valid: {valid_chains}"
        
        if animal not in valid_animals:
            return False, f"Invalid animal '{animal}'. Valid: {valid_animals}"
        
        if not sample:
            return False, "Sample/library code not found"
        
        return True, "Validation passed"
    
    def normalize_id(self, original_id):
        """Normalize ID with full logging"""
        self.stats["total"] += 1
        
        details = {
            "original": original_id,
            "timestamp": datetime.now().isoformat(),
            "steps": []
        }
        
        # Step 1: Regex parsing
        parsed, pattern_used, parse_success = self.parse_id_with_regex(original_id)
        
        details["steps"].append({
            "step": 1,
            "name": "Regex Parsing",
            "pattern": pattern_used,
            "success": parse_success,
            "result": parsed
        })
        
        if not parse_success:
            details["error"] = "Failed to match any regex pattern"
            self.stats["failed"] += 1
            return None, False, details
        
        self.stats["patterns_used"][pattern_used] = self.stats["patterns_used"].get(pattern_used, 0) + 1
        
        # Step 2: Validation
        is_valid, validation_msg = self.validate_parsed_id(parsed)
        
        details["steps"].append({
            "step": 2,
            "name": "Validation",
            "valid": is_valid,
            "message": validation_msg,
            "parsed_components": {
                "chain": parsed.get("chain"),
                "animal": parsed.get("animal"),
                "sample": parsed.get("sample")
            }
        })
        
        if not is_valid:
            if self.id_config.get("strict_mode", False):
                details["error"] = validation_msg
                self.stats["failed"] += 1
                return None, False, details
            else:
                details["warning"] = f"Validation failed ({validation_msg}), using original ID"
                details["steps"].append({
                    "step": 3,
                    "name": "Format Output",
                    "mode": "fallback",
                    "normalized": original_id
                })
                self.stats["success"] += 1
                self.stats["skipped"] += 1
                return original_id, False, details
        
        # Step 3: Format output
        output_format = self.id_config.get("output_format", "{chain}_raw_cDNA_{animal}_{sample}")
        normalized = output_format.format(**parsed)
        
        details["steps"].append({
            "step": 3,
            "name": "Format Output",
            "format_template": output_format,
            "normalized": normalized
        })
        
        details["success"] = True
        self.stats["success"] += 1
        
        return normalized, True, details
    
    def process_ab1_files(self, raw_data_dir, dry_run=False):
    
        ab1_files = sorted(raw_data_dir.glob("*.ab1"))
        
        if not ab1_files:
            print(f"No .ab1 files found matching {input_pattern}")
            return []
        
        print(f"\n{'='*100}")
        print(f"Processing and Renaming {len(ab1_files)} .ab1 Files")
        print(f"{'='*100}")
        print(f"Dry run: {dry_run}")
        print(f"Output directory: {ab1_files[0].parent}\n")
        print(f"{'#':<4} {'Original Name':<50} {'→':<2} {'New Name':<35} {'Status':<10}")
        print(f"{'-'*101}")
        
        for idx, ab1_file in enumerate(ab1_files, 1):
            original_id = ab1_file.stem
            normalized_id, success, details = self.normalize_id(original_id)
            
            status = "✓ OK" if success else "⚠ WARN"
            new_filename = f"{(normalized_id or 'FAILED')}.ab1" if normalized_id else "FAILED.ab1"
            
            print(f"{idx:<4} {ab1_file.name:<50} → {new_filename:<35} {status:<10}")
            
            if not success and "error" in details:
                print(f"     └─ ERROR: {details['error']}")
                self.stats["failed"] += 1
            elif not success and "warning" in details:
                print(f"     └─ WARNING: {details['warning']}")
            
            result = {
                "original_file": ab1_file,
                "original_id": original_id,
                "normalized_id": normalized_id,
                "new_filename": new_filename,
                "success": success,
                "details": details,
                "renamed": False
            }
            
            # Perform rename if not dry run
            if success and not dry_run:
                try:
                    new_filepath = ab1_file.parent / new_filename
                    shutil.move(str(ab1_file), str(new_filepath))
                    result["renamed"] = True
                    self.stats["renamed"] += 1
                except Exception as e:
                    print(f"     └─ RENAME ERROR: {e}")
                    result["rename_error"] = str(e)
            
            self.results.append(result)
        
        return self.results
    
    def write_mapping_log(self, log_file):
        """Write simple mapping log (original → normalized)"""
        with open(log_file, 'w') as f:
            f.write("VHH Pipeline - .ab1 File Renaming Log\n")
            f.write("=" * 90 + "\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n")
            f.write(f"Total processed: {self.stats['total']}\n")
            f.write(f"Successfully renamed: {self.stats['renamed']}\n")
            f.write(f"Skipped: {self.stats['skipped']}\n")
            f.write(f"Failed: {self.stats['failed']}\n")
            f.write("\n" + "=" * 90 + "\n\n")
            
            f.write(f"{'Original File Name':<55} → {'New File Name':<40}\n")
            f.write("-" * 97 + "\n")
            
            for result in self.results:
                original = result["original_file"].name
                new = result["new_filename"]
                status = "✓" if result["renamed"] else "✗"
                f.write(f"{status} {original:<53} → {new:<38}\n")
    
    def write_detailed_log(self, log_file):
        """Write detailed processing log"""
        with open(log_file, 'w') as f:
            f.write("VHH Pipeline - Detailed File Renaming Log\n")
            f.write("=" * 100 + "\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n")
            f.write(f"Configuration:\n")
            f.write(f"  Strict mode: {self.id_config.get('strict_mode', False)}\n")
            f.write(f"  Output format: {self.id_config.get('output_format', 'N/A')}\n")
            f.write(f"  Number of regex patterns: {len(self.id_config.get('regex_patterns', []))}\n")
            f.write("\n" + "=" * 100 + "\n\n")
            
            for idx, result in enumerate(self.results, 1):
                details = result["details"]
                
                f.write(f"Entry #{idx}\n")
                f.write(f"{'─' * 100}\n")
                f.write(f"Original file:  {result['original_file'].name}\n")
                f.write(f"Original ID:    {result['original_id']}\n")
                f.write(f"Normalized ID:  {result['normalized_id'] or 'FAILED'}\n")
                f.write(f"New file name:  {result['new_filename']}\n")
                f.write(f"Renamed:        {result['renamed']}\n")
                f.write(f"Timestamp:      {details['timestamp']}\n")
                f.write(f"\nProcessing Steps:\n")
                
                for step in details.get("steps", []):
                    f.write(f"\n  Step {step['step']}: {step['name']}\n")
                    
                    if "pattern" in step:
                        f.write(f"    Pattern used: {step['pattern']}\n")
                    
                    if "parsed_components" in step:
                        f.write(f"    Parsed components:\n")
                        for key, val in step["parsed_components"].items():
                            f.write(f"      - {key}: {val}\n")
                    
                    if "format_template" in step:
                        f.write(f"    Format template: {step['format_template']}\n")
                    
                    if "normalized" in step:
                        f.write(f"    Result: {step['normalized']}\n")
                
                if "error" in details:
                    f.write(f"\nError:\n  {details['error']}\n")
                
                if "warning" in details:
                    f.write(f"\nWarning:\n  {details['warning']}\n")
                
                f.write(f"\n{'─' * 100}\n\n")
    
    def write_summary_log(self, log_file):
        """Write summary report"""
        with open(log_file, 'w') as f:
            f.write("VHH Pipeline - .ab1 File Renaming Summary Report\n")
            f.write("=" * 100 + "\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n\n")
            
            f.write("Statistics:\n")
            f.write(f"  Total files processed:        {self.stats['total']}\n")
            f.write(f"  Successfully renamed:        {self.stats['renamed']}\n")
            f.write(f"  Skipped (validation warn):   {self.stats['skipped']}\n")
            f.write(f"  Failed:                      {self.stats['failed']}\n")
            f.write(f"  Success rate:                {(self.stats['renamed']/self.stats['total']*100):.1f}%\n")
            
            f.write(f"\nRegex Patterns Used:\n")
            for pattern_name, count in self.stats["patterns_used"].items():
                f.write(f"  {pattern_name:<30} {count:>3} files\n")
            
            f.write(f"\nConfiguration Used:\n")
            f.write(f"  Strict mode:                 {self.id_config.get('strict_mode', False)}\n")
            f.write(f"  Output format:               {self.id_config.get('output_format', 'N/A')}\n")
            f.write(f"  Valid chains:                {', '.join(self.id_config.get('valid_chains', []))}\n")
            f.write(f"  Valid animals:               {', '.join(self.id_config.get('valid_animals', []))}\n")
            
            f.write(f"\nProcessed Files:\n")
            f.write(f"{'-' * 100}\n")
            f.write(f"{'Original':<55} {'→':<2} {'New Name':<40}\n")
            f.write(f"{'-' * 100}\n")
            
            for result in self.results:
                status = "✓" if result["renamed"] else "✗"
                f.write(f"{status} {result['original_file'].name:<53} → {result['new_filename']:<38}\n")
    
    def write_errors_log(self, log_file):
        """Write errors and warnings log"""
        errors = [r for r in self.results if not r["renamed"]]
        
        if not errors:
            with open(log_file, 'w') as f:
                f.write("No errors or warnings found.\n")
            return
        
        with open(log_file, 'w') as f:
            f.write("VHH Pipeline - File Renaming Errors and Warnings\n")
            f.write("=" * 100 + "\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n")
            f.write(f"Total errors/warnings: {len(errors)}\n\n")
            
            for idx, result in enumerate(errors, 1):
                details = result["details"]
                f.write(f"Entry #{idx}: {result['original_file'].name}\n")
                f.write(f"{'─' * 100}\n")
                
                if "error" in details:
                    f.write(f"ERROR: {details['error']}\n")
                
                if "warning" in details:
                    f.write(f"WARNING: {details['warning']}\n")
                
                if "rename_error" in result:
                    f.write(f"RENAME ERROR: {result['rename_error']}\n")
                
                f.write(f"\nAttempted transformation:\n")
                f.write(f"  {result['original_file'].name} → {result['new_filename']}\n")
                f.write("\n")

def main(config_file="config/config.yaml", dry_run=False):
    config = load_config(config_file)
    dirs = get_data_dirs()
    
    log_dir = dirs['logs_dir']  
    ab1_files = dirs['raw_data_dir'].glob("*.ab1")
    
    log_files = setup_logging(log_dir)
    
    print(f"\n.ab1 File Normalization and Renaming")
    print(f"  Method: {config['id_processing'].get('method', 'hybrid')}")
    print(f"  Strict mode: {config['id_processing'].get('strict_mode', False)}")
    print(f"  Output format: {config['id_processing'].get('output_format', 'N/A')}")
    print(f"  Log directory: {log_dir}")
    
    processor = IDProcessor(config)
    processor.process_ab1_files(dirs['raw_data_dir'], dry_run)
    
    if config.get("id_processing", {}).get("log_transformations", False):
        print(f"\n{'='*100}")
        print(f"Writing logs...")
        print(f"{'='*100}")
        
        processor.write_mapping_log(log_files["mapping"])
        print(f" Mapping log: {log_files['mapping']}")
        
        processor.write_detailed_log(log_files["detailed"])
        print(f" Detailed log: {log_files['detailed']}")
        
        processor.write_summary_log(log_files["summary"])
        print(f" Summary log: {log_files['summary']}")
        
        processor.write_errors_log(log_files["errors"])
        print(f" Errors log: {log_files['errors']}")
        
        print(f"\n{'='*100}")
        print(f"Summary:")
        print(f"  Total processed: {processor.stats['total']}")
        print(f"  Renamed: {processor.stats['renamed']}")
        print(f"  Skipped: {processor.stats['skipped']}")
        print(f"  Failed: {processor.stats['failed']}")
        print(f"{'='*100}\n")
    
    return processor

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Normalize and rename .ab1 files to standard format"
    )
    parser.add_argument("--config", default="config/config.yaml",
                        help="Path to config file")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be renamed without actually renaming")
    
    args = parser.parse_args()
    
    processor = main(args.config, args.dry_run)
    exit(0)

