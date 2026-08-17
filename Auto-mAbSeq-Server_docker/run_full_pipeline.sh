#!/bin/bash

# run_full_pipeline.sh
# Запустить полный пайплайн обработки последовательностей антител

set -e  # Выход при ошибке

CONFIG_FILE="${1:-config/config.yaml}"

echo "=========================================="
echo "Auto-mAbSeq Pipeline Full Run"
echo "=========================================="
echo "Config: $CONFIG_FILE"
echo "DATA_DIR: $DATA_DIR"
echo "RAW_DATA_DIR: $RAW_DATA_DIR"
echo ""

# Проверить что конфиг существует
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file not found: $CONFIG_FILE"
    exit 1
fi

# Проверить что есть AB1 файлы
AB1_COUNT=$(find "$RAW_DATA_DIR" -name "*.ab1" 2>/dev/null | wc -l)
if [ $AB1_COUNT -eq 0 ]; then
    echo "Error: No .ab1 files found in $RAW_DATA_DIR"
    exit 1
fi

echo "Found $AB1_COUNT .ab1 files"
echo ""

# ===== STEP 1: Normalize IDs =====
echo "Step 1: Normalizing .ab1 file IDs..."
python workflow/scripts/normalize_ids.py --config "$CONFIG_FILE" || exit 1
echo "✓ Step 1 complete"
echo ""

# ===== STEP 2: Convert AB1 =====
echo "Step 2: Converting .ab1 to FASTA (multiframe)..."
python workflow/scripts/convert_ab1_multiframe.py --config "$CONFIG_FILE" || exit 1
echo "✓ Step 2 complete"
echo ""

# ===== STEP 3: Run ANARCI =====
echo "Step 3: Running ANARCI and filtering sequences..."
python workflow/scripts/run_anarci_filter.py --config "$CONFIG_FILE" || exit 1
echo "✓ Step 3 complete"
echo ""

# ===== STEP 4: Extract Verified Nucleotides =====
echo "Step 4: Extracting verified nucleotide sequences..."
python workflow/scripts/extract_verified_nucleotides.py --config "$CONFIG_FILE" || exit 1
echo "✓ Step 4 complete"
echo ""

# ===== STEP 5: Create DNA Files =====
echo "Step 5: Creating DNA files..."
python workflow/scripts/create_dna_files.py --config "$CONFIG_FILE" || exit 1
echo "✓ Step 5 complete"
echo ""

# ===== STEP 6: Normalize Protein IDs =====
echo "Step 6: Normalizing protein sequence IDs..."
python workflow/scripts/normalize_protein_ids.py --config "$CONFIG_FILE" || exit 1
echo "✓ Step 6 complete"
echo ""

# ===== STEP 7: Build MSA =====
echo "Step 7: Building multiple sequence alignments..."
python workflow/scripts/build_msa.py --config "$CONFIG_FILE" || exit 1
echo "✓ Step 7 complete"
echo ""

# ===== STEP 8: Create Final SnapGene Files =====
echo "Step 8: Creating final SnapGene files..."
python workflow/scripts/create_final_snapgene_files.py --config "$CONFIG_FILE" || exit 1
echo "✓ Step 8 complete"
echo ""

echo "=========================================="
echo "✓ PIPELINE COMPLETE!"
echo "=========================================="
echo ""
echo "Results location:"
echo "  Data: $DATA_DIR"
echo "  Results: $RESULTS_DIR"
echo "  Logs: $LOGS_DIR"
echo ""

exit 0
