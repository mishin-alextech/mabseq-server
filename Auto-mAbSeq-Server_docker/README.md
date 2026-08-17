ext
# VHH Pipeline

Automated pipeline для обработки Sanger-данных (.ab1) VHH антител с нормализацией ID, ANARCI аннотацией и множественным выравниванием.

## Обзор Pipeline

.ab1 files
↓ (Step 0)
normalize_ids.py → ID mapping log
↓ (Step 1)
convert_ab1_multiframe.py → raw_seqs.fasta, translated_frames.fasta
↓ (Step 2)
run_anarci_filter.py → vhh_valid_anarci.fasta
↓ (Step 3)
extract_verified_nucleotides.py → verify_frames.fasta
↓ (Step 4)
create_dna_files.py → individual .dna files
↓ (Step 5)
normalize_protein_ids.py → normalized protein IDs
↓ (Step 6)
build_msa.py → aligned sequences (FASTA + Clustal)
↓ (Step 7)
create_final_snapgene_files.py → final SnapGene files with traces & CDS


## Структура папок

vhh-pipeline/
├── config/
│   └── config.yaml          # Основной конфиг
├── workflow/
│   ├── scripts/
│   │   ├── normalize_ids.py                # Step 0: Нормализация ID
│   │   ├── convert_ab1_multiframe.py       # Step 1: Конвертация AB1
│   │   ├── run_anarci_filter.py            # Step 2: ANARCI аннотация
│   │   ├── extract_verified_nucleotides.py # Step 3: Верификация
│   │   ├── create_dna_files.py             # Step 4: Создание .dna
│   │   ├── normalize_protein_ids.py        # Step 5: Нормализация protein ID
│   │   ├── build_msa.py                    # Step 6: MSA выравнивание
│   │   └── create_final_snapgene_files.py  # Step 7: Финальные SnapGene файлы
│   └── envs/
│       └── base.yaml                       # Conda окружение
├── data/
│   ├── raw/sanger/         # Входные .ab1 файлы
│   ├── intermediate/       # Промежуточные файлы
│   ├── results/
│   │   ├── alignment/      # MSA результаты
│   │   ├── dna_files/      # Индивидуальные .dna файлы
│   │   ├── qc/
│   │   └── final_snapgene/ # Финальные SnapGene файлы
│   └── logs/               # Логирование
├── run_full_pipeline.sh
└── README.md


## Установка

### 1. Создание Conda окружения

conda env create -f workflow/envs/base.yaml
conda activate vhh-pipeline

### 2. Подготовка данных

Положи входные .ab1 файлы в `data/raw/sanger/`

Поддерживаемые форматы ID:
- `--VHH_raw_seq_cDNA_Cm_TN95089_3_R-1-A01-A96Well.ab1`
- `A1_VHeavy_raw_cDNA_Cm_TN95077_20251014_155315.ab1`
- `--VHH_raw_seq_cDNA_Cm_lib#1_1_3-1-A01-A96Well.ab1`

## Запуск

### Полный pipeline

./run_full_pipeline.sh config/config.yaml

### Dry run

./run_full_pipeline.sh config/config.yaml --dry-run

### Отдельные шаги

# Step 0: Нормализация ID
python workflow/scripts/normalize_ids.py --config config/config.yaml

# Step 1: AB1 → FASTA (3 фрейма)
python workflow/scripts/convert_ab1_multiframe.py --config config/config.yaml

# Step 2: ANARCI аннотация
python workflow/scripts/run_anarci_filter.py --config config/config.yaml

# Step 3: Экстракция верифицированных нуклеотидов
python workflow/scripts/extract_verified_nucleotides.py --config config/config.yaml

# Step 4: Создание .dna файлов
python workflow/scripts/create_dna_files.py --config config/config.yaml

# Step 5: Нормализация protein ID после ANARCI
python workflow/scripts/normalize_protein_ids.py --config config/config.yaml

# Step 6: Множественное выравнивание
python workflow/scripts/build_msa.py --config config/config.yaml

# Step 7: Финальные SnapGene файлы с трассами и CDS
python workflow/scripts/create_final_snapgene_files.py --config config/config.yaml


## Конфигурация

Все параметры находятся в `config/config.yaml`:
ID Normalization (Step 0)
id_processing.enabled — включить/отключить
id_processing.strict_mode — строгая проверка
id_processing.log_transformations — логировать изменения

##AB1 Conversion (Step 1)
ab1_conversion.forward_line — True (forward), False (reverse complement)

##ANARCI (Step 2)
anarci.scheme — imgt, kabat, chothia
anarci.min_aa — мин. длина VHH (80)
anarci.max_aa — макс. длина VHH (150)
anarci.max_stops — макс. стоп-кодонов (3)

## Pipeline control

pipeline:
  steps:                    # Вкл/выкл отдельные шаги
    - name: "normalize_ids"
      enabled: true
  stop_on_error: true       # Остановить при ошибке
  dry_run: false            # Тестовый режим
  
## MSA Alignment (Step 6)

alignment.tool — clustalo или mafft
alignment.num_threads — количество потоков
alignment.group_by_chain — группировка по типам цепей

##Logging

logging.level — DEBUG, INFO, WARNING, ERROR
logging.console_output — вывод в консоль

### ID Normalization (Step 0)
- `id_processing.enabled` — включить/отключить нормализацию
- `id_processing.strict_mode` — строгая проверка
- `id_processing.log_transformations` — логировать преобразования

### AB1 Conversion (Step 1)
- `ab1_conversion.forward_line` — True для forward strand, False для reverse complement

### ANARCI (Step 2)
- `anarci.scheme` — схема нумерации (imgt, kabat, chothia)
- `anarci.min_aa` — минимальная длина VHH
- `anarci.max_aa` — максимальная длина VHH
- `anarci.max_stops` — максимум стоп-кодонов

### MSA Alignment (Step 5)
- `alignment.tool` — clustalo или mafft
- `alignment.num_threads` — количество потоков

## Выход

### Step 0 Output
- `data/logs/id_mapping.txt` — лог преобразования ID

### Step 1 Output
- `data/intermediate/raw_seqs.fasta` — исходные нуклеотидные последовательности
- `data/intermediate/translated_frames.fasta` — белки (3 фрейма)

### Step 2 Output
- `data/intermediate/vhh_valid_anarci.fasta` — верифицированные белки

### Step 3 Output
- `data/intermediate/verify_frames.fasta` — верифицированные нуклеотиды

### Step 4 Output
- `data/results/dna_files/*.dna` — индивидуальные .dna файлы

### Step 5 Output
- `data/results/alignment/vhh_aligned.fasta` — выравненные последовательности (FASTA)
- `data/results/alignment/vhh_alignment.clustal` — выравнивание (Clustal формат)


| Шаг | Файлы                | Путь                                                                                  |
| --- | -------------------- | ------------------------------------------------------------------------------------- |
| 0   | ID mapping log       | data/logs/id_mapping.txt                                                              |
| 1   | Raw sequences        | data/intermediate/raw_seqs.fasta                                                      |
| 2   | ANARCI-validated     | data/intermediate/vhh_valid_anarci.fasta                                              |
| 3   | Verified nucleotides | data/intermediate/verify_frames.fasta                                                 |
| 4   | Individual .dna      | data/results/dna_files/*.dna                                                          |
| 6   | MSA alignment        | data/results/alignment/vhh_aligned.fasta
			       data/results/alignment/vhh_alignment.clustal |
| 7   | Final SnapGene       | data/results/final_snapgene/*.snap(с трассами + CDS)                                  |

## Требования

- Python 3.11+
- Biopython
- PyYAML
- ANARCI (устанавливается в окружение)
- Clustal Omega (или MAFFT)
- autosnapgene
- hmmer 3.4

## Требования к программам

ANARCI должен быть установлен в домашней директории:
$HOME/ANARCI-master/
Если расположение другое, обнови переменную `ANARCI_DIR` в `run_full_pipeline.sh`

## Логирование

### Step 0 Logs
- `data/logs/id_mapping.txt` — преобразование оригинальных ID в нормализованные
- `data/logs/id_processing_detailed_*.log` — подробная информация о парсинге
- `data/logs/id_processing_summary_*.log` — итоговый отчёт
- `data/logs/id_processing_errors_*.log` — только ошибки

## Примеры использования

### Только forward strand
ab1_conversion:
forward_line: True


### Только reverse complement
ab1_conversion:
forward_line: False


### Более строгая фильтрация
anarci:
min_aa: 100
max_aa: 130
max_stops: 1


### MAFFT вместо Clustal
alignment:
tool: "mafft"


##Troubleshooting
ANARCI / hmmer не найден:

Проверьте $HOME/ANARCI-master/

Активируйте conda: conda activate vhh-pipeline
conda install -c bioconda hmmer=3.3.2 -y
conda install -c bioconda anarci

##ID не парсятся:

bash
python workflow/scripts/normalize_ids.py --config config/config.yaml --dry-run
Проверьте data/logs/id_processing_detailed_*.log

Pipeline остановился:
Включите logging.level: "DEBUG" и pipeline.stop_on_error: false.​
