# Dockerfile для Auto-mAbSeq Pipeline с полными зависимостями

FROM python:3.11-slim

# Установить системные зависимости
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    clustalo \
    mafft \
    hmmer \
    autoconf \
    automake \
    libtool \
    pkg-config \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Установить рабочую директорию
WORKDIR /app

# ===== INSTALL io_lib (Staden package) =====
RUN mkdir -p /tmp/io_lib && \
    cd /tmp/io_lib && \
    git clone --recursive https://github.com/jkbonfield/io_lib.git . && \
    ./bootstrap && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    ldconfig && \
    rm -rf /tmp/io_lib

# ===== INSTALL ANARCI =====
COPY anarci-master /tmp/anarci
RUN cd /tmp/anarci && \
    python setup.py install && \
    rm -rf /tmp/anarci

# ===== INSTALL AUTOSNAPGENE =====
RUN pip install --no-cache-dir autosnapgene

# ===== INSTALL Python requirements =====
COPY server/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# ===== COPY PROJECT FILES =====
# Копируем конфиг
COPY config/ ./config/

# Копируем скрипты пайплайна
COPY workflow/scripts/ ./workflow/scripts/

# Копируем app.py и templates из server
COPY server/app.py ./
COPY server/templates ./templates

# Копируем run_full_pipeline.sh
COPY run_full_pipeline.sh ./
RUN chmod +x ./run_full_pipeline.sh

# ===== CREATE DATA DIRECTORIES =====
RUN mkdir -p /data/raw/sanger \
              /data/intermediate \
              /data/results/alignment \
              /data/results/dna_files \
              /data/results/qc \
              /data/results/final_snapgene \
              /data/logs \
              /data/jobs \
              /data/archives

# ===== ENVIRONMENT VARIABLES =====
ENV DATA_DIR=/data \
    RAW_DATA_DIR=/data/raw/sanger \
    INTERMEDIATE_DIR=/data/intermediate \
    RESULTS_DIR=/data/results \
    LOGS_DIR=/data/logs \
    FLASK_APP=app.py \
    FLASK_ENV=production \
    PYTHONUNBUFFERED=1 \
    LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# ===== EXPOSE PORT =====
EXPOSE 5000

# ===== HEALTH CHECK =====
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:5000/api/health || exit 1

# ===== START APPLICATION =====
CMD ["python", "app.py"]
