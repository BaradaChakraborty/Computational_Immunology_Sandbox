#!/bin/bash

# ==============================================================================
# Script: 01_Linux_Data_Retrieval.sh
# Purpose: Automated retrieval of Melanoma scRNA-seq matrices via Linux CLI
# Project: BioGrademy scRNA-seq Pipeline (Cert: BL25CFB15048)
# ==============================================================================

echo "Initializing BioGrademy Data Pipeline..."

# 1. Create secure data directories if they do not exist
mkdir -p data/raw
mkdir -p data/processed
echo "Directories verified."

# 2. Download raw scRNA-seq count matrices via GEO/curl
# (Using GSE115978 Melanoma ICB dataset as the baseline)
echo "Fetching raw single-cell transcriptomic matrices from server..."
curl -o data/raw/GSE115978_counts.csv.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_counts.csv.gz"

echo "Data retrieval complete. Ready for R-based Quality Control."