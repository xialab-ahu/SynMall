#!/bin/bash

# Base URL for downloads
ROOT_URL="https://bioinfo.ahu.edu.cn/synMall/common/download-static?filename=analysis/"

# Create directories if they don't exist
mkdir -p ./Result/whole_genome
mkdir -p ./Data/Analysis

echo "Starting download of SynMall data..."

# Whole genome synScore
echo "Downloading synScore.whole.genome.txt..."
wget -c "${ROOT_URL}synScore.whole.gnome.txt" -O ./Result/whole_genome/synScore.whole.gnome.txt

# Analysis files
FILES=(
    "ACMG_actions_synonymous_variants.txt"
    "gene_information_canonical.txt"
    "InterPro.tsv"
    "gnomad.v4.1.constraint_sSNVs.txt"
    "gnomad.v4.1.constraint_genes.txt"
    "COSMIC_v100.txt"
    "COSMIC_v100_top.txt"
    "COSMIC_v100_top_Genes.txt"
    "COSMIC_v100_top_Transcript_Wise.csv"
    "COSMIC_v100_KEGG.csv"
)

for file in "${FILES[@]}"; do
    echo "Downloading $file..."
    wget -c "${ROOT_URL}${file}" -O ./Data/Analysis/"$file"
done

echo "All downloads completed."
