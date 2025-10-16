#!/bin/bash

# Base URL for downloads
ROOT_URL="https://bioinfo.ahu.edu.cn/synMall/common/download-static?filename=analysis/"


echo "Starting download of SynMall data..."

# Visualization files
FILES=(
    "train-SpliceBERT-RNA.pth"
    "ERINE-RNA-last-hidden.pth"
    "train_with_biological.pt"
)

for file in "${FILES[@]}"; do
    echo "Downloading $file..."
    wget -c "${ROOT_URL}${file}" -O ./"$file"
done

echo "All downloads completed."
