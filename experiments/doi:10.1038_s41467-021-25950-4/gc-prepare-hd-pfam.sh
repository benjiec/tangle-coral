#!/bin/bash
set -e  # Exit immediately if any command fails

WORKSPACE="/data"
HMM_FILE="$WORKSPACE/Pfam-A.hmm.h3i"

mkdir -p "$WORKSPACE"
chmod 777 "$WORKSPACE"

if [ ! -f "$HMM_FILE" ]; then
  echo "HMM file missing. Downloading from GCS..."
  gsutil -m -o "GSUtil:check_hashes=never" cp gs://needle-files/pfam-downloads/Pfam-A.hmm.* "$WORKSPACE/"
else
  echo "HMM file already exists. Skipping download."
fi
