#!/bin/bash
set -e  # Exit immediately if any command fails

WORKSPACE="/data"

# Clean up any previous task inputs and outputs
rm -f $WORKSPACE/input_*.faa
rm -f $WORKSPACE/output_*.tsv

echo "Downloading input FASTA for task ${BATCH_TASK_INDEX}..."
gsutil -o "GSUtil:check_hashes=never" cp \
  "gs://needle-files/experiments/doi:10.1038_s41467-021-25950-4/proteins_${BATCH_TASK_INDEX}.faa" \
  "$WORKSPACE/input_${BATCH_TASK_INDEX}.faa"
