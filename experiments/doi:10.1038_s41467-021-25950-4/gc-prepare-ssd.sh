#!/bin/bash
set -e  # Exit immediately if any command fails

SSD_PATH="/dev/disk/by-id/google-local-nvme-ssd-0"
MOUNT_PATH="/mnt/disks/ssd"

echo "Formatting and mounting Local SSD..."
mkfs.ext4 -F "$SSD_PATH"
mkdir -p "$MOUNT_PATH"
mount "$SSD_PATH" "$MOUNT_PATH"
chmod a+w "$MOUNT_PATH"

echo "Downloading HMM database..."
gsutil -m -o "GSUtil:check_hashes=never" cp gs://needle-files/kegg-downloads/ko.hmm.* "$MOUNT_PATH/"
gsutil -m -o "GSUtil:check_hashes=never" cp gs://needle-files/kegg-downloads/ko_thresholds.tsv "$MOUNT_PATH/"

echo "Downloading input FASTA for task ${BATCH_TASK_INDEX}..."
gsutil -o "GSUtil:check_hashes=never" cp \
  "gs://needle-files/experiments/doi:10.1038_s41467-021-25950-4/proteins_${BATCH_TASK_INDEX}.faa" \
  "$MOUNT_PATH/input.faa"
