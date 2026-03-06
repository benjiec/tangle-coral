#!/bin/bash
set -e  # Exit immediately if any command fails

SSD_PATH="/dev/disk/by-id/google-local-nvme-ssd-0"
MOUNT_PATH="/mnt/disks/ssd"

echo "Formatting and mounting Local SSD..."
mkfs.ext4 -F "$SSD_PATH"
mkdir -p "$MOUNT_PATH"
mount "$SSD_PATH" "$MOUNT_PATH"
chmod 777 "$MOUNT_PATH"

echo "Downloading HMM database..."
gsutil -m -o "GSUtil:check_hashes=never" cp -n gs://needle-files/kegg-downloads/ko.hmm.* "$MOUNT_PATH/"
gsutil -m -o "GSUtil:check_hashes=never" cp -n gs://needle-files/kegg-downloads/ko_thresholds.tsv "$MOUNT_PATH/"
