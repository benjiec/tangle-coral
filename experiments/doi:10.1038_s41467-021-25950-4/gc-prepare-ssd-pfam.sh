#!/bin/bash
set -e  # Exit immediately if any command fails

SSD_PATH="/dev/disk/by-id/google-local-nvme-ssd-0"
MOUNT_PATH="/mnt/disks/ssd"
HMM_FILE="$MOUNT_PATH/Pfam-A.hmm.h3i"

if mountpoint -q "$MOUNT_PATH"; then
  echo "SSD already mounted. Skipping format."
else
  echo "Formatting and mounting Local SSD..."
  mkfs.ext4 -F "$SSD_PATH"
  mkdir -p "$MOUNT_PATH"
  mount "$SSD_PATH" "$MOUNT_PATH"
  chmod 777 "$MOUNT_PATH"
fi

if [ ! -f "$HMM_FILE" ]; then
  echo "HMM file missing. Downloading from GCS..."
  gsutil -m -o "GSUtil:check_hashes=never" cp gs://needle-files/pfam-downloads/Pfam-A.hmm.* "$MOUNT_PATH/"
else
  echo "HMM file already exists. Skipping download."
fi
