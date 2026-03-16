import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


import argparse

ap = argparse.ArgumentParser()
ap.add_argument("samfile")
ap.add_argument("pngfile")
args = ap.parse_args()


def final_robust_parse(file_path, transcript_len=1280):
    depth = np.zeros(transcript_len + 1)
    snps = {}
    clips = {}
    indels = {}
    
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split('\t')
            if len(parts) < 11: continue
            pos = int(parts[3])
            cigar = parts[5]
            seq = parts[9]
            
            curr_ref_pos = pos
            curr_read_pos = 0
            cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
            
            # Identify start clips
            if cigar_ops and cigar_ops[0][1] == 'S':
                if 1 <= pos <= transcript_len:
                    clips[pos] = clips.get(pos, 0) + 1
            
            for length, op in cigar_ops:
                length = int(length)
                if op in 'M=X':
                    for i in range(length):
                        if 1 <= curr_ref_pos <= transcript_len:
                            depth[curr_ref_pos] += 1
                            if seq != "*" and curr_read_pos < len(seq):
                                base = seq[curr_read_pos]
                                if curr_ref_pos not in snps: snps[curr_ref_pos] = {}
                                snps[curr_ref_pos][base] = snps[curr_ref_pos].get(base, 0) + 1
                        curr_ref_pos += 1
                        curr_read_pos += 1
                elif op == 'I':
                    if 1 <= curr_ref_pos <= transcript_len:
                        indels[curr_ref_pos] = indels.get(curr_ref_pos, 0) + 1
                    curr_read_pos += length
                elif op == 'D':
                    for i in range(length):
                        if 1 <= curr_ref_pos <= transcript_len:
                            indels[curr_ref_pos] = indels.get(curr_ref_pos, 0) + 1
                        curr_ref_pos += 1
                elif op == 'S':
                    if curr_read_pos > 0:
                        if 1 <= curr_ref_pos <= transcript_len:
                            clips[curr_ref_pos] = clips.get(curr_ref_pos, 0) + 1
                    curr_read_pos += length
                elif op == 'N':
                    curr_ref_pos += length
                    
    consensus = []
    for i in range(1, transcript_len + 1):
        if i in snps and snps[i]:
            consensus.append(max(snps[i], key=snps[i].get))
        else:
            consensus.append('N')
    consensus = "".join(consensus)
    
    snp_final, indel_final, clip_final = [], [], []
    for i in range(1, transcript_len + 1):
        d = depth[i]
        if d >= 1:
            if i in snps:
                con_base = consensus[i-1]
                mismatch = sum(v for k, v in snps[i].items() if k != con_base)
                freq = mismatch / d
                if freq >= 0.1: snp_final.append((i, freq))
            if i in indels:
                freq = indels[i] / d
                if freq >= 0.1: indel_final.append((i, freq))
            if i in clips:
                freq = clips[i] / d
                if freq >= 0.1: clip_final.append((i, freq))
            
    return depth, consensus, snp_final, indel_final, clip_final


depth, consensus, snp_data, indel_data, clip_data = final_robust_parse(args.samfile)

window = 50
gc = []
for i in range(len(consensus)):
    sub = consensus[max(0, i-25):min(len(consensus), i+25)]
    gc.append((sub.count('G')+sub.count('C'))/len(sub) if sub else 0)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios':[5, 1]}, sharex=True)
ax1.plot(range(len(depth)), depth, color='dodgerblue', lw=2, label='Read Depth')

def get_alpha(f):
    return np.clip(2*f, 0.2, 1.0)

msize = 120
# Legend tracking
snp_label, indel_label, clip_label = "SNP/MNP", "Indel", "Clip"

for i, (p, f) in enumerate(snp_data):
    ax1.scatter(p, depth[p] + 30, marker='v', color='red', s=msize, alpha=get_alpha(f), 
                edgecolors='none', label=snp_label if i == 0 else "")
for i, (p, f) in enumerate(indel_data):
    ax1.scatter(p, depth[p] + 70, marker='s', color='green', s=msize, alpha=get_alpha(f), 
                edgecolors='none', label=indel_label if i == 0 else "")
for i, (p, f) in enumerate(clip_data):
    y_off = depth[p] + 110 if depth[p] > 50 else np.max(depth)*0.15
    ax1.scatter(p, y_off, marker='x', color='purple', s=msize+50, alpha=get_alpha(f), 
                edgecolors='none', label=clip_label if i == 0 else "")

ax1.set_ylabel('Read Depth')
ax1.set_title('Pileup Visualization: Sample 1 (Divergence Spectrum)')
ax1.set_xlim(-80, 1360)
ax1.legend(loc='upper right')

cmap = mcolors.LinearSegmentedColormap.from_list("GC", ["blue", "yellow", "red"])
ax2.imshow(np.array(gc).reshape(1, -1), aspect='auto', cmap=cmap, norm=mcolors.Normalize(0.2, 0.8), extent=[0, 1280, 0, 1])
ax2.set_yticks([])
ax2.set_ylabel('GC %')
ax2.set_xlabel('Transcript Position (bp)')

plt.tight_layout()
plt.savefig(args.pngfile)
print("Final plot with legend saved.")
