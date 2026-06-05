import csv

files = [
  "sequence_3di_ko_3di_c0.8.tsv",
  "sequence_3di_ko_3di_c0.5.tsv",
  "sequence_3di_ko_3di_c0.tsv",
  "sequence_3di_ko_aa-3di_c0.8.tsv",
  "sequence_3di_ko_aa-3di_c0.tsv"
]

ko_file = "sequence_ko.tsv"
ko_threshold_file = "ko_thresholds.tsv"
tx_protein_file = "transcript_proteins.tsv"

hmmer_assignments = {}
fs_param_assignments = {}
ko_consensus_length = {}
protein_length = {}

print("load ko thresholds")
with open(ko_threshold_file, "r") as f:
  reader = csv.DictReader(f, delimiter='\t')
  for row in reader:
    ko_consensus_length[row["model"]] = int(row["mlen"]) if row["mlen"] else -1

print("load transcript proteins")
with open(tx_protein_file, "r") as f:
  reader = csv.DictReader(f, delimiter='\t')
  for row in reader:
    if row["target_type"] == "protein":
      protein_length[row["target_accession"]] = int(row["target_end"])-int(row["target_start"])+1

print("load hmmer")
with open(ko_file, "r") as f:
  reader = csv.DictReader(f, delimiter='\t')
  for row in reader:
    protein_id = row["query_accession"]
    ko = row["target_accession"]
    eval = float(row["evalue"])
    if protein_id not in hmmer_assignments:
      hmmer_assignments[protein_id] = []
    hmmer_assignments[protein_id].append((ko, eval))

for result_fn in files:
  parameter = result_fn[len("sequence_3di_ko_"):-len(".tsv")]
  print("load", parameter)
  fs_param_assignments[parameter] = {}
  adict = fs_param_assignments[parameter]

  with open(result_fn, "r") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
      protein_id = row["query_accession"]
      ko = row["target_accession"][0:6]
      eval = float(row["evalue"])
      qs = int(row["qstart"])
      qe = int(row["qend"])
      ts = int(row["tstart"])
      te = int(row["tend"])
      if protein_id not in adict:
        adict[protein_id] = []
      adict[protein_id].append((ko, eval, qs, qe, ts, te))

print("ranking hmmer")
for k,v in hmmer_assignments.items():
  hmmer_assignments[k] = sorted(v, key=lambda t: t[1])

for param, adict in fs_param_assignments.items():
  print(f"\nranking {param}")

  total_proteins = len(adict)
  not_in_hmmer = 0
  overlap_with_hmmer = 0
  agrees_with_hmmer_top3 = 0
  agrees_with_hmmer_top1 = 0
  new_hit_plen_100 = 0
  new_hit_plen_100_klen_100 = 0
  new_hit_plen_100_perc_80_klen_100 = 0
  new_hit_plen_100_perc_80_klen_100_perc_50 = 0
  same_hit_klen_100 = 0
  same_hit_p_perc_80_klen_100 = 0
  same_hit_p_perc_80_klen_100_perc_50 = 0

  for k,v in adict.items():
    adict[k] = sorted(v, key=lambda t: t[1])
    protein_id = k
    assignments = adict[k]
    
    if protein_id not in hmmer_assignments:
      not_in_hmmer += 1
      if protein_length[protein_id] >= 100:
        new_hit_plen_100 += 1

        for i,(fs_ko, fs_eval, qs, qe, ts, te) in enumerate(assignments):
          if fs_ko in ko_consensus_length and ko_consensus_length[fs_ko] > -1 and ko_consensus_length[fs_ko] >= 100:
            new_hit_plen_100_klen_100 += 1
            if (qe-qs+1)/protein_length[protein_id] >= 0.8:
              new_hit_plen_100_perc_80_klen_100 += 1
              if (te-ts+1)/ko_consensus_length[fs_ko] >= 0.5:
                new_hit_plen_100_perc_80_klen_100_perc_50 += 1
            break

    else:

      found_hmmer = False
      found_top1 = False
      found_top3 = False
      found_klen_100 = False
      found_p_perc_80 = False
      found_k_perc_50 = False

      for i,(fs_ko,fs_eval,qs,qe,ts,te) in enumerate(assignments):
        for j,(hm_ko,hm_eval) in enumerate(hmmer_assignments[protein_id]):
          if fs_ko == hm_ko:
            found_hmmer = True

            if fs_ko in ko_consensus_length and ko_consensus_length[fs_ko] > -1 and ko_consensus_length[fs_ko] >= 100:
              found_klen_100 = True
              if (qe-qs+1)/protein_length[protein_id] >= 0.8:
                found_p_perc_80 = True
                if (te-ts+1)/ko_consensus_length[fs_ko] >= 0.5:
                  found_k_perc_50 = True

            if i < 3 and j < 3:
              found_top3 = True
              if i < 1 and j < 1:
                found_top1 = True

      if found_hmmer:
        overlap_with_hmmer += 1

        if found_klen_100:
          same_hit_klen_100 += 1
          if found_p_perc_80:
            same_hit_p_perc_80_klen_100 += 1
            if found_k_perc_50:
              same_hit_p_perc_80_klen_100_perc_50 += 1

        if found_top3:
          agrees_with_hmmer_top3 += 1
          if found_top1:
            agrees_with_hmmer_top1 += 1

  print(f"total proteins: {total_proteins}")
  print(f"new hits not in hmm results: {not_in_hmmer}")

  in_hmmer = total_proteins-not_in_hmmer
  print(f"found at least one of hmmer hits: {overlap_with_hmmer}, {100*overlap_with_hmmer/in_hmmer:.0f}%")
  print(f"  and ko over 100 aa: {same_hit_klen_100}, {100*same_hit_klen_100/in_hmmer:.0f}%")
  print(f"    and % protein matched over 80%: {same_hit_p_perc_80_klen_100}, {100*same_hit_p_perc_80_klen_100/in_hmmer:.0f}%")
  print(f"      and % ko matched over 50%: {same_hit_p_perc_80_klen_100_perc_50}, {100*same_hit_p_perc_80_klen_100_perc_50/in_hmmer:.0f}%")
  print(f"  found top3 in both: {agrees_with_hmmer_top3}, {100*agrees_with_hmmer_top3/in_hmmer:.0f}%")
  print(f"    found top1 in both: {agrees_with_hmmer_top1}, {100*agrees_with_hmmer_top1/in_hmmer:.0f}%")

  print(f"new hit protein over 100 aa: {new_hit_plen_100}, {100*new_hit_plen_100/not_in_hmmer:.0f}%")
  print(f"  and ko over 100 aa: {new_hit_plen_100_klen_100}, {100*new_hit_plen_100_klen_100/not_in_hmmer:.0f}%")
  print(f"    and % protein matched over 80%: {new_hit_plen_100_perc_80_klen_100}, {100*new_hit_plen_100_perc_80_klen_100/not_in_hmmer:.0f}%")
  print(f"      and % ko matched over 50%: {new_hit_plen_100_perc_80_klen_100_perc_50}, {100*new_hit_plen_100_perc_80_klen_100_perc_50/not_in_hmmer:.0f}%")
