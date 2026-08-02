# Coral-symbiont SOD phylogeny (v2, 13 transcripts)

## Method
Independent phylogenetic test of the metal identity and grouping of 13 K04564
(Fe/Mn-SOD) transcript proteins recovered from *Symbiodiniaceae* symbionts of
four coral hosts (*Orbicella*, *Pseudodiploria*, *Siderastrea*) plus two
*Cladocopium goreaui* isoforms. The 13 queries were placed against a balanced
52-sequence reference panel (Fe- and Mn-SODs spanning cyanobacteria, land
plants, green/red/brown algae, diatoms, alveolates, enterobacteria, fungi,
mammals). Full-length ORFs (leaders included) were aligned with MAFFT L-INS-i;
trimAl `-gappyout` removed the divergent, gappy N-terminal leader columns,
leaving 185 columns of the conserved mature SOD domain; IQ-TREE inferred the ML
tree (ModelFinder selected WAG+R4; 1000 ultrafast bootstrap + 1000 SH-aLRT).
Because trimAl discards exactly the leader residues that encode localization,
the tree is independent of both the metal-residue fingerprint and the
DeepLoc/leader-based targeting analysis.

## Result
The 13 transcripts are monophyletic and resolve into the same two groups the
leader/DeepLoc analysis found. **Group A** (3 secreted transcripts;
SP-only leader) forms a tight, strongly-supported clade (SH-aLRT 99.9 / UFBoot
100). **Group B** (10 transcripts; bipartite SP+TP plastid-type leader,
including both new *C. goreaui* isoforms) is sister to two diatom Fe-SODs
(*Phaeodactylum*, *Eolimna*; that node UFBoot 100), within a stramenopile
Fe-SOD neighborhood with a brown-algal SOD (*Ectocarpus*) at its base. Critically,
Groups A and B sit inside one well-supported symbiont clade (UFBoot 93) whose
only annotated metal references are the two diatom **Fe**-SODs. Branch-length
(long-branch-attraction) risk is low: symbiont tip branches are short (the
transcripts are near-identical to one another). Internal resolution *within*
Group B is weak (several near-identical or identical sequences across coral
hosts and *Cladocopium*), but the node tying the whole group to the diatom
Fe-SODs is robust.

## Interpretation for annotation
The tree resolves whole-domain **ancestry**, not the metal each protein binds
today (Fe/Mn-SODs are structurally near-identical and switch metal among close
homologs), and it is mute on subcellular targeting (leaders trimmed). Both
groups fall in an Fe-SOD lineage; the A/B split is a real, strongly-supported
topological distinction but is **not** a metal distinction the tree can see — it
tracks the leader difference (SP-only vs bipartite), not a switch to Mn. A prior
"Group A nearest homolog = green-algal MnSOD" call rested on single-nearest-tip
patristic distance (Chlamydomonas Mn ~0.79-0.83 vs diatom Fe ~0.88-0.92, a
~0.09 margin) that the topology contradicts: Chlamydomonas Mn is *outside* the
symbiont clade and joins only at a weakly-supported deeper node. The only
evidence calling any of these transcripts "Mn" is the residue fingerprint; the
tree places all 13 in an Fe-SOD neighborhood. **Weakest link:** distinguishing
"same Fe-neighborhood, different sub-branch" from "genuinely different metal" is
what whole-domain topology cannot resolve, because the discriminating signal
lives in a handful of active-site residues. **Validating test:** score the
canonical Fe-vs-Mn discriminant (metal-coordinating and second-shell) residues
for Groups A and B against the diatom Fe and Chlamydomonas Mn references at
those columns specifically; definitively, express one protein per group and
determine bound metal directly (ICP-MS/EPR + H2O2/azide inhibitor assay).

## Files
- `sod_phylogeny_v2.png` — annotated ML tree; 13 transcripts labeled with full
  accessions; Fe refs red, Mn refs blue, Group A/B queries bold (square/circle);
  UFBoot >=80 shown.
- `symb_phylo_grouping_v2.tsv` — per-transcript group, leader architecture,
  fingerprint metal, patristic-nearest reference + metal, and topological
  placement/metal-neighborhood.
- `sod_ML_v2.treefile` — Newick tree with SH-aLRT/UFBoot support.
- `sod_alignment_trimmed_v2.faa` — 185-column trimAl-trimmed mature-domain
  alignment the tree was inferred from.
