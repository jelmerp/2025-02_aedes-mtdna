#!/bin/bash

# Several sequences are of poor quality - remove these:
seqkit grep -v -p 'Wooster_2024_14**' data/Clean_COI_Sequences.fa |
    seqkit grep -v -p 'Medina_2023_4' |
    seqkit grep -v -p 'Medina_2023_2' |
    seqkit grep -v -p 'Medina_2023_3' |
    seqkit grep -v -p 'Wooster_2024_22' |
    seqkit grep -v -p 'Wooster_2022_8' |
    seqkit grep -v -p 'Wooster_2022_5' |
    seqkit grep -v -p 'Wooster_2024_18**' |
    seqkit grep -v -p 'Wooster_2022_2' |
    seqkit grep -v -p 'Wooster_2024.10' \
    > data/COI.fa
grep -c ">" data/Clean_COI_Sequences.fa data/COI.fa

# Add reference sequences
ref_fa_init=/fs/ess/PAS0471/jelmer/assist/01_archive/2022-04_ferd_phylo/results/unalign_woutgroups.fa
seqkit grep -f data/ref_ids.txt "$ref_fa_init" > data/COI_refs.fa
cat data/COI.fa data/COI_refs.fa > data/COI_withrefs.fa

# Align the sequences 
bash mcic-scripts/align/align_fa.sh -i data/COI.fa -o results/aligned
bash mcic-scripts/align/align_fa.sh -i data/COI_withrefs.fa -o results/aligned

# Clip the alignment with ClipKit
sbatch mcic-scripts/align/clipkit.sh -i results/aligned/COI_aln.fa -o results/clipkit
sbatch mcic-scripts/align/clipkit.sh -i results/aligned/COI_withrefs_aln.fa -o results/clipkit

# Make a tree with IQ-Tree
sbatch -t30 mcic-scripts/trees/iqtree.sh \
    --infile results/clipkit/COI_withrefs_aln.fa \
    --outdir results/iqtree \
    --root "HQ944971.1" \
    --nboot 1000 \
    --opts "--msub mitochondrial"
# results/iqtree/COI_withrefs_aln.treefile
