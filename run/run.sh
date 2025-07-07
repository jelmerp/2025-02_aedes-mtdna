#!/bin/bash

# ==============================================================================
#                               COI
# ==============================================================================
# Define input files
ref_fa_init=/fs/ess/PAS0471/jelmer/assist/01_archive/2022-04_ferd_phylo/results/COIunalign_woutgroups.fa
fa_init=data/Clean_COI_Sequences.fa 
ref_ids=metadata/ref_ids.txt

# Define output files
outdir=results/COI && mkdir -p "$outdir"
fa="$outdir"/fasta/COI.fa

# Settings
outgroup=HQ944971.1

# Several sequences are of poor quality - remove these:
seqkit grep -v -p 'Wooster_2024_14**' "$fa_init" |
    seqkit grep -v -p 'Medina_2023_4' |
    seqkit grep -v -p 'Medina_2023_2' |
    seqkit grep -v -p 'Medina_2023_3' |
    seqkit grep -v -p 'Wooster_2024_22' |
    seqkit grep -v -p 'Wooster_2022_8' |
    seqkit grep -v -p 'Wooster_2022_5' |
    seqkit grep -v -p 'Wooster_2024_18**' |
    seqkit grep -v -p 'Wooster_2022_2' |
    seqkit grep -v -p 'Wooster_2024.10' \
    > "$fa"
grep -c ">" "$fa_init" "$fa"

# Add reference sequences
seqkit grep -f "$ref_ids" "$ref_fa_init" > "$outdir"/COI_refs.fa
cat "$fa" "$outdir"/COI_refs.fa > data/COI_withrefs.fa

# Align the sequences (with and without reference sequences)
bash mcic-scripts/align/align_fa.sh -i "$fa" -o "$outdir"/aligned
bash mcic-scripts/align/align_fa.sh -i data/COI_withrefs.fa -o "$outdir"/aligned

# Clip the alignment with ClipKit (with and without reference sequences)
sbatch mcic-scripts/align/clipkit.sh \
    -i "$outdir"/aligned/COI_aln.fa -o "$outdir"/clipkit
sbatch mcic-scripts/align/clipkit.sh \
    -i "$outdir"/aligned/COI_withrefs_aln.fa -o "$outdir"/clipkit

# Make a tree with IQ-Tree
sbatch mcic-scripts/trees/iqtree.sh \
    --infile "$outdir"/clipkit/COI_withrefs_aln.fa \
    --outdir "$outdir"/iqtree \
    --root "$outgroup" \
    --nboot 1000 \
    --more_opts "--msub mitochondrial"


# ==============================================================================
#                               NAD4
# ==============================================================================
# Define input files
fa_aln=data/WBNAD4Aligned.fasta

# Define output files
outdir=results/NAD4 && mkdir -p "$outdir"

# Clip the alignment with ClipKit
sbatch mcic-scripts/align/clipkit.sh -i "$fa_aln" -o "$outdir"/clipkit

# Make a tree with IQ-Tree
sbatch mcic-scripts/trees/iqtree.sh \
    --infile "$outdir"/clipkit/WBNAD4Aligned.fasta \
    --outdir "$outdir"/iqtree \
    --nboot 1000 \
    --more_opts "--msub mitochondrial"

#--root "$outgroup" \
