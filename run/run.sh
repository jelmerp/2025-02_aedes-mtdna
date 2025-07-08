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
fa_in=data/NAD4CleanedUnaligned.txt
ls -lh "$fa_in" && grep -c ">" "$fa_in"

# Define output files
outdir=results/NAD4 && mkdir -p "$outdir"/fasta
fa_withrefs="$outdir"/fasta/NAD4_withrefs.fa

# Add reference sequence
cat "$fa_in" results/NAD4/refdata/DQ176828.2.fna > "$fa_withrefs"
ls -lh "$fa_withrefs" && grep -c ">" "$fa_withrefs"

# Align the sequences (with reference sequence)
sbatch mcic-scripts/align/align_fa.sh -i "$fa_withrefs" -o "$outdir"/aligned
ls -lh "$outdir"/aligned/NAD4_withrefs_aln.fa

# Clip the alignment with ClipKit
sbatch mcic-scripts/align/clipkit.sh -i "$outdir"/aligned/NAD4_withrefs_aln.fa -o "$outdir"/clipkit
ls -lh "$outdir"/clipkit/NAD4_withrefs_aln.fa

# Make a tree with IQ-Tree
sbatch mcic-scripts/trees/iqtree.sh \
    --infile "$outdir"/clipkit/NAD4_withrefs_aln.fa \
    --outdir "$outdir"/iqtree \
    --nboot 1000 \
    --more_opts "--msub mitochondrial"
