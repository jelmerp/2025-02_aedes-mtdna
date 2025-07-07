# 2025-06-04, Jelmer Poelstra, R 4.4.2 on own computer

# SETUP ------------------------------------------------------------------------
# See https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0243927

# Load packages
library(tidyverse)
library(adegenet)
library(ape)
library(ips) # Not working on OSC
library(haplotypes)
library(pegas)
library(phytools)

# Settings
pop_names <- c("Cleveland", "Medina", "Wooster")
pop_cols_uniq <- c("purple3", "chartreuse4", "goldenrod3")
nuc_cols <- c("rosybrown", "sienna1", "lightgoldenrod1", "lightskyblue1", "grey")

# Other settings
options(ignore.negative.edge = TRUE)

# Define input files
fasta_file <- "results/COI/aligned/COI_aln.fa"
tree_file <- "results/COI/iqtree/COI_withrefs_aln.treefile"
hap_xy_file <- "results/COI/hapnetwork_coords.rds"

# Define output files
outdir <- "results/COI"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
network_file <- file.path(outdir, "haplotype_network.png")
fst_file <- file.path(outdir, "fst.tsv")


# DATA PREP --------------------------------------------------------------------
# Read the sequences
seqs <- fasta2DNAbin(fasta_file)

# Read the haplotype network coordinates
xy <- readRDS(hap_xy_file)

# Filter the alignment
# A) Remove ends with Ns and gaps:
seqs <- ips::trimEnds(seqs, min.n.seq = round(1 * nrow(seqs)))
# B) This seems to remove all gaps, while gap.max = 0 leaves no sequences:
seqs <- ips::deleteGaps(seqs, gap.max = 1)

# Check the alignment
# write.FASTA(seqs, file = "tmp.fa")
# DECIPHER::BrowseSeqs(Biostrings::readDNAStringSet("tmp.fa"))


# METADATA PREP ----------------------------------------------------------------
# Simplify sample labels in sequence object
rownames(seqs) <- rownames(seqs) |>
  gsub("\\*+", "", x = _) |>
  gsub("_([0-9]+$)", ".\\1", x = _) |>
  gsub("^ON[0-9]+_", "", x = _)

# Infer haplotypes
kh <- haplotypes::haplotype(as.dna(seqs))


# POPULATION-LEVEL HAPLOTYPE NETWORK -------------------------------------------
# After an initial plot,
# I used the following to interactively (!) change the location of one of the nodes:
# xy <- replot()
# saveRDS(xy, "results/hapnetwork_coords.rds")

# Create the network
h <- pegas::haplotype(seqs, strict = FALSE, trailingGapsAsN = TRUE)
hname <- paste("H", 1:nrow(h), sep = "")
rownames(h) <- paste(hname)
net <- haploNet(h, d = NULL, getProb = TRUE)
ind.hap <- with(
  utils::stack(setNames(attr(h, "index"), rownames(h))),
  table(hap = ind, individuals = rownames(seqs)[values])
)

# Set population colors and plot margins
pop_cols <- c(
  rep(pop_cols_uniq[1], 19),
  rep(pop_cols_uniq[2], 1),
  rep(pop_cols_uniq[3], 34)
  )

# Open a PNG to plot in
png("results/haplotype_network.png",
    width = 480 * (600/72),
    height = 480 * (600/72),
    units = "px",
    res = 600
    )

# Set the plot margins
par(mar = c(0.001, 0.001, 0.001, 0.001))

# Create the plot
plot(
  net,
  size = attr(net, "freq"),
  col = "grey",
  bg = pop_cols,
  xy = xy,
  scale.ratio = 2,
  cex = 0.7,
  labels = FALSE,
  pie = ind.hap,
  show.mutation = 1,
  font = 2,
  fast = TRUE,
  )

# Add a legend
legend(
  x = -20,
  y = 25,
  pop_names,
  fill = pop_cols_uniq,
  cex = 1.2,
  ncol = 1,
  bty = "n",
  x.intersp = 0.2
  )

# Close the graphics device
dev.off()


# FST --------------------------------------------------------------------------
# Convert sequences to a genind object and add population info
gi <- adegenet::DNAbin2genind(seqs)
pop(gi) <- factor(meta$Location[match(rownames(seqs), meta$ID)])

# Compute FST and its confidence interval
fst_results <- hierfstat::pairwise.WCfst(gi)
fst_CI <- hierfstat::boot.ppfst(gi)  # Get confidence intervals for FST values

# Prepare a results table
fst_ll <- data.frame(fst_CI$ll) |>
  rownames_to_column("pop1") |>
  pivot_longer(cols = !pop1, names_to = "pop2", values_to = "lower_limit")

fst_final <- data.frame(fst_results) |>
  rownames_to_column("pop1") |>
  pivot_longer(cols = !pop1, names_to = "pop2", values_to = "FST") |>
  left_join(fst_ll, by = c("pop1", "pop2")) |>
  drop_na() |>
  mutate(is_significant = ifelse(test = lower_limit > 0, yes = TRUE, no = FALSE))

write_tsv(fst_final, fst_file)
