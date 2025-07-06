# 2025-06-03 -- Haplotype network for Will Balding
# Jelmer Poelstra, R 4.4.2 on own computer

# SETUP ------------------------------------------------------------------------
# Load packages
library(tidyverse)
library(adegenet)
library(ape)
library(ips)          # Sequence trimming
library(pegas)        # Haplotype network inference
library(hierfstat)    # FST calculation

# Define input files
fasta_file <- "results/will_balding/WBNAD4Aligned.fasta"

# Define output files
outdir <- "results/will_balding"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
network_file <- file.path(outdir, "haplotype_network.png")
fst_file <- file.path(outdir, "fst.tsv")

# DATA PREP --------------------------------------------------------------------
# Read the sequences
seqs <- adegenet::fasta2DNAbin(fasta_file)

# Filter the alignment
# A) Remove ends with Ns and gaps:
seqs <- ips::trimEnds(seqs, min.n.seq = round(1 * nrow(seqs)))
# B) This seems to remove all gaps, while gap.max = 0 leaves no sequences:
seqs <- ips::deleteGaps(seqs, gap.max = 1)

# Visually check the alignment
# write.FASTA(seqs, file = "tmp.fa")
# DECIPHER::BrowseSeqs(Biostrings::readDNAStringSet("tmp.fa"))

# Populations and colors
pop_cols_uniq <- c(Cleveland = "purple3", Wooster = "chartreuse4") 
pops <- sub("WB_([[:alpha:]]+)_.*", "\\1", rownames(seqs))
pops_uniq <- sort(unique(pops))
pop_cols <- pop_cols_uniq[match(pops, names(pop_cols_uniq))]


# POPULATION-LEVEL HAPLOTYPE NETWORK -------------------------------------------
# Create the network
h <- pegas::haplotype(seqs, strict = FALSE, trailingGapsAsN = TRUE)
hname <- paste("H", 1:nrow(h), sep = "")
rownames(h) <- paste(hname)

net <- pegas::haploNet(h, d = NULL, getProb = TRUE)

ind.hap <- with(
  utils::stack(setNames(attr(h, "index"), rownames(h))),
  table(hap = ind, individuals = rownames(seqs)[values])
)

# Function to create the plot
make_plot <- function() {
  # Open a PNG to plot in
  png(
    network_file,
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
    x = -9,
    y = 17,
    pops_uniq,
    fill = pop_cols_uniq,
    cex = 1.2,
    ncol = 1,
    bty = "n",
    x.intersp = 0.2
    )
  
  # Close the graphics device
  dev.off()
}

# Create the plot
make_plot()


# FST --------------------------------------------------------------------------
# Convert sequences to a genind object and add population info
gi <- adegenet::DNAbin2genind(seqs)
pop(gi) <- factor(pops)

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
