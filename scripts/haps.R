# 2025-06-04, Jelmer Poelstra, R 4.4.2 on own computer

# SETUP ------------------------------------------------------------------------
# See https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0243927

# Load packages
library(tidyverse)
library(adegenet)
library(ape)
library(ggtree)
library(ips)
library(haplotypes)
library(pegas)
library(phytools)
library(ggnewscale)

# Settings
pop_names <- c("Cleveland", "Medina", "Wooster")
pop_cols_uniq <- c("purple3", "chartreuse4", "goldenrod3")
nuc_cols <- c("rosybrown", "sienna1", "lightgoldenrod1", "lightskyblue1", "grey")

# Other settings
options(ignore.negative.edge = TRUE)

# Define input files
fasta_file <- "results/COI_aln.fa"
tree_file <- "results/COI_withrefs_aln.treefile"
meta_refs_file <- "data/refs_meta.tsv"
hap_xy_file <- "results/hapnetwork_coords.rds"

# Define output files
outdir <- "results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
network_file <- file.path(outdir, "haplotype_network.png")
fst_file <- file.path(outdir, "fst.tsv")


# DATA PREP --------------------------------------------------------------------
# Read the sequences
seqs <- fasta2DNAbin(fasta_file)

# Read the metadata
meta_refs_init <- read_tsv(meta_refs_file, show_col_types = FALSE)

# Read the tree
iqtree <- read.tree(tree_file)
iqtree <- ape::root(iqtree, outgroup = "HQ944971.1", resolve.root = TRUE)

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

# Simplify sample labels in tree object
iqtree$tip.label <- iqtree$tip.label |>
  gsub("\\*+", "", x = _) |>
  gsub("_([0-9]+$)", ".\\1", x = _) |>
  gsub("^ON[0-9]+_", "", x = _) |>
  sub("HQ944971.1", "HQ944971.1_Aedes canadensis", x = _) |>
  sub("GQ254798.1", "GQ254798.1_Iowa_2008", x = _) |> 
  sub("GQ254793.1_Iowa_USA.2007", "GQ254793.1_Illinois_2007", x = _) |>
  sub("\\.(20[0-9]+)$", "_\\1", x = _)

# Infer haplotypes
kh <- haplotypes::haplotype(as.dna(seqs))

# Prep metadata for ref sequences
meta_refs <- meta_refs_init |>
  mutate(
    ID = sub("HQ944971.1", "HQ944971.1_Aedes canadensis", ID),
    ID = sub("GQ254798.1", "GQ254798.1_Iowa_2008", ID), 
    ID = sub("GQ254793.1_Iowa_USA.2007", "GQ254793.1_Illinois_2007", ID),
    ID = sub("MW475673.1_Belgium.2019", "MW475673.1_Belgium_2019", ID),
    Year = as.character(Year),
    focal = FALSE,
    Area = Location
  )     

# Prep metadata for own sequences and combine with ref
meta <- tibble(full_id = names(unlist(kh@hapind))) |>
  mutate(
    ID = sub("haplotype[0-9]+\\.", "", full_id),
    Haplotype = sub(".*haplotype([0-9]+).*", "\\1", full_id),
    Location = sub("_.*", "", ID),
    Area = "Ohio",
    Year = sub(".*_([0-9]+)\\..*", "\\1", ID),
    Species = "Aedes japonicus",
    focal = TRUE
  ) |>
  select(-full_id) |>
  bind_rows(meta_refs) 


# RECTANGULAR TREE PLOT --------------------------------------------------------
ggtree(iqtree, cex = 0.8) %<+% meta +
  geom_tiplab(
    aes(color = focal),
    align = TRUE,
    linesize = 0.25,
    size = 2.6,
    ) +
  geom_tiplab(
    geom = "label",
    aes(label = Area, fill = Area),
    offset = 0.3,
    align = TRUE,
    linesize = 0,
    size = 2.1,
    fontface = "bold"
  ) +
  scale_color_manual(values = c("grey70", "black"), na.translate = FALSE) +
  scale_fill_brewer(palette = "Pastel2") +
  guides(fill = "none", color = "none") +
  geom_treescale(y = - 5, color = "coral4", fontsize = 4) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "top",
    plot.margin = margin(0.5, 1, 0.05, 0.05, "cm"),
    )

ggsave("results/tree.png", width = 6, height = 8, dpi = 600)

#p +
#  new_scale_fill() +
#  geom_tiplab(
#    geom = "label",
#    aes(label = Year, fill = Year),
#    offset = 0.012,
#    align = TRUE,
#    linesize = 0,
#    size = 2.3,
#    fontface = "bold"
#  ) +
#  guides(fill = "none")


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
