# 2025-07-06, Jelmer Poelstra at OSC-Cardinal, R 4.4.0

# SETUP ------------------------------------------------------------------------
# Load packages
library(ape)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)

# Settings
outgroup <- "DQ176828.2"
focal_locs <- c("Cleveland", "Wooster", "other")
loc_cols <- c("purple3", "goldenrod3", "grey80")

# Define input files
tree_file <- "results/NAD4/iqtree/NAD4_withrefs_aln.treefile"

# Define output files
outdir <- "results/NAD4/plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
plotfile_png <- file.path(outdir, "NAD4-tree.png")
plotfile_svg <- file.path(outdir, "NAD4-tree.svg")

# Read the input files
tree <- read.tree(tree_file)

# Set the outgroup
tree <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE)
tree$edge.length[which.max(tree$edge.length)] <- 0.1


# METADATA PREP ----------------------------------------------------------------
meta <- tibble(ID = tree$tip.label) |>
  mutate(
    plotlab = sub("WB_", "", ID),
    location = sub("_.*", "", plotlab),
    plotlab = sub("_[0-9]+$", "", plotlab),
    plotlab = sub("_", " ", plotlab)
  ) |>
  mutate(
    location = ifelse(location == "DQ176828.2", "other", location),
    location = factor(location, levels = focal_locs),
    plotlab = ifelse(plotlab == "DQ176828.2", "(Ae. aegypti)", plotlab),
    )


# PLOT -------------------------------------------------------------------------
# Plot the tree - circular
ggtree(tree, size = 0.5, layout = "circular") %<+%
  meta +
  geom_tiplab(align = TRUE, linesize = 0.25, size = 0) +
  geom_tiplab(
    aes(label = plotlab, color = location),
    align = TRUE,
    linesize = 0,
    size = 3,
    offset = 0.01
  ) +
  geom_rootedge(rootedge = 0.005) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill = location),
    color = "black",
    pwidth = 0.01,
    width = 0.01
  ) +
  scale_fill_manual(values = loc_cols, name = "Location") +
  scale_color_manual(values = loc_cols, name = "Location") +
  guides(color = "none") +
  theme(
    plot.margin = margin(1, 1, 1, 0, "cm"),
    legend.title = element_text(face = "bold")
  )

ggsave(plotfile_png, width = 8, height = 8, dpi = 600)
ggsave(plotfile_svg, width = 14, height = 7)
