# 2025-07-06, Jelmer Poelstra at OSC-Cardinal, R 4.4.0

# SETUP ------------------------------------------------------------------------
# Load packages
library(ape)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)

# Settings
outgroup <- "HQ944971.1"
focal_locs <- c("Cleveland", "Medina", "Wooster", "other")
loc_cols <- c("purple3", "chartreuse4", "goldenrod3", "grey80")

# Define input files
tree_file <- "results/COI/iqtree/COI_withrefs_aln.treefile"
meta_refs_file <- "metadata/refs_meta.tsv"

# Define output files
outdir <- "results/COI/plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
plotfile_png <- file.path(outdir, "COI-tree.png")
plotfile_svg <- file.path(outdir, "COI-tree.svg")
plotfile_lab_png <- file.path(outdir, "COI-tree_labs.png")
plotfile_lab_svg <- file.path(outdir, "COI-tree_labs.svg")
plotfile_rect_png <- file.path(outdir, "COI-tree_rect.png")
plotfile_rect_svg <- file.path(outdir, "COI-tree_rect.svg")

# Read the input files
meta_refs_init <- read_tsv(meta_refs_file, show_col_types = FALSE)
tree <- read.tree(tree_file)

# Set the outgroup
tree <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE)

# Adjust edge length of outgroup because it is extremely long
tree$edge.length[which.max(tree$edge.length)] <- 0.2


# METADATA PREP ----------------------------------------------------------------
# Prep tip labels
tree$tip.label <- tree$tip.label |>
  gsub("\\*+", "", x = _) |>
  gsub("_([0-9]+$)", ".\\1", x = _) |>
  gsub("^ON[0-9]+_", "", x = _) |>
  sub("HQ944971.1", "HQ944971.1_Aedes canadensis", x = _) |>
  sub("GQ254798.1", "GQ254798.1_Iowa_2008", x = _) |> 
  sub("GQ254793.1_Iowa_USA.2007", "GQ254793.1_Illinois_2007", x = _) |>
  sub("\\.(20[0-9]+)$", "_\\1", x = _)

# Prep metadata for ref sequences
meta_refs <- meta_refs_init |>
  mutate(
    ID = sub("HQ944971.1", "HQ944971.1_Aedes canadensis", ID),
    ID = sub("GQ254798.1", "GQ254798.1_Iowa_2008", ID), 
    ID = sub("GQ254793.1_Iowa_USA.2007", "GQ254793.1_Illinois_2007", ID),
    plotlab = sub(".*\\.1_", "", ID),
    plotlab = sub("Aedes canadensis", "(Aedes canadensis)", plotlab),
    Year = as.character(Year),
    is_focal = FALSE,
    Area = Location
  ) 

# Prep metadata for own sequences
meta_focal <- tibble(ID = tree$tip.label) |>
  mutate(
    Location = sub("_.*", "", ID),
    Area = "Ohio",
    Year = sub(".*_([0-9]+)\\..*", "\\1", ID),
    Species = "Aedes japonicus",
    is_focal = TRUE,
    plotlab = ID,
    plotlab = sub("\\.[0-9]+$", "", plotlab)
  ) |>
  filter(Location %in% focal_locs)

# Combine  
meta <- bind_rows(meta_refs, meta_focal) |>
  mutate(
    loc_focal = case_when(
      grepl("Wooster", plotlab) ~ "Wooster",
      grepl("Cleveland", plotlab) ~ "Cleveland",
      grepl("Medina", plotlab) ~ "Medina",
      .default = "other"
    ),
    loc_focal = factor(loc_focal, levels = focal_locs),
    plotlab = sub("_", " ", plotlab)
  )


# PLOT -------------------------------------------------------------------------
# Plot the tree - circular
p <- ggtree(tree, size = 0.5, layout = "circular") %<+%
  meta +
  geom_tiplab(align = TRUE, linesize = 0.25, size = 0) +
  geom_rootedge(rootedge = 0.005) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill = loc_focal),
    color = "black",
    pwidth = 0.03,
    width = 0.03
  ) +
  scale_fill_manual(values = loc_cols, name = "Location") +
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),
    legend.title = element_text(face = "bold")
  )
p

ggsave(plotfile_png, width = 8, height = 7, dpi = 400)
ggsave(plotfile_svg, width = 8, height = 7)

# With labels
p +
  geom_tiplab(
    aes(label = plotlab, color = is_focal),
    align = TRUE,
    linesize = 0,
    size = 3,
    offset = 0.03
  ) +
  scale_color_manual(values = c("grey70", "black"), na.translate = FALSE) +
  guides(color = "none") +
  theme(
    legend.box.spacing = unit(50, "pt"),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm")
  )

ggsave(plotfile_lab_png, width = 8, height = 7, dpi = 400)
ggsave(plotfile_lab_svg, width = 8, height = 7)


# ALTERNATIVES -----------------------------------------------------------------
# Plot the tree - rectangular
ggtree(tree, size = 0.5, layout = "rectangular") %<+%
  meta +
  geom_tiplab(
    aes(label = plotlab, color = loc_focal),
    align = TRUE,
    linesize = 0.2,
    size = 3
  ) +
  geom_treescale(x = 0.01, y = length(tree$tip.label) - 3, color = "grey50") +
  geom_rootedge(rootedge = 0.005) +
  scale_color_manual(values = loc_cols) +
  coord_cartesian(clip = "off") +
  guides(color = "none") +
  theme(plot.margin = margin(0.2, 5, 0.2, 0.2, "cm"))

ggsave(plotfile_rect_png, width = 8, height = 7, dpi = 400)
ggsave(plotfile_rect_svg, width = 8, height = 7)

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

# Alternative
#ggtree(tree, cex = 0.8) %<+%
#  meta +
#  geom_tiplab(
#    aes(color = focal),
#    align = TRUE,
#    linesize = 0.25,
#    size = 2.6,
#  ) +
#  geom_tiplab(
#    geom = "label",
#    aes(label = Area, fill = Area),
#    offset = 0.3,
#    align = TRUE,
#    linesize = 0,
#    size = 2.1,
#    fontface = "bold"
#  ) +
#  scale_color_manual(values = c("grey70", "black"), na.translate = FALSE) +
#  scale_fill_brewer(palette = "Pastel2") +
#  guides(fill = "none", color = "none") +
#  geom_treescale(y = -5, color = "coral4", fontsize = 4) +
#  coord_cartesian(clip = "off") +
#  theme(
#    legend.position = "top",
#    plot.margin = margin(0.5, 1, 0.05, 0.05, "cm"),
#  )