setwd("")

library(ape)
library(phytools)
library(ggtree)
library(ggplot2)


#type I

# Read the tree
tr <- read.tree("type1_keratins_modified.fasta.treefile")

# Root with the outgroups, then remove them from the plotted tree
outgroups <- c("axolotl_krt80", "caecilian_krt80", "human_KRT80")

missing_outgroups <- setdiff(outgroups, tr$tip.label)
if (length(missing_outgroups) > 0) {
  stop("Outgroups not found: ", paste(missing_outgroups, collapse = ", "))
}

tr_plot <- tr |>
  unroot() |>
  root(outgroup = outgroups, resolve.root = TRUE) |>
  drop.tip(outgroups) |>
  ladderize(right = FALSE)

# Shorten all Leptobrachium tip labels
tr_plot$tip.label <- sub("^Leptobrachium_.*$", "Leptobrachium", tr_plot$tip.label)

# Transform branch lengths and make the tree ultrametric
tr_plot$edge.length <- log1p(tr_plot$edge.length)
tr_plot <- chronos(tr_plot)
tr_plot <- force.ultrametric(tr_plot, method = "extend")

# use UFBoot labels
get_ufboot <- function(label) {
  values <- strsplit(as.character(label), "/", fixed = TRUE)[[1]]
  suppressWarnings(as.numeric(if (length(values) >= 2) values[2] else values[1]))
}

support_data <- data.frame(
  node = Ntip(tr_plot) + seq_len(tr_plot$Nnode),
  bootstrap = vapply(tr_plot$node.label, get_ufboot, numeric(1))
)

# Circular ultrametric tree with black circles at bootstrap support >= 95.
p <- ggtree(
  tr_plot,
  layout = "circular",
  linewidth = 0.35,
  open.angle = 8
) %<+% support_data +
  geom_tiplab2(
    size = 3,
    offset = 0.002,
    linesize = 0.2,
    color = "grey15"
  ) +
  geom_point2(
    aes(subset = !isTip & bootstrap >= 95),
    size = 1.5,
    color = "black"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(20, 100, 20, 20),
    plot.background = element_rect(fill = "white", color = NA)
  )

p


#type II

# Read the tree
tr <- read.tree("type2_keratins_modified.fasta.treefile")

# Root with the outgroups, then remove them from the plotted tree
outgroups <- c("human_KRT23", "axolotl_krt23.1", "caecilian_krt23")

missing_outgroups <- setdiff(outgroups, tr$tip.label)
if (length(missing_outgroups) > 0) {
  stop("Outgroups not found: ", paste(missing_outgroups, collapse = ", "))
}

tr_plot <- tr |>
  unroot() |>
  root(outgroup = outgroups, resolve.root = TRUE) |>
  drop.tip(outgroups) |>
  ladderize(right = FALSE)

# Shorten all Leptobrachium tip labels
tr_plot$tip.label <- sub("^Leptobrachium_.*$", "Leptobrachium", tr_plot$tip.label)

# Make the tree ultrametric
tr_plot <- force.ultrametric(tr_plot, method = "extend")

# use UFBoot
get_ufboot <- function(label) {
  values <- strsplit(as.character(label), "/", fixed = TRUE)[[1]]
  suppressWarnings(as.numeric(if (length(values) >= 2) values[2] else values[1]))
}

support_data <- data.frame(
  node = Ntip(tr_plot) + seq_len(tr_plot$Nnode),
  bootstrap = vapply(tr_plot$node.label, get_ufboot, numeric(1))
)

# Circular ultrametric tree with black circles at bootstrap support >= 95.
p <- ggtree(
  tr_plot,
  layout = "circular",
  linewidth = 0.35,
  open.angle = 8
) %<+% support_data +
  geom_tiplab2(
    size = 3,
    offset = 0.002,
    linesize = 0.2,
    color = "grey15"
  ) +
  geom_point2(
    aes(subset = !isTip & bootstrap >= 95),
    size = 1.5,
    color = "black"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(20, 100, 20, 20),
    plot.background = element_rect(fill = "white", color = NA)
  )

p
