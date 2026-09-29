# stuff to run to start R & Rstudio in terminal
# export PATH=$PATH:/project/inbreh/conda_envs/gggenesR/bin/
# export RSTUDIO_WHICH_R=/project/inbreh/conda_envs/gggenesR/bin/R
# module load rstudio/2024.04.2
# rstudio


library(readxl)
library(dplyr)
library(ggplot2)
library(gggenes)


TYPE_I_XLSX  <- "/project/inbreh/frog_kaps/bombina_typeI_HMM.xlsx"
TYPE_II_XLSX <- "/project/inbreh/frog_kaps/bombina_typeII_HMM.xlsx"
OUT_DIR      <- "/project/inbreh/frog_kaps/gene_cluster_plots"
dir.create(OUT_DIR, showWarnings = FALSE)






# Cluster boundaries from GFF anchors
TYPE_I_CLUSTER  <- list(chrom = "NC_069499.1",
                        left_anchor  = list(name = "SMARCE1", start = 1128761400, end = 1128851600, strand = "-"),
                        right_anchor = list(name = "EIF1",    start = 1134391316, end = 1134394959, strand = "+"))

TYPE_II_CLUSTER <- list(chrom = "NC_069501.1",
                        left_anchor  = list(name = "TNS2",     start = 1062590657, end = 1062972746, strand = "-"),
                        right_anchor = list(name = "BCDIN3D",  start = 1066695348, end = 1066715532, strand = "-"))

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
prep_keratin_df <- function(xlsx_path, cluster, chrom_col, gene_col) {
  df <- read_excel(xlsx_path)
  # normalize column names
  df <- df %>%
    rename(chrom       = !!chrom_col,
           display_gene = !!gene_col) %>%
    filter(as.character(chrom) == sub("NC_", "", cluster$chrom) |   # in case chrom stored as just number
             as.character(chrom) == cluster$chrom |
             as.character(chrom) == "1" & cluster$chrom == "NC_069499.1" |
             as.character(chrom) == "3" & cluster$chrom == "NC_069501.1") %>%
    filter(start_1based >= cluster$left_anchor$start,
           end_1based   <= cluster$right_anchor$end) %>%
    mutate(strand_val = ifelse(strand == "+", 1, -1),
           category   = "keratin")
  df
}

add_anchor <- function(df, anchor, category = "anchor") {
  bind_rows(df,
            tibble(display_gene = anchor$name,
                   start_1based = anchor$start,
                   end_1based   = anchor$end,
                   strand       = anchor$strand,
                   strand_val   = ifelse(anchor$strand == "+", 1, -1),
                   category     = category,
                   keratin_tree_ID = NA_character_))
}

# --------------------------------------------------------------------------
# Prepare Type I
# --------------------------------------------------------------------------
type_i <- prep_keratin_df(TYPE_I_XLSX, TYPE_I_CLUSTER, "Chr", "gene")

# Highlight category based on keratin_tree_ID
type_i <- type_i %>%
  mutate(category = case_when(
    grepl("Krt34", keratin_tree_ID, ignore.case = TRUE) ~ "Krt34 paralog",
    TRUE                                                 ~ "keratin"
  ))

type_i <- type_i %>%
  add_anchor(TYPE_I_CLUSTER$left_anchor) %>%
  add_anchor(TYPE_I_CLUSTER$right_anchor)

# --------------------------------------------------------------------------
# Prepare Type II
# --------------------------------------------------------------------------
type_ii <- prep_keratin_df(TYPE_II_XLSX, TYPE_II_CLUSTER, "chromosome", "Gene")

type_ii <- type_ii %>%
  mutate(category = case_when(
    grepl("Krt59", keratin_tree_ID, ignore.case = TRUE) ~ "Krt59 paralog",
    TRUE                                                 ~ "keratin"
  ))

type_ii <- type_ii %>%
  add_anchor(TYPE_II_CLUSTER$left_anchor) %>%
  add_anchor(TYPE_II_CLUSTER$right_anchor)



# --------------------------------------------------------------------------
# Convert to evenly-spaced synthetic coordinates for readability
# --------------------------------------------------------------------------
make_synthetic_coords <- function(df, gene_width = 8, gene_gap = 4) {
  # Sort by real start, then assign evenly-spaced synthetic coordinates
  df <- df %>% arrange(start_1based)
  step <- gene_width + gene_gap
  df$syn_start <- seq(0, by = step, length.out = nrow(df))
  df$syn_end   <- df$syn_start + gene_width
  df
}

type_i  <- make_synthetic_coords(type_i)
type_ii <- make_synthetic_coords(type_ii)

# --------------------------------------------------------------------------
# Plotting function (uses synthetic coords)
# --------------------------------------------------------------------------
plot_cluster <- function(df, title, highlight_label, highlight_color) {
  fill_values <- c(
    "anchor"  = "#333333",
    "keratin" = "#B8B8B8"
  )
  fill_values[highlight_label] <- highlight_color
  
  ggplot(df, aes(xmin = syn_start, xmax = syn_end,
                 y = 1, forward = strand_val, fill = category,
                 label = display_gene)) +
    geom_gene_arrow(arrowhead_height  = grid::unit(8, "mm"),
                    arrowhead_width   = grid::unit(4, "mm"),
                    arrow_body_height = grid::unit(8, "mm")) +
    geom_gene_label(align = "centre", grow = FALSE,
                    height = grid::unit(6, "mm"),
                    padding.x = grid::unit(0.5, "mm")) +
    scale_fill_manual(values = fill_values, name = NULL) +
    labs(title    = title,
         subtitle = "Gene order preserved; sizes and spacing not to scale",
         x = NULL, y = NULL) +
    theme_genes() +
    theme(axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          axis.text.y     = element_blank(),
          axis.ticks.y    = element_blank(),
          legend.position = "top",
          plot.title      = element_text(face = "bold"))
}

# --------------------------------------------------------------------------
# Draw and save
# --------------------------------------------------------------------------
p1 <- plot_cluster(type_i,
                   title = sprintf("Bombina bombina — type I keratin cluster (%s)", TYPE_I_CLUSTER$chrom),
                   highlight_label = "Krt34 paralog",
                   highlight_color = "#C1272D")

p2 <- plot_cluster(type_ii,
                   title = sprintf("Bombina bombina — type II keratin cluster (%s)", TYPE_II_CLUSTER$chrom),
                   highlight_label = "Krt59 paralog",
                   highlight_color = "#2A6A9E")

n_i  <- nrow(type_i)
n_ii <- nrow(type_ii)

ggsave(file.path(OUT_DIR, "bombina_typeI_cluster.pdf"),  p1,
       width = max(10, n_i * 0.4),  height = 3.5, limitsize = FALSE)
ggsave(file.path(OUT_DIR, "bombina_typeII_cluster.pdf"), p2,
       width = max(10, n_ii * 0.4), height = 3.5, limitsize = FALSE)

ggsave(file.path(OUT_DIR, "bombina_typeI_cluster.png"),  p1,
       width = max(10, n_i * 0.4),  height = 3.5, dpi = 200, limitsize = FALSE)
ggsave(file.path(OUT_DIR, "bombina_typeII_cluster.png"), p2,
       width = max(10, n_ii * 0.4), height = 3.5, dpi = 200, limitsize = FALSE)

cat("Type I:  ", n_i,  "features (", sum(type_i$category  == "Krt34 paralog"), "Krt34 paralogs)\n")
cat("Type II: ", n_ii, "features (", sum(type_ii$category == "Krt59 paralog"), "Krt59 paralogs)\n")









# --------------------------------------------------------------------------
# Also draw versions with real coordinates preserved (spacing + size to scale)
# --------------------------------------------------------------------------
plot_cluster_scaled <- function(df, title, highlight_label, highlight_color) {
  fill_values <- c(
    "anchor"  = "#333333",
    "keratin" = "#B8B8B8"
  )
  fill_values[highlight_label] <- highlight_color
  
  ggplot(df, aes(xmin = start_1based, xmax = end_1based,
                 y = 1, forward = strand_val,
                 fill = category, color = category,
                 label = display_gene)) +
    geom_gene_arrow(arrowhead_height  = grid::unit(6, "mm"),
                    arrowhead_width   = grid::unit(3, "mm"),
                    arrow_body_height = grid::unit(6, "mm")) +
    geom_gene_label(align = "centre", grow = FALSE,
                    height = grid::unit(5, "mm"),
                    padding.x = grid::unit(0.5, "mm")) +
    scale_fill_manual(values = fill_values,  name = NULL) +
    scale_color_manual(values = fill_values, name = NULL, guide = "none") +
    scale_x_continuous(labels = function(x) paste0(round(x / 1e6, 2), " Mb")) +
    labs(title = title, x = NULL, y = NULL) +
    theme_genes() +
    theme(axis.text.y     = element_blank(),
          axis.ticks.y    = element_blank(),
          legend.position = "top",
          plot.title      = element_text(face = "bold"))
}

p1_scaled <- plot_cluster_scaled(type_i,
                                 title = sprintf("Bombina bombina — type I keratin cluster (%s, to scale)", TYPE_I_CLUSTER$chrom),
                                 highlight_label = "Krt34 paralog",
                                 highlight_color = "#C1272D")

p2_scaled <- plot_cluster_scaled(type_ii,
                                 title = sprintf("Bombina bombina — type II keratin cluster (%s, to scale)", TYPE_II_CLUSTER$chrom),
                                 highlight_label = "Krt59 paralog",
                                 highlight_color = "#2A6A9E")

ggsave(file.path(OUT_DIR, "bombina_typeI_cluster_toscale.pdf"),  p1_scaled,
       width = max(10, n_i * 0.4),  height = 3.5, limitsize = FALSE)
ggsave(file.path(OUT_DIR, "bombina_typeII_cluster_toscale.pdf"), p2_scaled,
       width = max(10, n_ii * 0.4), height = 3.5, limitsize = FALSE)
