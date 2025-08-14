
# === Script Description =======================================================
# Calculate % identity between groups and within groups
# ==============================================================================

# Clear global environment
rm(list = ls())

# Define a user-specific library directory
setwd("C:/Users/")
getwd()


# === User settings ============================================================
#input_csv     <- "pairwise distance matrix from MEGA on bW-HD1 AA alignment.csv"
input_csv     <- "Distance Data BE-HD2 pairwise from MEGA.csv"
output_prefix <- "hd2_identity"  # all outputs will start with this
locus_pattern <- "(BE[12]-HD2)" # regex to parse locus from sequence names
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
})

# ---- 1) Read MEGA distance matrix -------------------------------------------
# MEGA exports a triangular matrix; keep names exactly as in the file
dist_df <- read_csv(input_csv, show_col_types = FALSE)
# First column is row names; set them and drop the first column from data
rn <- dist_df[[1]]
dist_df <- dist_df[,-1]
dist_mat <- as.matrix(dist_df)
rownames(dist_mat) <- rn
colnames(dist_mat) <- colnames(dist_df)

# Coerce to numeric (MEGA writes blanks/NA in one triangle)
mode(dist_mat) <- "numeric"

# ---- 2) Symmetrize and set diagonal to 0 ------------------------------------
# Fill the upper triangle with the transpose of the lower (or vice versa)
dist_mat[upper.tri(dist_mat)] <- t(dist_mat)[upper.tri(dist_mat)]
diag(dist_mat) <- 0

# ---- 3) Convert distances to % identity -------------------------------------
id_mat <- 100 * (1 - dist_mat)   # % identity
id_mat <- round(id_mat, 1)

# ---- 4) Parse locus groups from sequence names ------------------------------
seq_names <- colnames(id_mat)
groups <- str_extract(seq_names, locus_pattern)
if (any(is.na(groups))) {
  warning("Some sequence names did not match locus_pattern. Update 'locus_pattern' if needed.")
}

# ---- 5) Build a long table of upper-triangle pairs --------------------------
ut <- which(upper.tri(id_mat), arr.ind = TRUE)
pairs_tbl <- tibble(
  seq1 = seq_names[ut[,1]],
  seq2 = seq_names[ut[,2]],
  group1 = groups[ut[,1]],
  group2 = groups[ut[,2]],
  pct_identity = id_mat[ut]
)

# ---- 6) Summaries: within-group and between-group ---------------------------
within_summary <- pairs_tbl %>%
  filter(group1 == group2) %>%
  group_by(group = group1) %>%
  summarise(
    n_pairs = n(),
    mean = round(mean(pct_identity), 1),
    sd   = round(sd(pct_identity), 1),
    min  = round(min(pct_identity), 1),
    max  = round(max(pct_identity), 1),
    .groups = "drop"
  )

between_summary <- pairs_tbl %>%
  filter(group1 != group2) %>%
  mutate(comparison = paste(group1, "vs", group2)) %>%
  group_by(comparison) %>%
  summarise(
    n_pairs = n(),
    mean = round(mean(pct_identity), 1),
    sd   = round(sd(pct_identity), 1),
    min  = round(min(pct_identity), 1),
    max  = round(max(pct_identity), 1),
    .groups = "drop"
  ) %>%
  arrange(comparison)

# ---- 7) Write outputs --------------------------------------------------------
write_csv(as.data.frame(id_mat) %>% mutate(seq = rownames(id_mat), .before = 1),
          paste0(output_prefix, "_identity_matrix.csv"))
write_csv(pairs_tbl,          paste0(output_prefix, "_pairs.csv"))
write_csv(within_summary,     paste0(output_prefix, "_within_summary.csv"))
write_csv(between_summary,    paste0(output_prefix, "_between_summary.csv"))

# ---- 8) Quick plots ----------------------------------------------
# Within-group boxplot
p1 =ggplot(pairs_tbl %>% filter(group1 == group2),
       aes(x = group1, y = pct_identity)) +
  geom_boxplot() + geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  labs(x = "Genes", y = "Pairwise % identity") +
  theme_minimal() +
  ggtitle("Within-group % identity") +
  theme(plot.title = element_text(hjust = 0.5))

# Between-group boxplot
p2 = ggplot(pairs_tbl %>% filter(group1 != group2),
       aes(x = paste(group1, "vs", group2), y = pct_identity)) +
  geom_boxplot() + geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  labs(x = "Comparison", y = "Pairwise % identity") +
  theme_minimal() +
  ggtitle("Between-group % identity") +
  theme(plot.title = element_text(hjust = 0.5))
  
grid.arrange(p1, p2, ncol=2)

# ---- 9) Print a concise summary to console ----------------------------------
cat("\nWithin-group % identity summary:\n")
print(within_summary)
cat("\nBetween-group % identity summary:\n")
print(between_summary)
