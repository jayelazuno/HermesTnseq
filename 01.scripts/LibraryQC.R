suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(purrr)
  library(gt)
  library(tidyverse)
})

# QC_library_stats)

# ============================================================
# Load stats.csv files


summary_dir <- here("05.summary")

stats_files <- list.files(
  summary_dir,
  pattern = "^stats\\.csv$",
  recursive = TRUE,
  full.names = TRUE
) %>%
  keep(~ dirname(.x) != summary_dir)

stopifnot(length(stats_files) > 0)

stats_tbl <- map_dfr(stats_files, function(f) {
  
  sample <- basename(dirname(f))
  
  read_csv(f, show_col_types = FALSE) %>%
    mutate(sample = sample)
  
})

#####################################################################
# Define condition 

stats_tbl <- stats_tbl %>%
  mutate(
    condition = if_else(str_detect(sample, "parent"), "parent", "H2O2-treated-facs"),
    condition = factor(condition, levels = c("parent", "H2O2-treated-facs"))
  )

# Clean numeric columns

percent_cols <- c(
  "% of hits in features",
  "% of intergenic hits",
  "% of features hit"
)

stats_tbl <- stats_tbl %>%
  mutate(
    across(
      all_of(percent_cols),
      ~ round(as.numeric(str_remove(.x, "%")), 1)
    ),
    `Mean Reads Per Hit` = round(as.numeric(`Mean Reads Per Hit`), 1),
    `Mean reads per hit in feature` =
      round(as.numeric(`Mean reads per hit in feature`), 1)
  )


# make a shared theme

qc_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      size = 12, color = "black", face = "bold"
    ),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title  = element_text(face = "bold", size = 14),
    strip.text  = element_text(face = "bold", size = 11),
    legend.title = element_blank(),
    legend.text  = element_text(size = 13, color = "black", face = "bold")
  )

# Plot 1: Library complexity


p_complexity <- ggplot(
  stats_tbl,
  aes(x = sample, `Total Hits`, fill = sample)
) +
  geom_col(width = 0.7) +
  facet_wrap(~ condition, scales = "free_x") +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = " ",
    x = " ",
    y = "Total unique insertion sites"
  ) +
  qc_theme +
  theme(legend.position = "none")

# Plot 2: % genomic features hit (FACETED → no sample colors)

p_features_hit <- ggplot(
  stats_tbl,
  aes(x = sample, y = `% of features hit`, fill = sample)
) +
  geom_col(width = 0.7) +
  facet_wrap(~ condition, scales = "free_x") +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title = " ",
    x = " ",
    y = "% features hit"
  ) +
  qc_theme +
  theme(legend.position = "none")


# Plot 3: Mean reads per insertion site


p_reads_per_hit <- ggplot(
  stats_tbl,
  aes(x = sample, y = `Mean Reads Per Hit`, fill = sample)
) +
  geom_col(width = 0.7) +
  scale_y_log10() +
  scale_fill_manual(values = sample_colors) +
  facet_wrap(~ condition, scales = "free_x") +
  labs(
    title = " ",
    x = " ",
    y = "Mean reads per hit (log10)"
  ) +
  qc_theme +
  theme(legend.position = "none")


# Plot 4: Feature vs intergenic distribution


stats_long <- stats_tbl %>%
  select(
    sample, condition,
    `% of hits in features`,
    `% of intergenic hits`
  ) %>%
  pivot_longer(
    cols = starts_with("%"),
    names_to = "category",
    values_to = "percent"
  )

p_feature_intergenic <- ggplot(
  stats_long,
  aes(x = sample, y = percent, fill = category)
) +
  geom_col(width = 0.7) +
  facet_wrap(~ condition, scales = "free_x") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title = " ",
    x = " ",
    y = "Percent of hits"
  ) +
  qc_theme


# Print plots

p_complexity
p_features_hit
p_reads_per_hit
p_feature_intergenic

# ggsave(file.path(plots_qc, "library-complexity.pdf"), p_complexity, width = 10, height = 10)
# 
# ggsave(file.path(plots_qc, "percent_features_hit.pdf"), p_features_hit, width = 10, height = 10)
# 
# ggsave(file.path(plots_qc, "percent_reads_per_hit.pdf"), p_reads_per_hit, width = 10, height = 10)
# 
# ggsave(file.path(plots_qc, "percent_feature_intergenic.pdf"), p_feature_intergenic, width = 10, height = 10)

############################################################
# Bin reads into 1k and 10kb bin size and plot hits and reads 

# Sample order

samples <- c(
  "yH298-parent-pool1",
  "yH298-parent-pool2",
  "yH299-parent-pool3",
  "yH299-parent-pool4",
  "yH298-H2O2-treated-facs-pool1",
  "yH298-H2O2-treated-facs-pool2",
  "yH299-H2O2-treated-facs-pool3",
  "yH299-H2O2-treated-facs-pool4"
)
# Locate *.all_hits.csv

all_hits_files <- list.files(
  summary_dir,
  pattern = "\\.all_hits\\.csv$",
  recursive = TRUE,
  full.names = TRUE
) %>%
  keep(~ dirname(.x) != summary_dir)

stopifnot(length(all_hits_files) > 0)

# Read hits

hits_tbl <- map_dfr(all_hits_files, function(f) {
  sample <- basename(dirname(f))
  read_csv(f, show_col_types = FALSE) %>%
    mutate(sample = sample)
})

stopifnot(all(c("Chromosome", "Position", "Reads") %in% colnames(hits_tbl)))

# Chromosome relabeling (ONCE)

chrom_map <- tibble(
  Chromosome = c(
    "CP048230.1","CP048231.1","CP048232.1","CP048233.1",
    "CP048234.1","CP048235.1","CP048236.1","CP048237.1",
    "CP048238.1","CP048239.1","CP048240.1","CP048241.1",
    "CP048242.1"
  ),
  ChromLabel = paste("Chr", LETTERS[1:13])
)

hits_tbl <- hits_tbl %>%
  left_join(chrom_map, by = "Chromosome") %>%
  mutate(
    sample     = factor(sample, levels = samples),
    Chromosome = factor(ChromLabel, levels = chrom_map$ChromLabel)
  ) %>%
  select(-ChromLabel)

# Bin sizes

bin_sizes <- c(1000, 10000)

# Bin hits

binned_tbl <- map_dfr(bin_sizes, function(bin_size) {
  hits_tbl %>%
    mutate(
      bin_size = bin_size,
      bin = floor(Position / bin_size)
    ) %>%
    group_by(sample, Chromosome, bin_size, bin) %>%
    summarise(
      hits  = n(),
      reads = sum(Reads),
      .groups = "drop"
    )
})


# Plot function 

plot_chrom_bins <- function(df, value_col, title_suffix) {
  ggplot(df, aes(x = bin, y = .data[[value_col]])) +
    geom_col(width = 1, fill = "steelblue") +
    facet_grid(
      Chromosome ~ sample,
      scales = "free_x",
      space  = "free_x"
    ) +
    scale_y_log10() +
    labs(
      title = paste0(title_suffix, " (raw, binned)"),
      x = "Genomic bin",
      y = paste0("log10(", value_col, " per bin)")
    ) +
    theme_bw() +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(face = "bold"),
      strip.text.x = element_text(angle = 90, face = "bold", size = 11),
      strip.text.y = element_text(face = "bold"),
      axis.title   = element_text(face = "bold", size = 14)
    )
}

# Plots


p_reads_1kb  <- binned_tbl %>% filter(bin_size == 1000)  %>% plot_chrom_bins("reads", "Reads per 1 kb bin")
p_reads_10kb <- binned_tbl %>% filter(bin_size == 10000) %>% plot_chrom_bins("reads", "Reads per 10 kb bin")

p_hits_1kb   <- binned_tbl %>% filter(bin_size == 1000)  %>% plot_chrom_bins("hits", "Hits per 1 kb bin")
p_hits_10kb  <- binned_tbl %>% filter(bin_size == 10000) %>% plot_chrom_bins("hits", "Hits per 10 kb bin")

p_reads_1kb
p_reads_10kb
p_hits_1kb
p_hits_10kb

# ggsave(file.path(plots_qc, "reads_binned_1kb.pdf"), p_reads_1kb, width = 20, height = 15)
# 
# ggsave(file.path(plots_qc, "reads_binned_10kb.pdf"), p_reads_10kb, width = 20, height = 15)
# 
# 
# ggsave(file.path(plots_qc, "hit_binned_1kb.pdf"), p_hits_1kb, width = 20, height = 15)
# 
# ggsave(file.path(plots_qc, "hit_binned_10kb.pdf"), p_hits_10kb, width = 20, height = 15)


################################################################# 

# plot pairwise correlation among parental pools and treated pools 


parent_samples  <- samples[grepl("parent", samples)]
treated_samples <- samples[grepl("treated", samples)]

metric_unnorm <- "reads"   # can switch to reads, NI, etc.

# Reader (same logic, as previous)

read_feature_table_unnorm <- function(f, sample) {
  read_csv(f, skip = 1, show_col_types = FALSE) %>%
    mutate(sample = sample)
}

feature_file_unnorm <- function(sample) {
  file.path(
    summary_dir,
    sample,
    paste0(sample, ".feature_table.RDF_1.csv")
  )
}

make_gene_vec_unnorm <- function(df, metric = "reads") {
  df %>%
    select(standard_name, value = all_of(metric)) %>%
    mutate(value = as.numeric(value)) %>%
    distinct(standard_name, .keep_all = TRUE)
}

# Load all unnormalized feature tables
feature_tbl_unnorm <- map_dfr(samples, function(s) {
  f <- feature_file_unnorm(s)
  stopifnot(file.exists(f))
  read_feature_table_unnorm(f, s)
})

# Pairwise correlation helper (unnormalized)

pairwise_corr_unnorm <- function(sample_vec, label) {
  
  combn(sample_vec, 2, simplify = FALSE) %>%
    map_dfr(function(pair) {
      
      a <- feature_tbl_unnorm %>%
        filter(sample == pair[1]) %>%
        make_gene_vec_unnorm(metric_unnorm)
      
      b <- feature_tbl_unnorm %>%
        filter(sample == pair[2]) %>%
        make_gene_vec_unnorm(metric_unnorm)
      
      joined <- inner_join(a, b, by = "standard_name",
                           suffix = c("_1", "_2"))
      
      tibble(
        group         = label,
        sample_1      = pair[1],
        sample_2      = pair[2],
        n_genes       = nrow(joined),
        spearman_rho  = suppressWarnings(
          cor(joined$value_1, joined$value_2, method = "spearman")
        ),
        pearson_r     = suppressWarnings(
          cor(joined$value_1, joined$value_2, method = "pearson")
        )
      )
    })
}

# Compute correlations (unnormalized)

corr_tbl_unnorm <- bind_rows(
  pairwise_corr_unnorm(parent_samples,  "parent"),
  pairwise_corr_unnorm(treated_samples, "H2O2-treated-facs")
)

corr_tbl_unnorm <- corr_tbl_unnorm %>%
  mutate(
    group = factor(group, levels = c("parent", "H2O2-treated-facs"))
  )

corr_tbl_unnorm

p_spearman_unnorm <- ggplot(
  corr_tbl_unnorm,
  aes(x = interaction(sample_1, sample_2), y = spearman_rho)
) +
  geom_point(size = 3) +
  facet_wrap(~ group, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = paste0(
      "Pairwise Spearman correlations (unnormalized, metric: ",
      metric_unnorm, ")"
    ),
    x = "Library pair",
    y = "Spearman ρ"
  ) +
  qc_theme

p_spearman_unnorm


p_pearson_unnorm <- ggplot(
  corr_tbl_unnorm,
  aes(x = interaction(sample_1, sample_2), y = pearson_r)
) +
  geom_point(size = 3) +
  facet_wrap(~ group, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = paste0(
      "Pairwise Pearson correlations (unnormalized, metric: ",
      metric_unnorm, ")"
    ),
    x = "Library pair",
    y = "Pearson r"
  ) +
  qc_theme

p_pearson_unnorm

ggsave(file.path(plots_qc, "pairwise_spearman_unnorm_reads.pdf"), p_spearman_unnorm, width = 10, height = 10)

ggsave(file.path(plots_qc, "pairwise_pearson_unnorm_reads.pdf"), p_pearson_unnorm, width = 10, height = 10)




###########################################################################################
# A  read in the the midlc estimate files to visualize when unique insertion reaches saturation point

qc_dir   <- here("07.QC")
plots_qc <- here("07.QC", "plots")
dir.create(plots_qc, recursive = TRUE, showWarnings = FALSE)

# sample dirs 
sample_dirs <- file.path(qc_dir, c(
  "yH298-H2O2-treated-facs-pool1",
  "yH298-H2O2-treated-facs-pool2",
  "yH298-parent-pool1",
  "yH298-parent-pool2",
  "yH299-H2O2-treated-facs-pool3",
  "yH299-H2O2-treated-facs-pool4",
  "yH299-parent-pool3",
  "yH299-parent-pool4"
))

missing <- sample_dirs[!dir.exists(sample_dirs)]
if (length(missing) > 0) {
  stop("These sample dirs do not exist:\n", paste(missing, collapse = "\n"))
}

midlc_files <- unlist(lapply(sample_dirs, function(d) {
  list.files(d, pattern = "\\.midlc\\.csv$", recursive = TRUE, full.names = TRUE)
}), use.names = FALSE)

if (length(midlc_files) == 0) {
  stop("No .midlc.csv files found under the specified sample dirs.")
}


read_one_midlc <- function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  names(df) <- str_replace_all(names(df), "\\s+", " ")
  if (!("Reads Sampled" %in% names(df))) {
    stop("Missing 'Reads Sampled' column in: ", f)
  }
  
  sample <- basename(sample_dirs[vapply(sample_dirs, function(d) startsWith(f, d), logical(1))][1])
  if (is.na(sample) || sample == "") {
    # fallback: immediate parent dir of the file
    sample <- basename(dirname(f))
  }
  
  df %>%
    pivot_longer(
      cols = -`Reads Sampled`,
      names_to = "trial",
      values_to = "unique_sites"
    ) %>%
    transmute(
      sample = sample,
      file = f,
      reads_sampled = as.numeric(`Reads Sampled`),
      trial = trial,
      unique_sites = as.numeric(unique_sites)
    )
}

midlc_long <- bind_rows(lapply(midlc_files, read_one_midlc))

# mean curve per sample
midlc_mean <- midlc_long %>%
  group_by(sample, reads_sampled) %>%
  summarize(unique_mean = mean(unique_sites, na.rm = TRUE), .groups = "drop")


# estimate MidLC_est = depth where mean unique reaches 50% of max
midlc_est_tbl <- midlc_mean %>%
  group_by(sample) %>%
  arrange(reads_sampled, .by_group = TRUE) %>%
  summarize(
    umax = max(unique_mean, na.rm = TRUE),
    midlc_est = {
      tgt <- umax / 2
      idx <- which(unique_mean >= tgt)[1]
      if (is.na(idx)) NA_real_ else reads_sampled[idx]
    },
    .groups = "drop"
  ) %>%
  mutate(
    condition = case_when(
      str_detect(sample, "parent") ~ "parent",
      str_detect(sample, "H2O2|treated|facs") ~ "H2O2-treated-facs",
      TRUE ~ "other"
    )
  )

midlc_mean2 <- midlc_mean2 %>%
  mutate(
    condition = factor(condition, levels = c("parent", "H2O2-treated-facs", "other")),
    pool = str_match(as.character(sample), "pool(\\d+)")[,2],
    pool = factor(pool, levels = c("1","2","3","4"))
  ) %>%
  arrange(sample, reads_sampled)

sample_colors <- c(
  "yH298-parent-pool1"            = "#F8766D",
  "yH298-H2O2-treated-facs-pool1" = "#F8766D",
  
  "yH298-parent-pool2"            = "#7CAE00",
  "yH298-H2O2-treated-facs-pool2" = "#7CAE00",
  
  "yH299-parent-pool3"            = "#00BFC4",
  "yH299-H2O2-treated-facs-pool3" = "#00BFC4",
  
  "yH299-parent-pool4"            = "#C77CFF",
  "yH299-H2O2-treated-facs-pool4" = "#C77CFF"
)

p_curve <- ggplot(
  midlc_mean2,
  aes(
    x = reads_sampled,
    y = unique_mean,
    group = sample,
    color = sample      
  )
) +
  geom_line(linewidth = 0.9, alpha = 0.95) +
  geom_point(size = 2.2, alpha = 0.95) +
  facet_wrap(~ condition, scales = "free_y", nrow = 1) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6),
    labels = scales::label_number(
      scale_cut = scales::cut_si(""),
      accuracy = 1
    )
  )+
  scale_color_manual(
    values = sample_colors,
    name = "Library"
  ) +
  labs(
    title = " ",
    x = "Reads sampled (log10 scale)",
    y = "Unique insertion sites"
  ) +
  theme_bw() +
  theme(
    axis.title   = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 14),
    strip.text   = element_text(face = "bold", size = 14),
    plot.title   = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(size = 13)
  )

p_curve_same_y <- p_curve + facet_wrap(~ condition, scales = "fixed", nrow = 1)


p_curve_parent  <- p_curve %+% dplyr::filter(midlc_mean2, condition == "parent") +
  facet_wrap(~ condition, scales = "free_y", nrow = 1)

p_curve_treated <- p_curve %+% dplyr::filter(midlc_mean2, condition == "H2O2-treated-facs") +
  facet_wrap(~ condition, scales = "free_y", nrow = 1)

# ggsave(file.path(plots_qc, "midlc_combined.pdf"),  p_curve_same_y,  width = 14, height = 6)

# ggsave(file.path(plots_qc, "midlc_parent.pdf"),  p_curve_parent,  width = 10, height = 6)

# ggsave(file.path(plots_qc, "midlc_treated.pdf"), p_curve_treated, width = 10, height = 6)



###########################################################################################
# Plot B: MidLC_est bars ... not really informative, I replaced with a simple tableqc_tbl_gt

###########################################################################################

p_est <- ggplot(midlc_est_tbl, aes(x = sample, y = midlc_est, fill = sample)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ condition, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(
    labels = scales::label_number(big.mark = ",")
  ) +
  labs(
    title = " ",
    x = NULL,
    y = "MidLC_est (reads)"
  ) +
  theme_bw() +
  theme(
    axis.title  = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 14, angle = 45, hjust = 1),
    strip.text  = element_text(face = "bold", size = 14),
    plot.title  = element_text(face = "bold", size = 14)
  )

# Separate calls (same y-scale implicitly; separate panels)
p_est_parent  <- p_est %+% dplyr::filter(midlc_est_tbl, condition == "parent") +
  facet_wrap(~ condition, scales = "free_x", nrow = 1)

p_est_treated <- p_est %+% dplyr::filter(midlc_est_tbl, condition == "H2O2-treated-facs") +
  facet_wrap(~ condition, scales = "free_x", nrow = 1)

# Call them:
p_est_parent
# p_est_treated

# Save only when you want:
#ggsave(file.path(plots_qc, "midlc_midlc_est_parent.pdf"),  p_est_parent,  width = 10, height = 4)
# ggsave(file.path(plots_qc, "midlc_midlc_est_treated.pdf"), p_est_treated, width = 10, height = 4)

# normalization of parental pool 

parent_norm_dir <- here("07.QC", "parent_norm")

summary_files <- list.files(
  parent_norm_dir,
  pattern = "library_diagnostics\\.summary\\.T\\d+\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

stopifnot(length(summary_files) > 0)

read_norm_summary <- function(f) {
  tgt <- str_match(basename(f), "T(\\d+)")[,2]
  
  read_csv(f, show_col_types = FALSE) %>%
    mutate(normalize_target = as.integer(tgt))
}

norm_tbl <- map_dfr(summary_files, read_norm_summary)

qc_summary_tbl <- norm_tbl %>%
  group_by(sample) %>%
  summarize(
    total_reads              = first(total_reads),
    unique_sites             = first(unique_sites),
    midlc_est                = first(midlc_est),
    depth_ratio_R_over_midlc = first(depth_ratio_R_over_midlc),
    .groups = "drop"
  ) %>%
  arrange(sample)

qc_tbl_gt <- qc_summary_tbl %>%
  mutate(sample = str_replace(sample, "-", " - ")) %>%
  rename(Library = sample) %>%
  gt() %>%
  tab_header(
    title = md("**Parent library QC summary**"),
    subtitle = md(" ")
  ) %>%
  cols_label(
    Library                  = "Library",
    total_reads              = "Total reads",
    unique_sites             = "Unique sites",
    midlc_est                = "MidLC",
    depth_ratio_R_over_midlc = "Reads / MidLC"
  ) %>%
  fmt_number(
    columns = c(total_reads, unique_sites, midlc_est),
    decimals = 0,
    sep_mark = ","
  ) %>%
  # 1) format the ratio to 1 decimal
  fmt_number(
    columns = depth_ratio_R_over_midlc,
    decimals = 1
  ) %>%
  # 2) append × to the already-formatted text
  text_transform(
    locations = cells_body(columns = depth_ratio_R_over_midlc),
    fn = function(x) paste0(x, "×")
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Library)
  ) %>%
  tab_options(
    table.border.top.style = "solid",
    table.border.top.width = px(8),
    table.border.bottom.style = "solid",
    table.border.bottom.width = px(8),
    column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.width = px(4),
    table.font.names = "Arial",
    table.font.size = px(16),
    column_labels.font.weight = "bold",
    data_row.padding = px(8)
  ) %>%
  opt_table_lines("none") %>%
  opt_row_striping()

qc_tbl_gt

#gtsave(qc_tbl_gt, filename = file.path(plots_qc, "parent_QC_summary_table.pdf"))

###########################################################################################
# C  so after observing the midlc, we can see that the libraries are sequenced to different 
# depth and will need to be normalized. For each sample, the normalization target at depth 𝑇
# was defined as Rtarget =T × MidLC, where MidLC is the estimated mid–library complexity; normalization proceeded only while sufficient reads were available to meet this target.
norm_tbl_plot <- norm_tbl %>%
  arrange(sample, normalize_target)

p_norm_reads <- ggplot(
  norm_tbl_plot,
  aes(
    x = normalize_target,
    y = reads_after_norm,
    group = sample,
    color = sample
  )
) +
  # lines stop naturally at NA (pool2 will stop at T60)
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  
  # points with a tiny x-nudge so overlapping curves become visible
  geom_point(
    size = 2.4,
    alpha = 0.95,
    position = position_dodge(width = 1.2),
    na.rm = TRUE
  ) +
  
  scale_x_continuous(
    breaks = c(20, 30, 40, 50, 60, 80, 100)
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      scale = 1e-6,
      suffix = "M"
    )

  ) +
  labs(
    title = " ",
    x = "Normalization target (T)",
    y = "Reads retained after normalization",
    color = "Library"
  ) +
  theme_bw() +
  theme(
    axis.title   = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 14),
    plot.title   = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(size = 13),
    
    # facet label (grey strip at top)
    strip.text.x = element_text(face = "bold", size = 14)
  )

p_norm_reads 

#p_norm_reads + facet_wrap(~ sample, nrow = 1, scales = "free_y")

p_norm_reads_facet <- p_norm_reads + facet_wrap(~ sample, ncol = 1, scales = "free_y")
#ggsave(file.path(plots_qc, "Reads_retained_after_normalization.pdf"), p_norm_reads_facet, width = 10, height = 15)

###########################################################################################

# Correlation gene-level read counts, insertion counts (hits), and derived features including the 
# neighborhood index (NI) between normalized libraries at increasing normalization targets (T) and the corresponding unnormalized data
###########################################################################################

# =========================

parent_samples <- c(
  "yH298-parent-pool1",
  "yH298-parent-pool2",
  "yH299-parent-pool3",
  "yH299-parent-pool4"
)

targets <- c(20, 30, 40, 50, 60, 80, 100)

unnorm_dir <- here("05.summary")
norm_dir   <- here("05.summary_normalized")

metric <- "hits"  # any gene level metric, reads, hits, NI, etc 

# =========================

read_feature_table <- function(f, sample, T = NA_integer_) {
  read_csv(f, skip = 1, show_col_types = FALSE) %>%
    mutate(
      sample = sample,
      normalize_target = T
    )
}

unnorm_file <- function(sample) {
  file.path(
    unnorm_dir,
    sample,
    paste0(sample, ".feature_table.RDF_1.csv")
  )
}

norm_file <- function(sample, T) {
  file.path(
    norm_dir,
    sprintf("T%03d", T),
    sample,
    paste0(sample, "_normalized.feature_table.RDF_1.csv")
  )
}

make_gene_vec <- function(df, metric = "hits") {
  df %>%
    select(standard_name, value = all_of(metric)) %>%
    mutate(value = as.numeric(value)) %>%
    distinct(standard_name, .keep_all = TRUE)
}

# =========================
# load unnormalized (baseline)

unnorm_tbl <- map_dfr(parent_samples, function(s) {
  f <- unnorm_file(s)
  stopifnot(file.exists(f))
  read_feature_table(f, sample = s, T = 0)
})

# =========================
# compute correlations per sample × T

corr_vs_unnorm <- map_dfr(parent_samples, function(s) {
  
  base <- unnorm_tbl %>%
    filter(sample == s, normalize_target == 0) %>%
    make_gene_vec(metric)
  
  map_dfr(targets, function(T) {
    
    f <- norm_file(s, T)
    
    if (!file.exists(f)) {
      return(NULL)  # <-- skip missing Ts per sample for samples 
    }
    
    cur <- read_feature_table(f, sample = s, T = T) %>%
      make_gene_vec(metric)
    
    joined <- inner_join(
      base, cur,
      by = "standard_name",
      suffix = c("_unnorm", "_norm")
    )
    
    tibble(
      sample = s,
      normalize_target = T,
      n_genes = nrow(joined),
      spearman_rho = suppressWarnings(
        cor(joined$value_unnorm, joined$value_norm, method = "spearman")
      ),
      pearson_r = suppressWarnings(
        cor(joined$value_unnorm, joined$value_norm, method = "pearson")
      )
    )
  })
})

corr_vs_unnorm

p_spearman <- ggplot(
  corr_vs_unnorm,
  aes(x = normalize_target, y = spearman_rho)
) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_point(size = 2.4, na.rm = TRUE) +
  facet_wrap(~ sample, ncol = 1) +
  scale_x_continuous(
    limits = c(min(targets), max(targets)),
    breaks = targets
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = paste0("Spearman correlation vs unnormalized (metric: ", metric, ")"),
    x = "Normalization target (T)",
    y = "Spearman ρ"
  ) +
  theme_bw() +
  theme(
    axis.title   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 16),
    plot.title   = element_text(face = "bold", size = 16),
    strip.text.x = element_text(face = "bold", size = 16)
  )

p_spearman


p_pearson <- ggplot(
  corr_vs_unnorm,
  aes(x = normalize_target, y = pearson_r)
) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_point(size = 2.4, na.rm = TRUE) +
  facet_wrap(~ sample, ncol = 1) +
  scale_x_continuous(
    limits = c(min(targets), max(targets)),
    breaks = targets
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = paste0("Pearson correlation vs unnormalized (metric: ", metric, ")"),
    x = "Normalization target (T)",
    y = "Pearson r"
  ) +
  theme_bw() +
  theme(
    axis.title   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 16),
    plot.title   = element_text(face = "bold", size = 18),
    strip.text.x = element_text(face = "bold", size = 16)
  )

p_pearson

#ggsave(file.path(plots_qc, "Spearman-correlation_vs_unnormalized_neighborhood_index_metric.pdf"), p_spearman, width = 10, height = 10)

#ggsave(file.path(plots_qc, "Pearson-correlation_vs_unnormalized_neighborhood_index_metric.pdf"), p_pearson, width = 10, height = 10)

###########################################################################################

# Insertion-site sequence bias was quantified by extracting reference genome bases at 
#fixed offsets (+2 and +7 bp) relative to each mapped insertion 
#coordinate and comparing observed base and dinucleotide frequencies to genomic background,

# Bias(b,k) = Pr(base b at offset k∣insertion) / Pr(base b at offset k∣genome) 


seqbias_tbl <- map_dfr(samples, function(s) {
  
  f <- list.files(
    file.path(qc_dir, s),
    pattern = "\\.seqbias_2_7\\.tsv$",
    full.names = TRUE
  )
  stopifnot(length(f) == 1)
  
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(sample = s)
})

seqbias_tbl <- seqbias_tbl %>%
  mutate(
    base_p2 = substr(pair, 1, 1),   # +2 base
    base_p7 = substr(pair, 2, 2),   # +7 base
    condition = if_else(str_detect(sample, "parent"), "parent", "treated")
  )

p_dinucl <- ggplot(
  seqbias_tbl,
  aes(x = pair, y = enrichment, fill = condition)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ sample, ncol = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  labs(
    title = " ",
    x = "Dinucleotide (+2/+7)",
    y = "Enrichment"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, color = "black", face = "bold"),
    axis.title  = element_text(face = "bold", size = 14),
    strip.text  = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 13, color = "black", face = "bold")
    
  )

p_dinucl

#ggsave(file.path(plots_qc, "Insertion-site-sequence-bias.pdf"), p_dinucl, width = 10, height = 10)




