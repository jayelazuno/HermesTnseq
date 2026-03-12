install.packages("slider")   # run once

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyverse)
  library(readxl)
  library(ComplexUpset)
  #library(ComplexHeatmap)
  library(circlize)  
  library(UpSetR)
  library(ggrepel)
  library(scales)
  library(slider)
  library(minpack.lm)
})

results_dir <- here("02.results")
############ Mapping stats
mapping_stats <- read_csv(here("02.results", "04.combined", "mapping_stats.csv"), show_col_types = FALSE)
sample_order <- c(
  "yH298-parent-pool1",
  "yH298-parent-pool2",
  "yH299-parent-pool3",
  "yH299-parent-pool4",
  "yH298-H2O2-treated-facs-pool1",
  "yH298-H2O2-treated-facs-pool2",
  "yH299-H2O2-treated-facs-pool3",
  "yH299-H2O2-treated-facs-pool4"
)
#################### 
mapping_stats <- mapping_stats %>%
  mutate(sample = factor(sample, levels = rev(sample_order)))
perc_long <- mapping_stats %>%
  select(sample, percent_mapped, percent_mapq_ge20) %>%
  pivot_longer(
    -sample,
    names_to = "metric",
    values_to = "percent"
  ) %>%
  mutate(
    metric = recode(
      metric,
      percent_mapped   = "% mapped",
      percent_mapq_ge20 = "% HQ (MAPQ≥20)"
    )
  )

scale_factor <- max(mapping_stats$total_records, na.rm = TRUE) / 100

p1.mapstats <- ggplot() +
  geom_col(data = mapping_stats, aes(x = sample, y = total_records / scale_factor), alpha = 0.4) +
  geom_point(data = perc_long, aes(x = sample, y = percent, color = metric), size = 4) +
  geom_line(data = perc_long, aes(x = sample, y = percent, color = metric, group = sample),
            linewidth = 0.6, alpha = 0.7) +
  coord_flip() +
  scale_y_continuous(
    name = "Percent",
    limits = c(0, 100),
    labels = label_percent(scale = 1),
    sec.axis = sec_axis(~ . * scale_factor,
                        name = "Total reads",
                        labels = label_number(scale_cut = cut_si("")))
  ) +
  labs(title = " ", x = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title.x     = element_text(face = "bold", size = 13),
    axis.title.y     = element_text(face = "bold", size = 13),
    axis.title.y.right = element_text(face = "bold", size = 13),  # secondary axis
    axis.text        = element_text(face = "bold")                # tick labels
  )
p1.mapstats


# ============================================================
# parent library QC 
# Plot 1: Library complexity
# Plot 2: % genomic features hit (FACETED → no sample colors)
# Plot 3: Mean reads per insertion site
# Plot 4: Feature vs intergenic distribution

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
    legend.text  = element_text(size = 13, color = "black", face = "bold"),
    
    axis.ticks.length = unit(3, "mm")
  )

# Plot 1: Library complexity

p_complexity <- ggplot(
  stats_tbl,
  aes(x = sample, y = `Total Hits`)
) +
  geom_col(
    width = 0.7,
    fill = "grey96",
    color = "black",
    linewidth = 1.1
  ) +
  facet_wrap(~ condition, scales = "free_x") +
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
  aes(x = sample, y = `% of features hit`)
) +
  geom_col(
    width = 0.7,
    fill = "grey96",
    color = "black",
    linewidth = 1.1
  ) +
  facet_wrap(~ condition, scales = "free_x") +
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
  aes(x = sample, y = `Mean Reads Per Hit`)
) +
  geom_col(
    width = 0.7,
    fill = "grey96",
    color = "black",
    linewidth = 1.1
  ) +
  scale_y_log10() +
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
  geom_col(
    width = 0.7,
    color = "black",     
    linewidth = 1.1      
  ) +
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

#  ggsave(file.path(plots_qc, "library-complexity.pdf"), p_complexity, width = 10, height = 10)
# # 
# ggsave(file.path(plots_qc, "percent_features_hit.pdf"), p_features_hit, width = 10, height = 10)
# # 
# ggsave(file.path(plots_qc, "percent_reads_per_hit.pdf"), p_reads_per_hit, width = 10, height = 10)
# # 
# ggsave(file.path(plots_qc, "percent_feature_intergenic.pdf"), p_feature_intergenic, width = 10, height = 10)

############################################################
# we want to look at the genome wide reads and insertion dynamics 

# Bin reads into 1k and 10kb bin size and plot hits and reads 

# p_reads_1kb
# p_reads_10kb
# p_hits_1kb
# p_hits_10kb

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

# Chromosome relabeling 

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
# next, we want look at the general concordance among the parentlibraries 

# plot pairwise correlation among parental pools and treated pools 
# p_spearman_unnorm
# p_pearsom_unnorm

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
# Next, we want to diagnose library complexiity, jackpots, and library saturations points 
# A  read in the the midlc estimate files to visualize when unique insertion reaches saturation point
# Plot MidLC_est 

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

p_curve <- ggplot(
  midlc_mean2,
  aes(
    x = reads_sampled,
    y = unique_mean,
    group = sample,
    shape = pool          # ← encode pool by shape
  )
) +
  geom_line(
    linewidth = 0.9,
    alpha = 0.95,
    color = "black"       # ← all lines black
  ) +
  geom_point(
    size = 3,
    alpha = 1.5,
    color = "black"       # ← all points black
  ) +
  facet_wrap(~ condition, scales = "free_y", nrow = 1) +
  scale_shape_manual(
    values = c("1" = 16,  # ●
               "2" = 17,  # ▲
               "3" = 15,  # ■
               "4" = 18), # ◆
    name = "Pool"
  ) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6),
    labels = scales::label_number(
      scale_cut = scales::cut_si(""),
      accuracy = 1
    )
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
    legend.text  = element_text(size = 15)
  )


p_curve_same_y <- p_curve + facet_wrap(~ condition, scales = "fixed", nrow = 1)

p_curve_same_y

p_curve_parent  <- p_curve %+% dplyr::filter(midlc_mean2, condition == "parent") +
  facet_wrap(~ condition, scales = "free_y", nrow = 1)

p_curve_treated <- p_curve %+% dplyr::filter(midlc_mean2, condition == "H2O2-treated-facs") +
  facet_wrap(~ condition, scales = "free_y", nrow = 1)

# ggsave(file.path(plots_qc, "midlc_combined.pdf"),  p_curve_same_y,  width = 14, height = 6)
# 
# ggsave(file.path(plots_qc, "midlc_parent.pdf"),  p_curve_parent,  width = 10, height = 6)
# 
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
# was defined as Rtarget =T × MidLC, where MidLC is the estimated mid–library complexity; 
#normalization proceeded only while sufficient reads were available to meet this target.
norm_tbl_plot <- norm_tbl %>%
  arrange(sample, normalize_target)

p_norm_reads <- ggplot(
  norm_tbl_plot,
  aes(
    x = normalize_target,
    y = reads_after_norm,
    group = sample,
    shape = sample      
  )
) +
  # lines stop naturally at NA
  geom_line(
    linewidth = 0.9,
    color = "black",   
    na.rm = TRUE
  ) +
  
  # points (cex ≈ size = 3)
  geom_point(
    size = 3,           
    alpha = 0.95,
    color = "black",   
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
    shape = "Library"   # ← legend now shows shapes
  ) +
  theme_bw() +
  theme(
    axis.title   = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 14),
    plot.title   = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(size = 13),
    strip.text.x = element_text(face = "bold", size = 14)
  )

p_norm_reads 

#p_norm_reads + facet_wrap(~ sample, nrow = 1, scales = "free_y")

p_norm_reads_facet <- p_norm_reads +
  facet_wrap(~ sample, ncol = 1, scales = "free_y")

#ggsave(file.path(plots_qc, "Reads_retained_after_normalization.pdf"), p_norm_reads_facet, width = 10, height = 15)

###########################################################################################

# Correlation gene-level read counts, insertion counts (hits), and derived features including the 
# neighborhood index (NI) between normalized libraries at increasing normalization 
# targets (T) and the corresponding unnormalized data
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
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 1.1) +
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


##################################################################################################
# after running the ML model we run pairwise correlation plots among the parental samples for unnormalized, normalized T50 and T60 


classifier_dir <- here("06.classifier")

paths <- list(
  unnormalized = here(
    "06.classifier", "unormalized",
    "scer_train__parents_only__fpr0p10__5566905__20260207_174345__677",
    "classification - Scer training - filtered (table mode)",
    "tables.xlsx"
  ),
  T050 = here(
    "06.classifier", "normalized", "T050",
    "scer_train__parents_only__fpr0p10__5566926__20260207_181145__13177",
    "classification - Scer training - filtered (table mode)",
    "tables.xlsx"
  ),
  T060 = here(
    "06.classifier", "normalized", "T060",
    "scer_train__parents_only__fpr0p10__5566923__20260207_180933__26597",
    "classification - Scer training - filtered (table mode)",
    "tables.xlsx"
  )
)

stopifnot(all(file.exists(unlist(paths))))

# Helper: read first 4 sheets (pool1–pool4)

read_classifier_tables <- function(xlsx, run_label, samples = parent_samples) {
  
  sheets <- excel_sheets(xlsx)[seq_along(samples)]
  stopifnot(length(sheets) == length(samples))
  
  map_dfr(seq_along(sheets), function(i) {
    
    read_excel(xlsx, sheet = sheets[i]) %>%
      transmute(
        standard_name = `Standard name`,
        common_name   = `Common name`,
        sc_ortholog   = `Sc ortholog`,
        rf_score      = as.numeric(`RF - G4`),
        rf_ess_0p10   = `RF - G4 - ess. for FPR 0.100`,
        pool          = samples[i],   # canonical pool name
        run           = run_label
      )
  })
}

read_classifier_tables

# Read all three runs (12 sheets total)

cls_unnorm <- read_classifier_tables(paths$unnormalized, "unnormalized")
cls_T050   <- read_classifier_tables(paths$T050, "T050")
cls_T060   <- read_classifier_tables(paths$T060, "T060")

classifier_tbl <- bind_rows(cls_unnorm, cls_T050, cls_T060) %>%
  mutate(
    pool = factor(pool, levels = parent_samples),
    run  = factor(run, levels = c("unnormalized", "T050", "T060"))
  )

classifier_tbl


#Prepare wide table for correlations


rf_wide <- classifier_tbl %>%
  select(standard_name, pool, run, rf_score) %>%
  pivot_wider(
    names_from  = run,
    values_from = rf_score
  )

rf_wide <- rf_wide %>%
  mutate(
    unnormalized = as.numeric(unnormalized),
    T050         = as.numeric(T050),
    T060         = as.numeric(T060)
  )


rf_wide



pairwise_pool_corr <- function(df, value_col, run_label) {
  
  pools <- unique(df$pool)
  
  combn(pools, 2, simplify = FALSE) %>%
    purrr::map_dfr(function(p) {
      
      a <- df %>% filter(pool == p[1]) %>% select(standard_name, value = all_of(value_col))
      b <- df %>% filter(pool == p[2]) %>% select(standard_name, value = all_of(value_col))
      
      joined <- inner_join(a, b, by = "standard_name",
                           suffix = c("_1", "_2"))
      
      tibble(
        run        = run_label,
        pool_1     = p[1],
        pool_2     = p[2],
        spearman   = cor(joined$value_1, joined$value_2,
                         method = "spearman", use = "pairwise.complete.obs"),
        pearson    = cor(joined$value_1, joined$value_2,
                         method = "pearson", use = "pairwise.complete.obs")
      )
    })
}

corr_pairwise <- bind_rows(
  pairwise_pool_corr(rf_wide, "unnormalized", "Unnormalized"),
  pairwise_pool_corr(rf_wide, "T050",         "T050"),
  pairwise_pool_corr(rf_wide, "T060",         "T060")
)

corr_long <- corr_pairwise %>%
  pivot_longer(
    cols = c(spearman, pearson),
    names_to = "method",
    values_to = "correlation"
  ) %>%
  mutate(
    method = str_to_title(method),
    run    = factor(run, levels = c("Unnormalized", "T050", "T060"))
  )

p_pairwise <- ggplot(
  corr_long,
  aes(x = interaction(pool_1, pool_2), y = correlation)
) +
  geom_point(size = 3) +
  facet_grid(method ~ run, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Pairwise RF essentiality score correlations among parental pools",
    x = " ",
    y = "Correlation coefficient"
  ) +
  qc_theme

p_pairwise


ess_T060 <- classifier_tbl %>%
  filter(run == "T060") %>%
  group_by(standard_name) %>%
  summarize(
    n_pools = sum(!is.na(rf_ess_0p10)),
    ess_calls = sum(rf_ess_0p10 == "yes", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # require gene to be scored in all 4 pools
  filter(n_pools == 4) %>%
  mutate(
    essential = ess_calls >= 2  # majority rule (≥2/4)
  )


ess_summary <- ess_T060 %>%
  summarize(
    total_genes = n(),
    essential_genes = sum(essential),
    percent_essential = round(100 * essential_genes / total_genes, 1)
  )

ess_summary

##################################################################################################

# next we will find the mean of the ess scores from our four libs, determine essentiality
# plot venn diagrams, run correlations ess library scores published 
# by Gale et al 2020 using the same ML 

lit_dir <- here("02.results", "05.literature_data")

stopifnot(dir.exists(lit_dir))

# Read Gale et al. Supplemental Table S1
gale_xlsx <- file.path(
  lit_dir,
  "20200820Supplemental_table_S1_Gale.xlsx"
)

stopifnot(file.exists(gale_xlsx))
# Inspect available sheets
#excel_sheets(gale_xlsx) we want the gene sheet, which is sheet 1

gale_tbl <- read_excel(gale_xlsx, sheet = 1)

#glimpse(gale_tbl)

# Gale et al used CBS138 strain with the CGD id, we will need to map them usinf this file 
map_file <- file.path(
  lit_dir,
  "qng_gwk_cagl_threeway_map_complete.csv"
)

stopifnot(file.exists(map_file))

id_map_tbl <- read_csv(map_file, show_col_types = FALSE)

#glimpse(id_map_tbl)

# Ensure IDs are character (defensive)
gale_tbl <- gale_tbl %>%
  mutate(`Cg-ORF` = as.character(`Cg-ORF`))

id_map_tbl <- id_map_tbl %>%
  mutate(cagl_id = as.character(cagl_id))

# Join: keep all Gale rows, add GWK/QNG info where available
gale_tbl <- gale_tbl %>%
  left_join(
    id_map_tbl %>%
      select(-gene_name),
    by = c("Cg-ORF" = "cagl_id")
  )

## find the mean ess scores from our library 
rf_mean_tbl <- classifier_tbl %>%
  filter(run == "T060") %>%
  group_by(standard_name) %>%
  summarize(
    avg_ess_score = mean(rf_score, na.rm = TRUE),
    .groups = "drop"
  )


# to determine essentiallity threshold , we use the training labels (SGD-essentiality: viable vs inviable) 
# to empirically determine an RF score cutoff that separates insertion-intolerant (“inviable”) from tolerant (“viable”) 
# genes, and then apply that numeric cutoff to your C. glabrata RF scores.

train_tbl <- gale_tbl %>%
  filter(!is.na(`SGD-essentiality`)) %>%
  mutate(
    `SGD-essentiality` = factor(
      `SGD-essentiality`,
      levels = c("viable", "inviable")
    )
  )

# are going to use 3 methods to try and estamate the sc_ess_cut off 

# Median of essential genes
median_sc_ess <- median(
  train_tbl$`Sc-Ess-Score`[train_tbl$`SGD-essentiality` == "inviable"],
  na.rm = TRUE
)

# Trimmed mean (robust average)

mean_sc_ess <- mean(
  train_tbl$`Sc-Ess-Score`[train_tbl$`SGD-essentiality` == "inviable"],
  trim = 0.10,
  na.rm = TRUE
)

# IQR-based filtering then mean
x <- train_tbl$`Sc-Ess-Score`[train_tbl$`SGD-essentiality` == "inviable"]

q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
iqr <- q[2] - q[1]

x_filt <- x[x >= (q[1] - 1.5*iqr) & x <= (q[2] + 1.5*iqr)]

IQR_sc_ess <- mean(x_filt, na.rm = TRUE)

# All three values fall in a very tight range ~0.90–0.93

sc_ess_cutoff <- 0.925

## apply score to mean table 
rf_mean_tbl <- rf_mean_tbl %>%
  mutate(
    essentia_0.925_BG2 = avg_ess_score >= sc_ess_cutoff
  )

# compute % essential genes 
per_ess_cg <- rf_mean_tbl %>%
  summarize(
    n_genes = n(),
    n_essential = sum(essentia_0.925_BG2, na.rm = TRUE),
    percent_essential = 100 * mean(essentia_0.925_BG2 , na.rm = TRUE)
  )

### apply the same cut off to gale et table 
gale_tbl <- gale_tbl %>%
  mutate(
    essential_0.925_CBS138 = `Cg_Ess-Score` >= sc_ess_cutoff
  )

###### applying this cuff off gave 16% of essential genes which means Gale et used 
# and even permissive cut off 
per_ess_cg_gale <- gale_tbl %>%
  summarize(
    n_genes = sum(!is.na(`Cg_Ess-Score`)),
    n_essential = sum(essential_0.925_CBS138, na.rm = TRUE),
    percent_essential = 100 * mean(essential_0.925_CBS138, na.rm = TRUE)
  )

# so we compute the cutoff that gives 25%

gale_cutoff_25pct <- quantile(
  gale_tbl$`Cg_Ess-Score`,
  probs = 0.75,
  na.rm = TRUE
)

gale_cutoff_25pct

## apply to the gale df 
gale_tbl <- gale_tbl %>%
  mutate(
    essential_25pct_CBS138 = `Cg_Ess-Score` >= gale_cutoff_25pct
  )

# apply to our existing df 
rf_mean_tbl <- rf_mean_tbl %>%
  mutate(
    essential_25pct_BG2 = avg_ess_score >= gale_cutoff_25pct
  )


### map the two data frames 
rf_vs_gale_tbl <- rf_mean_tbl %>%
  inner_join(
    gale_tbl,
    by = c("standard_name" = "gwk60_id")
  ) %>%
  filter(!is.na(`Cg_Ess-Score`))


################################## subset df for plotting 
upset_0925 <- rf_vs_gale_tbl %>%
  distinct(
    standard_name,
    essentia_0.925_BG2,
    essential_0.925_CBS138
  ) %>%
  transmute(
    gene_id = standard_name,
    ess_cutoff_0.925_BG2 = as.integer(essentia_0.925_BG2),
    ess_cutoff_0.925_CBS138 = as.integer(essential_0.925_CBS138)
  ) %>%
  as.data.frame()

upset_0535 <- rf_vs_gale_tbl %>%
  distinct(
    standard_name,
    essential_25pct_BG2,
    essential_25pct_CBS138
  ) %>%
  transmute(
    gene_id = standard_name,
    ess_galecutoff_0.535_BG2 = as.integer(essential_25pct_BG2),
    ess_galecutoff_0.535_CBS138 = as.integer(essential_25pct_CBS138)
  ) %>%
  as.data.frame()


op <- par(font = 2, cex = 4)   # turn EVERYTHING bold

upsetPlot_0.925 <- UpSetR::upset(
  upset_0925,
  sets = c("ess_cutoff_0.925_BG2", "ess_cutoff_0.925_CBS138"),
  order.by = "freq",
  mainbar.y.label = "Intersection size",
  sets.x.label = " ",
  show.numbers = "yes",
  mb.ratio = c(0.65, 0.35),
  text.scale = c(
    3.0,  # intersection size title
    3.0,  # set size title
    3.0,  # intersection numbers
    3.0,  # set size numbers
    3.0,  # set names
    3.0  # axis numbers
  )
)
upsetPlot_0.925
#par(op)  # restore previous graphics state

upsetPlot_0.535 <- UpSetR::upset(
  upset_0535,
  sets = c("ess_galecutoff_0.535_BG2", "ess_galecutoff_0.535_CBS138"),
  order.by = "freq",
  mainbar.y.label = "Intersection size",
  sets.x.label = " ",
  show.numbers = "yes",
  mb.ratio = c(0.65, 0.35),
  text.scale = c(
    3.0,  
    3.0,  
    3.0,  
    3.0,  
    3.0, 
    3.0  
  )
)
upsetPlot_0.535 



# RF vs C. glabrata essentiality
cor_rf_vs_gale_Cg <- tibble(
  spearman = cor(
    rf_vs_gale_tbl$avg_ess_score,
    rf_vs_gale_tbl$`Cg_Ess-Score`,
    method = "spearman",
    use = "pairwise.complete.obs"
  ),
  pearson = cor(
    rf_vs_gale_tbl$avg_ess_score,
    rf_vs_gale_tbl$`Cg_Ess-Score`,
    method = "pearson",
    use = "pairwise.complete.obs"
  ),
  n_genes = sum(!is.na(rf_vs_gale_tbl$`Cg_Ess-Score`))
)

cor_rf_vs_gale_Cg

# RF vs S. cerevisiae essentiality
rf_vs_gale_Sc_tbl <- rf_vs_gale_tbl %>%
  filter(!is.na(`Sc-Ess-Score`))


cor_rf_vs_gale_Sc <- tibble(
  spearman = cor(
    rf_vs_gale_Sc_tbl$avg_ess_score,
    rf_vs_gale_Sc_tbl$`Sc-Ess-Score`,
    method = "spearman",
    use = "pairwise.complete.obs"
  ),
  pearson = cor(
    rf_vs_gale_Sc_tbl$avg_ess_score,
    rf_vs_gale_Sc_tbl$`Sc-Ess-Score`,
    method = "pearson",
    use = "pairwise.complete.obs"
  ),
  n_genes = nrow(rf_vs_gale_Sc_tbl)
)

cor_rf_vs_gale_Sc

# RF vs Sc essentiality (ortholog-limited)
p_rf_vs_sc <- ggplot(
  rf_vs_gale_Sc_tbl,
  aes(x = `Sc-Ess-Score`, y = avg_ess_score)
) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title = "RF essentiality vs S. cerevisiae essentiality",
    subtitle = "Only genes with Sc orthologs",
    x = "S. cerevisiae essentiality score",
    y = "RF essentiality score (BG2, mean of pools)"
  ) +
  qc_theme

p_rf_vs_sc



#############################################################

# we use CPM normalizes for sequencing depth while preserving selection-driven 
# differences, making it the correct normalization for enrichment analysis in FACS-selected Tn-seq libraries.

unnorm_dir <- here("05.summary")


parent_samples <- c(
  "yH298-parent-pool1",
  "yH298-parent-pool2",
  "yH299-parent-pool3",
  "yH299-parent-pool4"
)

treated_samples <- c(
  "yH298-H2O2-treated-facs-pool1",
  "yH298-H2O2-treated-facs-pool2",
  "yH299-H2O2-treated-facs-pool3",
  "yH299-H2O2-treated-facs-pool4"
)

all_samples <- c(parent_samples, treated_samples)

read_feature_table <- function(f, sample) {
  read_csv(f, skip = 1, show_col_types = FALSE) %>%
    mutate(sample = sample)
}

unnorm_file <- function(sample) {
  file.path(
    unnorm_dir,
    sample,
    paste0(sample, ".feature_table.RDF_1.csv")
  )
}

feature_tbl <- purrr::map_dfr(all_samples, function(s) {
  
  f <- unnorm_file(s)
  stopifnot(file.exists(f))
  
  read_feature_table(f, sample = s)
})
feature_tbl <- purrr::map_dfr(all_samples, function(s) {
  
  f <- unnorm_file(s)
  stopifnot(file.exists(f))
  
  read_feature_table(f, sample = s)
})
#stopifnot("reads" %in% colnames(feature_tbl))
feature_tbl <- feature_tbl %>%
  mutate(
    condition = if_else(
      stringr::str_detect(sample, "parent"),
      "parent",
      "treated"
    )
  )

# Add pool + background
feature_tbl2 <- feature_tbl %>%
  mutate(
    pool = stringr::str_match(sample, "pool(\\d+)")[,2],
    pool = factor(pool, levels = c("1","2","3","4")),
    background = stringr::str_extract(sample, "H\\d+"),
    background = factor(background, levels = c("H298","H299"))
  )

stopifnot("reads" %in% colnames(feature_tbl2))


#Function to fit  curve for one background

#uses control-control log2 ratios (M_ctrl)

#computes 40-gene running SD vs running mean abundance

#fits sd = m1 + m2 * A^m3 using nls

#by default excludes very low abundance (min_A = 6 like Gale; adjust if CPM scale)
# Build a gene×sample read table (wide)

reads_wide <- feature_tbl2 %>%
  select(standard_name, sample, reads) %>%
  group_by(standard_name, sample) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%   # safety if duplicated rows
  pivot_wider(names_from = sample, values_from = reads, values_fill = 0)

normalize_pair <- function(x1, x2) {
  # "slight normalization" (size-factor scaling) for a replicate pair
  # works on raw reads
  g  <- sqrt((x1 + 1) * (x2 + 1))     # pseudo-geometric mean per gene
  r1 <- (x1 + 1) / g
  r2 <- (x2 + 1) / g
  s1 <- median(r1[is.finite(r1)], na.rm = TRUE)
  s2 <- median(r2[is.finite(r2)], na.rm = TRUE)
  
  list(
    x1n = x1 / s1,
    x2n = x2 / s2,
    s1 = s1, s2 = s2
  )
}

make_control_running <- function(reads_wide, ctrl1, ctrl2, w = 40, minA = 6) {
  
  norm <- normalize_pair(reads_wide[[ctrl1]], reads_wide[[ctrl2]])
  c1 <- norm$x1n
  c2 <- norm$x2n
  
  ctrl_tab <- tibble(
    standard_name = reads_wide$standard_name,
    A_ctrl = (c1 + c2) / 2,
    M_ctrl = log2((c1 + 1) / (c2 + 1))
  ) %>%
    arrange(desc(A_ctrl)) %>%
    mutate(
      run_A  = slide_dbl(A_ctrl, mean, .before = floor(w/2), .after = floor(w/2), .complete = TRUE),
      run_SD = slide_dbl(M_ctrl, sd,   .before = floor(w/2), .after = floor(w/2), .complete = TRUE)
    )
  
  fit_dat <- ctrl_tab %>%
    filter(!is.na(run_A), !is.na(run_SD), run_A >= minA, run_SD > 0)
  
  # ---- Data-driven starting values (key) ----
  # Fit log(run_SD) ~ log(run_A) to get a reasonable exponent start
  lm0 <- lm(log(run_SD) ~ log(run_A), data = fit_dat)
  m3_start <- as.numeric(coef(lm0)[2])          # slope
  # Clamp to bounds we allow
  m3_start <- max(min(m3_start, 0.3), -2.5)
  
  # Use low-A end to guess m1 (asymptote-ish) and overall range for m2
  m1_start <- as.numeric(quantile(fit_dat$run_SD, 0.02, na.rm = TRUE))
  # pick a representative x to scale m2
  x_ref <- as.numeric(median(fit_dat$run_A, na.rm = TRUE))
  y_ref <- as.numeric(median(fit_dat$run_SD, na.rm = TRUE))
  m2_start <- max((y_ref - m1_start) / (x_ref ^ m3_start), 1e-6)
  
  start <- list(m1 = m1_start, m2 = m2_start, m3 = m3_start)
  
  fit <- minpack.lm::nlsLM(
    run_SD ~ m1 + m2 * (run_A ^ m3),
    data = fit_dat,
    start = start,
    lower = c(m1 = 0,    m2 = 0,    m3 = -5),
    upper = c(m1 = Inf,  m2 = Inf,  m3 =  0.5),
    control = minpack.lm::nls.lm.control(maxiter = 500)
  )
  
  list(ctrl_tab = ctrl_tab, fit = fit, norm = norm, fit_dat = fit_dat, start = start)
}
#Fit H298 and H299 exactly

fit_H298 <- make_control_running(reads_wide, "yH298-parent-pool1", "yH298-parent-pool2")
fit_H299 <- make_control_running(reads_wide, "yH299-parent-pool3", "yH299-parent-pool4")

coef(fit_H298$fit)
coef(fit_H299$fit)
fit_H298$start
fit_H299$start

sd_from_fit <- function(fit_obj, A) {
  cf <- coef(fit_obj)
  m1 <- cf[["m1"]]; m2 <- cf[["m2"]]; m3 <- cf[["m3"]]
  pmax(m1 + m2 * (A ^ m3), 1e-8)
}

compute_pool_z <- function(reads_wide, treated, parent, fit_obj) {
  
  # slight normalization of exp vs control replicates
  norm <- normalize_pair(reads_wide[[treated]], reads_wide[[parent]])
  t <- norm$x1n
  p <- norm$x2n
  
  A_exp <- (t + p) / 2
  M_exp <- log2((t + 1) / (p + 1))              # <-- log2ratio experiment/control
  sd_local <- sd_from_fit(fit_obj, A_exp)       # <-- local SD from fitted power curve
  z <- M_exp / sd_local                         # <-- z-score 
  
  tibble(
    standard_name = reads_wide$standard_name,
    treated_sample = treated,
    parent_sample  = parent,
    A_exp = A_exp,
    M_exp = M_exp,
    sd_local = sd_local,
    z = z
  )
}

z_pool1 <- compute_pool_z(reads_wide,
                          "yH298-H2O2-treated-facs-pool1", "yH298-parent-pool1",
                          fit_H298$fit)

z_pool2 <- compute_pool_z(reads_wide,
                          "yH298-H2O2-treated-facs-pool2", "yH298-parent-pool2",
                          fit_H298$fit)

z_pool3 <- compute_pool_z(reads_wide,
                          "yH299-H2O2-treated-facs-pool3", "yH299-parent-pool3",
                          fit_H299$fit)

z_pool4 <- compute_pool_z(reads_wide,
                          "yH299-H2O2-treated-facs-pool4", "yH299-parent-pool4",
                          fit_H299$fit)

z_all <- bind_rows(z_pool1, z_pool2, z_pool3, z_pool4) %>%
  mutate(pool = stringr::str_match(treated_sample, "pool(\\d+)")[,2])


# Make a histogram of z_ctrl for each background (control-control):For H298 and H299
# check tials to see which replicate is more noisier 
m298 <- fit_H298$ctrl_tab$M_ctrl
m299 <- fit_H299$ctrl_tab$M_ctrl

bind_rows(
  tibble(background="H298", M=m298),
  tibble(background="H299", M=m299)
) %>%
  ggplot(aes(M)) +
  geom_histogram(bins=120) +
  facet_wrap(~background, scales="free_y") +
  theme_bw() +
  labs(x="control-control log2 ratio (M_ctrl)", y="count")

z_ctrl_H298 <- fit_H298$ctrl_tab$M_ctrl / sd_from_fit(fit_H298$fit, fit_H298$ctrl_tab$A_ctrl)
z_ctrl_H299 <- fit_H299$ctrl_tab$M_ctrl / sd_from_fit(fit_H299$fit, fit_H299$ctrl_tab$A_ctrl)

bind_rows(
  tibble(background="H298", z=z_ctrl_H298),
  tibble(background="H299", z=z_ctrl_H299)
) %>%
  ggplot(aes(x = z)) +
  geom_histogram(bins = 100) +
  facet_wrap(~ background, scales = "free_y") +
  theme_bw() +
  labs(x = "control-control z", y = "count")


Zthr_995_H298 <- quantile(abs(z_ctrl_H298), 0.995, na.rm=TRUE)
Zthr_999_H298 <- quantile(abs(z_ctrl_H298), 0.999, na.rm=TRUE)
# #> Zthr_995_H298 
# 99.5% 
# 4.03227 
# > Zthr_999_H298
# 99.9% 
# 4.67864 

Zthr_995_H299 <- quantile(abs(z_ctrl_H299), 0.995, na.rm=TRUE)
Zthr_999_H299 <- quantile(abs(z_ctrl_H299), 0.999, na.rm=TRUE)

# > Zthr_995_H299
# 99.5% 
# 2.1255 
# > Zthr_999_H299
# 99.9% 
# 2.312504 


# choose which stringency you want
Zthr_H298 <- Zthr_995_H298  # 4.03227
Zthr_H299 <- Zthr_995_H299  # 2.1255

z_ma <- z_all %>%
  mutate(
    background = if_else(pool %in% c("1","2"), "H298", "H299"),
    Zthr = if_else(background == "H298", Zthr_H298, Zthr_H299),
    call = case_when(
      z >=  Zthr ~ "enriched",
      z <= -Zthr ~ "depleted",
      TRUE ~ "none"
    )
  )

ggplot(z_ma, aes(x = A_exp + 1, y = M_exp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(color = call), alpha = 0.35, size = 1.2) +
  scale_x_log10() +
  facet_grid(background ~ pool) +
  theme_bw() +
  labs(
    x = "A = average reads (treated+parent)/2 (log10)",
    y = "M = log2(treated/parent)",
    color = "Call (z cutoff)"
  )


OSR_curated <- read_csv(here("02.results", "04.combined", "OSR_Scer_mapped.csv"))

OSR_curated_Mapped <- enriched_genes %>%
  inner_join(OSR_curated, by = c("Sc-ORF" = "gene_id")) %>%
  select(
    standard_name,
    `Cg-ORF`,
    `Cg-length`,
    avg_ess_score,
    mean_treated_cpm,
    `Sc-ORF`,
    `Sc-name`,
    `SGD-essentiality`,
    `SGD-description`,
    `Cg-to-Sc-relationship`,
    `pre-WGD-Ancestor`,
    n_pools, 
    pools_present
  ) %>%
  arrange(desc(mean_treated_cpm))


# ---------------------------
# 0) Use z0.99 cutoffs (per background)
# ---------------------------
Zthr_H298 <- quantile(abs(z_ctrl_H298), 0.99, na.rm = TRUE)
Zthr_H299 <- quantile(abs(z_ctrl_H299), 0.99, na.rm = TRUE)

# ---------------------------
# 1) z_all (per gene x pool) + background-specific calls
# ---------------------------
z_ma <- z_all %>%
  mutate(
    pool = as.character(pool),
    background = if_else(pool %in% c("1","2"), "H298", "H299"),
    Zthr = if_else(background == "H298", Zthr_H298, Zthr_H299),
    call = case_when(
      z >=  Zthr ~ "enriched",
      z <= -Zthr ~ "depleted",
      TRUE ~ "none"
    )
  )

# ---------------------------
# 2) Join z_ma to feature_tbl (raw read feature table)
#    Keep relevant annotation columns from feature_tbl
#    (edit selects if you want more/less)
# ---------------------------
# Keep ONLY the columns you actually want downstream (your "relevant cols")
feature_annot <- rf_vs_gale_tbl %>%
  select(
    standard_name,
    `Cg-ORF`,
    `Cg-length`,
    avg_ess_score,
    `Sc-ORF`,
    `Sc-name`,
    `SGD-essentiality`,
    `SGD-description`,
    `Cg_Ess-Score`,
    `Sc-Ess-Score`,
    `Cg-Sc-∆Ess-Score`,
    `Cg-to-Sc-relationship`,
    `pre-WGD-Ancestor`,
    `qng_id`
  ) %>%
  distinct(standard_name, .keep_all = TRUE)

# Join z table to the full annotation table
z_ma_annot <- z_ma %>%
  left_join(feature_annot, by = "standard_name")

# 3) add a global label column from the joined annotation table
# keep your factor levels
z_ma_annot2 <- z_ma_annot %>%
  mutate(
    label_name = coalesce(na_if(str_trim(`Sc-name`), ""), standard_name),
    call = factor(call, levels = c("depleted","none","enriched"))
  )

# then in your plot add:
scale_color_manual(
  values = c(
    depleted = "blue",
    none     = "grey70",
    enriched = "red"
  )
)
# 4) label set = ALL significant genes (global)
sig_tbl <- z_ma_annot2 %>%
  filter(call %in% c("enriched", "depleted"))

# (optional) limit labels to top hits per facet to reduce clutter
sig_tbl <- sig_tbl %>%
   group_by(background, pool, call) %>%
   slice_max(order_by = abs(z), n = 10) %>%
   ungroup()

base_ma <- function(df) {
  ggplot(df, aes(x = A_exp + 1, y = M_exp)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(aes(color = call), alpha = 0.35, size = 1.2) +
    scale_color_manual(values = c(depleted="blue", none="grey70", enriched="red")) +
    scale_x_log10() +
    guides(color = guide_legend(title = "Call (z0.99 cutoff)")) +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                  size = 10, color = "black", face = "bold"),
      axis.title   = element_text(face = "bold", size = 14),
      strip.text   = element_text(face = "bold", size = 10),
      axis.text.y  = element_text(size = 12, color = "black", face = "bold"),
      legend.title = element_text(size = 13, color = "black", face = "bold"),
      legend.text  = element_text(size = 13, color = "black", face = "bold")
    ) +
    labs(
      x = "A = average reads (treated+parent)/2 (log10)",
      y = "M = log2(treated/parent)"
    )
}

# 5) two separate plots (H298 and H299), label ALL significant genes
ma_plot_298 <- base_ma(z_ma_annot2 %>% filter(background == "H298")) +
  facet_wrap(~ pool, nrow = 1) +
  ggrepel::geom_text_repel(
    data = sig_tbl %>% filter(background == "H298"),
    aes(label = label_name),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.2,
    min.segment.length = 0,
    seed = 1
  ) +
  ggtitle("H298 MA plots (pools 1–2)")

ma_plot_299 <- base_ma(z_ma_annot2 %>% filter(background == "H299")) +
  facet_wrap(~ pool, nrow = 1) +
  ggrepel::geom_text_repel(
    data = sig_tbl %>% filter(background == "H299"),
    aes(label = label_name),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.2,
    min.segment.length = 0,
    seed = 1
  ) +
  ggtitle("H299 MA plots (pools 3–4)")

ma_plot_298
ma_plot_299

######################
OSR_curated <- read_csv(here("02.results", "04.combined", "OSR_Scer_mapped.csv"))

# Join OSR list onto your z/annotation table (key: Sc-ORF in z_ma_annot2 == gene_id in OSR_curated)
OSR_curated_Mapped <- z_ma_annot2 %>%
  inner_join(OSR_curated, by = c("Sc-ORF" = "gene_id")) %>%
  dplyr::select(-any_of(c("gene_id", "scer_gene_name", "functional_category",
                          "category_order", "reference", "comment", "mapped_by")))

# 1) OSR label table: use scer_gene_name from OSR_curated, fallback to standard_name
osr_labels <- OSR_curated_Mapped %>%
  transmute(
    standard_name,
    plot_name = coalesce(na_if(str_trim(`Sc-name`), ""), standard_name)
  ) %>%
  distinct()

# 2) Join OSR labels into z_ma_annot2 (this table already has call/background/pool/etc)
#    Use suffixes so plot_name can't silently rename if it already exists.
z_ma_annot_osr <- z_ma_annot2 %>%
  left_join(osr_labels, by = "standard_name", suffix = c("", "_osr")) %>%
  mutate(
    is_OSR = !is.na(plot_name),
    # if you prefer an explicit label column:
    osr_label = plot_name
  )

# 3) Split backgrounds FROM THE JOINED TABLE
z_ma_298_osr <- z_ma_annot_osr %>% filter(background == "H298")
z_ma_299_osr <- z_ma_annot_osr %>% filter(background == "H299")

# 4) Reuse base_ma() and add OSR labels
ma_plot_298_osr <- base_ma(z_ma_298_osr) +
  facet_wrap(~ pool, nrow = 1) +
  ggrepel::geom_text_repel(
    data = z_ma_298_osr %>% filter(is_OSR),
    aes(label = osr_label),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.2,
    min.segment.length = 0,
    seed = 1
  ) +
  ggtitle("H298 MA plots (pools 1–2)")

ma_plot_299_osr <- base_ma(z_ma_299_osr) +
  facet_wrap(~ pool, nrow = 1) +
  ggrepel::geom_text_repel(
    data = z_ma_299_osr %>% filter(is_OSR),
    aes(label = osr_label),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.2,
    min.segment.length = 0,
    seed = 1
  ) +
  ggtitle("H299 MA plots (pools 3–4)")

ma_plot_298_osr
ma_plot_299_osr



# ### QC check total reads in all libraries
# feature_tbl %>%
#   group_by(sample) %>%
#   summarize(
#     total_reads = sum(reads, na.rm = TRUE),
#     n_genes     = n()
#   )

# # CPM normalization (per sample, using reads)
# feature_cpm <- feature_tbl %>%
#   group_by(sample) %>%
#   mutate(
#     total_reads = sum(reads, na.rm = TRUE),
#     cpm = (reads / total_reads) * 1e6
#   ) %>%
#   ungroup()
#
#
# ######################
# # first we treat each sample separate to diagnose library concordance

# # Add condition + pool id
# feature_cpm2 <- feature_cpm %>%
#   mutate(
#     condition = if_else(str_detect(sample, "parent"), "parent", "treated"),
#     pool = str_match(sample, "pool(\\d+)")[, 2],
#     pool = factor(pool, levels = c("1","2","3","4"))
#   )
# # Add a “background” label
# feature_cpm2 <- feature_cpm2 %>%
#   mutate(
#     background = str_extract(sample, "H\\d+") #%>% replace_na("bg_unknown")
#   )
# 
# # Build per-gene parent replicate pairs within each background 
# parent_wide_bg <- parent_cpm %>%
#   left_join(feature_cpm2 %>% distinct(sample, background), by = c("parent_sample" = "sample")) %>%
#   select(standard_name, background, pool, parent_cpm) %>%
#   mutate(pool = as.character(pool)) %>%
#   # wide: one row per gene per background with two parent pools
#   pivot_wider(names_from = pool, values_from = parent_cpm, names_prefix = "pool") 



# fit_gale_sd_curve <- function(df_bg, poolA, poolB, w = 40, min_A = 6) {
#   # df_bg: rows are genes; must have columns pool{A}, pool{B}
#   a_col <- paste0("pool", poolA)
#   b_col <- paste0("pool", poolB)
#   
#   dat <- df_bg %>%
#     transmute(
#       standard_name,
#       A_ctrl = ((.data[[a_col]] %||% NA_real_) + (.data[[b_col]] %||% NA_real_)) / 2,
#       M_ctrl = log2((.data[[a_col]] + 1) / (.data[[b_col]] + 1))
#     ) %>%
#     filter(!is.na(A_ctrl), !is.na(M_ctrl)) %>%
#     arrange(desc(A_ctrl)) %>%
#     mutate(
#       run_sd = slide_dbl(M_ctrl, sd,
#                          .before = floor(w/2), .after = floor(w/2),
#                          .complete = TRUE),
#       run_A  = slide_dbl(A_ctrl, mean,
#                          .before = floor(w/2), .after = floor(w/2),
#                          .complete = TRUE)
#     ) %>%
#     filter(!is.na(run_sd), !is.na(run_A), run_A >= min_A, run_sd > 0)
#   
#   if (nrow(dat) < 50) stop("Not enough genes to fit SD curve (check thresholds / CPM scale).")
#   
#   nls(
#     run_sd ~ m1 + m2 * (run_A ^ m3),
#     data = dat,
#     start = list(m1 = 0.05, m2 = 1, m3 = -0.5),
#     control = nls.control(maxiter = 500, warnOnly = TRUE)
#   )
# }
# 
# #Fit SD curves per  genomic background
# bg_levels <- sort(unique(parent_wide_bg$background))
# 
# # Fit curves (edit pool pairs if your design differs)
# fit_list <- list()
# for (bg in bg_levels) {
#   df_bg <- parent_wide_bg %>% filter(background == bg)
#   
#   if (all(c("pool1","pool2") %in% names(df_bg)) && bg == "H298") {
#     fit_list[[bg]] <- fit_gale_sd_curve(df_bg, poolA = "1", poolB = "2")
#   } else if (all(c("pool3","pool4") %in% names(df_bg)) && bg == "H299") {
#     fit_list[[bg]] <- fit_gale_sd_curve(df_bg, poolA = "3", poolB = "4")
#   } else {
#     message("Skipping background ", bg, ": pool pair not defined/found.")
#   }
# }
# 
# # Split parent vs treated CPM tables
# parent_cpm <- feature_cpm2 %>%
#   filter(condition == "parent") %>%
#   select(standard_name, pool, parent_sample = sample, parent_cpm = cpm)
# 
# 
# treated_cpm <- feature_cpm2 %>%
#   filter(condition == "treated") %>%
#   select(standard_name, pool, treated_sample = sample, treated_cpm = cpm)

##############
# lets just create upset plots of parents and treated to diagnose replicated 
# concordance

# # Filter parent pools by CPM threshold
# parent_cpm_filtered <- parent_cpm %>%
#   filter(parent_cpm >= 1) %>%
#   select(standard_name, pool) %>%
#   mutate(present = 1) %>%
#   pivot_wider(names_from = pool, values_from = present, values_fill = 0) %>%
#   rename(pool1 = `1`, pool2 = `2`, pool3 = `3`, pool4 = `4`) %>%
#   as.data.frame()
# 
# # Create UpSet plot for parent pools
# par(font = 2, font.lab = 2, font.axis = 2, font.main = 2)
# 
# parent_cpm_upset <- upset(
#   parent_cpm_filtered,
#   sets = c("pool4", "pool3", "pool2", "pool1"),
#   order.by = "freq",
#   keep.order = TRUE,
#   mainbar.y.label = "Intersection size",
#   sets.x.label = "Parent pools (CPM ≥ 1)",
#   show.numbers = "yes",
#   text.scale = c(1.6, 1.2, 1.2, 0, 1.6, 1.6),
#   set_size.show = FALSE,
#   mb.ratio = c(0.6, 0.4),
#   sets.bar.color = "black",
#   main.bar.color = "black"
# )
# parent_cpm_upset
# 
# # Create UpSet plot for treated pools
# treated_cpm_filtered <- treated_cpm %>%
#   filter(treated_cpm >= 1) %>%
#   select(standard_name, pool) %>%
#   mutate(present = 1) %>%
#   pivot_wider(names_from = pool, values_from = present, values_fill = 0) %>%
#   rename(treated_pool1 = `1`, treated_pool2 = `2`, treated_pool3 = `3`, treated_pool4 = `4`) %>%
#   as.data.frame()
# 
# # # Create UpSet plot for parent pools
# treated_cpm_upset <- upset(
#   treated_cpm_filtered,
#   sets = c("treated_pool4", "treated_pool3", "treated_pool2", "treated_pool1"),
#   order.by = "freq",
#   keep.order = TRUE,
#   mainbar.y.label = "Intersection size",
#   sets.x.label = "H2O2-treated pools (CPM ≥ 1)",
#   show.numbers = "yes",
#   text.scale = c(1.6, 1.2, 1.2, 1.2, 1.6, 1.6),
#   set_size.show = FALSE,
#   mb.ratio = c(0.6, 0.4),
#   sets.bar.color = "black",
#   main.bar.color = "black"
# )
# treated_cpm_upset 
# 
# ####################################
# # Now we are going to filter for enriched genes from the treated FACS fool 
# # Primary hit list: ≥3/4 pools
# facs_enriched_primary <- treated_cpm %>%
#   filter(treated_cpm >= 1) %>%
#   group_by(standard_name) %>%
#   filter(n() >= 2) %>%  # At least 2 of 4 pools
#   summarise(
#     mean_treated_cpm = mean(treated_cpm),
#     n_pools = n(),
#     pools_present = paste(pool, collapse = ",")
#   ) %>%
#   arrange(desc(mean_treated_cpm))
# 
# # High-confidence subset: 4/4 pools
# facs_enriched_all4 <- facs_enriched_primary %>%
#   filter(n_pools == 4)
# 
# ######### Build annotated tables 
# 
# # >= 3 pools
# enriched_genes <- rf_vs_gale_tbl %>%
#   inner_join(facs_enriched_primary, by = "standard_name") %>%
#   select(
#     standard_name,
#     `Cg-ORF`,
#     `Cg-length`,
#     avg_ess_score,
#     mean_treated_cpm,
#     `Sc-ORF`,
#     `Sc-name`,
#     `SGD-essentiality`,
#     `SGD-description`,
#     `Cg-to-Sc-relationship`,
#     `pre-WGD-Ancestor`,
#     n_pools, 
#     pools_present
#   ) %>%
#   arrange(desc(mean_treated_cpm))
# 
# # subset curated OSR genes from Scer. load the Scer OSR curated genes from the 
# # 02.results 04.combined. OSR_Scer_mapped.csv
# 
# 
# # 4 pools 
# enriched_genes_4 <- rf_vs_gale_tbl %>%
#   inner_join(facs_enriched_all4, by = "standard_name") %>%
#   select(
#     standard_name,
#     `Cg-ORF`,
#     `Cg-length`,
#     avg_ess_score,
#     mean_treated_cpm,
#     `Sc-ORF`,
#     `Sc-name`,
#     `SGD-essentiality`,
#     `SGD-description`,
#     `Cg-to-Sc-relationship`,
#     `pre-WGD-Ancestor`,
#     n_pools, 
#     pools_present
#   ) %>%
#   arrange(desc(mean_treated_cpm))
# 
# OSR_curated_Mapped <- OSR_curated_Mapped %>%
#   mutate(
#     plot_name = coalesce(na_if(str_trim(`Sc-name`), ""), standard_name)
#   )
# 
# 
# OSR_curated_Mapped2 <- OSR_curated_Mapped %>%
#   mutate(
#     n_pools_int = as.integer(round(n_pools)),
#     SGD_ess_plot = coalesce(as.character(`SGD-essentiality`), "NA")
#   ) %>%
#   filter(n_pools_int %in% c(2,3,4)) %>%
#   mutate(n_pools_int = factor(n_pools_int, levels = c(2,3,4)))
# 
# osr_plot <- ggplot(OSR_curated_Mapped2,
#                    aes(x = avg_ess_score,
#                        y = mean_treated_cpm,
#                        color = SGD_ess_plot)) +
#   geom_point(aes(size = n_pools_int), alpha = 0.85) +
#   ggrepel::geom_text_repel(
#     aes(label = plot_name),     # ONLY label, no size mapping
#     size = 5,                   # constant font size for all labels
#     box.padding = 0.25,
#     point.padding = 0.2,
#     max.overlaps = Inf,
#     min.segment.length = 0,
#     seed = 1
#   ) +
#   scale_y_log10() +
#   scale_size_manual(values = c(`2` = 2.5, `3` = 3.5, `4` = 4.5)) +
#   labs(
#     x = "avg_ess_score",
#     y = "mean treated CPM (log10)",
#     size = "n_pools",
#     color = "SGD-essentiality"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title   = element_text(face = "bold", size = 14),
#     axis.text    = element_text(face = "bold", size = 14),
#     strip.text   = element_text(face = "bold", size = 14),
#     plot.title   = element_text(face = "bold", size = 16),
#     legend.title = element_text(face = "bold", size = 14),
#     legend.text  = element_text(size = 15)
#   )
# 
# osr_plot
# 
# # write_csv(
# #   enriched_genes,
# #   file = file.path(results_dir, "H2O2-treated-facs-enriched>=3pools100cpm.csv")
# # )
# # 
# # write_csv(
# #   enriched_genes_4,
# #   file = file.path(results_dir, "H2O2-treated-facs-enriched-4pools100cpm.csv")
# # )
# 
# ########################### Now let's make some exploratory plots 
# 
# # 1) Histogram (distribution shape) + density overlay
# p_hist <- ggplot(enriched_genes, aes(x = as.numeric(avg_ess_score))) +
#   geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.6) +
#   geom_density(linewidth = 1) +
#   labs(
#     title = "Essentiality score distribution (enriched genes)",
#     x = "avg_ess_score",
#     y = "Density"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title  = element_text(face = "bold", size = 14),
#     axis.text   = element_text(face = "bold"),
#     plot.title  = element_text(face = "bold", size = 16)
#   )
# 
# # 2) ECDF (good for comparing quantiles / thresholds)
# p_ecdf <- ggplot(enriched_genes, aes(x = as.numeric(avg_ess_score))) +
#   stat_ecdf(linewidth = 1) +
#   labs(
#     title = "ECDF of essentiality scores (enriched genes)",
#     x = "avg_ess_score",
#     y = "Cumulative fraction"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title  = element_text(face = "bold", size = 14),
#     axis.text   = element_text(face = "bold"),
#     plot.title  = element_text(face = "bold", size = 16)
#   )
# 
# # n values (drop NA scores just in case)
# n_enriched_3of4 <- enriched_genes %>%
#   filter(!is.na(avg_ess_score)) %>%
#   nrow()
# 
# # Violin + box with n in subtitle
# p_violin <- ggplot(enriched_genes, aes(x = "", y = as.numeric(avg_ess_score))) +
#   geom_violin(trim = FALSE, alpha = 0.6) +
#   geom_boxplot(width = 0.15, outlier.alpha = 0.4) +
#   labs(
#     title = "Essentiality score summary (enriched genes ≥3 pools)",
#     subtitle = paste0("n = ", n_enriched_3of4),
#     x = NULL,
#     y = "avg_ess_score"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title   = element_text(face = "bold", size = 14),
#     axis.text    = element_text(face = "bold", size = 14),
#     plot.title   = element_text(face = "bold", size = 16),
#     plot.subtitle= element_text(face = "bold", size = 13)
#   )
# 
# # Print
# p_hist
# 
# p_ecdf
# 
# p_violin
# 
# ############## all 4
# # 3) Violin + box (compact summary + tails)
# n_enriched_4of4 <- enriched_genes_4 %>%
#   filter(!is.na(avg_ess_score)) %>%
#   nrow()
# 
# p_violin_4 <- ggplot(enriched_genes_4, aes(x = "", y = as.numeric(avg_ess_score))) +
#   geom_violin(trim = FALSE, alpha = 0.6) +
#   geom_boxplot(width = 0.15, outlier.alpha = 0.4) +
#   labs(
#     title = "Essentiality score summary (enriched genes 4/4 pools)",
#     subtitle = paste0("n = ", n_enriched_4of4),
#     x = NULL,
#     y = "avg_ess_score"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title   = element_text(face = "bold", size = 14),
#     axis.text    = element_text(face = "bold", size = 14),
#     plot.title   = element_text(face = "bold", size = 16),
#     plot.subtitle= element_text(face = "bold", size = 13)
#   )
# # Print
# p_violin_4
# 
# 
# ############### now it's look at the parents 
# # 1) Identify "parent ~0 CPM" genes (across ALL parent pools)
# #    Here: require parent_cpm == 0 in >=3 pools (and separately 4/4).
# 
# parent_zero_tbl <- parent_cpm %>%
#   mutate(parent_cpm = as.numeric(parent_cpm)) %>%
#   group_by(standard_name) %>%
#   summarise(
#     n_pools_zero = sum(parent_cpm == 0, na.rm = TRUE),
#     n_pools_obs  = n(),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     parent_zero_3of4 = n_pools_zero >= 3,
#     parent_zero_4of4 = n_pools_zero == 4
#   )
# 
# enriched_genes_parent0_3of4 <- rf_vs_gale_tbl %>%
#   inner_join(parent_zero_tbl %>% filter(parent_zero_3of4), by = "standard_name")
# 
# enriched_genes4_parent0_4of4 <- rf_vs_gale_tbl %>%
#   inner_join(parent_zero_tbl %>% filter(parent_zero_4of4), by = "standard_name")
# 
# 
# n_3of4 <- enriched_genes_parent0_3of4 %>%
#   filter(!is.na(avg_ess_score)) %>%
#   nrow()
# 
# n_4of4 <- enriched_genes4_parent0_4of4 %>%
#   filter(!is.na(avg_ess_score)) %>%
#   nrow()
# 
# p_violin_parent0_3of4 <- ggplot(enriched_genes_parent0_3of4,
#                                 aes(x = "", y = as.numeric(avg_ess_score))) +
#   geom_violin(trim = FALSE, alpha = 0.6) +
#   geom_boxplot(width = 0.15, outlier.alpha = 0.4) +
#   labs(
#     title = paste0(
#       "Essentiality score (enriched genes ≥3 pools; parent CPM=0 in ≥3 pools)  n=", n_3of4
#     ),
#     x = NULL,
#     y = "avg_ess_score"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title  = element_text(face = "bold", size = 14),
#     axis.text   = element_text(face = "bold", size = 14),
#     plot.title  = element_text(face = "bold", size = 16)
#   )
# 
# p_violin_parent0_4of4 <- ggplot(enriched_genes4_parent0_4of4,
#                                 aes(x = "", y = as.numeric(avg_ess_score))) +
#   geom_violin(trim = FALSE, alpha = 0.6) +
#   geom_boxplot(width = 0.15, outlier.alpha = 0.4) +
#   labs(
#     title = paste0(
#       "Essentiality score (enriched genes 4/4 pools; parent CPM=0 in 4/4 pools)  n=", n_4of4
#     ),
#     x = NULL,
#     y = "avg_ess_score"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title  = element_text(face = "bold", size = 14),
#     axis.text   = element_text(face = "bold", size = 14),
#     plot.title  = element_text(face = "bold", size = 16)
#   )
# 
# p_violin_parent0_3of4
# 
# p_violin_parent0_4of4
# 


