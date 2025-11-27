# THIS FILE CONTAINS: 
# prune()
# run_pcoa_analysis()
# run_fava_minimal()


# prune  ------------------------------------------------------------------
prune <- function(
  df,
  pattern_asal   = "^ASAL",
  pattern_mig    = "^G",
  pattern = NULL,          # e.g. "1_6" to select rarefied
  keep_cols      = "clade_name",
  sd_threshold   = 0,              # SD cut across ALL selected samples
  prevalence_threshold = 0   # across all selected samples (0–1)
) {
# First, restrict to rarefied/unrarefied if pattern provided
  if (!is.null(pattern)) {
    df <- df %>% select(any_of(keep_cols), matches(pattern))
  }
# Identify ASAL and miG columns from the (possibly restricted) df
  asal_cols   <- grep(pattern_asal, names(df), value = TRUE)
  mig_cols    <- grep(pattern_mig,  names(df), value = TRUE)
  sample_cols <- c(asal_cols, mig_cols)
# 1) Keep only clades present in BOTH groups
  df_common <- df %>%
    rowwise() %>%
    mutate(
      present_asal = any(c_across(all_of(asal_cols)) > 0, na.rm = TRUE),
      present_mig  = any(c_across(all_of(mig_cols))  > 0, na.rm = TRUE)
    ) %>%
    filter(present_asal & present_mig) %>%
    ungroup() %>%
    select(-present_asal, -present_mig)
# 2) Apply SD + prevalence filters
  df_common %>%
    rowwise() %>%
    mutate(
      row_sd = sd(c_across(all_of(sample_cols)), na.rm = TRUE),
      prevalence_all  = mean(c_across(all_of(sample_cols)) > 0, na.rm = TRUE)
    ) %>%
    filter(
      row_sd > sd_threshold,
      prevalence_all  >= prevalence_threshold
    ) %>%
    ungroup() %>%
    select(-row_sd, -prevalence_all) %>% 
    column_to_rownames("clade_name") %>%
    as.matrix()
}
# run_pcoa_analysis  ------------------------------------------------------------------
# Define a function to run PCoA + PERMANOVA + BETADISPER
run_pcoa_analysis <- function(dist_method, data, group) {
    library(vegan)
  # 1. Distance matrix
  dist_mat <- vegdist(data, method = dist_method)

  # 2. PCoA (classical MDS)
  pcoa <- cmdscale(dist_mat, eig = TRUE, k = 2)  # keep top 2 axes

  # 3. Format for plotting
  coords <- as_tibble(pcoa$points) %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    mutate(group = group) %>%
    mutate(sample = rownames(pcoa$points), .before = PCoA1) %>%
    mutate(sample_type = ifelse(grepl("ASAL0(1[5-9]|20)", sample), "pooled", "individual"))

  # 4. PERMANOVA
  permanova <- adonis2(dist_mat ~ group, permutations = 999)

  # 5. BETADISPER
  bd <- betadisper(dist_mat, group)
  bd_test <- anova(bd)

  # 6. Variance explained
  eigvals <- pcoa$eig / sum(pcoa$eig)
  x_lab <- paste0("PCoA1 (", round(eigvals[1] * 100, 1), "%)")
  y_lab <- paste0("PCoA2 (", round(eigvals[2] * 100, 1), "%)")

  # 7. Plot
  centroids <- coords %>%
    group_by(group) %>%
    summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop")

  annot <- paste0("permanova \nR² = ", round(permanova$R2[1], 3),
                  " p = ", format.pval(permanova$`Pr(>F)`[1], digits = 3, eps = 0.001),
                  "\nbetadisper \np = ", format.pval(bd_test$`Pr(>F)`[1], digits = 3))

  p <- ggplot(coords, aes(x = PCoA1, y = PCoA2, color = group, shape = sample_type)) +
    geom_point(size = 5, alpha = 1) +
  stat_ellipse(
    aes(x = PCoA1, y = PCoA2, color = group, group = group),
    level = 0.95, linetype = "dashed", linewidth = 0.5,
    inherit.aes = FALSE
  ) +
    geom_point(data = centroids, aes(x = PCoA1, y = PCoA2),
               color = "red", shape = 18, size = 5, inherit.aes = FALSE) +
    geom_text(data = centroids, aes(x = PCoA1, y = PCoA2, label = group),
              vjust = -1, fontface = "bold", color = "black", inherit.aes = FALSE) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, label = annot,
             size = 4, fontface = "italic", color = "black") +
    scale_color_manual(values = c("ASAL" = "black", "miG" = "grey75")) +
    labs(
      #title = paste(dist_method, "distance"),
      x = x_lab,
      y = y_lab
    ) +
    ggpubr::theme_pubr() +
    theme(
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.position = "none"
    )

  return(list(df = coords, plot = p, permanova = permanova, betadisper = bd_test))
}
# run_fava_minimal ------------------------------------------------------------------
run_fava_minimal <- function(
  df,
  levels = 1:7,
  subsets = c("rarefied", "unrarefied"),
  bootstrap_reps = 500,
  miG_n = 14,
  ASAL_n = 14,
  mix_each = 7,
  seed = NULL
) {
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(tibble)

  if (!is.null(seed)) set.seed(seed)

  # helper: build FAVA-ready matrix for one level/subset
  prep_fava_matrix <- function(tbl, level_n, subset_label) {
    sylph_n <- tbl %>% dplyr::filter(stringr::str_count(clade_name, "\\|") == level_n)
    sylph_sub <- if (subset_label == "rarefied") {
      sylph_n %>% dplyr::select(clade_name, dplyr::matches("1_6"))
    } else {
      sylph_n %>% dplyr::select(clade_name, !dplyr::matches("1_"))
    }
    if (ncol(sylph_sub) <= 1) return(NULL)

    data_fava <- sylph_sub %>%
      tidyr::pivot_longer(-clade_name, names_to = "sample", values_to = "abundance") %>%
      tidyr::pivot_wider(names_from = clade_name, values_from = abundance)

    if (nrow(data_fava) == 0 || ncol(data_fava) <= 1) return(NULL)

    M <- data_fava %>% tibble::column_to_rownames("sample") %>% as.matrix()
    rs <- rowSums(M)
    keep <- rs > 0
    if (!any(keep)) return(NULL)
    M <- sweep(M[keep, , drop = FALSE], 1, rowSums(M[keep, , drop = FALSE]), "/")

    df_ready <- as.data.frame(M) %>%
      tibble::rownames_to_column("sample") %>%
      dplyr::mutate(group = ifelse(stringr::str_detect(sample, "ASAL"), "ASAL", "miG"),
                    .before = sample)

    list(df = df_ready, K = ncol(df_ready) - 2)
  }

  out <- list()

  for (lvl in levels) {
    for (subset_label in subsets) {
      prep <- prep_fava_matrix(df, lvl, subset_label)
      if (is.null(prep)) next

      dat <- prep$df
      K <- prep$K
      if (K <= 0) next

      asal_data <- dat %>% dplyr::filter(group == "ASAL")
      mig_data  <- dat %>% dplyr::filter(group == "miG")

      # (2) ASAL all together (n = 14), no bootstrap
      if (nrow(asal_data) > 0) {
        asal_all <- asal_data %>% dplyr::select(-group)
        fava_asal_all <- FAVA::fava(asal_all, group = NULL, K = K, normalized = FALSE)
        out[[length(out) + 1]] <- tibble::tibble(
          scenario = "ASAL_all14",
          level = lvl, subset = subset_label,
          bootstrap = NA_integer_, FAVA = fava_asal_all
        )
      }

      # (1) miG bootstrapped n = 14
      if (nrow(mig_data) > 0) {
        for (b in seq_len(bootstrap_reps)) {
          mig_boot <- mig_data %>% dplyr::slice_sample(n = miG_n, replace = TRUE) %>%
            dplyr::select(-group)
          fava_mig <- FAVA::fava(mig_boot, group = NULL, K = K, normalized = FALSE)
          out[[length(out) + 1]] <- tibble::tibble(
            scenario = "miG_boot14",
            level = lvl, subset = subset_label,
            bootstrap = b, FAVA = fava_mig
          )
        }
      }

      # (3) mixed bootstrapped 7 + 7
      if (nrow(asal_data) > 0 && nrow(mig_data) > 0) {
        for (b in seq_len(bootstrap_reps)) {
          asal_boot <- asal_data %>% dplyr::slice_sample(n = mix_each, replace = TRUE)
          mig_boot  <- mig_data  %>% dplyr::slice_sample(n = mix_each, replace = TRUE)
          mix_7_7   <- dplyr::bind_rows(asal_boot, mig_boot) %>% dplyr::select(-group)
          fava_mix  <- FAVA::fava(mix_7_7, group = NULL, K = K, normalized = FALSE)
          out[[length(out) + 1]] <- tibble::tibble(
            scenario = "mix_boot7_7",
            level = lvl, subset = subset_label,
            bootstrap = b, FAVA = fava_mix
          )
        }
      }
    }
  }

  dplyr::bind_rows(out)
}
