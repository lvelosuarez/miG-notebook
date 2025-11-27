# helper_functions.R
# ------------------------------------------------------------------------------
# This file provides three utilities for salivary metagenomics analyses:
#   - prune():    filter a feature table by presence in both groups, variability,
#                 and prevalence, optionally on a rarefied/unrarefied subset.
#   - run_pcoa_analysis(): compute PCoA + PERMANOVA + BETADISPER and return a
#                 publication-ready ggplot alongside stats.
#   - run_fava_minimal(): build level- and subset-specific matrices and run FAVA
#                 across three scenarios (miG bootstraps, ASAL all, mixed 7+7).
# ------------------------------------------------------------------------------

# NOTE: This script avoids library() calls in the global environment.
# It uses pkg::fun() where possible. Make sure these packages are installed:
# dplyr, tidyr, tibble, stringr, ggplot2, ggpubr, vegan, FAVA

# ------------------------------------------------------------------------------
# prune()
# ------------------------------------------------------------------------------

#' Prune a feature (taxon) table by group presence, variability, and prevalence
#'
#' @param df A data.frame/tibble with one feature-identifying column and
#'   sample-abundance columns (non-negative numeric). Typically contains both
#'   ASAL and miG columns.
#' @param pattern_asal Regex to identify ASAL columns. Default: "^ASAL".
#' @param pattern_mig  Regex to identify miG columns.  Default: "^G".
#' @param pattern Optional regex to restrict columns first (e.g. "1_6" for rarefied).
#'   If supplied, only `keep_cols` + columns matching `pattern` are kept.
#' @param keep_cols Name(s) of non-sample columns to retain. First element must
#'   be the feature name column (defaults to "clade_name").
#' @param sd_threshold Keep features whose SD across *all selected sample columns*
#'   exceeds this value. Default: 0.
#' @param prevalence_threshold Keep features whose prevalence (fraction of samples
#'   with value > 0) across *all selected sample columns* is >= this value (0–1).
#'   Default: 0.
#' @return A numeric matrix with features as rownames and selected samples as
#'   columns. If nothing passes filters, returns a 0x0 matrix.
#' @examples
#' # M <- prune(df, pattern = "1_6", sd_threshold = 0.001, prevalence_threshold = 0.2)
prune <- function(
  df,
  pattern_asal   = "^ASAL",
  pattern_mig    = "^G",
  pattern        = NULL,     # e.g., "1_6" to select rarefied
  keep_cols      = "clade_name",
  sd_threshold   = 0,
  prevalence_threshold = 0
) {
  # Basic checks
  stopifnot(is.data.frame(df))
  stopifnot(length(keep_cols) >= 1)
  feature_col <- keep_cols[1]
  if (!feature_col %in% names(df)) {
    stop(sprintf("Feature column '%s' not found in df.", feature_col))
  }

  # If a subset pattern is provided, keep only those columns (plus keep_cols)
  if (!is.null(pattern)) {
    df <- dplyr::select(df, dplyr::any_of(keep_cols), dplyr::matches(pattern))
  }

  # Identify ASAL and miG columns (after any restriction above)
  asal_cols   <- grep(pattern_asal, names(df), value = TRUE)
  mig_cols    <- grep(pattern_mig,  names(df), value = TRUE)
  sample_cols <- c(asal_cols, mig_cols)

  if (length(asal_cols) == 0 || length(mig_cols) == 0) {
    warning("No ASAL or no miG columns detected after filtering—returning empty matrix.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  # Ensure numeric for sample columns
  df[sample_cols] <- lapply(df[sample_cols], function(x) as.numeric(x))

  # 1) keep only features present (>0) in BOTH groups at least once
  df_common <- df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      present_asal = any(dplyr::c_across(dplyr::all_of(asal_cols)) > 0, na.rm = TRUE),
      present_mig  = any(dplyr::c_across(dplyr::all_of(mig_cols))  > 0, na.rm = TRUE)
    ) |>
    dplyr::filter(present_asal & present_mig) |>
    dplyr::ungroup() |>
    dplyr::select(-present_asal, -present_mig)

  if (nrow(df_common) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  # 2) apply SD + prevalence filters across ALL selected sample columns
  res <- df_common |>
    dplyr::rowwise() |>
    dplyr::mutate(
      row_sd = stats::sd(dplyr::c_across(dplyr::all_of(sample_cols)), na.rm = TRUE),
      prevalence_all = mean(dplyr::c_across(dplyr::all_of(sample_cols)) > 0, na.rm = TRUE)
    ) |>
    dplyr::filter(
      row_sd > sd_threshold,
      prevalence_all >= prevalence_threshold
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-row_sd, -prevalence_all)

  if (nrow(res) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  # Construct matrix with feature names as rownames
  if (!feature_col %in% names(res)) {
    stop(sprintf("Feature column '%s' not found after filtering.", feature_col))
  }

  res_mat <- res |>
    tibble::column_to_rownames(feature_col) |>
    as.matrix()

  return(res_mat)
}

# ------------------------------------------------------------------------------
# run_pcoa_analysis()
# ------------------------------------------------------------------------------

#' PCoA + PERMANOVA + BETADISPER with a publication-ready plot
#'
#' @param dist_method Distance metric for vegan::vegdist (e.g., "bray", "euclidean").
#' @param data Numeric matrix/data.frame (samples x features) with non-negative
#'   values. Rows are samples; columns are features.
#' @param group A factor/character vector of group labels (length = nrow(data)).
#' @param pooled_regex Regex to tag pooled samples for shape aesthetics.
#'   Default matches "ASAL015-020".
#' @param seed Optional integer seed for reproducibility of permutations.
#' @return A list with:
#'   \item{df}{data.frame with PCoA1, PCoA2, sample, group, sample_type}
#'   \item{plot}{ggplot object}
#'   \item{permanova}{adonis2 result}
#'   \item{betadisper}{ANOVA table on betadisper}
#' @examples
#' # out <- run_pcoa_analysis("bray", M, group)
run_pcoa_analysis <- function(
  dist_method,
  data,
  group,
  pooled_regex = "ASAL0(1[5-9]|20)",
  seed = NULL
) {
  stopifnot(nrow(data) == length(group))

  if (!is.null(seed)) set.seed(seed)

  # 1. Distance matrix
  dist_mat <- vegan::vegdist(data, method = dist_method)

  # 2. PCoA
  pcoa <- stats::cmdscale(dist_mat, eig = TRUE, k = 2)  # top 2 axes

  # 3. Coordinates for plotting
  coords <- tibble::as_tibble(pcoa$points, .name_repair = "minimal") |>
    dplyr::rename(PCoA1 = 1, PCoA2 = 2) |>
    dplyr::mutate(
      group = group,
      sample = rownames(pcoa$points),
      .before = PCoA1
    ) |>
    dplyr::mutate(
      sample_type = ifelse(stringr::str_detect(sample, pooled_regex), "pooled", "individual")
    )

  # 4. PERMANOVA
  permanova <- vegan::adonis2(dist_mat ~ group, permutations = 999)

  # 5. BETADISPER
  bd <- vegan::betadisper(dist_mat, group)
  bd_test <- stats::anova(bd)

  # 6. Variance explained labels
  eigvals <- pcoa$eig / sum(pcoa$eig)
  x_lab <- paste0("PCoA1 (", round(eigvals[1] * 100, 1), "%)")
  y_lab <- paste0("PCoA2 (", round(eigvals[2] * 100, 1), "%)")

  # 7. Group centroids for labeling
  centroids <- coords |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      PCoA1 = mean(PCoA1),
      PCoA2 = mean(PCoA2),
      .groups = "drop"
    )

  annot <- paste0(
    "permanova \nR² = ", round(permanova$R2[1], 3),
    "  p = ", stats::format.pval(permanova$`Pr(>F)`[1], digits = 3, eps = 0.001),
    "\nbetadisper \n p = ", stats::format.pval(bd_test$`Pr(>F)`[1], digits = 3)
  )

  p <- ggplot2::ggplot(coords, ggplot2::aes(x = PCoA1, y = PCoA2, color = group, shape = sample_type)) +
    ggplot2::geom_point(size = 5, alpha = 1) +
    ggplot2::stat_ellipse(
      ggplot2::aes(x = PCoA1, y = PCoA2, color = group, group = group),
      level = 0.95, linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE
    ) +
    ggplot2::geom_point(data = centroids, ggplot2::aes(x = PCoA1, y = PCoA2),
                        color = "red", shape = 18, size = 5, inherit.aes = FALSE) +
    ggplot2::geom_text(data = centroids, ggplot2::aes(x = PCoA1, y = PCoA2, label = group),
                       vjust = -1, fontface = "bold", color = "black", inherit.aes = FALSE) +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, label = annot,
                      size = 4, fontface = "italic", color = "black") +
    ggplot2::scale_color_manual(values = c("ASAL" = "black", "miG" = "grey75")) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 18, face = "bold"),
      axis.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_blank(),
      legend.position = "none"
    )

  list(df = coords, plot = p, permanova = permanova, betadisper = bd_test)
}

# ------------------------------------------------------------------------------
# run_fava_minimal()
# ------------------------------------------------------------------------------

#' Minimal FAVA workflow across levels and subsets with bootstrap scenarios
#'
#' @param df Long-ish table with a `clade_name` column and sample columns
#'   (wide) containing abundances/frequencies. Column names must let you
#'   distinguish ASAL vs miG (via "ASAL" substring).
#' @param levels Integer vector of taxonomic depths to evaluate. The function
#'   counts '|' in clade_name and keeps rows with exactly `levels`.
#' @param subsets Character vector among c("rarefied", "unrarefied").
#'   "rarefied" selects columns matching "1_6"; "unrarefied" selects columns
#'   NOT matching "1_".
#' @param bootstrap_reps Number of bootstrap replicates for miG and mix scenarios.
#' @param miG_n Bootstrap sample size for miG-only scenario.
#' @param ASAL_n Unused in current version (kept for API stability).
#' @param mix_each Number drawn with replacement from each group in mix scenario.
#' @param seed Optional integer seed for reproducibility.
#' @return A tibble with columns: scenario, level, subset, bootstrap, FAVA.
#'   Each row holds the scalar FAVA value for that run.
#' @examples
#' # fava_tbl <- run_fava_minimal(df, levels = 1:7, subsets = c("rarefied","unrarefied"))
run_fava_minimal <- function(
  df,
  levels = 1:7,
  subsets = c("rarefied", "unrarefied"),
  bootstrap_reps = 500,
  miG_n = 14,
  ASAL_n = 14,      # kept for API stability; not used directly
  mix_each = 7,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # helper: build one FAVA-ready matrix for a level/subset
  prep_fava_matrix <- function(tbl, level_n, subset_label) {
    stopifnot("clade_name" %in% names(tbl))

    # keep rows with exactly 'level_n' pipes in clade_name
    sylph_n <- dplyr::filter(tbl, stringr::str_count(.data$clade_name, "\\|") == level_n)

    sylph_sub <- if (subset_label == "rarefied") {
      dplyr::select(sylph_n, .data$clade_name, dplyr::matches("1_6"))
    } else {
      dplyr::select(sylph_n, .data$clade_name, !dplyr::matches("1_"))
    }

    if (ncol(sylph_sub) <= 1L) return(NULL)

    data_fava <- sylph_sub |>
      tidyr::pivot_longer(-.data$clade_name, names_to = "sample", values_to = "abundance") |>
      tidyr::pivot_wider(names_from = .data$clade_name, values_from = .data$abundance)

    if (nrow(data_fava) == 0L || ncol(data_fava) <= 1L) return(NULL)

    # rows=samples, cols=features; normalize rows to sum=1 after removing all-zero rows
    M <- data_fava |>
      tibble::column_to_rownames("sample") |>
      as.matrix()

    rs <- rowSums(M, na.rm = TRUE)
    keep <- rs > 0
    if (!any(keep)) return(NULL)

    M <- sweep(M[keep, , drop = FALSE], 1, rowSums(M[keep, , drop = FALSE]), "/")

    df_ready <- M |>
      as.data.frame() |>
      tibble::rownames_to_column("sample") |>
      dplyr::mutate(group = ifelse(stringr::str_detect(.data$sample, "ASAL"), "ASAL", "miG"),
                    .before = "sample")

    list(df = df_ready, K = ncol(df_ready) - 2L)
  }

  out <- list()

  for (lvl in levels) {
    for (subset_label in subsets) {
      prep <- prep_fava_matrix(df, lvl, subset_label)
      if (is.null(prep)) next

      dat <- prep$df
      K <- prep$K
      if (K <= 0) next

      asal_data <- dplyr::filter(dat, .data$group == "ASAL")
      mig_data  <- dplyr::filter(dat, .data$group == "miG")

      # (2) ASAL all together (n = all ASAL rows), no bootstrap
      if (nrow(asal_data) > 0) {
        asal_all <- dplyr::select(asal_data, - .data$group)
        fava_asal_all <- FAVA::fava(asal_all, group = NULL, K = K, normalized = FALSE)
        out[[length(out) + 1L]] <- tibble::tibble(
          scenario  = "ASAL_all",
          level     = lvl,
          subset    = subset_label,
          bootstrap = NA_integer_,
          FAVA      = fava_asal_all
        )
      }

      # (1) miG bootstrapped n = miG_n
      if (nrow(mig_data) > 0) {
        for (b in seq_len(bootstrap_reps)) {
          mig_boot <- dplyr::slice_sample(mig_data, n = miG_n, replace = TRUE) |>
            dplyr::select(- .data$group)
          fava_mig <- FAVA::fava(mig_boot, group = NULL, K = K, normalized = FALSE)
          out[[length(out) + 1L]] <- tibble::tibble(
            scenario  = "miG_boot",
            level     = lvl,
            subset    = subset_label,
            bootstrap = b,
            FAVA      = fava_mig
          )
        }
      }

      # (3) mixed bootstrapped 7 + 7
      if (nrow(asal_data) > 0 && nrow(mig_data) > 0) {
        for (b in seq_len(bootstrap_reps)) {
          asal_boot <- dplyr::slice_sample(asal_data, n = mix_each, replace = TRUE)
          mig_boot  <- dplyr::slice_sample(mig_data,  n = mix_each, replace = TRUE)
          mix_7_7   <- dplyr::bind_rows(asal_boot, mig_boot) |>
            dplyr::select(- .data$group)
          fava_mix  <- FAVA::fava(mix_7_7, group = NULL, K = K, normalized = FALSE)
          out[[length(out) + 1L]] <- tibble::tibble(
            scenario  = "mix_boot_7_7",
            level     = lvl,
            subset    = subset_label,
            bootstrap = b,
            FAVA      = fava_mix
          )
        }
      }
    }
  }

  if (length(out) == 0) {
    return(tibble::tibble(
      scenario = character(), level = integer(), subset = character(),
      bootstrap = integer(), FAVA = numeric()
    ))
  }

  dplyr::bind_rows(out)
}
