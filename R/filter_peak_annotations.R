#' @title Filter promoter peak-gene annotations with distance-adaptive thresholds
#'
#' @description Applies distance-adaptive chromatin accessibility quantile thresholds to refine promoter peak-gene annotations,
#' integrates gene expression quality adjustment, and supports an optional per-gene fallback rule.
#' Built-in disk caching accelerates repeated runs.
#'
#' @param results_identify_TGs Output list from \code{\link{identify_TGs}}.
#' @param scRNAseq Seurat object with a \code{counts} layer. Used for gene quality score calculation.
#' @param scATACseq Seurat object with a \code{counts} layer. Used for peak accessibility feature calculation.
#' @param forced_feature_calculation Logical.
#' If \code{TRUE}, force recomputation of gene/peak features and overwrite cache files.
#' Defaults to \code{FALSE}: cached RDS files are loaded when available to skip the computationally expensive feature calculation step.
#' @param distance_points Numeric vector of genomic distances (bp) to TSS,
#' used as anchor points for threshold curve fitting.
#' Must be strictly increasing and match the length of \code{min_quantile_points}.
#' Defaults to \code{c(200, 5000, 500000)}.
#' @param min_quantile_points Positive numeric vector of minimum quantile thresholds corresponding to each distance anchor.
#' Must be strictly increasing. Defaults to \code{c(0, 0.6, 1)}.
#' @param distance_quantile_fit_start Named list (\code{a}, \code{b}, \code{c}) of starting values for constrained nonlinear least squares curve fitting.
#' Default values are optimized for the default anchor points.
#' @param one_peak_at_least_per_gene Logical.
#' If \code{TRUE}, guarantee at least one retained peak per target gene:
#' single-peak genes are kept unconditionally;
#' multi-peak genes with zero passing peaks retain the one with highest overall accessibility.
#' Defaults to \code{TRUE}.
#' @param gene_quality_factor_min Numeric in \[0, 1\].
#' Minimum scaling factor for high-quality genes (smaller value = more relaxed thresholds). Defaults to \code{0.8}.
#' @param weight_gene_nonzero_freq Numeric in \[0, 1\].
#' Weight of non-zero cell frequency in the composite gene quality score.
#' @param weight_gene_nonzero_mean Numeric in \[0, 1\].
#' Weight of non-zero mean expression in the composite gene quality score.
#' @param weight_gene_overall_mean Numeric ranging from 0 to 1.
#' Weight of overall mean expression in the composite gene quality score.
#' The three weights must sum to 1. Defaults: \code{0.2}, \code{0.2}, \code{0.6}.
#' @param min_peak_nonzero_freq_quantile_scale Positive numeric.
#' Scaling factor for the peak non-zero frequency quantile threshold (smaller value = more relaxed thresholds). Default: \code{0.7}.
#' @param min_peak_nonzero_mean_scale Positive numeric.
#' Scaling factor for the peak non-zero mean accessibility quantile threshold (smaller value = more relaxed thresholds). Default: \code{0.7}.
#' @param min_peak_overall_mean_quantile_scale Positive numeric.
#' Scaling factor for the peak overall mean accessibility quantile threshold (smaller value = more relaxed thresholds). Default: \code{1}.
#'
#' @returns Updated \code{results_identify_TGs} list with the following elements:
#' \itemize{
#'   \item \code{peak_anno_promoter}: Filtered peak-gene annotation table with quality metrics.
#'   \item \code{TGs}: Character vector of retained target gene symbols.
#'   \item \code{peak_anno_promoter_all}: Original unfiltered annotation as backup.
#'   \item \code{params}: Full parameter record for reproducibility.
#' }
#'
#' @details
#' The function implements a 6-step filtering pipeline:
#' \enumerate{
#'   \item \strong{Feature quantification with cache}
#'   Compute three metrics (non-zero frequency, non-zero mean, overall mean) and their genome-wide percentile ranks for all genes and peaks.
#'   Results are cached to disk for reuse.
#'
#'   \item \strong{Annotation expansion}
#'   Expand each peak to one row per associated gene (primary + flanking genes), and retain only genes present in the scRNA-seq dataset.
#'   Flanking genes use the average distance across all flanking loci.
#'
#'   \item \strong{Distance-threshold curve fitting}
#'   Fit an inverse-distance model \eqn{min\_quantile = -1/(a \cdot distance + b) + c} via constrained nonlinear least squares,
#'   where required stringency increases with distance from TSS.
#'
#'   \item \strong{Gene quality adjustment}
#'   Calculate a weighted composite gene quality score from expression metrics.
#'   High-quality genes receive a lower scaling factor (relaxed thresholds), while low-quality genes remain under stricter filtering.
#'
#'   \item \strong{Per-metric threshold scaling}
#'   Apply independent scaling factors to the base threshold for three peak accessibility metrics. All final thresholds are clamped to the \[0, 1\] range.
#'
#'   \item \strong{Per-gene peak filtering}
#'   Retain peaks that pass all three metric thresholds. When enabled, the fallback rule ensures every gene retains at least its strongest peak.
#' }
#'
#' @note
#' \itemize{
#'   \item Distances outside the anchor point range are clamped to boundary values to avoid extrapolation artifacts.
#'   \item Feature cache files are stored under \code{./1.2 Data Preprocessing - Identify TGs By Annotation/}.
#'   \item Input Seurat objects must contain a \code{counts} layer for raw metric calculation.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results_identify_TGs = base::readRDS(
#'   "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds"
#' )
#' scRNAseq = readRDS("./scRNAseq.rds")
#' scATACseq = base::readRDS("./scATACseq.rds")
#' results_filtered = filter_peak_annotations(
#'   results_identify_TGs = results_identify_TGs,
#'   scRNAseq = scRNAseq,
#'   scATACseq = scATACseq,
#'   distance_points = c(200, 3000, 500000),
#'   min_quantile_points = c(0, 0.3, 1),
#'   one_peak_at_least_per_gene = TRUE
#' )
#' base::saveRDS(
#'   results_identify_TGs,
#'   file = "./1.2 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds"
#' )
#' }
filter_peak_annotations = function(
    results_identify_TGs,
    scRNAseq,
    scATACseq,
    forced_feature_calculation = FALSE,
    distance_points = base::c(200, 5000, 500000),
    min_quantile_points = base::c(0, 0.6, 1),
    distance_quantile_fit_start = base::list(a = 0.1648084, b = 1.2657021, c = 0.7791173),
    one_peak_at_least_per_gene = TRUE,
    gene_quality_factor_min = 0.8,
    weight_gene_nonzero_freq = 0.2,
    weight_gene_nonzero_mean = 0.2,
    weight_gene_overall_mean = 0.6,
    min_peak_nonzero_freq_quantile_scale = 0.7,
    min_peak_nonzero_mean_scale = 0.7,
    min_peak_overall_mean_quantile_scale = 1
) {
  ### Input check
  # Validate distance anchor points
  if (base::length(distance_points) != base::length(min_quantile_points)) {
    stop("'distance_points' and 'min_quantile_points' must have the same length.")
  }
  if (base::any(distance_points < 0)) {
    warning("It is recommended that 'distance_points' contain only non-negative genomic distances.")
  }
  if (base::is.unsorted(distance_points, strictly = TRUE)) {
    stop("'distance_points' must be sorted in strictly increasing order.")
  }

  # Validate quantile anchor points
  if (base::any(min_quantile_points < 0) || base::any(min_quantile_points > 1)) {
    warning("It is recommended that 'min_quantile_points' contain values within the range [0, 1].")
  }
  if (base::is.unsorted(min_quantile_points, strictly = TRUE)) {
    stop("'min_quantile_points' must be sorted in strictly increasing order.")
  }

  # Validate curve fitting starting parameters
  if (!base::is.list(distance_quantile_fit_start) ||
      !base::all(c("a", "b", "c") %in% base::names(distance_quantile_fit_start))) {
    stop("'distance_quantile_fit_start' must be a named list with elements 'a', 'b', 'c'.")
  }
  if (distance_quantile_fit_start$a <= 0) {
    stop("'distance_quantile_fit_start$a' must be a positive numeric value.")
  }

  # Validate fallback rule flag
  if (!base::is.logical(one_peak_at_least_per_gene) ||
      base::length(one_peak_at_least_per_gene) != 1 ||
      base::is.na(one_peak_at_least_per_gene)) {
    stop("'one_peak_at_least_per_gene' must be a single non-NA logical value (TRUE or FALSE).")
  }

  # Validate gene quality parameters
  if (gene_quality_factor_min < 0) {
    stop("'gene_quality_factor_min' must be a non-negative numeric value.")
  }
  if (gene_quality_factor_min > 1) {
    warning("It is recommended that 'gene_quality_factor_min' is not greater than 1.")
  }

  if (weight_gene_nonzero_freq < 0 || weight_gene_nonzero_freq > 1) {
    stop("'weight_gene_nonzero_freq' must be a numeric value between 0 and 1.")
  }
  if (weight_gene_nonzero_mean < 0 || weight_gene_nonzero_mean > 1) {
    stop("'weight_gene_nonzero_mean' must be a numeric value between 0 and 1.")
  }
  if (weight_gene_overall_mean < 0 || weight_gene_overall_mean > 1) {
    stop("'weight_gene_overall_mean' must be a numeric value between 0 and 1.")
  }
  if (base::abs(
    weight_gene_nonzero_freq + weight_gene_nonzero_mean + weight_gene_overall_mean - 1
  ) > 1e-6) {
    stop("The sum of the three gene quality weights must equal 1.")
  }

  # Validate peak threshold scaling factors
  if (min_peak_nonzero_freq_quantile_scale < 0) {
    stop("'min_peak_nonzero_freq_quantile_scale' must be a positive numeric value.")
  }
  if (min_peak_nonzero_freq_quantile_scale > 1) {
    warning("'min_peak_nonzero_freq_quantile_scale' is greater than 1.")
  }
  if (min_peak_nonzero_mean_scale < 0) {
    stop("'min_peak_nonzero_mean_scale' must be a positive numeric value.")
  }
  if (min_peak_nonzero_mean_scale > 1) {
    warning("'min_peak_nonzero_mean_scale' is greater than 1.")
  }
  if (min_peak_overall_mean_quantile_scale < 0) {
    stop("'min_peak_overall_mean_quantile_scale' must be a positive numeric value.")
  }
  # if (min_peak_overall_mean_quantile_scale > 1) {
  #   warning("'min_peak_overall_mean_quantile_scale' is greater than 1.")
  # }
  # if (min_peak_nonzero_freq_quantile_scale > min_peak_overall_mean_quantile_scale) {
  #   warning("'min_peak_nonzero_freq_quantile_scale' is greater than 'min_peak_overall_mean_quantile_scale'.")
  # }
  # if (min_peak_nonzero_mean_scale > min_peak_overall_mean_quantile_scale) {
  #   warning("'min_peak_nonzero_mean_scale' is greater than 'min_peak_overall_mean_quantile_scale'.")
  # }

  ### Start.
  t_start = base::Sys.time()
  message("Run: Filtering peak annotation results. Time: ", t_start, ".")
  cache_dir = "./1.2 Data Preprocessing - Identify TGs By Annotation/"
  if (!base::dir.exists(cache_dir)) base::dir.create(cache_dir, recursive = TRUE)

  ### Compute or load gene and peak features (with disk cache).
  if (
    forced_feature_calculation ||
    !base::file.exists("./1.2 Data Preprocessing - Identify TGs By Annotation/gene_features_all.rds") ||
    !base::file.exists("./1.2 Data Preprocessing - Identify TGs By Annotation/peak_features_all.rds")
  ) {
    # Extract raw count matrices from Seurat objects.
    message("Extract raw gene-by-cell count matrix.")
    gene_counts_all = Seurat::GetAssayData(scRNAseq, layer = "counts")

    message("Extract raw peak-by-cell count matrix.")
    peak_counts_all = Seurat::GetAssayData(scATACseq, layer = "counts")

    # Compute three per-gene statistics: detection rate, mean in non-zero cells, overall mean.
    message("Compute per-gene statistics: non-zero frequency, mean of non-zero values, overall mean.")
    n_cells_rna = base::ncol(gene_counts_all)
    gene_row_sums = Matrix::rowSums(gene_counts_all)
    gene_n_nonzero = Matrix::rowSums(gene_counts_all != 0)
    gene_features_all = base::data.frame(
      nonzero_freq = gene_n_nonzero / n_cells_rna,
      nonzero_mean = gene_row_sums / gene_n_nonzero,
      overall_mean = gene_row_sums / n_cells_rna,
      row.names = base::rownames(gene_counts_all),
      check.names = FALSE
    )
    base::rm(gene_counts_all, gene_row_sums, gene_n_nonzero, n_cells_rna)

    # Compute three per-peak statistics: detection rate, mean in non-zero cells, overall mean.
    message("Compute per-peak statistics: non-zero frequency, mean of non-zero values, overall mean.")
    n_cells_atac = base::ncol(peak_counts_all)
    peak_row_sums = Matrix::rowSums(peak_counts_all)
    peak_n_nonzero = Matrix::rowSums(peak_counts_all != 0)
    peak_features_all = base::data.frame(
      nonzero_freq = peak_n_nonzero / n_cells_atac,
      nonzero_mean = peak_row_sums / peak_n_nonzero,
      overall_mean = peak_row_sums / n_cells_atac,
      row.names = base::rownames(peak_counts_all),
      check.names = FALSE
    )
    base::rm(peak_counts_all, peak_row_sums, peak_n_nonzero, n_cells_atac)

    # Convert raw statistics to percentile ranks and cache to disk.
    message("Add percentile ranks for each gene feature (0 = lowest, 1 = highest).")
    gene_features_all[base::paste0(base::names(gene_features_all), "_quantile")] =
      base::lapply(gene_features_all, function(x) base::rank(x) / base::length(x))
    base::colnames(gene_features_all) = base::paste0("gene_", base::colnames(gene_features_all))
    base::saveRDS(gene_features_all, "./1.2 Data Preprocessing - Identify TGs By Annotation/gene_features_all.rds")

    message("Add percentile ranks for each peak feature (0 = lowest, 1 = highest).")
    peak_features_all[base::paste0(base::names(peak_features_all), "_quantile")] =
      base::lapply(peak_features_all, function(x) base::rank(x) / base::length(x))
    base::colnames(peak_features_all) = base::paste0("peak_", base::colnames(peak_features_all))
    base::saveRDS(peak_features_all, "./1.2 Data Preprocessing - Identify TGs By Annotation/peak_features_all.rds")
  } else {
    # Load pre-computed features from cache to skip re-calculation.
    message("Read feature data of genes and peaks.")
    gene_features_all = base::readRDS("./1.2 Data Preprocessing - Identify TGs By Annotation/gene_features_all.rds")
    peak_features_all = base::readRDS("./1.2 Data Preprocessing - Identify TGs By Annotation/peak_features_all.rds")
  }

  ### Expand peak annotations with flanking genes.
  message("Expand peak annotation to one row per peak-gene pair (primary + flanking genes).")
  anno = results_identify_TGs[["peak_anno_promoter"]]
  anno$peak = base::rownames(anno)
  base::rownames(anno) = NULL
  anno = anno[, c("peak", base::colnames(anno)[base::colnames(anno) != "peak"])]
  anno[base::is.na(anno)] = ""

  # Expand one peak row into multiple rows, one per flanking target gene.
  message("Expand rows for each unique flanking gene.")
  anno_all = base::do.call(base::rbind, base::lapply(base::seq_len(base::nrow(anno)), function(i) {
    cur = anno[i, ]

    # Parse flanking genes, exclude primary gene (SYMBOL), keep non-empty.
    genes = base::unlist(base::strsplit(cur$flank_geneSymbols, ";"))
    genes = genes[genes != ""]
    genes_remain = genes[genes != cur$SYMBOL]
    n = base::length(genes_remain)

    if (n == 0) return(cur) # Return original row if no additional genes

    # Compute average distance across all flanking genes.
    avg_dist = base::round(base::mean(base::as.numeric(base::unlist(base::strsplit(cur$flank_gene_distances, ";")))))

    # Duplicate row for each remaining flanking gene.
    new_rows = cur[base::rep(1, n), ]
    new_rows$SYMBOL = genes_remain
    new_rows$distanceToTSS = avg_dist
    new_rows[, c(
      "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand", "geneId",
      "transcriptId", "ENSEMBL", "GENENAME",
      "flank_txIds", "flank_geneIds", "flank_gene_distances", "flank_geneSymbols"
    )] = ""

    # Combine original and expanded rows.
    base::rbind(cur, new_rows)
  }))

  ### Merge feature data and prepare for filtering.
  message("Prepare annotation: compute absolute TSS distance, merge gene and peak features.")
  # Keep only genes present in the scRNA-seq dataset.
  anno_all = anno_all[anno_all$SYMBOL %in% base::rownames(scRNAseq), ]
  # Use absolute distance to TSS for threshold evaluation.
  anno_all$absDistanceToTSS = base::abs(anno_all$distanceToTSS)
  # Attach per-gene and per-peak quality features.
  anno_all = base::cbind(anno_all, gene_features_all[anno_all$SYMBOL, ]); rm(gene_features_all)
  anno_all = base::cbind(anno_all, peak_features_all[anno_all$peak, ]); rm(peak_features_all)

  ### Fit inverse-distance threshold curve.
  message("Fit inverse-distance model: min_quantile = -1 / (distance * a + b) + c, a>0.")
  fit = minpack.lm::nlsLM(
    min_quantile_points ~ -1 / (distance_points * a + b) + c,
    start = distance_quantile_fit_start,
    # algorithm = "port",
    lower = c(0, -Inf, -Inf),
    upper = c(Inf, Inf, Inf)
  )

  message("Extract fitted coefficients.")
  a = base::as.numeric(stats::coef(fit)["a"])
  b = base::as.numeric(stats::coef(fit)["b"])
  c = base::as.numeric(stats::coef(fit)["c"]); rm(fit)
  message("a = ", base::round(a, digits = 6),
          ", b = ", base::round(b, digits = 6),
          ", c = ", base::round(c, digits = 6))

  ### Filter peaks per target gene.
  message("Initialize vector to store retained peak IDs.")
  peaks_tg_keep = base::character(0)

  message("For each target gene, filter peaks by distance-adaptive quantile threshold.")
  for (tg in results_identify_TGs[["TGs"]]) {
    # Subset annotation for current target gene.
    anno_tg = anno_all[anno_all$SYMBOL == tg, ]

    if (one_peak_at_least_per_gene) {
      # --- Mode 1: with at-least-one-peak fallback ---
      if (base::nrow(anno_tg) > 1) {
        # Clamp distance to anchor range to avoid extrapolation artifacts
        anno_tg$fixed_absDistanceToTSS = base::pmin(
          base::pmax(anno_tg$absDistanceToTSS, base::min(distance_points) + 0.01),
          base::max(distance_points) - 0.01
        )

        # Compute base distance-dependent quantile threshold
        anno_tg$min_quantile = -1 / (anno_tg$fixed_absDistanceToTSS * a + b) + c

        # Compute composite gene quality factor (higher quality → smaller factor → relaxed threshold).
        gene_nonzero_freq_quantile = base::mean(anno_tg$gene_nonzero_freq_quantile)
        gene_nonzero_mean_quantile = base::mean(anno_tg$gene_nonzero_mean_quantile)
        gene_overall_mean_quantile = base::mean(anno_tg$gene_overall_mean_quantile)
        gene_factor = gene_quality_factor_min + (1 - gene_quality_factor_min) * (1 - (
          gene_nonzero_freq_quantile * weight_gene_nonzero_freq +
            gene_nonzero_mean_quantile * weight_gene_nonzero_mean +
            gene_overall_mean_quantile * weight_gene_overall_mean
        ))

        # Apply per-metric scaling factors to derive three independent thresholds.
        anno_tg$min_peak_nonzero_freq_quantile = anno_tg$min_quantile * gene_factor * min_peak_nonzero_freq_quantile_scale
        anno_tg$min_peak_nonzero_mean_quantile = anno_tg$min_quantile * gene_factor * min_peak_nonzero_mean_scale
        anno_tg$min_peak_overall_mean_quantile = anno_tg$min_quantile * gene_factor * min_peak_overall_mean_quantile_scale
        # Clamp all thresholds to valid [0, 1] range.
        anno_tg$min_peak_nonzero_freq_quantile[anno_tg$min_peak_nonzero_freq_quantile > 1] = 1
        anno_tg$min_peak_nonzero_freq_quantile[anno_tg$min_peak_nonzero_freq_quantile < 0] = 0
        anno_tg$min_peak_nonzero_mean_quantile[anno_tg$min_peak_nonzero_mean_quantile > 1] = 1
        anno_tg$min_peak_nonzero_mean_quantile[anno_tg$min_peak_nonzero_mean_quantile < 0] = 0
        anno_tg$min_peak_overall_mean_quantile[anno_tg$min_peak_overall_mean_quantile > 1] = 1
        anno_tg$min_peak_overall_mean_quantile[anno_tg$min_peak_overall_mean_quantile < 0] = 0

        # Keep peaks where all 3 metrics meet their respective thresholds.
        anno_tg$keep = base::with(
          anno_tg,
          peak_nonzero_freq_quantile >= min_peak_nonzero_freq_quantile &
            peak_nonzero_mean_quantile >= min_peak_nonzero_mean_quantile &
            peak_overall_mean_quantile >= min_peak_overall_mean_quantile
        )

        # Fallback: if no peak passes, retain the one with highest overall mean accessibility.
        if (base::sum(anno_tg$keep) == 0) {
          anno_tg[anno_tg$peak_overall_mean_quantile == base::max(anno_tg$peak_overall_mean_quantile), "keep"] = TRUE
        }

        # Append retained peaks.
        peaks_tg_keep = base::c(peaks_tg_keep, anno_tg[anno_tg$keep, ]$peak)
      } else {
        # Single-peak gene: keep unconditionally.
        peaks_tg_keep = base::c(peaks_tg_keep, anno_tg$peak)
      }
    } else {
      # --- Mode 2: strict filtering without fallback ---
      # Clamp distance to anchor range to avoid extrapolation artifacts.
      anno_tg$fixed_absDistanceToTSS = base::pmin(
        base::pmax(anno_tg$absDistanceToTSS, base::min(distance_points) + 0.01),
        base::max(distance_points) - 0.01
      )

      # Compute base distance-dependent quantile threshold.
      anno_tg$min_quantile = -1 / (anno_tg$fixed_absDistanceToTSS * a + b) + c

      # Compute composite gene quality factor.
      gene_nonzero_freq_quantile = base::mean(anno_tg$gene_nonzero_freq_quantile)
      gene_nonzero_mean_quantile = base::mean(anno_tg$gene_nonzero_mean_quantile)
      gene_overall_mean_quantile = base::mean(anno_tg$gene_overall_mean_quantile)
      gene_factor = gene_quality_factor_min + (1 - gene_quality_factor_min) * (1 - (
        gene_nonzero_freq_quantile * weight_gene_nonzero_freq +
          gene_nonzero_mean_quantile * weight_gene_nonzero_mean +
          gene_overall_mean_quantile * weight_gene_overall_mean
      ))

      # Apply per-metric scaling factors and clamp to [0, 1].
      anno_tg$min_peak_nonzero_freq_quantile = anno_tg$min_quantile * gene_factor * min_peak_nonzero_freq_quantile_scale
      anno_tg$min_peak_nonzero_mean_quantile = anno_tg$min_quantile * gene_factor * min_peak_nonzero_mean_scale
      anno_tg$min_peak_overall_mean_quantile = anno_tg$min_quantile * gene_factor * min_peak_overall_mean_quantile_scale
      anno_tg$min_peak_nonzero_freq_quantile[anno_tg$min_peak_nonzero_freq_quantile > 1] = 1
      anno_tg$min_peak_nonzero_freq_quantile[anno_tg$min_peak_nonzero_freq_quantile < 0] = 0
      anno_tg$min_peak_nonzero_mean_quantile[anno_tg$min_peak_nonzero_mean_quantile > 1] = 1
      anno_tg$min_peak_nonzero_mean_quantile[anno_tg$min_peak_nonzero_mean_quantile < 0] = 0
      anno_tg$min_peak_overall_mean_quantile[anno_tg$min_peak_overall_mean_quantile > 1] = 1
      anno_tg$min_peak_overall_mean_quantile[anno_tg$min_peak_overall_mean_quantile < 0] = 0

      # Keep peaks where all 3 metrics exceed their respective thresholds.
      anno_tg$keep = base::with(
        anno_tg,
        peak_nonzero_freq_quantile >= min_peak_nonzero_freq_quantile &
          peak_nonzero_mean_quantile >= min_peak_nonzero_mean_quantile &
          peak_overall_mean_quantile >= min_peak_overall_mean_quantile
      )

      # Append retained peaks (gene is dropped entirely if no peak qualifies).
      if (base::sum(anno_tg$keep) > 0) {
        peaks_tg_keep = base::c(peaks_tg_keep, anno_tg[anno_tg$keep, ]$peak)
      }
    }
  }

  ### Assemble final output.
  message("Backup original full annotation and replace with filtered set.")
  results_identify_TGs[["peak_anno_promoter_all"]] = results_identify_TGs[["peak_anno_promoter"]]
  results_identify_TGs[["peak_anno_promoter"]] = anno_all[anno_all$peak %in% peaks_tg_keep, ]
  results_identify_TGs[["TGs"]] = unique(results_identify_TGs[["peak_anno_promoter"]][, "SYMBOL"])

  # Record all parameter values for reproducibility.
  results_identify_TGs[["params"]] = base::list()
  results_identify_TGs[["params"]][["forced_feature_calculation"]] = forced_feature_calculation
  results_identify_TGs[["params"]][["distance_points"]] = distance_points
  results_identify_TGs[["params"]][["min_quantile_points"]] = min_quantile_points
  results_identify_TGs[["params"]][["one_peak_at_least_per_gene"]] = one_peak_at_least_per_gene
  results_identify_TGs[["params"]][["gene_quality_factor_min"]] = gene_quality_factor_min
  results_identify_TGs[["params"]][["weight_gene_nonzero_freq"]] = weight_gene_nonzero_freq
  results_identify_TGs[["params"]][["weight_gene_nonzero_mean"]] = weight_gene_nonzero_mean
  results_identify_TGs[["params"]][["weight_gene_overall_mean"]] = weight_gene_overall_mean
  results_identify_TGs[["params"]][["min_peak_nonzero_freq_quantile_scale"]] = min_peak_nonzero_freq_quantile_scale
  results_identify_TGs[["params"]][["min_peak_nonzero_mean_scale"]] = min_peak_nonzero_mean_scale
  results_identify_TGs[["params"]][["min_peak_overall_mean_quantile_scale"]] = min_peak_overall_mean_quantile_scale

  ### End.
  t_end = base::Sys.time()
  message("Finish: Filtering peak annotation results. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(results_identify_TGs)
}
