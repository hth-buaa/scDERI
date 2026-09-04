#' @title Filter Enhancer-Driven Gene Regulatory Networks via Multi-Factor Adjusted Quantile Thresholds
#'
#' @description
#' Filters cell-type-specific predictive enhancer-driven gene regulatory networks (peGRNs) using an adjusted quantile threshold strategy.
#' A base threshold curve is fitted as an inverse function of predicted TF-TG regulatory strength (\code{FG-theta_p}),
#' then adjusted by 11 feature-specific factors (penalty or relaxation). Two sequential rescue mechanisms with quality floors
#' ensure minimum network connectivity: minimum retained TFs per peak-TG pair, and minimum retained peak-TG pairs per TF.
#'
#' @param interest_cell_type_peGRN Output list from \code{\link{infer_eGRN}}.
#' @param theta_p_points Numeric vector. Three anchor points of \code{FG-theta_p} (range 0–1) for threshold curve fitting.
#'   Must be sorted ascending and match the length of \code{initial_quantile_points}. Default: \code{c(0, 0.5, 1)}.
#' @param initial_quantile_points Numeric vector. Corresponding initial quantile thresholds at each \code{theta_p} anchor.
#'   Must match the length of \code{theta_p_points}. Default: \code{c(1.5, 0.9, 0.6)}.
#' @param theta_p_quantile_fit_start Named list with initial values \code{a}, \code{b}, \code{c} for inverse curve fitting
#'   via \code{\link[minpack.lm]{nlsLM}}. Parameter \code{a} is constrained to be positive. Defaults to pre-estimated values.
#' @param rescue_peakTG Logical. Whether to perform peak-TG pair rescue to guarantee minimum retained TFs per pair.
#'   Default: \code{TRUE}.
#' @param peakTG_rescue_quantile_threshold Numeric in \[0, 1\]. Quality floor quantile of \code{FP-S-R} across peak-TG rescue candidates.
#'   Only candidates above this threshold are eligible for rescue (0 for unconditional rescue; 1 for effectively disable). Defaults to \code{0.1}.
#' @param min_TF_per_peakTG Positive integer. Minimum number of TFs to retain per peak-TG pair after rescue.
#'   Only effective when \code{rescue_peakTG = TRUE}. Default: \code{1}.
#' @param rescue_TF Logical. Whether to perform TF rescue to guarantee minimum retained peak-TG pairs per TF.
#'   Default: \code{TRUE}.
#' @param TF_rescue_quantile_threshold Numeric in \[0, 1\].
#' Quality floor quantile of \code{FP-S-R} across TF rescue candidates.
#'   Only candidates above this threshold are eligible for rescue (0 for unconditional rescue; 1 for effectively disable). Defaults to \code{0}.
#' @param min_peakTG_per_TF Positive integer. Minimum number of peak-TG pairs to retain per TF after rescue.
#'   Only effective when \code{rescue_TF = TRUE}. Default: \code{2}.
#' @param factor_balance Named list. Balance center values for each feature-specific adjustment factor, controlling the neutral position of adjustment.
#'   Names must match \code{factor_range} and \code{factor_attribute}. Defaults to 0.2 for all 11 features.
#' @param factor_range Named list. Adjustment amplitude for each feature. Larger values produce stronger penalty/relaxation effects.
#'   Defaults are feature-specific.
#' @param factor_attribute Named list. Adjustment direction for each feature:
#'   \code{"penalty"} raises the threshold (stricter filtering); \code{"relaxation"} lowers the threshold (more lenient).
#'
#' @returns
#' The input \code{interest_cell_type_peGRN} list with two additional elements per stringency level (\code{loose}, \code{strict}) within each cell type:
#' \itemize{
#'   \item \code{<level>_filtered}: Full edge table with all original columns plus computed fields:
#'     \code{peak-TG}, \code{IQ}, 11 feature adjustment factors (with \code{-F} suffix), \code{AQ},
#'     three mechanism-specific retention flags (\code{filter_keep}, \code{rescue_peakTG_keep}, \code{rescue_TF_keep}),
#'     and the final retention flag \code{keep}.
#'   \item \code{<level>_filtered_simplified}: Simplified three-column table (\code{TF}, \code{peak}, \code{TG}) for direct downstream network analysis.
#' }
#'
#' @details
#' The function executes in three main steps: curve fitting runs once globally, while filtering and rescue are applied independently per cell type and per stringency level.
#'
#' **1. Global inverse threshold curve fitting**
#' Fits the model \deqn{initial\_quantile = \dfrac{1}{a \cdot \theta_p + b} + c} with constraint \eqn{a > 0}
#' using the Levenberg–Marquardt algorithm (\code{minpack.lm::nlsLM}). Higher \code{FG-theta_p} corresponds to a lower base quantile threshold.
#'
#' **2. Multi-factor adjusted thresholding**
#' For each edge in each cell type / stringency level:
#' - Compute the initial quantile threshold (\code{IQ}) from the fitted curve based on the edge’s \code{FG-theta_p}.
#' - Calculate 11 feature-specific adjustment factors. Penalty factors increase the final threshold; relaxation factors decrease it.
#' - Derive the adjusted quantile threshold (\code{AQ}) by multiplying \code{IQ} with all adjustment factors, then clamp values to the valid range \[0, 1\].
#' - Primary filtering rule: retain edges where \code{FP-S-R ≥ AQ}.
#'
#' **3. Two-level sequential rescue**
#' All candidate edges are pre-ranked by priority (descending unless otherwise noted):
#' \code{FP-S-R} > \code{P-AA-R} > \code{P-NZAR-R} > \code{FG-theta_i-R} > \code{FG-theta_p} > \code{FG-NL} (ascending).
#'
#' - **Peak-TG rescue**: For each peak-TG pair with fewer retained edges than \code{min_TF_per_peakTG},
#'   restore top-ranked TF candidates up to the minimum requirement, subject to the quality floor threshold.
#' - **TF rescue**: For each TF with fewer retained edges than \code{min_peakTG_per_TF} (after peak-TG rescue),
#'   restore top-ranked peak-TG candidates up to the minimum requirement, subject to the quality floor threshold.
#'
#' The final retention status is the union of the primary filter, peak-TG rescue and TF rescue results.
#'
#' @note
#' \itemize{
#'   \item The threshold curve is fitted once globally and shared across all cell types and stringency levels.
#'   \item Adjusted thresholds (\code{AQ}) are clamped to \[0, 1\] to avoid invalid values from extreme factor combinations.
#'   \item Rescue steps are applied sequentially: TF rescue runs on the combined result of primary filtering and peak-TG rescue.
#'   \item Minimum connectivity parameters define the lower bound; actual retained edges may be higher if the primary filter already meets the requirement.
#'   \item Retention flags allow tracing which mechanism retained each individual edge.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_peGRN = readRDS("./5 Infer GRN/interest_cell_type_peGRN.rds")
#' interest_cell_type_peGRN = filter_eGRN(
#'   interest_cell_type_peGRN = interest_cell_type_peGRN,
#'   theta_p_points = c(0, 0.5, 1),
#'   initial_quantile_points = c(1.5, 0.9, 0.6),
#'   theta_p_quantile_fit_start = list(a = 0.333333, b = 0.333333, c = -0.5),
#'   rescue_peakTG = TRUE,
#'   peakTG_rescue_quantile_threshold = 0.1,
#'   min_TF_per_peakTG = 1,
#'   rescue_TF = TRUE,
#'   TF_rescue_quantile_threshold = 0.0,
#'   min_peakTG_per_TF = 2,
#'   factor_range = list(
#'     `FG-NL` = 0.25,
#'     `FG-theta_i-R` = 0.25,
#'     `F-NZER-R` = 0.05, `F-NZEA-R` = 0.05, `F-EA-R` = 0.05,
#'     `P-NZAR-R` = 0.25, `P-NZAA-R` = 0.25, `P-AA-R` = 0.25,
#'     `G-NZER-R` = 0.15, `G-NZEA-R` = 0.15, `G-EA-R` = 0.15
#'   )
#' )
#' saveRDS(interest_cell_type_peGRN, file = "./5 Infer GRN/interest_cell_type_peGRN.rds")
#' }
filter_eGRN = function(
    interest_cell_type_peGRN,
    theta_p_points = base::c(0, 0.5, 1), initial_quantile_points = base::c(1.5, 0.9, 0.6),
    theta_p_quantile_fit_start = base::list(a = 0.333333, b = 0.333333, c = -0.5),
    rescue_peakTG = TRUE, peakTG_rescue_quantile_threshold = 0.1, min_TF_per_peakTG = 1,
    rescue_TF = TRUE, TF_rescue_quantile_threshold = 0.0, min_peakTG_per_TF = 7,
    factor_balance = base::list(
      `FG-NL` = 0.2, `FG-theta_i-R` = 0.2,
      `F-NZER-R` = 0.2, `F-NZEA-R` = 0.2, `F-EA-R` = 0.2,
      `P-NZAR-R` = 0.2, `P-NZAA-R` = 0.2, `P-AA-R` = 0.2,
      `G-NZER-R` = 0.2, `G-NZEA-R` = 0.2, `G-EA-R` = 0.2
    ),
    factor_range = base::list(
      `FG-NL` = 0.25,
      `FG-theta_i-R` = 0.25,
      `F-NZER-R` = 0.05, `F-NZEA-R` = 0.05, `F-EA-R` = 0.05,
      `P-NZAR-R` = 0.25, `P-NZAA-R` = 0.25, `P-AA-R` = 0.25,
      `G-NZER-R` = 0.15, `G-NZEA-R` = 0.15, `G-EA-R` = 0.15
    ),
    factor_attribute = base::list(
      `FG-NL` = "penalty",
      `FG-theta_i-R` = "relaxation",
      `F-NZER-R` = "relaxation", `F-NZEA-R` = "relaxation", `F-EA-R` = "relaxation",
      `P-NZAR-R` = "relaxation", `P-NZAA-R` = "relaxation", `P-AA-R` = "relaxation",
      `G-NZER-R` = "relaxation", `G-NZEA-R` = "relaxation", `G-EA-R` = "relaxation"
    )
) {
  ### Input validation
  if (length(theta_p_points) != length(initial_quantile_points)) {
    stop("theta_p_points and initial_quantile_points must have the same length")
  }

  feature_names = c(
    "FG-NL", "FG-theta_i-R",
    "F-NZER-R", "F-NZEA-R", "F-EA-R",
    "P-NZAR-R", "P-NZAA-R", "P-AA-R",
    "G-NZER-R", "G-NZEA-R", "G-EA-R"
  )
  if (!all(feature_names %in% names(factor_balance)) ||
      !identical(sort(names(factor_balance)), sort(names(factor_range))) ||
      !identical(sort(names(factor_balance)), sort(names(factor_attribute)))) {
    stop("factor_balance, factor_range, factor_attribute must have identical feature names")
  }

  if (peakTG_rescue_quantile_threshold < 0 | peakTG_rescue_quantile_threshold > 1) {
    stop("peakTG_rescue_quantile_threshold must be between 0 and 1")
  }
  if (TF_rescue_quantile_threshold < 0 | TF_rescue_quantile_threshold > 1) {
    stop("TF_rescue_quantile_threshold must be between 0 and 1")
  }

  if (rescue_peakTG && (min_TF_per_peakTG < 1 || !is.numeric(min_TF_per_peakTG) || min_TF_per_peakTG != as.integer(min_TF_per_peakTG))) {
    stop("min_TF_per_peakTG must be a positive integer")
  }
  if (rescue_TF && (min_peakTG_per_TF < 1 || !is.numeric(min_peakTG_per_TF) || min_peakTG_per_TF != as.integer(min_peakTG_per_TF))) {
    stop("min_peakTG_per_TF must be a positive integer")
  }

  ### Initialize timing
  t_start = base::Sys.time()
  message("[Start] Filtering peGRN network at ", t_start)

  ### Fit inverse quantile threshold curve
  # Model form: initial_quantile = 1 / (theta_p * a + b) + c, with a > 0
  message("Fitting inverse threshold model: initial_quantile ~ 1/(theta_p * a + b) + c")
  fit = minpack.lm::nlsLM(
    initial_quantile_points ~ 1 / (theta_p_points * a + b) + c,
    start = theta_p_quantile_fit_start,
    # algorithm = "port",
    lower = c(1e-9, -Inf, -Inf),
    upper = c(Inf, Inf, Inf)
  )

  # Extract fitted model coefficients
  message("Extracting fitted model coefficients")
  a = base::as.numeric(stats::coef(fit)["a"])
  b = base::as.numeric(stats::coef(fit)["b"])
  c = base::as.numeric(stats::coef(fit)["c"]); rm(fit)
  message("Fitted params: a=", base::round(a, 6), ", b=", base::round(b, 6), ", c=", base::round(c, 6))

  ### Filter eGRN per cell type and stringency level
  for (cell_type in base::names(interest_cell_type_peGRN)) {
    for (level in c("loose", "strict")) {
      eGRN = interest_cell_type_peGRN[[cell_type]][[level]]
      eGRN[, "peak-TG"] = base::paste0(eGRN$peak, "_to_", eGRN$TG)
      base::rownames(eGRN) = 1:base::nrow(eGRN)

      peakTG_all = base::unique(eGRN[, "peak-TG"])
      TF_all = base::unique(eGRN[, "TF"])

      # Calculate initial quantile (IQ) threshold
      message("[", cell_type, " - ", level, "] Calculating initial quantile (IQ) threshold")
      eGRN[, "IQ"] = 1 / (eGRN[, "FG-theta_p"] * a + b) + c

      # Compute feature-specific adjustment factors
      # Formula:
      #   relaxation (fa=1): factor = 1 + fr * (fb - feature_value)
      #   penalty    (fa=-1): factor = 1 + fr * (feature_value - fb)
      # message("[", cell_type, " - ", level, "] Computing ", length(factor_features), " feature-specific adjustment factors")
      factor_features = base::names(factor_balance)
      for (ff in factor_features) {
        message("[", cell_type, " - ", level, "] Computing ", factor_attribute[[ff]], " factor for ", ff)
        fb = factor_balance[[ff]]
        fr = factor_range[[ff]]
        fa = (-1)^base::as.numeric(factor_attribute[[ff]] == "penalty")
        eGRN[, base::paste0(ff, "-F")] = (1 + fa * fb * fr) - fa * eGRN[, ff] * fr
      }

      # Calculate adjusted quantile (AQ) via vectorized row product (optimized)
      message("[", cell_type, " - ", level, "] Calculating adjusted quantile (AQ) threshold")
      aq_cols = c("IQ", base::paste0(factor_features, "-F"))
      # Clamp values to positive range to avoid log(NaN)
      aq_matrix = base::as.matrix(eGRN[, aq_cols])
      aq_matrix[aq_matrix <= 0] = .Machine$double.eps
      eGRN[, "AQ"] = base::exp(base::rowSums(base::log(aq_matrix)))
      # Clamp AQ to valid [0, 1] range
      eGRN$AQ = base::pmin(base::pmax(eGRN$AQ, 0), 1)

      # Primary filtering by AQ threshold
      message("[", cell_type, " - ", level, "] Performing primary threshold filtering")
      eGRN$filter_keep = eGRN[, "FP-S-R"] >= eGRN[, "AQ"]

      # Initialize rescue flags
      eGRN$rescue_peakTG_keep = FALSE
      eGRN$rescue_TF_keep = FALSE
      eGRN$keep = eGRN$filter_keep

      # Pre-sort once by priority; reuse for both rescue steps
      ord = base::order(
        -eGRN$`FP-S-R`, -eGRN$`P-AA-R`, -eGRN$`P-NZAR-R`,
        -eGRN$`FG-theta_i-R`, -eGRN$`FG-theta_p`, eGRN$`FG-NL`
      )
      eGRN_ordered = eGRN[ord, ]

      ### Rescue step 1: peak-TG pair rescue
      if (rescue_peakTG) {
        message("[", cell_type, " - ", level, "] Rescuing peak-TG pairs (min ", min_TF_per_peakTG, " TF per pair)")

        # Rank entries within each peak-TG group (follows pre-sorted priority order)
        eGRN_ordered$peakTG_rank = stats::ave(
          base::seq_len(base::nrow(eGRN_ordered)),
          eGRN_ordered$`peak-TG`,
          FUN = base::seq_along
        )

        # Select top N candidates per peak-TG
        peakTG_candidates = eGRN_ordered[eGRN_ordered$peakTG_rank <= min_TF_per_peakTG, ]

        # Quality filter by quantile threshold
        q_thresh = stats::quantile(peakTG_candidates$`FP-S-R`, peakTG_rescue_quantile_threshold)
        peakTG_candidates = peakTG_candidates[peakTG_candidates$`FP-S-R` >= q_thresh, ]

        # Count currently kept entries per peak-TG
        keep_count = base::tapply(eGRN$keep, eGRN$`peak-TG`, sum)
        candidate_keep_count = keep_count[peakTG_candidates$`peak-TG`]

        # Vectorized rescue decision: fill connectivity gap up to minimum requirement
        need_rescue = peakTG_candidates$peakTG_rank <= (min_TF_per_peakTG - candidate_keep_count) &
          candidate_keep_count < min_TF_per_peakTG

        # Mark rescued rows by rowname for precise indexing
        rescue_rows = base::rownames(peakTG_candidates)[need_rescue]
        eGRN[rescue_rows, "rescue_peakTG_keep"] = TRUE

        rm(peakTG_candidates, q_thresh, keep_count, candidate_keep_count, need_rescue, rescue_rows)
      }

      # Update keep status before TF rescue
      eGRN$keep = eGRN$filter_keep | eGRN$rescue_peakTG_keep

      ### Rescue step 2: TF rescue
      if (rescue_TF) {
        message("[", cell_type, " - ", level, "] Rescuing TFs (min ", min_peakTG_per_TF, " peak-TG per TF)")

        # Rank entries within each TF group (reuses pre-sorted order)
        eGRN_ordered$TF_rank = stats::ave(
          base::seq_len(base::nrow(eGRN_ordered)),
          eGRN_ordered$TF,
          FUN = base::seq_along
        )

        # Select top N candidates per TF
        TF_candidates = eGRN_ordered[eGRN_ordered$TF_rank <= min_peakTG_per_TF, ]

        # Quality filter by quantile threshold
        q_thresh = stats::quantile(TF_candidates$`FP-S-R`, TF_rescue_quantile_threshold)
        TF_candidates = TF_candidates[TF_candidates$`FP-S-R` >= q_thresh, ]

        # Count currently kept entries per TF (after peak-TG rescue)
        keep_count = base::tapply(eGRN$keep, eGRN$TF, sum)
        candidate_keep_count = keep_count[TF_candidates$TF]

        # Vectorized rescue decision
        need_rescue = TF_candidates$TF_rank <= (min_peakTG_per_TF - candidate_keep_count) &
          candidate_keep_count < min_peakTG_per_TF

        # Mark rescued rows by rowname
        rescue_rows = base::rownames(TF_candidates)[need_rescue]
        eGRN[rescue_rows, "rescue_TF_keep"] = TRUE

        rm(TF_candidates, q_thresh, keep_count, candidate_keep_count, need_rescue, rescue_rows)
      }

      rm(ord, eGRN_ordered)

      # Assemble final filtered outputs
      eGRN$keep = eGRN$filter_keep | eGRN$rescue_peakTG_keep | eGRN$rescue_TF_keep
      filtered_eGRN = eGRN[eGRN$keep, ]
      simplified_filtered_eGRN = filtered_eGRN[, c("TF", "peak", "TG")]

      interest_cell_type_peGRN[[cell_type]][[base::paste0(level, "_filtered")]] = filtered_eGRN
      interest_cell_type_peGRN[[cell_type]][[base::paste0(level, "_filtered_simplified")]] = simplified_filtered_eGRN
    }
  }

  ### Finalize and report runtime
  t_end = base::Sys.time()
  message("[Finish] Filtering complete at ", t_end)
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_peGRN)
}
