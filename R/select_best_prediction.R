#' @title Select Best GRN Predictions with Multi-Level Filtering
#'
#' @description
#' Post-processes inferred gene regulatory network results to select high-confidence regulatory predictions.
#' For each TF-TG pair per branch, it selects the optimal prediction from multiple trained models via a four-level hierarchical filtering strategy, combining quantile thresholds and absolute bounds.
#' It further integrates branch-level results into cell-type-level summary using cell-count weighted averaging.
#'
#' @param interest_cell_type_pGRN List. Output from \code{\link{infer_GRN}}.
#' @param interest_cell_type_tau Named numeric vector. Output from \code{\link{get_cell_type_tau}}, cell-type-specific binarization thresholds.
#' @param interest_cell_type_genes_pseudotime_info List. Output from \code{\link{get_genes_pseudotime_info}}, containing cell counts per branch for weighted integration.
#' @param fit_loss_quantile_threshold Numeric between 0 and 1. Quantile threshold for filtering predictions based on loss. Defaults to 0.2.
#' @param theta_p_quantile_threshold Numeric between 0 and 1. Quantile threshold for regulatory strength filtering. Defaults to 0.8.
#' @param fit_loss_upper Positive numeric. Absolute upper bound of acceptable fitting loss. Defaults to 999 (no filtering by default). Intended for development and debugging only; quantile thresholds are recommended for routine use.
#' @param theta_p_lower Numeric between 0 and 1. Absolute lower bound of acceptable regulatory strength. Defaults to 0 (no filtering by default). Intended for development and debugging only.
#' @param each Positive integer. Maximum number of TF-TG pairs per chunk for parallel processing. Defaults to 1,000,000.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#'
#' @return
#' The input list \code{interest_cell_type_pGRN} with added elements, organized by cell type:
#' \itemize{
#'   \item For each branch: a \code{best_result} data frame is added, containing optimal regulatory strength (\code{theta_p}), corresponding fitting loss, and the selected model name for each TF-TG pair.
#'   \item A cell-type-level \code{pred_summary} data frame, which integrates results across all branches, including branch-specific values and cell-count-weighted averages of \code{theta_p} and \code{fit_loss}.
#' }
#'
#' @details
#' The function executes the following workflow independently for each cell type:
#' \enumerate{
#'   \item \strong{Threshold calculation}: Computes branch-specific filtering thresholds by fusing quantile statistics and absolute bounds, generating two loss thresholds (tight and standard) and one minimum regulatory strength threshold.
#'   \item \strong{Data partitioning}: Splits TF-TG pairs into chunks to control memory usage and enable parallel processing.
#'   \item \strong{Four-level hierarchical selection}: For each TF-TG pair, selects the optimal prediction in priority order with fallback logic to ensure complete coverage:
#'     \enumerate{
#'       \item Retain predictions below the tight loss threshold, select the one with maximum \code{theta_p}.
#'       \item Fall back to the standard loss threshold if no result meets the first criterion, select the one with maximum \code{theta_p}.
#'       \item Fall back to predictions above the minimum strength threshold, select the one with minimum loss.
#'       \item Fall back to the prediction with minimum loss directly if all thresholds fail.
#'     }
#'   \item \strong{Post-filtering}: Removes invalid boundary values (\code{theta_p} < 0.011) and infinite loss values.
#'   \item \strong{Cross-branch integration}: Aggregates branch-level optimal results, and computes cell-count-weighted averages to generate the cell-type-level final GRN.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and temporarily switches to the \code{5 Infer GRN} directory, and restores the original working directory on exit.
#'   \item The \code{each} parameter controls peak memory usage; reduce its value if out-of-memory errors occur.
#'   \item The function appends new elements to the input list; the original input object is not modified in place.
#'   \item A fixed threshold of 0.011 is applied to filter out boundary-invalid \code{theta_p} values after selection.
#'   \item Parallel computation relies on \code{\link[parallel]{mclapply}}; set \code{ncores = 1} on Windows systems.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
#' interest_cell_type_genes_pseudotime_info = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_pGRN = select_best_prediction(
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   interest_cell_type_tau = interest_cell_type_tau,
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   fit_loss_quantile_threshold = 0.2,
#'   theta_p_quantile_threshold = 0.8,
#'   fit_loss_upper = 999,
#'   theta_p_lower = 0,
#'   each = 1000,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
select_best_prediction = function(
    interest_cell_type_pGRN,
    interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info,
    fit_loss_quantile_threshold = 0.2,
    theta_p_quantile_threshold = 0.8,
    fit_loss_upper = 999,
    theta_p_lower = 0,
    each = 1000000,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Filtering GRN. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "5 Infer GRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Process.
  for (cell_type in base::names(interest_cell_type_pGRN)) {
    # Initialize cell-type-level integrated result container.
    interest_cell_type_pGRN[[cell_type]][["pred_summary"]] = list()

    # Process each branch independently.
    for (n in base::seq_along(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]])) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      re_pred = base::rownames(interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["theta_p"]])
      N_pred_all = base::length(re_pred)
      message("There are ", N_pred_all, " potential regulatory relationships in branch ", n, " of cell type ", cell_type, ".")

      # Split TF-TG pairs into chunks for parallel processing.
      message("Partitioning data.")
      n_part = base::ceiling(N_pred_all / each)
      # each_last = N_pred_all - each * (n_part - 1)

      re_pred_list = base::list()
      for (s in 1:n_part) {
        if (s != n_part) {
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):(each * s)]
        } else (
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):N_pred_all]
        )
      }; gc()

      # Initialize result containers for each chunk.
      pGRN_list = base::list()
      for (s in 1:n_part) {
        pGRN_list[[s]] = base::list()
        pGRN_list[[s]][["theta_p"]] = interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["theta_p"]][re_pred_list[[s]], ]
        pGRN_list[[s]][["loss"]] = interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["loss"]][re_pred_list[[s]], ]
        pGRN_list[[s]][["best_result"]] = base::data.frame(
          theta_p = base::rep(0, base::length(re_pred_list[[s]])),
          loss = base::rep(Inf, base::length(re_pred_list[[s]])),
          fit = base::rep(NA, base::length(re_pred_list[[s]])),
          row.names = re_pred_list[[s]]
        )
      }; gc()

      # Calculate branch-level dual loss thresholds and theta minimum.
      message("Selecting best theta for each TF-TG of branch ", n, " of cell type ", cell_type, ".")
      # Priority 1 loss threshold: tighter quantile (quantile / 2)
      loss_max_s1 = base::min(
        stats::quantile(
          base::unlist(interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]][["loss"]]),
          probs = fit_loss_quantile_threshold / 2,
          na.rm = TRUE
        ),
        fit_loss_upper
      ); gc()
      # Priority 2 loss threshold: standard quantile
      loss_max_s2 = base::min(
        stats::quantile(
          base::unlist(interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]][["loss"]]),
          probs = fit_loss_quantile_threshold,
          na.rm = TRUE
        ),
        fit_loss_upper
      ); gc()
      # Minimum valid theta: max of quantile, tau, and absolute lower bound
      theta_p_min = base::max(
        stats::quantile(
          base::unlist(interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]][["theta_p"]]),
          probs = theta_p_quantile_threshold,
          na.rm = TRUE
        ),
        interest_cell_type_tau[cell_type],
        theta_p_lower
      ); gc()

      # Core worker function: select best prediction for each pair in a chunk.
      process_data_slice = function(pGRN = pGRN_list[[s]],
                                    loss_max_s1 = loss_max_s1,
                                    loss_max_s2 = loss_max_s2,
                                    theta_p_min = theta_p_min) {
        for (re in base::rownames(pGRN[[1]])) {
          theta_p = pGRN[["theta_p"]][re, , drop = FALSE]
          loss = pGRN[["loss"]][re, , drop = FALSE]
          best_result = base::data.frame(base::matrix(NA, nrow = 1, ncol = 3, dimnames = base::list(re, c("theta_p", "loss", "fit"))))
          best_result[re, 1:2] = c(0, Inf)

          # 4-level hierarchical selection: fall back sequentially if higher priority has no valid results
          # Priority 1: loss below tight threshold, select max theta_p
          theta_p_selected1 = base::max(theta_p[, loss <= loss_max_s1, drop = FALSE])
          # Priority 2: loss below standard threshold, select max theta_p
          theta_p_selected2 = base::max(theta_p[, loss <= loss_max_s2, drop = FALSE])
          # Priority 3: theta_p above minimum threshold, select min loss
          theta_p_selected3 = theta_p[, theta_p >= theta_p_min, drop = FALSE]

          if (theta_p_selected1 != -Inf) {
            best_result[re, "theta_p"] = theta_p_selected1
            best_result[re, "fit"] = base::names(loss)[which(theta_p == theta_p_selected1)][1]
            best_result[re, "loss"] = loss[re, best_result[re, "fit"]]
          } else if (theta_p_selected2 != -Inf) {
            best_result[re, "theta_p"] = theta_p_selected2
            best_result[re, "fit"] = base::names(loss)[which(theta_p == theta_p_selected2)][1]
            best_result[re, "loss"] = loss[re, best_result[re, "fit"]]
          } else if (base::ncol(theta_p_selected3) > 0) {
            best_result[re, "loss"] = base::min(loss[base::names(theta_p_selected3)])
            best_result[re, "fit"] = base::names(loss)[base::which(loss == best_result[re, "loss"])][1]
            best_result[re, "theta_p"] = theta_p[re, best_result[re, "fit"]]
          } else {
            # Priority 4 (fallback): select min loss directly
            best_result[re, "loss"] = base::min(loss)
            best_result[re, "fit"] = base::names(loss)[base::which(loss == best_result[re, "loss"])][1]
            best_result[re, "theta_p"] = theta_p[re, best_result[re, "fit"]]
          }

          pGRN[[re]] = base::list(best_result = best_result)
        }

        # Merge pair-level results into a unified data frame.
        # The first 3 elements of pGRN are fixed: theta_p, loss, best_result; subsequent elements are per-pair results.
        pGRN[["best_result"]] = base::as.data.frame(data.table::rbindlist(base::lapply(
          pGRN[4:base::length(pGRN)], function(x) x[["best_result"]]
        )))
        base::rownames(pGRN[["best_result"]]) = base::names(pGRN[4:base::length(pGRN)])

        # Post-filter: remove invalid boundary values (theta_p < 0.011) and infinite loss
        pGRN[["best_result"]] = pGRN[["best_result"]][
          pGRN[["best_result"]]$theta_p >= 0.011 &
            # pGRN[["best_result"]]$theta_p <= 0.989 &
            !base::is.infinite(pGRN[["best_result"]]$loss), , drop = FALSE]

        return(pGRN[["best_result"]])
      }

      prediction_list = parallel::mclapply(
        X = 1:n_part,
        FUN = function(s) {process_data_slice(
          pGRN = pGRN_list[[s]],
          loss_max_s1 = loss_max_s1,
          loss_max_s2 = loss_max_s2,
          theta_p_min = theta_p_min
        )},
        mc.cores = ncores
      ); gc()

      # Merge chunk results into complete branch-level best_result.
      message("Summarizing best_result of this branch.")
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]] =
        base::do.call("rbind", prediction_list)
      base::colnames(interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]]) = base::paste0(
        base::colnames(interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]]), "_", "branch", n
      )
    }

    # Integrate results across all branches and compute weighted averages.
    message("Summarizing the predictive results of each branch.")
    # Collect all unique TF-TG pairs across branches.
    rows = base::sort(base::unique(base::Reduce(base::c, base::lapply(
      interest_cell_type_pGRN[[cell_type]][base::paste0("branch", base::seq_along(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]]))],
      function(x) base::rownames(x[["best_result"]])
    ))))
    # Collect all column names across branches.
    cols = base::Reduce(base::c, base::lapply(
      interest_cell_type_pGRN[[cell_type]][base::paste0("branch", base::seq_along(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]]))],
      function(x) base::colnames(x[["best_result"]])
    ))
    # Initialize integrated summary table.
    pred_summary = base::data.frame(base::matrix(NA, nrow = base::length(rows), ncol = base::length(cols), dimnames = base::list(rows, cols)))
    # Fill in branch-specific results.
    for (n in base::seq_along(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]])) {
      pred_summary[
        base::rownames(interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]]),
        base::colnames(interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]])
      ] = interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]]
    }

    n_branch = base::length(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]])
    if (n_branch > 1) {
      # Weighted by cell count of each branch.
      weights = base::unlist(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]])
      pred_summary$theta_p = base::apply(pred_summary[, base::grep("^theta_p_branch([0-9]{1,2})$", base::names(pred_summary), value = TRUE)], 1, function(x) {
        non_zero_idx = which(!base::is.na(x))
        if (base::length(non_zero_idx) == 0) 0 else {
          base::sum(x[non_zero_idx] * weights[non_zero_idx]) / base::sum(weights[non_zero_idx])
        }
      })
      pred_summary$fit_loss = base::apply(pred_summary[, base::grep("^loss_branch([0-9]{1,2})$", base::names(pred_summary), value = TRUE)], 1, function(x) {
        non_zero_idx = which(!base::is.na(x))
        if (base::length(non_zero_idx) == 0) Inf else {
          base::sum(x[non_zero_idx] * weights[non_zero_idx]) / base::sum(weights[non_zero_idx])
        }
      })
    } else if (n_branch == 1) {
      pred_summary$theta_p = pred_summary$theta_p_branch1
      pred_summary$fit_loss = pred_summary$loss_branch1
    } else {stop("The number of branches is less than 1.")}

    interest_cell_type_pGRN[[cell_type]][["pred_summary"]] = pred_summary
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Filtering GRN. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_pGRN)
}
