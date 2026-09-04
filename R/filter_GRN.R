#' @title Filter GRN Results with Dual Stringency Levels
#'
#' @description
#' Filters inferred gene regulatory network results at two stringency levels using dual thresholds combining quantile statistics and absolute bounds.
#' It retains high-confidence regulatory interactions based on fitting loss and predicted regulatory strength, and outputs both loose and strict filtered sets.
#'
#' @param interest_cell_type_pGRN List. Output from \code{\link{select_best_prediction}}.
#' @param fit_loss_quantile_threshold Numeric between 0 and 1. Quantile threshold for fitting loss. The final cutoff is the minimum of this quantile and \code{fit_loss_upper}. Lower values mean stricter filtering. Defaults to 0.2.
#' @param theta_p_quantile_threshold Numeric between 0 and 1. Quantile threshold for predicted regulatory strength. The final cutoff is the maximum of this quantile and \code{theta_p_lower}. Higher values mean stricter filtering. Defaults to 0.8.
#' @param fit_loss_upper Non-negative numeric. Absolute upper bound of acceptable fitting loss, serving as a safety cutoff to avoid overly lenient thresholds. Defaults to 0.2.
#' @param theta_p_lower Non-negative numeric. Absolute lower bound of acceptable regulatory strength, serving as a safety cutoff to avoid overly lenient thresholds. Defaults to 0.2.
#'
#' @return
#' The input nested list \code{interest_cell_type_pGRN} with two additional data frames per cell type:
#' \itemize{
#'   \item \code{pred_summary_loose}: Interactions passing the dual-threshold filtering (loss below cutoff and strength above cutoff).
#'   \item \code{pred_summary_strict}: High-confidence subset of the loose set, with all rows containing missing values removed.
#' }
#'
#' @details
#' The function performs filtering independently on the \code{pred_summary} data frame of each cell type:
#' \enumerate{
#'   \item \strong{Loose filtering}: Retains interactions satisfying both conditions with strict inequality:
#'     \itemize{
#'       \item Fitting loss is lower than the minimum of \code{fit_loss_upper} and the corresponding quantile value.
#'       \item Predicted regulatory strength is higher than the maximum of \code{theta_p_lower} and the corresponding quantile value.
#'     }
#'   \item \strong{Strict filtering}: Derived from the loose set by removing all rows with any missing values, yielding a complete-case high-confidence interaction set.
#' }
#' The dual-threshold design balances adaptability to different data distributions and guaranteed minimum quality standards.
#'
#' @note
#' \itemize{
#'   \item The function creates and temporarily switches to the \code{5 Infer GRN} directory, and restores the original working directory on exit.
#'   \item The loose set balances sensitivity and specificity, suitable for exploratory analysis and global network construction.
#'   \item The strict set provides conservative high-confidence interactions, recommended for functional validation and core regulatory network analysis.
#'   \item You can adjust filtering stringency via quantile parameters: lower \code{fit_loss_quantile_threshold} or higher \code{theta_p_quantile_threshold} increases stringency.
#'   \item The function appends new elements to the input list; the original \code{pred_summary} data remains fully preserved.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' fit_loss_quantile_threshold2 = 0.2
#' theta_p_quantile_threshold2 = 0.8
#' fit_loss_upper2 = 0.2
#' theta_p_lower2 = 0.2
#' interest_cell_type_pGRN = filter_GRN(
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   fit_loss_quantile_threshold = fit_loss_quantile_threshold2,
#'   theta_p_quantile_threshold = theta_p_quantile_threshold2,
#'   fit_loss_upper = fit_loss_upper2,
#'   theta_p_lower = theta_p_lower2
#' )
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
filter_GRN = function(
    interest_cell_type_pGRN,
    fit_loss_quantile_threshold = 0.2,
    theta_p_quantile_threshold = 0.8,
    fit_loss_upper = 0.2,
    theta_p_lower = 0.2
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
    message("Filtering GRN of cell type ", cell_type, ".")

    # Step 1: Loose filtering with dual thresholds.
    message("Generating loose-level GRN of cell type ", cell_type, ".")
    pred_summary_loose = interest_cell_type_pGRN[[cell_type]][["pred_summary"]]
    pred_summary_loose = pred_summary_loose[
      pred_summary_loose$fit_loss < base::min(fit_loss_upper, stats::quantile(pred_summary_loose$fit_loss, probs = fit_loss_quantile_threshold, na.rm = TRUE)) &
        pred_summary_loose$theta_p > base::max(theta_p_lower, stats::quantile(pred_summary_loose$theta_p, probs = theta_p_quantile_threshold, na.rm = TRUE)),
    ]
    interest_cell_type_pGRN[[cell_type]][["pred_summary_loose"]] = pred_summary_loose

    # Step 2: Strict filtering: complete-case filtering based on loose set.
    message("Generating strict-level GRN of cell type ", cell_type, ".")
    interest_cell_type_pGRN[[cell_type]][["pred_summary_strict"]] = stats::na.omit(pred_summary_loose)
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Filtering GRN. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_pGRN)
}
