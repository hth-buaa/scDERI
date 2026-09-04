#' @title Compute Cell-Type-Specific GRN Binarization Thresholds
#'
#' @description Computes a cell-type-specific binarization threshold (\eqn{\tau}) for each initial gene regulatory network.
#' The threshold is derived from the minimum non-zero interaction strength, with a fixed lower bound to avoid over-sensitivity.
#'
#' @param interest_cell_type_iGRN Named list. Output from \code{\link{get_iGRN_by_TFBS_pwm_by_JASPAR2024}}
#'   or \code{\link{filter_iGRN_by_ChIP_Atlas_priors}}, with each element being a TG × TF association matrix.
#'
#' @returns A named numeric vector.
#' Names correspond to cell types, and values are the computed binarization thresholds (\eqn{\tau}).
#'
#' @details
#' The threshold is calculated independently for each cell type following the rule:
#' \deqn{\tau = \max(\min(\theta > 0) - 0.005,\ 0.005)}
#' where \eqn{\theta} represents all non-zero interaction strengths in the iGRN matrix.
#'
#' Specifically:
#' \enumerate{
#'   \item Extracts all non-zero interaction strength values from the matrix (missing values automatically ignored).
#'   \item Takes the minimum non-zero value and subtracts a small offset of 0.005.
#'   \item Enforces a lower bound of 0.005 to prevent excessively small thresholds.
#' }
#'
#' This threshold is typically used to convert continuous regulatory strengths into binary interaction calls
#' (presence / absence of regulatory relationships) for downstream discrete modeling.
#'
#' @note
#' \itemize{
#'   \item Thresholds are cell-type-specific, adapting to the distribution of interaction strengths in each network.
#'   \item The 0.005 lower bound prevents overly sensitive thresholds when minimum non-zero values are extremely small.
#'   \item The function creates (if missing) and temporarily switches to the \code{3 get iGRN/} directory during execution,
#'     and restores the original working directory on completion.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
#' interest_cell_type_tau = get_cell_type_tau(interest_cell_type_iGRN)
#' base::saveRDS(interest_cell_type_tau, file = "./3 get iGRN/interest_cell_type_tau.rds")
#' }
get_cell_type_tau = function(interest_cell_type_iGRN)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Extracting tau for each cell type. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Process.
  message("Initializing threshold output vector.")
  interest_cell_type_tau = base::rep(999, base::length(interest_cell_type_iGRN))
  base::names(interest_cell_type_tau) = base::names(interest_cell_type_iGRN)

  message("Calculating binarization thresholds for each cell type.")
  # Compute tau: subtract 0.005 from minimum non-zero theta, enforce 0.005 lower bound
  for (cell_type in base::names(interest_cell_type_iGRN)) {
    interest_cell_type_tau[cell_type] = base::max(
      base::min(interest_cell_type_iGRN[[cell_type]][interest_cell_type_iGRN[[cell_type]] > 0], na.rm = TRUE) - 0.005,
      0.005
    )
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Extracting tau for each cell type. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_tau)
}
