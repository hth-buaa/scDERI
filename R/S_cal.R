#' @title Calculate S-type Regulatory Rate Correction Term
#'
#' @description Computes the S-type correction term used in gene expression dynamic equations.
#' It calculates a nonlinear regulatory rate coefficient based on current and steady-state Hill binding levels.
#'
#' @param s Positive numeric. Shape parameter controlling the nonlinearity of the S function.
#' @param U Numeric between 0 and 1. Hill function output at the current time point,
#'  representing the fractional binding level.
#' @param U_tilde Numeric between 0 and 1. Steady-state Hill function output,
#' representing the equilibrium fractional binding level.
#'
#' @return Numeric. S-type correction term used to adjust regulatory efficiency in dynamic equations.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' S_val <- S_cal(
#'   s = s1,
#'   U = U_current,
#'   U_tilde = U_steady
#' )
#' }
S_cal = function(s, U, U_tilde)
{
  Um = 1 - U
  Um_tilde = 1 - U_tilde
  Numerator = s * (U * Um)^(s-1) * (U_tilde - U) * ((1 - U) * (1 - U_tilde) + U_tilde * U)
  Denominator = (U_tilde * Um_tilde)^(s+1) * exp(s * (2 * U_tilde - 1) * (U_tilde - U) / (U_tilde * Um_tilde))
  S = Numerator / Denominator
  return(S)
}
