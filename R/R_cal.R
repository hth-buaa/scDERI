#' @title Calculate Normalized Regulatory Coefficient R
#'
#' @description Computes the activation regulatory coefficient used in gene expression dynamic equations.
#' Combines a linear term and a double-sigmoid saturation term, and normalizes by the value at theta = 1.
#'
#' @param theta_TF_TG Numeric vector. Regulatory strength of transcription factor on target gene.
#' @param tau Numeric between 0 and 1. Regulatory binarization threshold, also used to scale the linear term.
#' @param r1 Positive numeric. Hyperparameter controlling the scaling of the linear term. Generally greater than 100. Default is 200.
#' @param r2 Positive numeric. Hyperparameter controlling the overall steepness of the sigmoid saturation term.
#' @param r3 Positive numeric. Hyperparameter representing the inflection point position of the sigmoid curve. Default is 1.
#' @param r4 Positive numeric. Hyperparameter controlling the steepness of the second sigmoid segment. Generally greater than 100. Default is 1000.
#' @param r5 Positive numeric. Hyperparameter controlling the maximum amplitude of the saturation term. Default is 1.
#'
#' @return Numeric vector. Activation regulatory coefficient for scaling regulatory rates in dynamic equations.
#' With default parameters and theta in \[0, 1\], output values range from 0 to 1.
#'
#' @details The function follows a two-step calculation:
#' \enumerate{
#'   \item Computes the unnormalized value as the sum of a linear term and a multiplicative double-sigmoid term.
#'   \item Normalizes by the function value at theta = 1 to ensure consistent scaling.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' R_val <- R_cal(
#'   theta_TF_TG = c(0.2, 0.5, 0.8),
#'   tau = 0.3,
#'   r1 = 200,
#'   r2 = 1,
#'   r3 = 1,
#'   r4 = 1000,
#'   r5 = 1
#' )
#' }
R_cal = function(theta_TF_TG, tau, r1 = 200, r2, r3 = 1, r4 = 1000, r5 = 1)
{
  R_unnormalize = theta_TF_TG / (r1 * tau) + r5 / (
    (1 + 10^(r2 * (r3 - theta_TF_TG))) * (1 + 10^(r4 * r2 * (r3 - theta_TF_TG)))
  )
  r6 = 1 / (r1 * tau) + r5 / (
    (1 + 10^(r2 * (r3 - 1))) * (1 + 10^(r4 * r2 * (r3 - 1)))
  )
  R = R_unnormalize / r6
  return(R)
}
