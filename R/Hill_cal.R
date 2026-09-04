#' @title Calculate Hill Function Value
#'
#' @description Computes the Hill function widely used to describe molecular binding and gene regulatory kinetics.
#' The function returns the fractional saturation level corresponding to the input regulatory molecule concentration.
#'
#' @param x Positive numeric vector. Concentration or expression level of the regulatory molecule.
#' @param Dissociation_Constant Positive numeric.
#' The dissociation constant representing the concentration at which half-maximal binding occurs.
#' Smaller values indicate higher binding affinity.
#' @param Hill_Coefficient Positive numeric. The Hill coefficient representing the degree of cooperativity in binding.
#' Values greater than 1 indicate positive cooperativity, values equal to 1 indicate non-cooperative binding, and values less than 1 indicate negative cooperativity.
#'
#' @return Numeric vector with values between 0 and 1,
#' representing the fractional occupancy or binding probability of the regulatory molecule.
#' 0 indicates no binding and 1 indicates full saturation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' U <- Hill_cal(
#'   x = TGE_baseline,
#'   Dissociation_Constant = u11,
#'   Hill_Coefficient = u12
#' )
#' }
Hill_cal = function(x, Dissociation_Constant, Hill_Coefficient)
{
  x = base::as.numeric(x)
  Dissociation_Constant = base::as.numeric(Dissociation_Constant)
  Hill_Coefficient = base::as.numeric(Hill_Coefficient)
  U_scale = x ^ Hill_Coefficient
  U = U_scale / (Dissociation_Constant + U_scale)
  return(U)
}
