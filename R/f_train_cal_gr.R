#' @title Numerical Gradient Calculation for Model Training Objective Function
#'
#' @description Computes the numerical gradient of the scDERI training objective function
#' using forward finite difference approximation. Supports parallel computation of partial
#' derivatives for each model parameter, serving gradient-based optimization algorithms.
#'
#' @param params Numeric vector. Flattened model parameter vector for a single TF-TG pair.
#'   Order must match \code{params_names}.

#' @param h_max Tiny positive numeric. Maximum step size for finite difference approximation. Defaults to 0.01.
#' @param h_min_percent Tiny positive numeric. Minimum step size as a fraction of parameter range. Defaults to 0.01.
#' @param ncores Positive integer. Number of cores for parallel computation of partial derivatives. Defaults to 1.
#' @param params_names Character vector. Names corresponding to each parameter element.
#' @param train_data Data frame. Training data for one regulatory pair,
#' containing regulatory strength \code{theta_i} and time-series TFE, TGA, TGE features.
#' @param tau Positive numeric. Regulatory binarization threshold.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths.
#' @param params_upper Numeric vector. Upper bounds of parameters for step size calculation.
#' @param params_lower Numeric vector. Lower bounds of parameters for step size calculation.
#' @param params_init Numeric vector. Initial parameter values for regularization.
#'
#' @return Numeric vector of the same length as \code{params}. Partial derivatives of the objective
#' function with respect to each parameter, approximated by forward finite difference.
#' NA and infinite values are replaced with 0.
#'
#' @details The function uses an adaptive step size strategy:
#' \enumerate{
#'   \item For each parameter, the actual step size is the minimum of \code{h_max}
#'     and \code{h_min_percent} times the parameter range.
#'   \item Computes the baseline objective value once and reuses it for all partial derivatives.
#'   \item Computes each partial derivative independently in parallel via \code{parallel::mclapply}.
#'   \item Replaces NA and infinite gradient values with 0 for numerical stability.
#' }
#'
#' @note
#' \itemize{
#'   \item Internally calls \code{\link{f_train_cal}} to compute objective function values.
#'   \item Forward difference is computationally cheaper but less accurate than central difference.
#'   \item Parallel computation relies on \code{parallel::mclapply};
#'     set \code{ncores = 1} on Windows systems.
#' }
#'
#' @seealso \code{\link{f_train_cal}}, \code{\link{gradient_ascent_train}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Compute gradient at given parameter values
#' grad <- f_train_cal_gr(
#'   params = param_vec,
#'   h_max = 0.01,
#'   h_min_percent = 0.01,
#'   ncores = 1,
#'   params_names = param_names,
#'   train_data = train_data,
#'   tau = tau,
#'   Tslot_K = Tslot_K,
#'   params_upper = param_upper,
#'   params_lower = param_lower,
#'   params_init = param_init
#' )
#' }
f_train_cal_gr = function (params,
                           h_max = 0.01,
                           h_min_percent = 0.01,
                           ncores = 1,
                           params_names = params_names,
                           train_data = train_data,
                           tau = tau,
                           Tslot_K = Tslot_K,
                           params_upper = params_upper,
                           params_lower = params_lower,
                           params_init = params_init)
{
  vec = base::as.numeric(params)
  n = base::length(vec)
  f0 = scDERI::f_train_cal(
    params = vec,
    params_names = params_names,
    train_data = train_data,
    tau = tau,
    Tslot_K = Tslot_K,
    params_init = params_init,
    params_upper = params_upper,
    params_lower = params_lower
  )
  h_vec = base::pmin(base::rep(h_max, n), (params_upper - params_lower) * h_min_percent)
  grad = base::unlist(parallel::mclapply(1:n, function(i) {
    par_perturb = vec
    par_perturb[i] = vec[i] + h_vec[i]
    gr = (scDERI::f_train_cal(
      params = par_perturb,
      params_names = params_names,
      train_data = train_data,
      tau = tau,
      Tslot_K = Tslot_K,
      params_init = params_init,
      params_upper = params_upper,
      params_lower = params_lower
    ) - f0)/h_vec[i]
    rm(list = base::setdiff(base::ls(), "gr"))
    return(gr)
  }, mc.cores = ncores, mc.silent = TRUE, mc.allow.recursive = FALSE))
  grad[base::is.na(grad) | base::is.infinite(grad)] = 0
  return(grad)
}
