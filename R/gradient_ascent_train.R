#' @title Gradient Ascent for Local Parameter Refinement
#'
#' @description Performs gradient ascent optimization to maximize the scDERI training objective.
#' Implements numerical gradient computation, adaptive learning rate selection via Brent's method,
#' parameter boundary clipping and early stopping for local parameter refinement.
#'
#' @param iterations Positive integer. Maximum number of gradient ascent iterations.
#' @param alpha_lower Tiny positive numeric. Lower bound of learning rate search range. Defaults to 1e-5.
#' @param alpha_upper Tiny positive numeric. Upper bound of learning rate search range. Defaults to 1e-3.
#' @param alpha_guess Numeric between \code{alpha_lower} and \code{alpha_upper}. Initial guess for learning rate optimization. Defaults to 1e-4.
#' @param init_theta Numeric vector in \[0, 1\]. Initial parameter values to start optimization from.
#' @param h_max Tiny positive numeric. Maximum step size for numerical gradient approximation. Defaults to 0.01.
#' @param h_min_percent Tiny positive numeric in \[0, 1\].
#' Minimum step size for numerical gradient, as a proportion of the parameter value range. Defaults to 0.01.
#' @param ncores Positive integer. Number of cores for parallel gradient computation via \code{\link[parallel]{mclapply}}.
#'   Defaults to \code{1}. Parallel processing is only supported on Unix-like systems;
#'   set to \code{1} when running on Windows.
#' @param params_lower Numeric vector. Lower bounds of model parameters.
#' @param params_upper Numeric vector. Upper bounds of model parameters.
#' @param monitor Logical. Whether to print iteration progress. Defaults to \code{FALSE}.
#' @param params_names Character vector. Names corresponding to each parameter element.
#' @param train_data Data frame. Training data for a single TF-TG pair.
#' @param tau Positive numeric. Regulatory binarization threshold.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths.
#' @param params_init Numeric vector. Initial parameter values for regularization.
#'
#' @return List with two elements:
#' \itemize{
#'   \item \code{parameters}: Optimized parameter vector that maximizes the objective function.
#'   \item \code{objective_value}: Maximum objective function value achieved.
#' }
#' Optimization stops immediately if the objective value does not improve in one iteration.
#'
#' @details Each iteration follows the workflow below:
#' \enumerate{
#'   \item Computes numerical gradients for all parameters in parallel via \code{\link{f_train_cal_gr}}.
#'   \item Skips parameters with zero gradient to avoid redundant updates.
#'   \item Searches for the optimal learning rate within the given range using Brent's one-dimensional method via \code{\link[stats]{optim}}.
#'   \item Updates parameters and clips values to stay within predefined bounds.
#'   \item Terminates early if the objective value does not increase.
#' }
#'
#' @note
#' \itemize{
#'   \item Uses forward finite difference approximation, which is computationally cheaper
#'     but less accurate than central difference.
#'   \item Designed for local refinement after global search (e.g. genetic algorithm),
#'     not for full training from random initialization.
#'   \item Parallel gradient computation relies on \code{\link[parallel]{mclapply}};
#'     set \code{ncores = 1} on Windows systems.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Refine parameters from GA solution
#' grad_result <- gradient_ascent_train(
#'   iterations = 30,
#'   init_theta = as.numeric(ga_result@solution[1, ]),
#'   params_lower = param_lower,
#'   params_upper = param_upper,
#'   params_names = param_names,
#'   train_data = train_data,
#'   tau = tau,
#'   Tslot_K = Tslot_K,
#'   params_init = param_init,
#'   ncores = 1
#' )
#' }
gradient_ascent_train = function (iterations,
                                  alpha_lower = 1e-05,
                                  alpha_upper = 1e-03,
                                  alpha_guess = 1e-04,
                                  init_theta,
                                  h_max = 0.01,
                                  h_min_percent = 0.01,
                                  ncores = 1,
                                  params_lower = params_lower,
                                  params_upper = params_upper,
                                  monitor = FALSE,
                                  params_names = params_names,
                                  train_data = train_data,
                                  tau = tau,
                                  Tslot_K = Tslot_K,
                                  params_init = params_init)
{
  theta = base::as.numeric(init_theta)
  f0 = scDERI::f_train_cal(
    params = theta,
    params_names = params_names,
    train_data = train_data,
    tau = tau,
    Tslot_K = Tslot_K,
    params_init = params_init,
    params_upper = params_upper,
    params_lower = params_lower
  )
  for (i in 1:iterations) {
    grad = scDERI::f_train_cal_gr(
      params = theta,
      h_max = h_max,
      h_min_percent = h_min_percent,
      ncores = ncores,
      params_names = params_names,
      train_data = train_data,
      tau = tau,
      Tslot_K = Tslot_K,
      params_upper = params_upper,
      params_lower = params_lower,
      params_init = params_init
    )
    grad_need_cal = grad != 0
    best_alpha = stats::optim(par = alpha_guess, fn = function(al) {
      theta_new = theta
      theta_new[grad_need_cal] = base::pmin(
        base::pmax(theta[grad_need_cal] + al * grad[grad_need_cal], params_lower[grad_need_cal]),
        params_upper[grad_need_cal]
      )
      -scDERI::f_train_cal(
        params = theta_new,
        params_names = params_names,
        train_data = train_data,
        tau = tau,
        Tslot_K = Tslot_K,
        params_init = params_init,
        params_upper = params_upper,
        params_lower = params_lower
      )
    }, lower = alpha_lower, upper = alpha_upper, method = "Brent")
    theta[grad_need_cal] = base::pmin(
      base::pmax(theta[grad_need_cal] + best_alpha$par * grad[grad_need_cal], params_lower[grad_need_cal]),
      params_upper[grad_need_cal]
    )
    current_f = -best_alpha$value
    if (f0 >= current_f) {
      if (monitor) {
        message("Grad | iter = ", i, " | Best = ", base::round(f0, 6), " | learning_rate = ", base::round(best_alpha$par, 6), " | ncores = ", ncores)
      }
      return(base::list(parameters = theta, objective_value = f0))
    } else {
      if (monitor) {
        message("Grad | iter = ", i, " | Best = ", base::round(current_f, 6), " | learning_rate = ", base::round(best_alpha$par, 6), " | ncores = ", ncores)
      }
      f0 = current_f
    }
  }
  base::list(parameters = theta, objective_value = f0)
}
