#' @title Calculate Weighted Training Objective Function
#'
#' @description Computes the weighted objective function for scDERI model training, integrating
#' prediction error of expression dynamics, L2 regularization and boundary barrier penalty.
#' Returns negated and scaled values for compatibility with maximization-based optimizers such as genetic algorithms.
#'
#' @param params Numeric vector. Flattened model parameter vector for a single TF-TG pair.
#'   Order must match \code{params_names}.
#' @param params_names Character vector. Names corresponding to each element in \code{params}.
#' @param train_data Data frame. Training data for one regulatory pair, containing regulatory strength
#'   \code{theta_i} and time-series columns prefixed with \code{TFE_}, \code{TGA_}, \code{TGE_}.
#' @param tau Positive numeric. Regulatory binarization threshold, used to calculate regulatory coefficient R.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths between consecutive time points,
#'   used for TGE dynamics calculation.
#' @param params_init Numeric vector. Initial parameter values for regularization penalty. Order must match \code{params_names}.
#' @param params_upper Numeric vector. Upper bounds of model parameters. Order must match \code{params_names}.
#' @param params_lower Numeric vector. Lower bounds of model parameters. Order must match \code{params_names}.
#'
#' @return Single numeric value. Negated and scaled weighted objective function.
#' Higher values correspond to lower actual loss and better fitting performance.
#' The objective consists of:
#' \itemize{
#'   \item Weighted mean squared error of TFE, TGA and TGE dynamics (TFE/TGA weight = 2, TGE weight = 1)
#'   \item L2 regularization for beta coefficients and baseline expression parameters
#'   \item Boundary barrier penalty for all constrained parameters
#' }
#'
#' @details For a single TF-TG pair, the function executes the following steps:
#' \enumerate{
#'   \item Assigns names to the input parameter vector and extracts individual parameters.
#'   \item Computes normalized regulatory coefficient R via \code{\link{R_cal}}.
#'   \item Computes steady-state and time-step Hill function values via \code{\link{Hill_cal}}.
#'   \item Computes regulatory efficiency values via \code{\link{S_cal}} for TFE and TGA dynamics.
#'   \item Computes weighted squared prediction error for each time step.
#'   \item Adds regularization terms and boundary barrier penalties.
#'   \item Negates and scales the result by 10000 for optimizer compatibility.
#' }
#'
#' @note
#' \itemize{
#'   \item The return value is negated because the genetic algorithm performs maximization,
#'     while the actual loss is minimized.
#'   \item This function is called repeatedly during optimization and is computationally intensive.
#'   \item Regularization and boundary penalties improve generalization and keep parameters
#'     within biologically reasonable ranges.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Evaluate objective at given parameters
#' obj_value <- f_train_cal(
#'   params = param_vec,
#'   params_names = param_names,
#'   train_data = train_data,
#'   tau = tau,
#'   Tslot_K = Tslot_K,
#'   params_init = param_init,
#'   params_upper = param_upper,
#'   params_lower = param_lower
#' )
#' }
f_train_cal = function (params,
                        params_names = params_names,
                        train_data = train_data,
                        tau = tau,
                        Tslot_K = Tslot_K,
                        params_init = params_init,
                        params_upper = params_upper,
                        params_lower = params_lower)
{
  base::names(params) = params_names
  base::names(params_init) = params_names
  base::names(params_upper) = params_names
  base::names(params_lower) = params_names

  TFE = base::as.numeric(train_data[, base::grep("TFE_", base::colnames(train_data), value = TRUE)])
  TGA = base::as.numeric(train_data[, base::grep("TGA_", base::colnames(train_data), value = TRUE)])
  TGE = base::as.numeric(train_data[, base::grep("TGE_", base::colnames(train_data), value = TRUE)])

  Length = base::length(TFE)
  # seqLength = base::seq_len(Length - 1)

  r2 = params["r.r2"]
  r2_i = params_init["r.r2"]; r2_u = params_upper["r.r2"]; r2_l = params_lower["r.r2"]
  R = base::as.numeric(scDERI::R_cal(
    theta_TF_TG = base::as.numeric(train_data[, "theta_i"]),
    tau = tau, r1 = 200, r2 = r2, r3 = 1, r4 = 1000, r5 = 1
  ))

  a1 = params["alpha.alpha1"]
  a2 = params["alpha.alpha2"]
  a_i = params_init["alpha.alpha1"]; a_u = params_upper["alpha.alpha1"]; a_l = params_lower["alpha.alpha1"]

  b1 = params["beta.beta1"]
  b2 = params["beta.beta2"]
  b3 = params["beta.beta3"]

  s1 = params["s.s1"]
  s2 = params["s.s2"]
  s_i = params_init["s.s1"]; s_u = params_upper["s.s1"]; s_l = params_lower["s.s1"]

  u11 = params["u.u11"]
  u21 = params["u.u21"]
  u311 = params["u.u311"]
  u321 = params["u.u321"]
  u331 = params["u.u331"]
  u1_i = params_init["u.u11"]; u1_u = params_upper["u.u11"]; u1_l = params_lower["u.u11"]

  u12 = params["u.u12"]
  u22 = params["u.u22"]
  u312 = params["u.u312"]
  u322 = params["u.u322"]
  u332 = params["u.u332"]
  u2_i = params_init["u.u12"]; u2_u = params_upper["u.u12"]; u2_l = params_lower["u.u12"]

  TGE_tilde = params["E~.E~_TG"]
  TFE_tilde = params["E~.E~_TF"]

  U1_tilde = scDERI::Hill_cal(x = TGE_tilde, Dissociation_Constant = u11, Hill_Coefficient = u12)
  U2_tilde = scDERI::Hill_cal(x = TFE_tilde, Dissociation_Constant = u21, Hill_Coefficient = u22)
  U1 = scDERI::Hill_cal(x = TGE[-Length], Dissociation_Constant = u11, Hill_Coefficient = u12)
  U2 = scDERI::Hill_cal(x = TFE[-Length], Dissociation_Constant = u21, Hill_Coefficient = u22)
  U31 = scDERI::Hill_cal(x = TGA[-Length], Dissociation_Constant = u311, Hill_Coefficient = u312)
  U32 = scDERI::Hill_cal(x = TGE[-Length], Dissociation_Constant = u321, Hill_Coefficient = u322)
  U33 = scDERI::Hill_cal(x = TGE_tilde, Dissociation_Constant = u331, Hill_Coefficient = u332)

  v1 = params["v.v1"]
  v2 = params["v.v2"]
  v3 = params["v.v3"]
  v_i = params_init["v.v1"]; v_u = params_upper["v.v1"]; v_l = params_lower["v.v1"]

  a_c = a_l - a_u
  s_c = s_l - s_u
  u1_c = u1_l - u1_u
  u2_c = u2_l - u2_u
  v_c = v_l - v_u

  -10000 * (
    base::mean(base::sapply(2:Length, function(K) {
      prev_idx = K - 1
      2 * (TFE[prev_idx] + a1 * R * scDERI::S_cal(s = s1, U = U1[prev_idx], U_tilde = U1_tilde) + b1 - TFE[K])^2 +
        2 * (TGA[prev_idx] + a2 * R * scDERI::S_cal(s = s2, U = U2[prev_idx], U_tilde = U2_tilde) + b2 - TGA[K])^2 +
        (TGE[prev_idx] + R * (v1 * U31[prev_idx] - v2 * U32[prev_idx] - v3 * U33) * Tslot_K[K] + b3 - TGE[K])^2
    })) +
      b1^2 + b2^2 + b3^2 +
      0.5 * ((TGE_tilde - params_init["E~.E~_TG"])^2 + (TFE_tilde - params_init["E~.E~_TF"])^2) +
      # 0.000091 * (
      #   (r2_i - r2_l)^2 / (r2 - r2_l) + (r2_u - r2_i)^2 / (r2_u - r2) + r2_l - r2_u
      # ) +
      0.000100 * (
        (a_i - a_l)^2 / (a1 - a_l) + (a_u - a_i)^2 / (a_u - a1) + a_c +
          (a_i - a_l)^2 / (a2 - a_l) + (a_u - a_i)^2 / (a_u - a2) + a_c
      ) +
      0.000112 * (
        (s_i - s_l)^2 / (s1 - s_l) + (s_u - s_i)^2 / (s_u - s1) + s_c +
          (s_i - s_l)^2 / (s2 - s_l) + (s_u - s_i)^2 / (s_u - s2) + s_c
      ) +
      0.000110 * (
        (u1_i - u1_l)^2 / (u11 - u1_l) + (u1_u - u1_i)^2 / (u1_u - u11) + u1_c +
          (u1_i - u1_l)^2 / (u21 - u1_l) + (u1_u - u1_i)^2 / (u1_u - u21) + u1_c +
          (u1_i - u1_l)^2 / (u311 - u1_l) + (u1_u - u1_i)^2 / (u1_u - u311) + u1_c +
          (u1_i - u1_l)^2 / (u321 - u1_l) + (u1_u - u1_i)^2 / (u1_u - u321) + u1_c +
          (u1_i - u1_l)^2 / (u331 - u1_l) + (u1_u - u1_i)^2 / (u1_u - u331) + u1_c
      ) +
      0.000118 * (
        (u2_i - u2_l)^2 / (u12 - u2_l) + (u2_u - u2_i)^2 / (u2_u - u12) + u2_c +
          (u2_i - u2_l)^2 / (u22 - u2_l) + (u2_u - u2_i)^2 / (u2_u - u22) + u2_c +
          (u2_i - u2_l)^2 / (u312 - u2_l) + (u2_u - u2_i)^2 / (u2_u - u312) + u2_c +
          (u2_i - u2_l)^2 / (u322 - u2_l) + (u2_u - u2_i)^2 / (u2_u - u322) + u2_c +
          (u2_i - u2_l)^2 / (u332 - u2_l) + (u2_u - u2_i)^2 / (u2_u - u332) + u2_c
      ) +
      0.000116 * (
        (v_i - v_l)^2 / (v1 - v_l) + (v_u - v_i)^2 / (v_u - v1) + v_c +
          (v_i - v_l)^2 / (v2 - v_l) + (v_u - v_i)^2 / (v_u - v2) + v_c +
          (v_i - v_l)^2 / (v3 - v_l) + (v_u - v_i)^2 / (v_u - v3) + v_c
      )
  )
}
