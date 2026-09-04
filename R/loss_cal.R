#' @title Calculate Pure Prediction Loss for a Single TF-TG Pair
#'
#' @description Computes weighted mean squared prediction error for a single TF-TG pair
#' using optimized model parameters. Contains only prediction error without regularization
#' or boundary penalty, for model evaluation purposes.
#'
#' @param data Data frame. Input data for one TF-TG pair, containing \code{theta_i}
#'   and time-series columns prefixed with \code{TFE_}, \code{TGA_}, \code{TGE_}.
#' @param params Numeric vector. Optimized model parameter vector.
#' @param params_names Character vector. Names corresponding to each parameter element.
#' @param tau Positive numeric. Regulatory binarization threshold.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths.
#'
#' @return Single numeric value. Weighted mean squared prediction error.
#' TFE and TGA errors have a weight of 2, TGE error has a weight of 1.
#' Lower values indicate better prediction performance.
#'
#' @details The function calculates pure prediction loss following these steps:
#' \enumerate{
#'   \item Extracts TFE, TGA and TGE time-series values from input data.
#'   \item Computes normalized regulatory coefficient R via \code{\link{R_cal}}.
#'   \item Computes steady-state and time-step Hill function values via \code{\link{Hill_cal}}.
#'   \item Computes regulatory efficiency values via \code{\link{S_cal}} for TFE and TGA dynamics.
#'   \item Calculates weighted squared prediction error for each time step and returns the mean.
#' }
#'
#' @note
#' \itemize{
#'   \item Differs from \code{\link{f_train_cal}}: this function excludes regularization
#'     and boundary penalty terms, and is used only for performance evaluation.
#'   \item Input data must contain columns with prefixes \code{TFE_}, \code{TGA_}, \code{TGE_}
#'     and a column named \code{theta_i}.
#'   \item Can be applied to training, validation or test datasets.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Evaluate fitting performance on training data
#' train_loss <- loss_cal(
#'   data = train_data,
#'   params = best_params,
#'   params_names = param_names,
#'   tau = tau,
#'   Tslot_K = Tslot_K
#' )
#' }
loss_cal = function (
    data = train_data,
    params = result_gradient$parameters,
    params_names = params_names,
    tau = tau,
    Tslot_K = Tslot_K
)
{
  base::names(params) = params_names

  TFE = base::as.numeric(data[, base::grep("TFE_", base::colnames(data), value = TRUE)])
  TGA = base::as.numeric(data[, base::grep("TGA_", base::colnames(data), value = TRUE)])
  TGE = base::as.numeric(data[, base::grep("TGE_", base::colnames(data), value = TRUE)])

  Length = base::length(TFE)
  seqLength = base::seq_len(Length - 1)

  R = base::as.numeric(scDERI::R_cal(
    theta_TF_TG = base::as.numeric(data[, "theta_i"]),
    tau = tau,
    r1 = 200,
    r2 = params["r.r2"],
    r3 = 1,
    r4 = 1000,
    r5 = 1
  ))

  alpha1 = params["alpha.alpha1"]
  alpha2 = params["alpha.alpha2"]

  beta1 = params["beta.beta1"]
  beta2 = params["beta.beta2"]
  beta3 = params["beta.beta3"]

  s1 = params["s.s1"]
  s2 = params["s.s2"]

  u11 = params["u.u11"]
  u21 = params["u.u21"]
  u311 = params["u.u311"]
  u321 = params["u.u321"]
  u331 = params["u.u331"]

  u12 = params["u.u12"]
  u22 = params["u.u22"]
  u312 = params["u.u312"]
  u322 = params["u.u322"]
  u332 = params["u.u332"]

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

  base::mean(base::sapply(2:Length, function(K) {
    prev_idx = K - 1
    2 * (TFE[prev_idx] + alpha1 * R * scDERI::S_cal(s = s1, U = U1[prev_idx], U_tilde = U1_tilde) + beta1 - TFE[K])^2 +
      2 * (TGA[prev_idx] + alpha2 * R * scDERI::S_cal(s = s2, U = U2[prev_idx], U_tilde = U2_tilde) + beta2 - TGA[K])^2 +
      (TGE[prev_idx] + R * (v1 * U31[prev_idx] - v2 * U32[prev_idx] - v3 * U33) * Tslot_K[K] + beta3 - TGE[K])^2
  }))
}
